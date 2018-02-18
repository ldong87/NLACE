c**********************************************************
      subroutine elem32  (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     filters the displacement data to enforce the incompressibility
c     projecting displacements while min. divergence
c     can work on a structured non-regular grid (cf the datn section of the input file)
c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem
      parameter (mnode_elem = 27,mquad_elem=64)
      integer ielem, itask, nnode
      integer i, j, l, ninte, iinte, inode, jnode
      integer ievab, jevab, jdofn, idofn
      double precision xjaco, wtjac, tetfac
      double precision xforc 
      double precision alpha,dd(3),um(3),ident(3,3)
      double precision vect1(3), vect2(3)
      double precision shap(4,mnode_elem)
      double precision sg(mquad_elem),tg(mquad_elem),zg(mquad_elem)
      double precision wg(mquad_elem)
      double precision pelem(*)
      double precision estif(mevab,mevab)
      double precision estif2(mevab,mevab)
      double precision eforc(mevab)
c      double precision eforc2(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,*)
      double precision elemdata_nodal(nset_nodal,*)
      double precision elemdata_nodal2(3,3,nnode)
      double precision uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
      double precision uelem_meas2(mdofn,nnode)
      integer ndim != ndime (local copy)
c----------------------------------------------------------

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask
      
c     -----------------------------
c     initialize element parameters
c     -----------------------------
 1    elemvec_ndofn(1:nnode) = 3
      buildSymmetricMatrices=.false.
      return


c     --------------------------------------------------------
c     input material properties
c     --------------------------------------------------------
c     l     : no. of gauss pts/direction
 2    read(iread,*) l,alpha,dd(1),dd(2),dd(3)
      write(iwrit,205) l,alpha,dd(1),dd(2),dd(3)
 205  format(4x,'DISP. PROJ. W. INCOMP (3D,3DOF)'/
     +       4x,'gauss pts/dir                 ',i10  ,/
     +       4x,'alpha                         ',e15.6,/
     +       4x,'dd(1)                         ',e15.6,/
     +       4x,'dd(2)                         ',e15.6,/
     +       4x,'dd(3)                         ',e15.6)
      pelem(1) = dble(l)
      pelem(2) =alpha 
      pelem(3) = dd(1)
      pelem(4) = dd(2)
      pelem(5) = dd(3)
      return


c     ------------------------
c     form element stiffness
c     ------------------------
 3    l      = int(pelem(1))
      alpha  = pelem(2)
      dd(1)  = pelem(3)
      dd(2)  = pelem(4)
      dd(3)  = pelem(5)
      ndim = ndime

      ident(1:3,1:3) = 0.0d0! is actually a kronecker symbol (delta_{ij}) for convenience
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(3,3) = 1.0d0
      estif(:,:) = 0.0d0 ! elemental consistent tangent stiffness
      eforc(:) = 0.0d0! elemental RHS/residual
      
      if (nnode.eq.4) then
        call gausstet(l,ninte,sg,tg,zg,wg)
        tetfac= 0.166666666666667d0
      else
        call gauss3(l,ninte,sg,tg,zg,wg)
        tetfac = 1.0d0
      endif
     
c      print*,"elemdata_nodal: "
c      print*,elemdata_nodal(1,1),", ",elemdata_nodal(2,1),", ",
c     $       elemdata_nodal(3,1)
c      print*,elemdata_nodal(4,1),", ",elemdata_nodal(5,1),", ",
c     $       elemdata_nodal(6,1)
c      print*,elemdata_nodal(7,1),", ",elemdata_nodal(8,1),", ",
c     $       elemdata_nodal(9,1)
c     elemdata_nodal contains the coordinates of the curvilinear vectors in the cartesian grid
      ! JFD figure out the rotation is this one or the inverse...
      do idofn=1,3
        do jdofn=1,3
          do inode=1,nnode
            elemdata_nodal2(idofn,jdofn,inode)=
     &             elemdata_nodal((idofn-1)*3+jdofn,inode)
          enddo
        enddo
      enddo

      do iinte = 1, ninte
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shap,
     $               xjaco,.false.,nnode,ndim,elnods,xelem)
        wtjac = tetfac*wg(iinte)*xjaco

c       construct the stiffness matrix 
        ievab = 0
        do inode = 1, nnode
          do idofn = 1,elemvec_ndofn(inode)
c            Q(1,1,inode)=elemdata_nodal(1,inode)
c            Q(1,2,inode)=elemdata_nodal(2,inode)
c            Q(1,3,inode)=elemdata_nodal(3,inode)
c            Q(2,1,inode)=elemdata_nodal(4,inode)
c            Q(2,2,inode)=elemdata_nodal(5,inode)
c            Q(2,3,inode)=elemdata_nodal(6,inode)
c            Q(3,1,inode)=elemdata_nodal(7,inode)
c            Q(3,2,inode)=elemdata_nodal(8,inode)
c            Q(3,3,inode)=elemdata_nodal(9,inode)
            ievab = ievab +1 
            jevab = 0
            do  jnode = 1, nnode
              do jdofn = 1,elemvec_ndofn(jnode)
                jevab = jevab+1
c                estif(ievab,jevab) = estif(ievab,jevab)+
c     $                      (dd(idofn)*dd(idofn)*ident(idofn,jdofn)*
c     $                       shap(4,inode)*shap(4,jnode)+
c     $                       alpha*shap(idofn,inode)*
c     $                       shap(jdofn,jnode))*wtjac
                do i=1,3
                  vect1(i)=dd(i)*elemdata_nodal2(i,idofn,inode)*
     $                       shap(4,inode)
                  vect2(i)=dd(i)*elemdata_nodal2(i,jdofn,jnode)*
     $                       shap(4,jnode)
                 enddo
                 estif(ievab,jevab)=estif(ievab,jevab)+wtjac*
     $                   (vect1(1)*vect2(1)+vect1(2)*vect2(2)+
     $                    vect1(3)*vect2(3))
                
c                !do i=1,elemvec_ndofn(inode)
c                !  do j=1,elemvec_ndofn(jnode)
c                i=idofn
c                j=jdofn
c                    estif(ievab,jevab) = estif(ievab,jevab)+wtjac*
c     $                      (dd(i)*elemdata_nodal((i-1)*3+idofn,inode)*
c     $                       ident(i,j)*
c     $                       dd(j)*elemdata_nodal((j-1)*3+jdofn,jnode)*
c     $                       shap(4,inode)*shap(4,jnode))
c                !  enddo
c                !enddo
                estif(ievab,jevab) = estif(ievab,jevab) + wtjac*
     $              alpha*shap(idofn,inode)*shap(jdofn,jnode)! component coming from the div(u)^2 term
              enddo!jdofn
            enddo!jnode
          enddo!idofn
        enddo!inode
      enddo!iinte
c      estif2(:,:)=0.0d0
c      do inode=1,nnode
c        do jnode=1,nnode
c          do idofn=1,3
c            do jdofn=1,3
c              do i=1,3
c                estif2(idofn+(inode-1)*3,jdofn+(jnode-1)*3)=
c     &             estif2(idofn+(inode-1)*3,jdofn+(jnode-1)*3)+
c     &              estif(idofn+(inode-1)*3,i+(jnode-1)*3)*
c     &              elemdata_nodal2(i,jdofn,inode)! notranspose
c              enddo
c            enddo
c          enddo
c        enddo
c      enddo
c      estif(:,:)=estif2(:,:)
c      estif2(:,:)=0.0d0
c      do inode=1,nnode
c        do jnode=1,nnode
c          do idofn=1,3
c            do jdofn=1,3
c              do i=1,3
c                estif2(idofn+(inode-1)*3,jdofn+(jnode-1)*3)=
c     &             estif2(idofn+(inode-1)*3,jdofn+(jnode-1)*3)+
c     &              elemdata_nodal2(i,idofn,jnode)*! transpose
c     &              estif(i+(inode-1)*3,jdofn+(jnode-1)*3)
c              enddo
c            enddo
c          enddo
c        enddo
c      enddo
c      estif(:,:)=estif2(:,:)
c      print*,estif(1,1)," ",estif(1,2)," ",estif(1,3)," ",estif(1,4)
c      print*,estif(2,1)," ",estif(2,2)," ",estif(2,3)," ",estif(2,4)
c      print*,estif(3,1)," ",estif(3,2)," ",estif(3,3)," ",estif(3,4)
c      print*,estif(4,1)," ",estif(4,2)," ",estif(4,3)," ",estif(4,4)

c     compute element internal force vector (we use a non-linear solver and therefore want to compute only the increment - therefore we put K*U_{previous} in the rhs)
      ievab = 0
      do inode = 1,nnode
        do idofn = 1,elemvec_ndofn(inode)
          ievab = ievab+1
          jevab = 0
          do jnode = 1,nnode
            do jdofn = 1,elemvec_ndofn(jnode)
              jevab = jevab + 1
              xforc = estif(ievab,jevab)*uelem(jdofn,jnode)
              eforc(ievab) = eforc(ievab) + xforc
            enddo
          enddo
        enddo
      enddo 
      return

      
c     -------------------------------------
c     compute element volume force vector
c     -------------------------------------
 5    l      = int(pelem(1))
      dd(1)  = pelem(3)
      dd(2)  = pelem(4)
      dd(3)  = pelem(5)
      ndim = ndime

      eforc(:) = 0.0d0! elemental RHS/residual

      do idofn=1,3
        do jdofn=1,3
          do inode=1,nnode
            elemdata_nodal2(idofn,jdofn,inode)=
     &             elemdata_nodal((idofn-1)*3+jdofn,inode)
          enddo
        enddo
      enddo
c      uelem_meas2(:,1:nnode)=0.0d0
c      do inode=1,nnode
c        do idofn=1,3
c          do i=1,3
c            uelem_meas2(idofn,inode)=uelem_meas2(idofn,inode)+
cc     &           elemdata_nodal2(idofn,i,inode)*
c     &           elemdata_nodal2(i,idofn,inode)*! transpose
c     &           uelem_meas(i,inode)
c          enddo
c        enddo
c      enddo
c      uelem_meas(:,1:nnode)=uelem_meas2(:,1:nnode)

      if (nnode.eq.4) then
        call gausstet(l,ninte,sg,tg,zg,wg)
        tetfac= 0.166666666666667d0
      else
        call gauss3(l,ninte,sg,tg,zg,wg)
        tetfac = 1.0d0
      endif
     
      do iinte = 1, ninte
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shap,
     $               xjaco,.false.,nnode,ndim,elnods,xelem)
        wtjac = tetfac*wg(iinte)*xjaco
c       calculate the local value of um in which coordinate/axis system?
        um(1:3) = 0.0d0
        do inode = 1,nnode
          do idofn = 1,elemvec_ndofn(inode)
            um(idofn) = um(idofn)+
     $                  shap(4,inode)*uelem_meas(idofn,inode)! data at the node are not necessarily in the same coordinate system
          enddo
        enddo
         
        ievab = 0
        do inode = 1, nnode
          do idofn = 1,elemvec_ndofn(inode)
            ievab = ievab +1
            do i=1,3
              vect1(i)=dd(i)*elemdata_nodal2(i,idofn,inode)*
     $                       shap(4,inode)
              vect2(i)=dd(i)*um(i)
            enddo
            eforc(ievab) = eforc(ievab)+ wtjac*
     $           (vect1(1)*vect2(1)+vect1(2)*vect2(2)+
     $            vect1(3)*vect2(3))
c            eforc(ievab) = eforc(ievab)+ 
c     $           dd(idofn)*dd(idofn)*um(idofn)*
c     $           shap(4,inode)*wtjac
          enddo!idofn
        enddo!inode
      enddo!iinte
      return



 4    return
 6    return
 7    return
 8    return
 9    return
 10   return
 11   return
 12   return
 13   return
 14   return
 15   return
 16   return
 17   return
      end
