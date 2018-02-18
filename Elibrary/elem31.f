c**********************************************************
      subroutine elem31  (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     filters the displacement data to enforce the incompressibility
c     projecting displacements while min. divergence
c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem
      parameter (mnode_elem = 27,mquad_elem=64)
      integer ielem, itask, nnode
      integer l, ninte, iinte, inode, jnode
      integer ievab, jevab, jdofn, idofn
      double precision xjaco, wtjac, tetfac
      double precision xforc 
      double precision alpha,dd(3),um(3),ident(3,3)
      double precision shap(4,mnode_elem)
      double precision sg(mquad_elem),tg(mquad_elem),zg(mquad_elem)
      double precision wg(mquad_elem)
      double precision pelem(*)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,*)
      double precision elemdata_nodal(nset_nodal,*)
      double precision uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
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

      ident(1:3,1:3) = 0.0d0
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
     
      do iinte = 1, ninte
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shap,
     $               xjaco,.false.,nnode,ndim,elnods,xelem)
        wtjac = tetfac*wg(iinte)*xjaco

c       construct the stiffness matrix 
        ievab = 0
        do inode = 1, nnode
          do idofn = 1,elemvec_ndofn(inode)
            ievab = ievab +1 
            jevab = 0
            do  jnode = 1, nnode
              do jdofn = 1,elemvec_ndofn(jnode)
                jevab = jevab+1
                estif(ievab,jevab) = estif(ievab,jevab)+
     $                      (dd(idofn)*dd(idofn)*ident(idofn,jdofn)*
     $                       shap(4,inode)*shap(4,jnode)+
     $                       alpha*shap(idofn,inode)*
     $                       shap(jdofn,jnode))*wtjac
              enddo!jdofn
            enddo!jnode
          enddo!idofn
        enddo!inode
      enddo!iinte

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

      if (nnode.eq.4) then
        call gausstet(l,ninte,sg,tg,zg,wg)
        tetfac= 0.166666666666667d0
      else
        call gauss3(l,ninte,sg,tg,zg,wg)
        tetfac = 1.0d0
      endif
c     
      do iinte = 1, ninte
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shap,
     $               xjaco,.false.,nnode,ndim,elnods,xelem)
        wtjac = tetfac*wg(iinte)*xjaco
c       calculate the local value of um 
        um(1:3) = 0.0d0
        do inode = 1,nnode
          do idofn = 1,elemvec_ndofn(inode)
            um(idofn) = um(idofn)+
     $                  shap(4,inode)*uelem_meas(idofn,inode)
          enddo
        enddo
         
        ievab = 0
        do inode = 1, nnode
          do idofn = 1,elemvec_ndofn(inode)
            ievab = ievab +1 
            eforc(ievab) = eforc(ievab)+ 
     $           dd(idofn)*dd(idofn)*um(idofn)*
     $           shap(4,inode)*wtjac
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
