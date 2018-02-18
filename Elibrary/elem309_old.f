c**********************************************************************
      subroutine elem309 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c*********************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ielem, itask, nnode, nno, ndim, iinte, ninte
      integer i, j, k, l, inode, jnode, inter, iset, knode,linter
      integer ievab, idofn, jevab,jdofn,mnode_elem, mquad_elem
      integer logflag, inode2d, jnode2d
      parameter (mnode_elem = 4,mquad_elem=20)
      integer elnods(mnode) 
      integer elnods2d(3), node_num_2d(4)
      double precision pelem(*), eforc(mevab), xelem(ndime,*)
      double precision tmpVec(3)
      double precision xCent(3),n2(3),n3(3),tmp2,dotProd ,xT(3,3)
      double precision shap(4,mnode_elem),shap2d(3,3),xel(2,3)
      double precision elemdata_nodal(nset_nodal,*), uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
      double precision tmpCoor(3,3)
      double precision sg(mquad_elem),tg(mquad_elem), wg(mquad_elem)
      double precision ident(3,3), Finv(3,3), Fdef(3,3)
      double precision estif(mevab,mevab)
      double precision diffX(3,2),Unorm(3),tmp6
      double precision xjaco,wtjac
      double precision pres,Fdet
      double precision fpress
c----------------------------------------------------------------------      

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask 

c     -----------------------------------------------------------------
c     Initialize the elemental parameters
c     -----------------------------------------------------------------
 1    elemvec_ndofn(1:nnode) = 4
      buildSymmetricMatrices=.false.
      return
      
c     -----------------------------------------------------------------
c     Read and store the material properties
c     -----------------------------------------------------------------
c     l     : no. of gauss pts/direction
c     inter : interior node
c     fpress : fluid pressure on each patch
c     logflag: flag for log approach(0/1:off,on)
 2    read(iread,*) l, inter, fpress, logflag
      write(iwrit,205) l, inter, fpress, logflag
 205  format(/' finite elasticity (3D,4DOF) '/
     +       ' gauss pts/dir .........................',i12,/
     +       ' interior node .........................',i12,/
     +       ' fluid pressure  .......................',1p,e16.4,/
     +       ' logflag ...............................',i12,/)
      pelem(1) = dble(l)
      pelem(2) = dble(inter)
      pelem(3) = fpress 
      pelem(4) = dble(logflag)
      return
      
c     -----------------------------------------------------------------
c     Build the elemental lhs and rhs arising from pressure term 
c     -----------------------------------------------------------------
 3    l     = int(pelem(1))
      inter = int(pelem(2))
      fpress = pelem(3)
      logflag = int(pelem(4))
      estif(:,:)=0.0d0
      eforc(:) = 0.0d0
      xCent(:) = 0.0d0
      elnods2d(:)=0 
      tmpCoor(:,:)=0.0d0
      node_num_2d(:)=0 
      ident(:,:) = 0.0d0
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(3,3) = 1.0d0
      Fdef(:,:)=0.0d0
      Finv(:,:)=0.0d0
      diffX(:,:)=0.0d0
      xel(:,:) =0.0d0
c----------------------------------------------------------     
c associating flag to each node of the tet , flag is 0 , if it is an interior node
c and 1,2,3 value if it is a triangular face node
      
      i=0
      do inode =1,nnode
       if (elnods(inode).eq.inter) then
        node_num_2d(inode)=0
                           ! interior node is flagged 0
       else
        i=i+1
        node_num_2d(inode)=i
       endif
      enddo

c    Compute deformation gradient for tet and use that later on for 
c    integration over surface. This is justified, because the deformation
c    gradient is a constant
c
       call shape3(ielem,0.0d0,0.0d0,0.0d0,shap,
     $               xjaco,.false.,nnode,ndime,elnods,xelem)
c    Compute the deformation gradient at Gauss Point
        Fdef(:,:) = ident(:,:)
        do inode = 1,nnode
          do j = 1,ndime
            do i = 1,ndime
              Fdef(i,j)=Fdef(i,j)+uelem(i,inode)*shap(j,inode)
            end do
          end do
        end do

c       Fdet is the Jacobian, determinant of the deformation gradient
       Fdet = Fdef(1,1)*Fdef(2,2)*Fdef(3,3) +
     $        Fdef(1,2)*Fdef(2,3)*Fdef(3,1) +
     $        Fdef(1,3)*Fdef(2,1)*Fdef(3,2) -
     $        Fdef(1,3)*Fdef(2,2)*Fdef(3,1) -
     $        Fdef(1,2)*Fdef(2,1)*Fdef(3,3) -
     $        Fdef(1,1)*Fdef(2,3)*Fdef(3,2)

c       Finv is inverse of deformation gradient at Gauss Point
        Finv(1,1) = Fdef(2,2)*Fdef(3,3)-Fdef(2,3)*Fdef(3,2)
        Finv(1,2) = Fdef(1,3)*Fdef(3,2)-Fdef(1,2)*Fdef(3,3)
        Finv(1,3) = Fdef(1,2)*Fdef(2,3)-Fdef(1,3)*Fdef(2,2)
        Finv(2,1) = Fdef(2,3)*Fdef(3,1)-Fdef(2,1)*Fdef(3,3)
        Finv(2,2) = Fdef(1,1)*Fdef(3,3)-Fdef(1,3)*Fdef(3,1)
        Finv(2,3) = Fdef(1,3)*Fdef(2,1)-Fdef(1,1)*Fdef(2,3)
        Finv(3,1) = Fdef(2,1)*Fdef(3,2)-Fdef(2,2)*Fdef(3,1)
        Finv(3,2) = Fdef(1,2)*Fdef(3,1)-Fdef(1,1)*Fdef(3,2)
        Finv(3,3) = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)
        Finv(:,:) = (1.0d0/Fdet)*Finv(:,:)

       nno = 3 
               ! reduced number of nodes for triangular elements
       ndim = 2
               ! reduce dimension for triangular elements

       do inode =1,nnode 
        inode2d = node_num_2d(inode)
            if(inode2d.ne.0)then 
                                 ! not an interior node,
                elnods2d(inode2d)=elnods(inode) 
                                             ! global node numbers of the surface nodes
                tmpCoor(1:ndime,inode2d)=xelem(1:ndime,inode)
                                                   ! global cordinates of the surface nodes 
            endif
       enddo

c Compute the unit outward normal Unorm with crossproduct rule
       diffX(1:3,1)=tmpCoor(1:3,2)-tmpCoor(1:3,1)
       diffX(1:3,2)=tmpCoor(1:3,3)-tmpCoor(1:3,1)
       Unorm(1)=diffX(2,1)*diffX(3,2)-diffX(2,2)*diffX(3,1)
       Unorm(2)=diffX(3,1)*diffX(1,2)-diffX(3,2)*diffX(1,1)
       Unorm(3)=diffX(1,1)*diffX(2,2)-diffX(1,2)*diffX(2,1)

c Norming to unit length one
       tmp6=sqrt(Unorm(1)**2.0d0+Unorm(2)**2.0d0+Unorm(3)**2.0d0)
       Unorm(:)=(1.0d0/tmp6)*Unorm(:) 

c More general way to check if the normal points outward
c Calculate the centroid of the triangle 
       xCent(:)=0.0d0
       xCent(:) = (1.0d0/3.0d0)*(tmpCoor(:,1)+tmpCoor(:,2)+tmpCoor(:,3))

c Calculate the vector pointing from centroid to the interior node
       tmpVec(:) = xelem(:,linter)-xCent(:)

c Compute dot product of the tmpVec with the Unorm vector
      dotProd = Unorm(1)*tmpVec(1)+Unorm(2)*tmpVec(2)+Unorm(3)*tmpVec(3)

c Check the direction of Unorm
       if(dotProd.gt.0.0d0) then ! pointing inwards
       Unorm(:) = -Unorm(:) 
      endif

c Compute other 2 inplane bases for the csys attached with the face (csys2);
c  whose one basis is Unorm, and the centroid of the face (xcent) is the origin    
      n2(:)=tmpCoor(:,1) - xCent(:)

c Normalize the basis vector
      tmp2= sqrt(n2(1)**2.0d0+n2(2)**2.0d0+n2(3)**2.0d0)
      n2(:)=(1.0d0/tmp2)*n2(:)

c Compute the 3rd basis of csys2, using cross product of two unit bases
      n3(1)= Unorm(2)*n2(3)-Unorm(3)*n2(2)
      n3(2)= Unorm(3)*n2(1)-Unorm(1)*n2(3)
      n3(3)= Unorm(1)*n2(2)-Unorm(2)*n2(1)

c Translate reference csys1 to the centroid of the triangle.
c Cordinates of the 3 face nodes w.r.t the translated csys1.
      do i=1,3
         xT(:,i) =tmpCoor(:,i)-xCent(:)
      enddo

c Compute cordinates of the 3 face nodes w.r.t csys2
c Unorm is out of plane basis for csys2, and n2 and n3 are inplane bases
      xel(:,:) =0.0d0
      do i = 1,nno
         do j = 1,ndime
            xel(1,i) =xel(1,i)+xT(j,i)*n2(j)
            xel(2,i) =xel(2,i)+xT(j,i)*n3(j)
         enddo
      enddo
      
      call gausstri(l,ninte,sg,tg,wg)

       do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap2d,xjaco,
     $              .false.,nno,elnods2d,ndim,xel)
        wtjac = wg(iinte)*xjaco
c compute the deformation gradient at Gauss Point
        ievab = 0
c Construct element estiffness matrix ( lhs) for primal problem        
        do inode =1,nnode
            do idofn = 1,4
               ievab =ievab +1 
               jevab =0
                  do jnode =1,nnode
                     do jdofn =1,4
                        jevab =jevab+1
                        inode2d = node_num_2d(inode)
                           if (inode2d.ne.0) then ! it's a surface node
                              if (idofn.ne.4 .AND.jdofn.ne.4) then
c P_f*J*w_{a}*Finv_{b,c}*Finv_{d,a}*[grad(Du)]_{c,b}*N_{d} -----> surface integral
                              do k = 1,3
                                do l = 1,3
                                  estif(ievab,jevab)
     $                               =  estif(ievab,jevab)
     $                              +fpress*Fdet*shap2d(3,inode2d)
     $                              *Finv(l,jdofn)*shap(l,jnode)
     $                              *Finv(k,idofn)*Unorm(k)*wtjac
                                enddo
                              ! l
                              enddo 
                            ! k
cc -P_f*J*w_{a}*Finv_{b,a}*Finv_{c,d}*[grad(Du)]_{d,b}*N_{c} ---> surface integral
                              do k = 1,3
                                do l = 1,3
                                 estif(ievab,jevab)
     $                               =  estif(ievab,jevab)
     $                              -fpress*Fdet*shap2d(3,inode2d)
     $                              *Finv(l,idofn)*shap(l,jnode)
     $                              *Finv(k,jdofn)*Unorm(k)*wtjac
                                enddo ! l
                              enddo
                           ! k
                            endif
                       ! if(idofn.eq.4....)
                       endif 
                       ! if ( inode2d.ne.0...
                    enddo
                 ! jdofn
                enddo
            ! jnode
            enddo
           ! idofn
        enddo
        ! inode
c construct element force vector( rhs ) for primal problem
c P_f*J*w_{a}*Finv_{b,a}*N_{b}  ---> surface integral
      ievab =0

       do inode =1,nnode
         do idofn = 1,4
            ievab = ievab +1
            inode2d = node_num_2d(inode)
               if (inode2d.ne.0) then
                  if(idofn.ne.4) then 
                     do k = 1,3
                        eforc(ievab) = eforc(ievab)
     $                  +fpress*Fdet*shap2d(3,inode2d)*Unorm(k)
     $                  *Finv(k,idofn)*wtjac
                     enddo
                  ! k
                  endif
             ! idofn
               endif
           ! (inode_2d ...)
            enddo
           ! idofn
         enddo
         ! inod
      enddo
      ! iinte
      return
 4    return
 5    return
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
