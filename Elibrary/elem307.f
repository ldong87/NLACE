c*********************************************************************
      subroutine elem307 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     Material model by Sevan Goenezen, written by Sevan Goenezen     
c     This element is a postprocessing element. It computes for the 
c     unknowns using force and moment information.
c     It is being used with elem305.f.
c*********************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ielem, itask, nnode, iinte, ninte
      integer i, j, k, q, r, l, inode, iset, knode
      integer idofn, mnode_elem, mquad_elem
      parameter (mnode_elem = 4,mquad_elem=30)
      integer elnods(mnode), node(4)
      double precision pelem(*), xelem(ndime,*), nw(4)
      double precision shape(4,mnode_elem),estif(mevab,mevab)
      double precision elemdata_nodal(nset_nodal,*), uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
      double precision sg(mquad_elem),tg(mquad_elem),zg(mquad_elem)
      double precision wg(mquad_elem)
      double precision Fdef(3,3), Ctens(3,3), dWdC(3,3), dJdC(3,3)
      double precision Cinv(3,3), SecPK(3,3), ident(3,3), Finv(3,3)
      double precision K2(3,3),FS(3,3),eforc(mevab)
      double precision xjaco,wtjac,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
      double precision tmp7,pres,Cdet,Fdet
      double precision gamm,mu,K1,Inv1,w(3)
      
c----------------------------------------------------------------------      

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask 

c     -----------------------------------------------------------------
c     Initialize the elemental parameters
c     -----------------------------------------------------------------
 1    elemvec_ndofn(1:nnode) = 4
      buildSymmetricMatrices = .false.
      return
      
c     -----------------------------------------------------------------
c     Read and store the material properties
c     -----------------------------------------------------------------
c     l         : no. of gauss pts/direction
c     node(1:4) : global node numbers, if 0 then interior node
c     ForceTotal     : the total force
c     MomentTotal    : the total moment
 2    read(iread,*) l,node(1),node(2),node(3),node(4),
     +             ForceTotal,MomentTotal
c      write(iwrit,205) l,node(1),node(2),node(3),node(4),
c     +             ForceTotal,MomentTotal
 205  format(/' finite elasticity (3D,4DOF) '/
     +       ' gauss pts/dir .......................',i12,/
     +       ' global node .........................',i12,/
     +       ' global node .........................',i12,/
     +       ' global node .........................',i12,/
     +       ' global node .........................',i12,/
     +       ' total force .........................',1p,e16.4,/
     +       ' total moment ........................',1p,e16.4,/)
     
      pelem(1) = dble(l)
      pelem(2) = dble(node(1))
      pelem(3) = dble(node(2))
      pelem(4) = dble(node(3))
      pelem(5) = dble(node(4))
      pelem(6) = ForceTotal
      pelem(7) = MomentTotal
      return
  
c     --------------------------------------------------------------------     
 3     eforc(:)   = 0.0d0
       estif(:,:) = 0.0d0
       return
c     --------------------------------------------------------------------      

 4     l           = int(pelem(1))
       node(1)     = int(pelem(2))
       node(2)     = int(pelem(3))
       node(3)     = int(pelem(4))
       node(4)     = int(pelem(5))
       ForceTotal  = pelem(6)
       MomentTotal = pelem(7)

c    Compute deformation gradient for tet at only one point. It is 
c    constant over the elemnt, so no need to recompute at every 
c    Gaus point!

c     identity matrix ident
      ident(:,:) = 0.0d0
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(3,3) = 1.0d0 

        call shape3(ielem,0.0d0,0.0d0,0.0d0,shape,
     $               xjaco,.false.,nnode,ndime,elnods,xelem)

c       compute the deformation gradient at Gauss Point
        Fdef(:,:) = ident(:,:)
        do inode = 1,nnode
          do j = 1,ndime
            do i = 1,ndime
              Fdef(i,j)=Fdef(i,j)+uelem(i,inode)*shape(j,inode)
            end do
          end do
        end do 
      
c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet = Fdef(1,1)*Fdef(2,2)*Fdef(3,3) +
     $         Fdef(1,2)*Fdef(2,3)*Fdef(3,1) +
     $         Fdef(1,3)*Fdef(2,1)*Fdef(3,2) -
     $         Fdef(1,3)*Fdef(2,2)*Fdef(3,1) -
     $         Fdef(1,2)*Fdef(2,1)*Fdef(3,3) -
     $         Fdef(1,1)*Fdef(2,3)*Fdef(3,2)

c       compute the Cauchy tensor and its inverse
        Ctens(1:3,1:3) = 0.0d0
        do j = 1,ndime
          do i = j,ndime! Ctens is symmetric
            do r = 1,ndime
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do
c       enforce the symmetries
        Ctens(1,2)=Ctens(2,1)
        Ctens(1,3)=Ctens(3,1)
        Ctens(2,3)=Ctens(3,2)

        Cdet = Ctens(1,1)*Ctens(2,2)*Ctens(3,3) +
     $         Ctens(1,2)*Ctens(2,3)*Ctens(3,1) +
     $         Ctens(1,3)*Ctens(2,1)*Ctens(3,2) -
     $         Ctens(1,3)*Ctens(2,2)*Ctens(3,1) -
     $         Ctens(1,2)*Ctens(2,1)*Ctens(3,3) -
     $         Ctens(1,1)*Ctens(2,3)*Ctens(3,2)

        Cinv(1,1) = Ctens(2,2)*Ctens(3,3)-Ctens(2,3)*Ctens(3,2)
        Cinv(1,2) = Ctens(1,3)*Ctens(3,2)-Ctens(1,2)*Ctens(3,3)
        Cinv(1,3) = Ctens(1,2)*Ctens(2,3)-Ctens(1,3)*Ctens(2,2)
        Cinv(2,1) = Cinv(1,2)! = Ctens(2,3)*Ctens(3,1)-Ctens(2,1)*Ctens(3,3)
        Cinv(2,2) = Ctens(1,1)*Ctens(3,3)-Ctens(1,3)*Ctens(3,1)
        Cinv(2,3) = Ctens(1,3)*Ctens(2,1)-Ctens(1,1)*Ctens(2,3)
        Cinv(3,1) = Cinv(1,3)! = Ctens(2,1)*Ctens(3,2)-Ctens(2,2)*Ctens(3,1)
        Cinv(3,2) = Cinv(2,3)! = Ctens(1,2)*Ctens(3,1)-Ctens(1,1)*Ctens(3,2)
        Cinv(3,3) = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)
        Cinv(:,:) = (1.0d0/Cdet)*Cinv(:,:)

c       principal invariant Inv1
        Inv1 = Ctens(1,1)+Ctens(2,2)+Ctens(3,3)
      
c       dJdC is the derivative of the Jacobian with respect to 
c       the Cauchy Green tensor
        dJdC(1:3,1:3) = 0.5d0*Fdet*Cinv(1:3,1:3)

        tmp5 = Fdet**(-2.0d0/3.0d0)
        K1 = tmp5*Inv1-3.0d0
        tmp2 = 1.0d0/3.0d0 
        K2(1:3,1:3) = (ident(1:3,1:3) - 
     $                tmp2*Cinv(1:3,1:3)*Inv1)*tmp5
        tmp5=tmp5/3.0d0
      
      call gausstet(l,ninte,sg,tg,zg,wg)

      do iinte = 1,ninte! for all Gauss integration points
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shape,
     $               xjaco,.false.,nnode,ndime,elnods,xelem)
        wtjac = wg(iinte)*xjaco

c       compute the pressure at the Gauss Point
        pres = 0.0d0
        do inode = 1, nnode
          pres = pres + uelem(4,inode)*shape(4,inode)
        end do

        mu = 0.0d0! value of mu at the current Gauss point
        gamm = 0.0d0! value of gamma at the current Gauss point
        do inode = 1,nnode
          mu = mu + shape(4,inode)*elemdata_nodal(2,inode)
          gamm = gamm + shape(4,inode)*elemdata_nodal(1,inode)
        enddo

c       compute the second Piola Kirchhoff stress SecPK 
        tmp4=0.5d0*mu*exp(K1*gamm)

c        Print*,"Fdef(1,1)",Fdef(1,1)
c       Print*,"Fdef(2,2)",Fdef(2,2)
c       Print*,"Fdef(3,3)",Fdef(3,3)
c       Print*,"Fdef(1,2)",Fdef(1,2)
c       Print*,"Fdef(1,3)",Fdef(1,3)
c       Print*,"Fdef(2,1)",Fdef(2,1)
c       Print*,"Fdef(2,3)",Fdef(2,3)
c       Print*,"Fdef(3,1)",Fdef(3,1)
c       Print*,"Fdef(3,2)",Fdef(3,2)
c       stop


        dWdC(1:3,1:3) = tmp4*K2(1:3,1:3)

        SecPK(1:3,1:3) = 2.0d0*(dWdC(1:3,1:3)-pres*dJdC(1:3,1:3))

        FS(:,:) = 0.0d0
        do i = 1, ndime
          do j = 1, ndime
            do r = 1, ndime
               FS(i,j) = FS(i,j)+Fdef(i,r)*SecPK(r,j)
            end do
          end do
        end do

c     To compute the integral, test functions will be set to zero at inner nodes,
        w(1:3) = 0.0d0     !derivatives 1:3

        do inode=1, nnode
          if (node(inode) .eq. 0) then
            nw(inode)=0.0d0
          else
            nw(inode)=1.0d0
          endif
        enddo

        do i=1, 3
          do inode=1, nnode
            w(i) = w(i) + shape(i,inode)*nw(inode)
          enddo 
        enddo
c        Print*, pforce
c        Print*,FS(3,1)
c        Print*,FS(3,2)
c        Print*,FS(3,3)
        do i=1, 3
c          pforce = pforce + w(i)*FS(3,i)*wtjac! for a force along the z-axis
          pforce = pforce + w(i)*FS(2,i)*wtjac! for a force along the y-axis
        enddo

c       Print*, pforce
c       stop
      end do! iinte
      

      return
c     -------------------------------------------------------------------------------
 5    eforc(:)   = 0.0d0
c     no body forces
      return
      
c     -------------------------------------------------------------------------------    
 6    eforc(:)   = 0.0d0
      estif(:,:) = 0.0d0
      return
      
c     --------------------------------------------------------------------------------
 7    eforc(:)   = 0.0d0
      estif(:,:) = 0.0d0
      return

c     --------------------------------------------------------------------------------        


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






