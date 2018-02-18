c*********************************************************************
      subroutine elem306 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     Material model by Sevan Goenezen, written by Sevan Goenezen     
c     This element accounts for an additional force matching term in 
c     the objective function.
c     It is being used with elem305.f and adds additional contributions
c     arising from the additional penalty force matching term
c*********************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ielem, itask, nnode, nno, ndim, iinte, ninte
      integer i, j, k, q, r, l, inode, inter, iset, knode,linter
      integer ievab, idofn, mnode_elem, mquad_elem
      parameter (mnode_elem = 4,mquad_elem=20)
      integer elnods(mnode),elnods2d(3)
      double precision pelem(*), eforc(mevab), xelem(ndime,*)
      double precision shap(4,mnode_elem),shap2d(3,3),xel(2,3)
      double precision elemdata_nodal(nset_nodal,*), uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
      double precision elemdata_nodal2d(2,3),uelem2d(4,3),tmpCoor(3,3)
      double precision sg(mquad_elem),tg(mquad_elem), wg(mquad_elem)
      double precision Fdef(3,3), Ctens(3,3), dWdC(3,3), dJdC(3,3)
      double precision Cinv(3,3), SecPK(3,3), ident(3,3), Finv(3,3)
      double precision Ctang(3,3,3,3),d2JdC(3,3,3,3),estif(mevab,mevab)
      double precision diffX(3,2),Unorm(3),K2(3,3),egrad2d(2,3)
      double precision dCinvdC(3,3,3,3),dWdC_grad(3,3),SecPK_grad(3,3)
      double precision prop_grad(2,4),var_shape_prop(2,3)      
      double precision xjaco,wtjac,tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
      double precision tmp7,pres,dCinvdCt,Cdet,Fdet
      double precision gamm,mu,K1,Inv1,xi,force
      
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
c     xi    : penalty factor for force matching term
c     force : the total force to be matched
 2    read(iread,*) l, inter, xi, force
      write(iwrit,205) l, inter, xi, force
 205  format(/' finite elasticity (3D,4DOF) '/
     +       ' gauss pts/dir .........................',i12,/
     +       ' interior node .........................',i12,/
     +       ' penalty factor for force matching .....',1p,e16.4,/
     +       ' total force ...........................',1p,e16.4)
     
      pelem(1) = l
      pelem(2) = inter
      pelem(3) = xi
      pelem(4) = force
      return
      
 3    eforc(:)   = 0.0d0
      estif(:,:) = 0.0d0
c     primal solve does not change
      return
      
 4    l     = pelem(1)
      inter = pelem(2)
      
c     identity matrix ident
      ident(:,:) = 0.0d0
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(3,3) = 1.0d0     
      

cc   Computation of force matching needed for the residual dual
      
      
c    Compute deformation gradient for tet and use that later on for 
c    integration over surface. This is justified, because the deformation
c    gradient is a constant

        call shape3(ielem,0.0d0,0.0d0,0.0d0,shap,
     $               xjaco,.false.,nnode,ndime,elnods,xelem)
c       compute the deformation gradient at Gauss Point
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
c JFD: note Cdet=Fdet*Fdet...
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

      nno = 3 ! reduced number of nodes for triangular elements
      ndim = 2  ! reduce dimension for triangular elements
      
c     modify connectivity and element coordinates to 
c     reduce for triangular elements
      do i=1,4
        r=4
        if(elnods(i).eq.inter) then
          linter=i
          do j=1,(i-1)
            r=r-1
            elnods2d(r)=elnods(j)
            tmpCoor(1:3,r)=xelem(1:3,j)
            uelem2d(:,r)=uelem(:,j)
            elemdata_nodal2d(:,r)=elemdata_nodal(:,j)
          enddo
          do j=(i+1),4
            r=r-1
            elnods2d(r)=elnods(j)
            tmpCoor(1:3,r)=xelem(1:3,j)
            uelem2d(:,r)=uelem(:,j)
            elemdata_nodal2d(:,r)=elemdata_nodal(:,j)
          enddo
        go to 300
        endif
      enddo
      
 300  continue 
      
c     Compute the unit outward normal Unorm with crossproduct rule
      diffX(1:3,1)=tmpCoor(1:3,2)-tmpCoor(1:3,1)
      diffX(1:3,2)=tmpCoor(1:3,3)-tmpCoor(1:3,1)
      Unorm(1)=diffX(2,1)*diffX(3,2)-diffX(2,2)*diffX(3,1)
      Unorm(2)=diffX(3,1)*diffX(1,2)-diffX(3,2)*diffX(1,1)
      Unorm(3)=diffX(1,1)*diffX(2,2)-diffX(1,2)*diffX(2,1)
c     Check if this vector shows outwards 
      diffX(1:3,1)=xelem(1:3,linter)-tmpCoor(1:3,1)  ! i here was determined above
c      tmp6=Unorm(1)*diffX(1,1)+Unorm(2)*diffX(2,1)+Unorm(3)*diffX(3,1)
c      if(tmp6 .gt. 0.0d0) then
c        Unorm(:)=-Unorm(:)
c      endif
c     Norming to unit length one
        tmp6=sqrt(Unorm(1)**2.0d0+Unorm(2)**2.0d0+Unorm(3)**2.0d0)
        Unorm(:)=(1.0d0/tmp6)*Unorm(:)      
     
c     check along which axis compression is applied
c     and reduce coordinates to 2d xel appropriately
      if (Unorm(1).ne.0.0d0) then
        xel(1:2,:)=tmpCoor(2:3,:)
      elseif (Unorm(2).ne.0.0d0) then
        xel(1,:)=tmpCoor(1,:)
        xel(2,:)=tmpCoor(3,:)
      elseif (Unorm(3).ne.0.0d0) then
        xel(1:2,:)=tmpCoor(1:2,:)
      endif

      tmp6=0.0d0  !will be reused later on
      do i=1,3
        do j=1,3
          tmp6=tmp6+Unorm(i)*Cinv(i,j)*Unorm(j)
        end do
      end do
      tmp6=1.0d0/sqrt(tmp6)
      call gausstri(l,ninte,sg,tg,wg)
      
      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap2d,xjaco,
     $              .false.,nno,elnods2d,ndim,xel)
        wtjac = wg(iinte)*xjaco

        mu = 0.0d0! value of mu at the current Gauss point
        gamm = 0.0d0! value of gamma at the current Gauss point
        do inode = 1,nno
          mu = mu + shap2d(3,inode)*elemdata_nodal2d(2,inode)
          gamm = gamm + shap2d(3,inode)*elemdata_nodal2d(1,inode)
        enddo
c       compute the spatial derivatives of the pressure and the pressure at the Gauss Point
c        (shap2d contain the derivative of the shape functions and its value)
        pres = 0.0d0
          do inode = 1, nno
            pres = pres + uelem2d(4,inode)*shap2d(3,inode)
          end do
c     Compute the second Piola Kirchhoff stress SecPK
c     and the material tangent Ctang
        tmp5 = Fdet**(-2.0d0/3.0d0)! tmp5 is re-used here to avoid computing the power multiple times
        K1 = tmp5*Inv1-3.0d0
        tmp2 = 1.0d0/3.0d0 ! tmp2 is reused to avoid computing the fraction multiple times
        K2(1:3,1:3) = (ident(1:3,1:3) - 
     $                tmp2*Cinv(1:3,1:3)*Inv1)*tmp5
        tmp5=tmp5/3.0d0
        tmp4=0.5d0*mu*exp(K1*gamm)! tmp4 is re-used here to avoid computing the exp multiple times
        tmp3=0.25d0*Fdet! tmp3 is reused to avoid computing the product multiple times
        dWdC(1:3,1:3) = tmp4*K2(1:3,1:3)
        
        SecPK(1:3,1:3) = 2.0d0*(dWdC(1:3,1:3)-pres*dJdC(1:3,1:3))     
      
c    Computing an integral prefactor used after "enddo !iinte"
       do i=1,3
       do j=1,3
         forceM=forceM+tmp6*Unorm(i)*Unorm(j)*SecPK(i,j)*wtjac
       enddo
       enddo

      enddo !iinte
      
      return

 5    eforc(:)   = 0.0d0
c     no body forces
      return

c     -----------------------------------------------------------------
c     Build the elemental RHS arising from the penalty matching 
c     term for the dual problem
c     -----------------------------------------------------------------
 6    l     = pelem(1)
      inter = pelem(2)
      xi    = pelem(3)
      force = pelem(4)
      
      eforc(:) = 0.0d0
      
c     identity matrix ident
      ident(:,:) = 0.0d0
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(3,3) = 1.0d0     
        

c    Compute deformation gradient for tet and use that later on for 
c    integration over surface. This is justified, because the deformation
c    gradient is a constant

        call shape3(ielem,0.0d0,0.0d0,0.0d0,shap,
     $               xjaco,.false.,nnode,ndime,elnods,xelem)
c       compute the deformation gradient at Gauss Point
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
     $         Fdef(1,2)*Fdef(2,3)*Fdef(3,1) +
     $         Fdef(1,3)*Fdef(2,1)*Fdef(3,2) -
     $         Fdef(1,3)*Fdef(2,2)*Fdef(3,1) -
     $         Fdef(1,2)*Fdef(2,1)*Fdef(3,3) -
     $         Fdef(1,1)*Fdef(2,3)*Fdef(3,2)

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
c JFD: note Cdet=Fdet*Fdet...
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

        do i = 1, 3
          do j = 1, 3
            do q = 1, 3
              do r = 1, 3
                dCinvdC(i,j,q,r) = -0.5d0*(Cinv(i,q)*Cinv(j,r)+
     $                                   Cinv(i,r)*Cinv(j,q)) 
              end do
            end do
          end do
        end do

      nno = 3 ! reduced number of nodes for triangular elements
      ndim = 2  ! reduce dimension for triangular elements
      
c     modify connectivity and element coordinates to 
c     reduce for triangular elements
      do i=1,4
        r=4
        if(elnods(i).eq.inter) then
          linter=i
          do j=1,(i-1)
            r=r-1
            elnods2d(r)=elnods(j)
            tmpCoor(1:3,r)=xelem(1:3,j)
            uelem2d(:,r)=uelem(:,j)
            elemdata_nodal2d(:,r)=elemdata_nodal(:,j)
          enddo
          do j=(i+1),4
            r=r-1
            elnods2d(r)=elnods(j)
            tmpCoor(1:3,r)=xelem(1:3,j)
            uelem2d(:,r)=uelem(:,j)
            elemdata_nodal2d(:,r)=elemdata_nodal(:,j)
          enddo
          go to 100
        endif
      enddo
      
 100  continue 
      
c     Compute the unit outward normal Unorm with crossproduct rule
      diffX(1:3,1)=tmpCoor(1:3,2)-tmpCoor(1:3,1)
      diffX(1:3,2)=tmpCoor(1:3,3)-tmpCoor(1:3,1)
      Unorm(1)=diffX(2,1)*diffX(3,2)-diffX(2,2)*diffX(3,1)
      Unorm(2)=diffX(3,1)*diffX(1,2)-diffX(3,2)*diffX(1,1)
      Unorm(3)=diffX(1,1)*diffX(2,2)-diffX(1,2)*diffX(2,1)
c     Check if this vector shows outwards
      diffX(1:3,1)=xelem(1:3,linter)-tmpCoor(1:3,1)  ! i here was determined above
c      tmp6=Unorm(1)*diffX(1,1)+Unorm(2)*diffX(2,1)+Unorm(3)*diffX(3,1)
c      if(tmp6 .gt. 0.0d0) then
c        Unorm(:)=-Unorm(:)
c      endif
c     Norming to unit length one
        tmp6=sqrt(Unorm(1)**2.0d0+Unorm(2)**2.0d0+Unorm(3)**2.0d0)
        Unorm(:)=(1.0d0/tmp6)*Unorm(:) 

c     check along which axis compression is applied
c     and reduce coordinates to 2d xel appropriately
      if (Unorm(1).ne.0.0d0) then
        xel(1:2,:)=tmpCoor(2:3,:)
      elseif (Unorm(2).ne.0.0d0) then
        xel(1,:)=tmpCoor(1,:)
        xel(2,:)=tmpCoor(3,:)
      elseif (Unorm(3).ne.0.0d0) then
        xel(1:2,:)=tmpCoor(1:2,:)
      endif

      tmp6=0.0d0  !will be reused later on
      do i=1,3
        do j=1,3
          tmp6=tmp6+Unorm(i)*Cinv(i,j)*Unorm(j)
        end do
      end do
      tmp6=1.0d0/sqrt(tmp6)
      tmp7=tmp6**3.0d0

      call gausstri(l,ninte,sg,tg,wg)
      
      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap2d,xjaco,
     $              .false.,nno,elnods2d,ndim,xel)
        wtjac = wg(iinte)*xjaco

        mu = 0.0d0! value of mu at the current Gauss point
        gamm = 0.0d0! value of gamma at the current Gauss point
        do inode = 1,nno
          mu = mu + shap2d(3,inode)*elemdata_nodal2d(2,inode)
          gamm = gamm + shap2d(3,inode)*elemdata_nodal2d(1,inode)
        enddo

c       compute the spatial derivatives of the pressure and the pressure at the Gauss Point
c        (shap2d contain the derivative of the shape functions and its value)
        pres = 0.0d0
          do inode = 1, nno
            pres = pres + uelem2d(4,inode)*shap2d(3,inode)
          end do
  
c     Compute the second Piola Kirchhoff stress SecPK
c     and the material tangent Ctang
        tmp5 = Fdet**(-2.0d0/3.0d0)! tmp5 is re-used here to avoid computing the power multiple times
        K1 = tmp5*Inv1-3.0d0
        tmp2 = 1.0d0/3.0d0 ! tmp2 is reused to avoid computing the fraction multiple times
        K2(1:3,1:3) = (ident(1:3,1:3) - 
     $                tmp2*Cinv(1:3,1:3)*Inv1)*tmp5
        tmp5=tmp5/3.0d0
        tmp4=0.5d0*mu*exp(K1*gamm)! tmp4 is re-used here to avoid computing the exp multiple times
        tmp3=0.25d0*Fdet! tmp3 is reused to avoid computing the product multiple times
        dWdC(1:3,1:3) = tmp4*K2(1:3,1:3)
        
        SecPK(1:3,1:3) = 2.0d0*(dWdC(1:3,1:3)-pres*dJdC(1:3,1:3))

c        if (.false.) then
c        do i = 1, ndime
c          do j = 1, ndime
c            do q = 1, ndime
c              do r = 1, ndime
c                dCinvdC(i,j,q,r) = -0.5d0*(Cinv(i,q)*Cinv(j,r)+
c     $                                   Cinv(i,r)*Cinv(j,q)) 
c                Ctang(i,j,q,r) = tmp4*((gamm*K2(i,j)*K2(q,r))-
c     $                          (1.0d0/3.0d0)*K2(i,j)*Cinv(q,r)-
c     $                          tmp5*(dCinvdC(i,j,q,r)*Inv1
c     $                          +Cinv(i,j)*ident(q,r)))
c                d2JdC(i,j,q,r) = 0.25d0*Fdet*(Cinv(q,r)*Cinv(i,j)
c     $                           +2.0d0*dCinvdC(i,j,q,r))
c              end do
c            end do
c          end do
c        end do
c        else
cc        Ctangt(:,:,:,:)=Ctang(:,:,:,:)
cc        d2JdCt(:,:,:,:)=d2JdC(:,:,:,:)

c        Ctang and d2JdC have minor and major symmetries
c        that is each one has 21 independent components
c       i=1,j=1,q=1,r=1
        dCinvdCt=-Cinv(1,1)*Cinv(1,1)
        Ctang(1,1,1,1) = tmp4*(gamm*K2(1,1)*K2(1,1) - tmp2*
     $            K2(1,1)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(1,1)))
        d2JdC(1,1,1,1)=tmp3*dCinvdCt
c       i=2,j=1,q=1,r=1
        dCinvdCt=-Cinv(2,1)*Cinv(1,1)
        Ctang(2,1,1,1) = tmp4*(gamm*K2(2,1)*K2(1,1) - tmp2*
     $            K2(2,1)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(2,1)))
        d2JdC(2,1,1,1)=tmp3*dCinvdCt
        Ctang(1,2,1,1)=Ctang(2,1,1,1)
        Ctang(1,1,2,1)=Ctang(2,1,1,1)
        d2JdC(1,2,1,1)=d2JdC(2,1,1,1)
        d2JdC(1,1,2,1)=d2JdC(2,1,1,1)
c       i=3,j=1,q=1,r=1
        dCinvdCt=-Cinv(3,1)*Cinv(1,1)
        Ctang(3,1,1,1) = tmp4*(gamm*K2(3,1)*K2(1,1) - tmp2*
     $            K2(3,1)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(3,1)))
        d2JdC(3,1,1,1)=tmp3*dCinvdCt
        Ctang(1,3,1,1)=Ctang(3,1,1,1)
        Ctang(1,1,3,1)=Ctang(3,1,1,1)
        d2JdC(1,3,1,1)=d2JdC(3,1,1,1)
        d2JdC(1,1,3,1)=d2JdC(3,1,1,1)
c       i=2,j=2,q=1,r=1
        dCinvdCt=-Cinv(2,1)*Cinv(2,1)
        Ctang(2,2,1,1) = tmp4*(gamm*K2(2,2)*K2(1,1) - tmp2*
     $            K2(2,2)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(2,2)))
        d2JdC(2,2,1,1)=tmp3*(Cinv(1,1)*Cinv(2,2)+2.0d0*dCinvdCt)
        Ctang(1,1,2,2)=Ctang(2,2,1,1)
        d2JdC(1,1,2,2)=d2JdC(2,2,1,1)
c       i=3,j=2,q=1,r=1
        dCinvdCt=-Cinv(3,1)*Cinv(2,1)
        Ctang(3,2,1,1) = tmp4*(gamm*K2(3,2)*K2(1,1) - tmp2*
     $            K2(3,2)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(3,2)))
        d2JdC(3,2,1,1)=tmp3*(Cinv(1,1)*Cinv(3,2)+2.0d0*dCinvdCt)
        Ctang(2,3,1,1)=Ctang(3,2,1,1)
        Ctang(1,1,3,2)=Ctang(3,2,1,1)
        d2JdC(2,3,1,1)=d2JdC(3,2,1,1)
        d2JdC(1,1,3,2)=d2JdC(3,2,1,1)
c       i=3,j=3,q=1,r=1
        dCinvdCt=-Cinv(3,1)*Cinv(3,1)
        Ctang(3,3,1,1) = tmp4*(gamm*K2(3,3)*K2(1,1) - tmp2*
     $            K2(3,3)*Cinv(1,1) - tmp5*(dCinvdCt*Inv1+Cinv(3,3)))
        d2JdC(3,3,1,1)=tmp3*(Cinv(1,1)*Cinv(3,3)+2.0d0*dCinvdCt)
        Ctang(1,1,3,3)=Ctang(3,3,1,1)
        d2JdC(1,1,3,3)=d2JdC(3,3,1,1)
c       i=2,j=1,q=2,r=1
        dCinvdCt=-0.5d0*(Cinv(2,2)*Cinv(1,1)+Cinv(2,1)*Cinv(1,2))
        Ctang(2,1,2,1) = tmp4*(gamm*K2(2,1)*K2(2,1) - tmp2*
     $            K2(2,1)*Cinv(2,1) - tmp5*dCinvdCt*Inv1)
        d2JdC(2,1,2,1)=-tmp3*Cinv(2,2)*Cinv(1,1)
        Ctang(1,2,2,1)=Ctang(2,1,2,1)
        d2JdC(1,2,2,1)=d2JdC(2,1,2,1)
c       i=3,j=1,q=2,r=1
        dCinvdCt=-0.5d0*(Cinv(3,2)*Cinv(1,1)+Cinv(3,1)*Cinv(1,2))
        Ctang(3,1,2,1) = tmp4*(gamm*K2(3,1)*K2(2,1) - tmp2*
     $            K2(3,1)*Cinv(2,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,1,2,1)=-tmp3*Cinv(1,1)*Cinv(3,2)
        Ctang(1,3,2,1)=Ctang(3,1,2,1)
        Ctang(2,1,3,1)=Ctang(3,1,2,1)
        Ctang(1,2,3,1)=Ctang(3,1,2,1)
        d2JdC(1,3,2,1)=d2JdC(3,1,2,1)
        d2JdC(2,1,3,1)=d2JdC(3,1,2,1)
        d2JdC(1,2,3,1)=d2JdC(3,1,2,1)
c       i=2,j=2,q=2,r=1
        dCinvdCt=-Cinv(2,2)*Cinv(2,1)
        Ctang(2,2,2,1) = tmp4*(gamm*K2(2,2)*K2(2,1) - tmp2*
     $            K2(2,2)*Cinv(2,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(2,2,2,1)=tmp3*dCinvdCt
        Ctang(2,1,2,2)=Ctang(2,2,2,1)
        Ctang(1,2,2,2)=Ctang(2,2,2,1)
        d2JdC(2,1,2,2)=d2JdC(2,2,2,1)
        d2JdC(1,2,2,2)=d2JdC(2,2,2,1)
c       i=3,j=2,q=2,r=1
        dCinvdCt=-0.5d0*(Cinv(3,2)*Cinv(2,1)+Cinv(3,1)*Cinv(2,2))
        Ctang(3,2,2,1) = tmp4*(gamm*K2(3,2)*K2(2,1) - tmp2*
     $            K2(3,2)*Cinv(2,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,2,2,1)=-tmp3*Cinv(3,1)*Cinv(2,2)
        Ctang(2,3,2,1)=Ctang(3,2,2,1)
        Ctang(2,1,3,2)=Ctang(3,2,2,1)
        Ctang(1,2,3,2)=Ctang(3,2,2,1)
        d2JdC(2,3,2,1)=d2JdC(3,2,2,1)
        d2JdC(2,1,3,2)=d2JdC(3,2,2,1)
        d2JdC(1,2,3,2)=d2JdC(3,2,2,1)
c       i=3,j=3,q=2,r=1
        dCinvdCt=-Cinv(3,2)*Cinv(3,1)
        Ctang(3,3,2,1) = tmp4*(gamm*K2(3,3)*K2(2,1) - tmp2*
     $            K2(3,3)*Cinv(2,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,3,2,1)=tmp3*(Cinv(2,1)*Cinv(3,3)+2.0d0*dCinvdCt)
        Ctang(2,1,3,3)=Ctang(3,3,2,1)
        Ctang(1,2,3,3)=Ctang(3,3,2,1)
        d2JdC(2,1,3,3)=d2JdC(3,3,2,1)
        d2JdC(1,2,3,3)=d2JdC(3,3,2,1)
c       i=3,j=1,q=3,r=1
        dCinvdCt=-0.5d0*(Cinv(3,3)*Cinv(1,1)+Cinv(3,1)*Cinv(1,3))
        Ctang(3,1,3,1) = tmp4*(gamm*K2(3,1)*K2(3,1) - tmp2*
     $            K2(3,1)*Cinv(3,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,1,3,1)=-tmp3*Cinv(3,3)*Cinv(1,1)
        Ctang(1,3,3,1)=Ctang(3,1,3,1)
        d2JdC(1,3,3,1)=d2JdC(3,1,3,1)
c       i=2,j=2,q=3,r=1
        dCinvdCt=-Cinv(2,1)*Cinv(2,3)
        Ctang(2,2,3,1) = tmp4*(gamm*K2(2,2)*K2(3,1) - tmp2*
     $            K2(2,2)*Cinv(3,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(2,2,3,1)=tmp3*(Cinv(3,1)*Cinv(2,2)+2.0d0*dCinvdCt)
        Ctang(3,1,2,2)=Ctang(2,2,3,1)
        Ctang(1,3,2,2)=Ctang(2,2,3,1)
        d2JdC(3,1,2,2)=d2JdC(2,2,3,1)
        d2JdC(1,3,2,2)=d2JdC(2,2,3,1)
c       i=3,j=2,q=3,r=1
        dCinvdCt=-0.5d0*(Cinv(3,3)*Cinv(2,1)+Cinv(3,1)*Cinv(2,3))
        Ctang(3,2,3,1) = tmp4*(gamm*K2(3,2)*K2(3,1) - tmp2*
     $            K2(3,2)*Cinv(3,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,2,3,1)=-tmp3*Cinv(3,3)*Cinv(2,1)
        Ctang(2,3,3,1)=Ctang(3,2,3,1)
        Ctang(3,1,3,2)=Ctang(3,2,3,1)
        Ctang(1,3,3,2)=Ctang(3,2,3,1)
        d2JdC(2,3,3,1)=d2JdC(3,2,3,1)
        d2JdC(3,1,3,2)=d2JdC(3,2,3,1)
        d2JdC(1,3,3,2)=d2JdC(3,2,3,1)
c       i=3,j=3,q=3,r=1
        dCinvdCt=-Cinv(3,3)*Cinv(3,1)
        Ctang(3,3,3,1) = tmp4*(gamm*K2(3,3)*K2(3,1) - tmp2*
     $            K2(3,3)*Cinv(3,1)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,3,3,1)=tmp3*dCinvdCt
        Ctang(3,1,3,3)=Ctang(3,3,3,1)
        Ctang(1,3,3,3)=Ctang(3,3,3,1)
        d2JdC(3,1,3,3)=d2JdC(3,3,3,1)
        d2JdC(1,3,3,3)=d2JdC(3,3,3,1)
c       i=2,j=2,q=2,r=2
        dCinvdCt=-Cinv(2,2)*Cinv(2,2)
        Ctang(2,2,2,2) = tmp4*(gamm*K2(2,2)*K2(2,2) - tmp2*
     $            K2(2,2)*Cinv(2,2)-tmp5*(dCinvdCt*Inv1+Cinv(2,2)))
        d2JdC(2,2,2,2)=tmp3*dCinvdCt
c       i=3,j=2,q=2,r=2
        dCinvdCt=-Cinv(3,2)*Cinv(2,2)
        Ctang(3,2,2,2) = tmp4*(gamm*K2(3,2)*K2(2,2) - tmp2*
     $            K2(3,2)*Cinv(2,2)-tmp5*(dCinvdCt*Inv1+Cinv(3,2)))
        d2JdC(3,2,2,2)=tmp3*dCinvdCt
        Ctang(2,3,2,2)=Ctang(3,2,2,2)
        Ctang(2,2,3,2)=Ctang(3,2,2,2)
        d2JdC(2,3,2,2)=d2JdC(3,2,2,2)
        d2JdC(2,2,3,2)=d2JdC(3,2,2,2)
c       i=3,j=3,q=2,r=2
        dCinvdCt=-Cinv(3,2)*Cinv(3,2)
        Ctang(3,3,2,2) = tmp4*(gamm*K2(3,3)*K2(2,2) - tmp2*
     $            K2(3,3)*Cinv(2,2)-tmp5*(dCinvdCt*Inv1+Cinv(3,3)))
        d2JdC(3,3,2,2)=tmp3*(Cinv(2,2)*Cinv(3,3)+2.0d0*dCinvdCt)
        Ctang(2,2,3,3)=Ctang(3,3,2,2)
        d2JdC(2,2,3,3)=d2JdC(3,3,2,2)
c       i=3,j=2,q=3,r=2
        dCinvdCt=-0.5d0*(Cinv(3,3)*Cinv(2,2)+Cinv(3,2)*Cinv(3,2))
        Ctang(3,2,3,2) = tmp4*(gamm*K2(3,2)*K2(3,2) - tmp2*
     $            K2(3,2)*Cinv(3,2)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,2,3,2)=-tmp3*Cinv(3,3)*Cinv(2,2)
        Ctang(2,3,3,2)=Ctang(3,2,3,2)
        d2JdC(2,3,3,2)=d2JdC(3,2,3,2)
c       i=3,j=3,q=3,r=2
        dCinvdCt=-Cinv(3,3)*Cinv(3,2)
        Ctang(3,3,3,2) = tmp4*(gamm*K2(3,3)*K2(3,2) - tmp2*
     $            K2(3,3)*Cinv(3,2)-tmp5*dCinvdCt*Inv1)
        d2JdC(3,3,3,2)=tmp3*dCinvdCt
        Ctang(3,2,3,3)=Ctang(3,3,3,2)
        Ctang(2,3,3,3)=Ctang(3,3,3,2)
        d2JdC(3,2,3,3)=d2JdC(3,3,3,2)
        d2JdC(2,3,3,3)=d2JdC(3,3,3,2)
c       i=3,j=3,q=3,r=3
        dCinvdCt=-Cinv(3,3)*Cinv(3,3)
        Ctang(3,3,3,3) = tmp4*(gamm*K2(3,3)*K2(3,3) - tmp2*
     $            K2(3,3)*Cinv(3,3)-tmp5*(dCinvdCt*Inv1+Cinv(3,3)))
        d2JdC(3,3,3,3)=tmp3*dCinvdCt
c       enforce the symmetries
        Ctang(:,:,1,2)=Ctang(:,:,2,1)
        Ctang(:,:,1,3)=Ctang(:,:,3,1)
        Ctang(:,:,2,3)=Ctang(:,:,3,2)
        d2JdC(:,:,1,2)=d2JdC(:,:,2,1)
        d2JdC(:,:,1,3)=d2JdC(:,:,3,1)
        d2JdC(:,:,2,3)=d2JdC(:,:,3,2)

        ! build Ctang=Ctang-pres*d2JdC now
        Ctang(:,:,:,:)=4.0d0*(Ctang(:,:,:,:)-pres*d2JdC(:,:,:,:))

        ievab = 0
        r=4
        do inode = 1,4  !here all 4 nodes are considered
          do idofn = 1,4
            ievab = ievab+1
            if (idofn.ne.4) then
              do i=1,3
              do j=1,3
              do k=1,3
              do q=1,3
                eforc(ievab)=eforc(ievab)-tmp6*shap(q,inode)*
     $          Unorm(i)*Unorm(j)*Ctang(i,j,k,q)*Fdef(idofn,k)
     $          *wtjac  
              enddo
              enddo
              enddo
              enddo
            elseif (elnods(inode).ne.inter) then
              r=r-1                   !!!Check if assembly is ok
              do i=1,3
                do j=1,3
                  eforc(ievab)=eforc(ievab)+tmp6*shap2d(3,r)*
     $            Fdet*Unorm(i)*Unorm(j)*Cinv(i,j)*wtjac
                enddo
              enddo
            endif
          enddo
        enddo

        ievab = 0
        do inode = 1,4  !here all 4 nodes are considered
          do idofn = 1,4
            ievab = ievab+1
            if (idofn.ne.4) then
              do i=1,3
              do j=1,3
              do k=1,3
              do q=1,3
              do r=1,3
              do l=1,3
              eforc(ievab)=eforc(ievab)+tmp7*shap(r,inode)*
     $        Unorm(i)*Unorm(j)*SecPK(i,j)*Unorm(k)*Unorm(q)*
     $        dCinvdC(k,q,r,l)*Fdef(idofn,l)
     $        *wtjac
              enddo
            enddo
            enddo
            enddo
            enddo
            enddo
          endif
        enddo
      enddo

c    Computing an integral prefactor used after "enddo !iinte"

      enddo !iinte  

      eforc(:)=xi*(forceM-force)*eforc(:)

      return
      
 7    l     = pelem(1)
      inter = pelem(2)
      xi    = pelem(3)
      force = pelem(4)

      egrad(:,:) = 0.0d0! initialize egrad
      egrad2d(:,:) =0.0d0

c     identity matrix ident
      ident(:,:) = 0.0d0
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(3,3) = 1.0d0

        call shape3(ielem,0.0d0,0.0d0,0.0d0,shap,
     $               xjaco,.false.,nnode,ndime,elnods,xelem)

c       compute the deformation gradient at Gauss Point
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
     $         Fdef(1,2)*Fdef(2,3)*Fdef(3,1) +
     $         Fdef(1,3)*Fdef(2,1)*Fdef(3,2) -
     $         Fdef(1,3)*Fdef(2,2)*Fdef(3,1) -
     $         Fdef(1,2)*Fdef(2,1)*Fdef(3,3) -
     $         Fdef(1,1)*Fdef(2,3)*Fdef(3,2)

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
c JFD: note Cdet=Fdet*Fdet...
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

c for code check
c	if (ielem.eq. 7) then
c	Print*,"Fdef(1,1)",Fdef(1,1)
c	Print*,"Fdef(1,2)",Fdef(1,2)
c	Print*,"Fdef(1,3)",Fdef(1,3)
c	Print*,"Fdef(2,1)",Fdef(2,1)
c	Print*,"Fdef(2,2)",Fdef(2,2)
c	Print*,"Fdef(2,3)",Fdef(2,3)
c	Print*,"Fdef(3,1)",Fdef(3,1)
c	Print*,"Fdef(3,2)",Fdef(3,2)
c	Print*,"Fdef(3,3)",Fdef(3,3)
c	Print*,"Fdet     ",Fdet
c        Print*,"Ctens(1,1)",Ctens(1,1)
c	Print*,"Ctens(1,2)",Ctens(1,2)
c	Print*,"Ctens(1,3)",Ctens(1,3)
c	Print*,"Ctens(2,1)",Ctens(2,1)
c	Print*,"Ctens(2,2)",Ctens(2,2)
c	Print*,"Ctens(2,3)",Ctens(2,3)
c	Print*,"Ctens(3,1)",Ctens(3,1)
c	Print*,"Ctens(3,2)",Ctens(3,2)
c	Print*,"Ctens(3,3)",Ctens(3,3)
c	Print*,"Cinv(1,1)",Cinv(1,1)
c	Print*,"Cinv(1,2)",Cinv(1,2)
c	Print*,"Cinv(1,3)",Cinv(1,3)
c	Print*,"Cinv(2,1)",Cinv(2,1)
c	Print*,"Cinv(2,2)",Cinv(2,2)
c	Print*,"Cinv(2,3)",Cinv(2,3)
c	Print*,"Cinv(3,1)",Cinv(3,1)
c	Print*,"Cinv(3,2)",Cinv(3,2)
c	Print*,"Cinv(3,3)",Cinv(3,3)
c	endif
c  end code check

c       principal invariant Inv1
        Inv1 = Ctens(1,1)+Ctens(2,2)+Ctens(3,3)      
      
c       dJdC is the derivative of the Jacobian with respect to 
c       the Cauchy Green tensor
        dJdC(1:3,1:3) = 0.5d0*Fdet*Cinv(1:3,1:3)
      
c       create prop_grad (contains prop gradient. The value will be computed using
c       triangular shapefunctions later on) 
        prop_grad(:,:)  = 0.0d0
        do iset = 1,2
          do idofn = 1,3      !three gradients
            do inode = 1,nnode
               prop_grad(iset,idofn) = prop_grad(iset,idofn)+
     $                            shap(idofn,inode)*
     $                      elemdata_nodal(iset,inode)
            enddo
          enddo
        end do      

      nno = 3 ! reduced number of nodes for triangular elements
      ndim = 2  ! reduce dimension for triangular elements
      
c     modify connectivity and element coordinates to 
c     reduce for triangular elements
      do i=1,4
        r=4
        if(elnods(i).eq.inter) then
          linter=i
          do j=1,(i-1)
            r=r-1
            elnods2d(r)=elnods(j)
            tmpCoor(1:3,r)=xelem(1:3,j)
            uelem2d(:,r)=uelem(:,j)
            elemdata_nodal2d(:,r)=elemdata_nodal(:,j)
          enddo
          do j=(i+1),4
            r=r-1
            elnods2d(r)=elnods(j)
            tmpCoor(1:3,r)=xelem(1:3,j)
            uelem2d(:,r)=uelem(:,j)
            elemdata_nodal2d(:,r)=elemdata_nodal(:,j)
          enddo
          go to 200
        endif
      enddo
 200  continue
 
c  for code check
c       if (ielem.eq.6) then
c       Print*,"ielem",ielem
c       Print*,"elnods2d",elnods2d
c       Print*,"elnods",elnods(1:4)
c       Print*,"elemdata_nodal2d(2,:)",elemdata_nodal2d(2,:)
c       Print*,"elemdata_nodal(2,1:4)",elemdata_nodal(2,1:4)
c       endif
c       if (ielem.eq.7) then
c       Print*,"ielem",ielem
c       Print*,"elnods2d",elnods2d
c       Print*,"elnods",elnods(1:4)
c       Print*,"elemdata_nodal2d(2,:)",elemdata_nodal2d(2,:)
c       Print*,"elemdata_nodal(2,1:4)",elemdata_nodal(2,1:4)
c       stop
c       endif
c  end code check 
      
c     Compute the unit outward normal Unorm with crossproduct rule
      diffX(1:3,1)=tmpCoor(1:3,2)-tmpCoor(1:3,1)
      diffX(1:3,2)=tmpCoor(1:3,3)-tmpCoor(1:3,1)
      Unorm(1)=diffX(2,1)*diffX(3,2)-diffX(2,2)*diffX(3,1)
      Unorm(2)=diffX(3,1)*diffX(1,2)-diffX(3,2)*diffX(1,1)
      Unorm(3)=diffX(1,1)*diffX(2,2)-diffX(1,2)*diffX(2,1)
c     Check if this vector shows outwards 
      diffX(1:3,1)=xelem(1:3,linter)-tmpCoor(1:3,1)
c      tmp6=Unorm(1)*diffX(1,1)+Unorm(2)*diffX(2,1)+Unorm(3)*diffX(3,1)
c      if(tmp6 .gt.0.0d0) then
c        Unorm(:)=-Unorm(:)
c      endif     
c     Norming to unit length one
      tmp6=sqrt(Unorm(1)**2.0d0+Unorm(2)**2.0d0+Unorm(3)**2.0d0)
      Unorm(:)=(1.0d0/tmp6)*Unorm(:)    

c     this part checks along which axis compression is applied
c     and reduce coordinates to 2d xel appropriately
      if (Unorm(1).ne.0.0d0) then
        xel(1:2,:)=tmpCoor(2:3,:)
      elseif (Unorm(2).ne.0.0d0) then
        xel(1,:)=tmpCoor(1,:)
        xel(2,:)=tmpCoor(3,:)
      elseif (Unorm(3).ne.0.0d0) then
        xel(1:2,:)=tmpCoor(1:2,:)
      endif

      tmp6=0.0d0  !will be reused later on
      do i=1,3
        do j=1,3
          tmp6=tmp6+Unorm(i)*Cinv(i,j)*Unorm(j)
        end do
      end do
      tmp6=1.0d0/sqrt(tmp6)

      call gausstri(l,ninte,sg,tg,wg)
      
      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap2d,xjaco,
     $              .false.,nno,elnods2d,ndim,xel)
     
c  code check
c        if (ielem.eq.7) then
c        Print*,"elnods2d",elnods2d
c	Print*,"ndim    ",ndim
c	Print*,"nno     ",nno
c	Print*,"xjaco   ",xjaco
c	Print*,"ielem    ",ielem
c        Print*,"xel(:,1)",xel(:,1)
c	Print*,"xel(:,2)",xel(:,2)
c	Print*,"xel(:,3)",xel(:,3)
c	stop
c	endif
c  end code check  
     
        wtjac = wg(iinte)*xjaco

c       create prop_grad for the functional value 
        prop_grad(:,4)=0.0d0
        do iset = 1,2
          do inode = 1,nno
               prop_grad(iset,4) = prop_grad(iset,4)+
     $                            shap2d(3,inode)*
     $                      elemdata_nodal2d(iset,inode)
          enddo
        enddo 

c       compute the spatial derivatives of the pressure and the pressure at the Gauss Point
c        (shap2d contain the derivative of the shape functions and its value)
        pres = 0.0d0
          do inode = 1, nno
            pres = pres + uelem2d(4,inode)*shap2d(3,inode)
          end do
  
c     Compute the second Piola Kirchhoff stress SecPK
c     and the material tangent Ctang
        tmp5 = Fdet**(-2.0d0/3.0d0)! tmp5 is re-used here to avoid computing the power multiple times
        K1 = tmp5*Inv1-3.0d0
        tmp2 = 1.0d0/3.0d0 ! tmp2 is reused to avoid computing the fraction multiple times
        K2(1:3,1:3) = (ident(1:3,1:3) - 
     $                tmp2*Cinv(1:3,1:3)*Inv1)*tmp5
        tmp4=0.5d0*prop_grad(2,4)*exp(K1*prop_grad(1,4))! tmp4 is re-used here to avoid computing the exp multiple times
        dWdC(1:3,1:3) = tmp4*K2(1:3,1:3)
        
        SecPK(1:3,1:3) = 2.0d0*(dWdC(1:3,1:3)-pres*dJdC(1:3,1:3))

c for checking purpose
c        if (ielem.eq.7) then
c	Print*,"SecPK(1,1)",SecPK(1,1)
c	Print*,"SecPK(1,2)",SecPK(1,2)
c	Print*,"SecPK(1,3)",SecPK(1,3)
c	Print*,"SecPK(2,1)",SecPK(2,1)
c	Print*,"SecPK(2,2)",SecPK(2,2)
c	Print*,"SecPK(2,3)",SecPK(2,3)
c	Print*,"SecPK(3,1)",SecPK(3,1)
c	Print*,"SecPK(3,2)",SecPK(3,2)
c	Print*,"SecPK(3,3)",SecPK(3,3)
c	Print*,"2.0d0*dWdC(1,1)",2.0d0*dWdC(1,1)
c	Print*,"2.0d0*dWdC(1,2)",2.0d0*dWdC(1,2)
c	Print*,"2.0d0*dWdC(1,3)",2.0d0*dWdC(1,3)
c	Print*,"2.0d0*dWdC(2,1)",2.0d0*dWdC(2,1)
c	Print*,"2.0d0*dWdC(2,2)",2.0d0*dWdC(2,2)
c	Print*,"2.0d0*dWdC(2,3)",2.0d0*dWdC(2,3)
c	Print*,"2.0d0*dWdC(3,1)",2.0d0*dWdC(3,1)
c	Print*,"2.0d0*dWdC(3,2)",2.0d0*dWdC(3,2)
c	Print*,"2.0d0*dWdC(3,3)",2.0d0*dWdC(3,3)
c	Print*,"pres     ", pres
c	stop
c        endif
c  ending checking 
  
        do iset = 1,2
          do knode = 1,nno
c           set variations in the shape function for the property 
            var_shape_prop(:,:) = 0.0d0
            if (ielem_nod_grad(iset,knode).eq.1) then
              do idofn = 1,3
                var_shape_prop(iset,idofn) = shap2d(idofn,knode)
              enddo
            endif

c           compute derivative of second Piola Kirchhoff stress SecPK 
c           with respect to the material parameters
         
            K1 = (Fdet**(-2.0d0/3.0d0))*Inv1-3.0d0
            K2(1:3,1:3) = (ident(1:3,1:3) - 
     $                   (1.0d0/3.0d0)*Cinv(1:3,1:3)*Inv1)
     $                   *(Fdet**(-2.0d0/3.0d0))

            if (iset.eq.1) then
              dWdC_grad(1:3,1:3) = 0.5d0*prop_grad(2,4)*
     $                         K1*var_shape_prop(iset,3)*
     $                         exp(K1*prop_grad(iset,4))*K2(1:3,1:3)
            else if (iset.eq.2) then
              dWdC_grad(1:3,1:3) = 0.5d0*var_shape_prop(iset,3)
     $                         *exp(K1*prop_grad(1,4))*K2(1:3,1:3)
            end if

            SecPK_grad(1:3,1:3) = 2.0d0*dWdC_grad(1:3,1:3)
        
          do i=1,3
            do j =1,3
              egrad2d(iset,knode)=egrad2d(iset,knode)+tmp6*Unorm(i)*
     $           Unorm(j)*SecPK_grad(i,j)*wtjac
            enddo
          enddo
  
        enddo  !knode
      enddo    !iset


c  for code check
c           Print*,"iinte",iinte
c           Print*,"wtjac",wtjac
c	   Print*,"xjaco",xjaco
c end code check

      enddo  !iinte
      
c     Rearrange egrad2d into egrad
      do i=1,4
        r=4
        if(elnods(i).eq.inter) then
          linter=i
          do j=1,(i-1)
            r=r-1
            egrad(:,j)=egrad2d(:,r)
          enddo
          do j=(i+1),4
            r=r-1
            egrad(:,j)=egrad2d(:,r)
          enddo
          go to 250
        endif
      enddo
 250  continue     
      
      egrad(:,:) = xi*(forceM-force)*egrad(:,:)


c     Important to note that the input file has to be organized in such a way
c      that elem306 is the last element being read!!!
c      if (ielem.eq.nelem) then
c        forceMatch = forceMatch+0.5d0*xi*(forceM-force)**2.0d0
c        Print*,"forceMatch",forceMatch
c        stop
c      endif
      
c     Gxi and Gforce are the global variable to xi and force to be processed in lmain.f      
      Gxi = xi
      Gforce = force
      
      
      return

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
