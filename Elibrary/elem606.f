c*************************************************************
      subroutine elem606 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     Written by Sevan Goenezen, and optimized by Jean-Francois     
c     2D (plane stress), quadrilaterals, finite elasticity, incompressible
c     The strain energy density function is Veronda-Westman
c*************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem
      parameter (mnode_elem = 9,mquad_elem=25)
      integer ielem, itask, ireg, iset
      integer nGauss, ninte, iinte, inode, jnode, knode
      integer ievab, jevab, jdofn, idofn, i, j, q, r
      double precision xjaco, wtjac
      double precision powe, deno, alpha(2), beta(2), Treg(2)! regularization variables
      double precision temp_dual(mnode_elem*2),prop_grad(2,3)
      double precision gamm, mu, Cdet, Inv1, Inv2, K1
      double precision CdetI, expVal! temporary variables to save time
      double precision shap(3,mnode_elem)
      double precision sg(mquad_elem),tg(mquad_elem),wg(mquad_elem)
      double precision Fdef(2,2), Ctens(2,2)
      double precision var_shape_prop(2,3)! 2 = nset
      double precision SecPK_grad(2,2)
      double precision Cinv(2,2), SecPK(2,2), ident(2,2)
      double precision dK1(2,2),dK2(2,2),Ctang(2,2,2,2)
      double precision d2K1t,d2K2t,Bt
      double precision Igeo(4,4),Lmat(2,2,2,2), Dgeomat(4,4)
      double precision T(4), FS(2,2), bb(mnode_elem,4,2),udiff(2)
      double precision regulv
c      double precision K2, xforc! unused variables
c      double precision d2K1(2,2,2,2),d2K2(2,2,2,2),B(2,2,2,2)! old variables
c     Localized information
      integer nnode
      double precision pelem(*)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,*)
      double precision elemdata_nodal(nset_nodal,*)
      double precision uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
      integer ndim! = ndime (local copy)
c      integer TID,OMP_GET_THREAD_NUM! uncomment for openmp
c-------------------------------------------------------------

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask
      
c     --------------------------------------------------------
c     Initialize the elemental parameters
c     --------------------------------------------------------
 1    elemvec_ndofn(1:nnode) = 2
c      buildSymmetricMatrices = .true.
      return


c     --------------------------------------------------------
c     Read and store the material properties
c     --------------------------------------------------------
c     nGauss     : no. of gauss pts/direction
c     ireg       : type of regularization (0/1/2/3/4:none/H1/TVD/newTV/power)
c     alpha(1)   : regularization weight for the non-linear parameter
c     alpha(2)   : regularization weight for the shear modulus
c     beta(1)    : extra parameter for TVD - JFD ask Sevan for a better description
c     beta(2)    : extra parameter for TVD - JFD idem
c     Treg(1)    : extra parameter for new TVD - JFD idem
c     Treg(2)    : extra parameter for new TVD - JFD idem
 2    read(iread,*) nGauss,ireg,alpha(1),alpha(2),beta(1),beta(2),
     +    Treg(1), Treg(2)
      write(iwrit,205) nGauss,ireg,alpha(1),alpha(2),beta(1),beta(2),
     +    Treg(1), Treg(2)
 205  format(/' finite elasticity (2D,plane stress approx.,2DOFs)'/
     +       ' gauss pts/dir .........................',i12,/
     +       ' reg type(0/1/2/4:none/H1/TVD/log, etc).',i12,/
     +       ' regularization parameter for gamma ....',1p,e16.4,/
     +       ' regularization parameter for mu .......',1p,e16.4,/
     +       ' extra parameter (TVD) .................',1p,e16.4,/! modify description JFD
     +       ' extra parameter (TVD) .................',1p,e16.4,/! idem JFD
     +       ' second extra parameter (TVD) ..........',1p,e16.4,/! modify description JFD
     +       ' second extra parameter (TVD) ..........',1p,e16.4)! idem JFD
c      nprop = 8
c     check for error in input file
      if (alpha(1).lt.0.0d0) then
        print*,'elem606.f: alpha(1) is less than 0.0d0: exiting'
        stop
      elseif (alpha(2).lt.0.0d0) then
        print*,'elem606.f: alpha(2) is less than 0.0d0: exiting'
        stop
      elseif (beta(1).le.0.0d0) then
        print*,'elem606.f: beta(1) is less or equal to 0.0d0: exiting'
        stop
      elseif (beta(2).le.0.0d0) then
        print*,'elem606.f: beta(2) is less or equal to 0.0d0: exiting'
        stop
      elseif (Treg(1).le.0.0d0) then
        print*,'elem606.f: Treg(1) is less or equal to 0.0d0: exiting'
        stop
      elseif (Treg(2).le.0.0d0) then
        print*,'elem606.f: Treg(2) is less or equal to 0.0d0: exiting'
        stop
      endif
c     note: ireg is checked in itask.eq.7 for keeping things simple
      pelem(1) = dble(nGauss)
      pelem(2) = dble(ireg)
      pelem(3) = alpha(1)
      pelem(4) = alpha(2)
      pelem(5) = beta(1)
      pelem(6) = beta(2)
      pelem(7) = Treg(1)
      pelem(8) = Treg(2)
      return


c     --------------------------------------------------------
c     Build the elemental consistent tangent stiffness matrix (estif)
c       and the elemental RHS (eforc)
c     --------------------------------------------------------
 3    nGauss = int(pelem(1))
c     initialize some variables
      estif(:,:) = 0.0d0! elemental tangent stiffness
      eforc(:) = 0.0d0! elemental RHS/residual
      ident(1,1) = 1.0d0! identity matrix
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0
      ndim = ndime
      if (nnode.eq.3) then
         call gausstri(nGauss,ninte,sg,tg,wg)
      else
         call gauss2(nGauss,ninte,sg,tg,wg)
      endif

cc      The following lines test openmp
c!$OMP PARALLEL PRIVATE(TID)
c      TID = OMP_GET_THREAD_NUM()
c      print*,'hello from thread ',TID
c!$OMP END PARALLEL

cc      note: this parallelization is rather slow but allows to test openmp 
cc            without allocatable variables if needed (last tested 2010-03-20;
cc            some variables may have been modified/added/removed since this date)
c!$OMP PARALLEL DO
c!$OMP& PRIVATE(iinte,shap,xjaco,wtjac,inode,mu,gamm,i,j,r,
c!$OMP&   Fdef,Ctens,Cdet,CdetI,Cinv,Inv1,Inv2,K1,dK1,dK2,expVal,
c!$OMP&   SecPK,Igeo,d2K1t,d2K2t,Ctang,Bt,q,Lmat,Dgeomat,bb,ievab,
c!$OMP&   idofn,jevab,jnode,jdofn,FS,T,TID)
c!$OMP& SHARED(ninte,ielem,sg,tg,nnode,elnods,ndim,xelem,wg,
c!$OMP&   elemdata_nodal,ident,uelem)
c!$OMP% DEFAULT(none)
cc!$OMP& SCHEDULE(STATIC,2)
c!$OMP& REDUCTION(+:estif,eforc)
      do iinte = 1,ninte! for all Gauss integration points
cc       The following 2 lines test openmp
c        TID = OMP_GET_THREAD_NUM()
c        print*,'hello from thread ',TID
        call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.false.,nnode,
     $      elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco
 
        mu = 0.0d0! value of mu at the Gauss point
        gamm = 0.0d0! value of gamma at the Gauss point
        do inode = 1,nnode
          gamm = gamm + shap(3,inode)*elemdata_nodal(1,inode)
          mu = mu + shap(3,inode)*elemdata_nodal(2,inode)
        enddo
         
c       compute the deformation gradient
        Fdef(:,:) = ident(:,:)
        do inode = 1,nnode
          do j = 1,ndim
            do i = 1,ndim
              Fdef(i,j)=Fdef(i,j)+uelem(i,inode)*shap(j,inode)
            end do
          end do
        end do

c       compute the Cauchy tensor and its inverse
        Ctens(1:2,1:2) = 0.0d0
        do j = 1,ndim
          do i = j,ndim! Ctens is symmetric
            do r = 1,ndim
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do
        Ctens(1,2)=Ctens(2,1)! Ctens is symmetric

        Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)
        CdetI = 1.0d0/Cdet

        Cinv(1,1) = CdetI*Ctens(2,2)
        Cinv(2,2) = CdetI*Ctens(1,1)
        Cinv(1,2) = -CdetI*Ctens(1,2)
        Cinv(2,1) = Cinv(1,2)! Cinv is symmetric

c       principal invariants Inv1 and Inv2
        Inv1 = Ctens(1,1)+Ctens(2,2)
        Inv2 = Cdet

c       compute the second Piola Kirchhoff stress SecPK and the material tangent Ctang
        K1 = Inv1+CdetI-3.0d0

        dK1(1:2,1:2) = ident(1:2,1:2)-(Cinv(1:2,1:2)*CdetI)
        dK2(1:2,1:2) = (ident(1:2,1:2)*CdetI) + 
     $                 (Inv2-(Inv1*CdetI))*Cinv(1:2,1:2)

        expVal=2.0d0*exp(K1*gamm)
        SecPK(1:2,1:2) = mu*(expVal*dK1(1:2,1:2)-dK2(1:2,1:2))

        Igeo(1,1) = SecPK(1,1)
        Igeo(1,2) = 0.0d0
        Igeo(1,3) = 0.5d0*SecPK(1,2)
        Igeo(1,4) = Igeo(1,3)! = 0.5d0*SecPK(1,2)
        Igeo(2,2) = SecPK(2,2)
        Igeo(2,3) = 0.5d0*SecPK(2,1)
        Igeo(2,4) = -Igeo(2,3)! = -0.5d0*SecPK(2,1)
        Igeo(3,3) = 0.25d0*(SecPK(2,2)+SecPK(1,1))
        Igeo(3,4) = 0.25d0*(SecPK(2,2)-SecPK(1,1))
        Igeo(4,4) = Igeo(3,3)! = 0.25d0*(SecPK(2,2)+SecPK(1,1))
c       enforce the symmetries
        Igeo(2,1) = Igeo(1,2)
        Igeo(3,1) = Igeo(1,3)
        Igeo(4,1) = Igeo(1,4)
        Igeo(3,2) = Igeo(2,3)
        Igeo(4,2) = Igeo(2,4)
        Igeo(4,3) = Igeo(3,4)

c       The following commented lines are the unoptimized code for computing Ctang but
c        Ctang has both minor and major symmetries (that is 6 independent components)
cc        do r = 1, 2
cc          do q = 1, 2
cc            do j = 1, 2
cc              do i = 1, 2
cc               B(i,j,q,r) = -0.5d0*(Cinv(i,q)*Cinv(j,r) + 
cc     $                      Cinv(i,r)*Cinv(j,q)) 
cc               d2K1(i,j,q,r) = (1.0d0/Inv2)*(Cinv(i,j)*Cinv(q,r)-
cc     $                          B(i,j,q,r))
cc               d2K2(i,j,q,r) = (Inv2+(Inv1/Inv2))*Cinv(i,j)*
cc     $                           Cinv(q,r)-(1.0d0/Inv2)*(Cinv(i,j)*
cc     $                           ident(q,r)+Cinv(q,r)*ident(i,j))+
cc     $                           (Inv2-(Inv1/Inv2))*B(i,j,q,r)
cc               Ctang(i,j,q,r) = 2.0d0*mu*(2.0d0*exp(gamm*K1)* 
cc     $                      (gamm*dK1(i,j)*dK1(q,r)+
cc     $                      d2K1(i,j,q,r))-d2K2(i,j,q,r))
cc              end do
cc            end do
cc          end do
cc        end do
c
c        if (.false.) then! use only the minor symmetries
c        do r = 1, 2
c          do q = r, 2! use minor symmetry
c            do j = 1, 2
c              do i = j, 2! use minor symmetry
c               Bt = -0.5d0*(Cinv(i,q)*Cinv(j,r) + 
c     $                      Cinv(i,r)*Cinv(j,q)) 
c               d2K1t = (1.0d0/Inv2)*(Cinv(i,j)*Cinv(q,r)-Bt)
c               d2K2t = (Inv2+(Inv1/Inv2))*Cinv(i,j)*
c     $                           Cinv(q,r)-(1.0d0/Inv2)*(Cinv(i,j)*
c     $                           ident(q,r)+Cinv(q,r)*ident(i,j))+
c     $                           (Inv2-(Inv1/Inv2))*Bt
c               Ctang(i,j,q,r) = 2.0d0*mu*(2.0d0*exp(gamm*K1)* 
c     $                      (gamm*dK1(i,j)*dK1(q,r)+
c     $                      d2K1t)-d2K2t)
c              end do
c            end do
c          end do
c        end do
c       enforce the minor symmetries
c        Ctang(1,2,1,1)=Ctang(2,1,1,1)
c        Ctang(1,2,2,1)=Ctang(2,1,2,1)
cc        Ctang(1,2,1,2)=Ctang(2,1,1,2)
c        Ctang(1,2,2,2)=Ctang(2,1,2,2)
c        Ctang(1,1,1,2)=Ctang(1,1,2,1)
c        Ctang(2,1,1,2)=Ctang(2,1,2,1)
c        Ctang(1,2,1,2)=Ctang(2,1,2,1)
c        Ctang(2,2,1,2)=Ctang(2,2,2,1)
c        else

c       r=1,q=1,j=1,i=1
        d2K1t = 2.0d0*CdetI*Cinv(1,1)*Cinv(1,1)
        d2K2t = Inv1*d2K1t-2.0d0*CdetI*Cinv(1,1)
        Ctang(1,1,1,1) = 2.0d0*mu*(expVal*(gamm*dK1(1,1)*dK1(1,1)+
     $               d2K1t)-d2K2t)
c       r=1,q=1,j=1,i=2
        d2K1t = 2.0d0*CdetI*Cinv(2,1)*Cinv(1,1)
        d2K2t = Inv1*d2K1t-CdetI*Cinv(2,1)
        Ctang(2,1,1,1) = 2.0d0*mu*(expVal*(gamm*dK1(2,1)*dK1(1,1)+
     $               d2K1t)-d2K2t)
        Ctang(1,2,1,1) = Ctang(2,1,1,1)
        Ctang(1,1,2,1) = Ctang(2,1,1,1)
        Ctang(1,1,1,2) = Ctang(2,1,1,1)
c       r=1,q=1,j=2,i=2
        Bt = -Cinv(2,1)*Cinv(2,1)
        d2K1t = CdetI*(Cinv(2,2)*Cinv(1,1)-Bt)
        d2K2t = (Inv2+(Inv1*CdetI))*Cinv(2,2)*Cinv(1,1)-
     $           CdetI*(Cinv(2,2)+Cinv(1,1))+
     $           (Inv2-(Inv1*CdetI))*Bt
        Ctang(2,2,1,1) = 2.0d0*mu*(expVal*(gamm*dK1(2,2)*dK1(1,1)+
     $               d2K1t)-d2K2t)
        Ctang(1,1,2,2) = Ctang(2,2,1,1)
c       i=2,q=2,j=1,r=1
        Bt = -0.5d0*(Cinv(2,2)*Cinv(1,1)+Cinv(2,1)*Cinv(1,2)) 
        d2K1t = CdetI*(Cinv(2,1)*Cinv(2,1)-Bt)
        d2K2t = (Inv2+(Inv1*CdetI))*Cinv(2,1)*Cinv(2,1)+
     $           (Inv2-(Inv1*CdetI))*Bt
        Ctang(2,1,2,1) = 2.0d0*mu*(expVal*(gamm*dK1(2,1)*dK1(2,1)+
     $               d2K1t)-d2K2t)
        Ctang(2,1,1,2) = Ctang(2,1,2,1)
        Ctang(1,2,1,2) = Ctang(2,1,2,1)
        Ctang(1,2,2,1) = Ctang(2,1,2,1)
c       i=2,q=2,j=1,r=2
        d2K1t = 2.0d0*CdetI*Cinv(2,1)*Cinv(2,2)
        d2K2t = Inv1*d2K1t-CdetI*Cinv(2,1)
        Ctang(2,1,2,2) = 2.0d0*mu*(expVal*(gamm*dK1(2,1)*dK1(2,2)+
     $               d2K1t)-d2K2t)
        Ctang(1,2,2,2) = Ctang(2,1,2,2)
        Ctang(2,2,1,2) = Ctang(2,1,2,2)
        Ctang(2,2,2,1) = Ctang(2,1,2,2)
c       i=2,q=2,j=2,r=2
        d2K1t = 2.0d0*CdetI*Cinv(2,2)*Cinv(2,2)
        d2K2t = Inv1*d2K1t-2.0d0*CdetI*Cinv(2,2)
        Ctang(2,2,2,2) = 2.0d0*mu*(expVal*(gamm*dK1(2,2)*dK1(2,2)+
     $               d2K1t)-d2K2t)
c        endif
        
c       Lmat has major symmetry
        do r = 1, ndim
          do q = 1, ndim
            do j = r, ndim! use r instead of 1
              do i = q, ndim! use q instead of 1
                Lmat(i,j,q,r) = Fdef(i,1)*Fdef(q,1)*Ctang(j,1,r,1)+
     $                          Fdef(i,1)*Fdef(q,2)*Ctang(j,1,r,2)+
     $                          Fdef(i,2)*Fdef(q,1)*Ctang(j,2,r,1)+
     $                          Fdef(i,2)*Fdef(q,2)*Ctang(j,2,r,2)
              end do
            end do
          end do
        end do
c       one term is missing: i=1, j=2, q=2, r=1
        Lmat(1,2,2,1) = Fdef(1,1)*Fdef(2,1)*Ctang(2,1,1,1)+
     $                  Fdef(1,1)*Fdef(2,2)*Ctang(2,1,1,2)+
     $                  Fdef(1,2)*Fdef(2,1)*Ctang(2,2,1,1)+
     $                  Fdef(1,2)*Fdef(2,2)*Ctang(2,2,1,2)
c       enforce the symmetries - not needed here (and therefore the lines are left commented)
c        Lmat(1,1,2,1) = Lmat(2,1,1,1)
c        Lmat(1,1,1,2) = Lmat(1,2,1,1)
c        Lmat(2,1,1,2) = Lmat(1,2,2,1)
c        Lmat(1,1,2,2) = Lmat(2,2,1,1)
c        Lmat(2,1,2,2) = Lmat(2,2,2,1)
c        Lmat(1,2,2,2) = Lmat(2,2,1,2)

        Dgeomat(1,1) = Lmat(1,1,1,1)
        Dgeomat(1,2) = Lmat(2,2,1,1)! = Lmat(1,1,2,2)
        Dgeomat(1,3) = 0.5d0*(Lmat(1,2,1,1)+Lmat(2,1,1,1))! = 0.5d0*(Lmat(1,1,1,2)+Lmat(1,1,2,1))
        Dgeomat(1,4) = 0.5d0*(Lmat(1,2,1,1)-Lmat(2,1,1,1))! = 0.5d0*(Lmat(1,1,1,2)-Lmat(1,1,2,1))
        Dgeomat(2,2) = Lmat(2,2,2,2)
        Dgeomat(2,3) = 0.5d0*(Lmat(2,2,1,2)+Lmat(2,2,2,1))
        Dgeomat(2,4) = 0.5d0*(Lmat(2,2,1,2)-Lmat(2,2,2,1))
        Dgeomat(3,3) = 0.25d0*Lmat(1,2,1,2)+0.5d0*Lmat(1,2,2,1)
     $                 +0.25d0*Lmat(2,1,2,1)
        Dgeomat(3,4) = 0.25d0*(Lmat(1,2,1,2)-Lmat(2,1,2,1))
        Dgeomat(4,4) = 0.25d0*Lmat(1,2,1,2)+0.25d0*Lmat(2,1,2,1)
     $                 -0.5d0*Lmat(1,2,2,1)
c       enforce the symmetries
        Dgeomat(2,1) = Dgeomat(1,2)
        Dgeomat(3,1) = Dgeomat(1,3)
        Dgeomat(3,2) = Dgeomat(2,3)
        Dgeomat(4,1) = Dgeomat(1,4)
        Dgeomat(4,2) = Dgeomat(2,4)
        Dgeomat(4,3) = Dgeomat(3,4)

        Dgeomat(1:4,1:4) = Dgeomat(1:4,1:4)+Igeo(1:4,1:4)

c       create the b-matrix
        do inode = 1, nnode
          bb(inode,1,1) = shap(1,inode)
          bb(inode,2,1) = 0.0d0
          bb(inode,3,1) = shap(2,inode)
          bb(inode,4,1) = shap(2,inode)
          bb(inode,1,2) = 0.0d0
          bb(inode,2,2) = shap(2,inode)
          bb(inode,3,2) = shap(1,inode)
          bb(inode,4,2) = -shap(1,inode)
        end do  

c       build the tangent stiffness matrix: estif is symmetric
        ievab = 0
        do inode = 1, nnode
          do idofn = 1,2! 2 = elemvec_ndofn(inode)
            ievab = ievab+1
            jevab = 2*(inode-1)! jevab=0 when not using the symmetries
            do jnode = inode, nnode! inode is used instead of 1
              do jdofn = 1,2! 2 = elemvec_ndofn(jnode)
                jevab = jevab+1
                do j = 1, 4
                  do i = 1, 4
                    estif(ievab,jevab)=estif(ievab,jevab)+
     $                   bb(inode,i,idofn)*Dgeomat(i,j)*
     $                   bb(jnode,j,jdofn)*wtjac
                  end do
                end do
              end do! jdofn
            end do! jnode
          end do! idofn
        end do! inode
c       note: the symmetric part of estif is filled after the Gauss point loop ends

c       create element residual for right hand side
        do inode = 1, nnode
c          bb(inode,1,1) = shap(1,inode)! was previously initialized so
c          bb(inode,2,1) = 0.0d0
c          bb(inode,3,1) = shap(2,inode)
          bb(inode,4,1) = 0.0d0
c          bb(inode,1,2) = 0.0d0
c          bb(inode,2,2) = shap(2,inode)
          bb(inode,3,2) = 0.0d0
          bb(inode,4,2) = shap(1,inode)
        end do

        FS(:,:) = 0.0d0
        do j = 1, ndim
          do r = 1, ndim
            do i = 1, ndim
              FS(i,j) = FS(i,j)+Fdef(i,r)*SecPK(r,j)
            end do
          end do
        end do

        T(1) = FS(1,1)! is FS usefull? maybe create T immediately JFD
        T(2) = FS(2,2)
        T(3) = FS(1,2)
        T(4) = FS(2,1)
         
        ievab = 0
        do inode = 1, nnode
          do i = 1,2! 2 = elemvec_ndofn(inode)
            ievab = ievab+1
            do j = 1,4
              eforc(ievab)=eforc(ievab)+bb(inode,j,i)*T(j)*wtjac
            end do
          end do
        end do! inode
     
      end do! iinte
c!$OMP END PARALLEL DO
      
c     symmetrise estif
      do i=1,2*nnode
        do j=i,2*nnode
          estif(j,i)=estif(i,j)
        enddo
      enddo
c      estif(3,1)=estif(1,3)
c      estif(4,1)=estif(1,4)
c      estif(5,1)=estif(1,5)
c      estif(6,1)=estif(1,6)
c      estif(7,1)=estif(1,7)
c      estif(8,1)=estif(1,8)
c      estif(3,2)=estif(2,3)
c      estif(4,2)=estif(2,4)
c      estif(5,2)=estif(2,5)
c      estif(6,2)=estif(2,6)
c      estif(7,2)=estif(2,7)
c      estif(8,2)=estif(2,8)
c      estif(5,3)=estif(3,5)
c      estif(6,3)=estif(3,6)
c      estif(7,3)=estif(3,7)
c      estif(8,3)=estif(3,8)
c      estif(5,4)=estif(4,5)
c      estif(6,4)=estif(4,6)
c      estif(7,4)=estif(4,7)
c      estif(8,4)=estif(4,8)
c      estif(7,5)=estif(5,7)
c      estif(8,5)=estif(5,8)
c      estif(7,6)=estif(6,7)
c      estif(8,6)=estif(6,8)
      
      return


c     --------------------------------------------------------
c     Evaluate the elemental RHS/residual (eforc) from a source term
c     --------------------------------------------------------
 5    continue
      eforc(:) = 0.0d0! elemental RHS/residual
c     turned off for now
      return


c     --------------------------------------------------------
c     Build the elemental RHS/residual (eforc) for the dual problem
c     --------------------------------------------------------
 6    nGauss = int(pelem(1))
      ndim = ndime

      eforc(:) = 0.0d0! elemental RHS/residual

      if (nnode.eq.3) then
         call gausstri(nGauss,ninte,sg,tg,wg)
      else
        call gauss2(nGauss,ninte,sg,tg,wg)
      endif

      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.false.,nnode,
     $      elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco

        ievab = 0
        do inode = 1, nnode
          do idofn = 1,2 ! 2 = elemvec_ndofn(inode)
            ievab = ievab+1
            do jnode = 1,nnode
c              do jdofn = 1,2 ! 2 = elemvec_ndofn(jnode)! useless do-loop
c                if (idofn.eq.jdofn) then
c                  eforc(ievab) = eforc(ievab)-
c     $                       shap(3,inode)*shap(3,jnode)
c     $                       *uelem_diff(jdofn,jnode)*wtjac
c                endif
c              enddo! jdofn
              jdofn = idofn
              eforc(ievab) = eforc(ievab) - shap(3,inode)*shap(3,jnode)
     $                       *uelem_diff(jdofn,jnode)*wtjac
            enddo! jnode
          enddo! idofn
        enddo! inode
      enddo! iinte

      return 
      

c     --------------------------------------------------------
c     Compute the objective function (dataMatch+regularization)
c      and its gradient (egrad) on the element
c     --------------------------------------------------------
 7    nGauss    = int(pelem(1))
      ireg      = int(pelem(2))
      alpha(1)  = pelem(3)
      alpha(2)  = pelem(4)
      beta(1)   = pelem(5)
      beta(2)   = pelem(6)
      Treg(1)   = pelem(7)
      Treg(2)   = pelem(8)

      powe = 0.5d0! power value for the power regularization (ireg.eq.4) - JFD should be in the input file?

      elemDataMatch(ielem) = 0.0d0
      elemRegul(ielem) = 0.0d0
      l2grad1(ielem) = 0.0d0
      l2grad2(ielem) = 0.0d0
      egrad(:,:) = 0.0d0! initialize egrad

      temp_dual(:) = 0.0d0
      ident(1,1)=1.0d0
      ident(2,2)=1.0d0
      ident(1,2)=0.0d0
      ident(2,1)=0.0d0
      ndim = ndime
      
      ievab=0
      do inode = 1,nnode
        do idofn = 1,2! 2 = elemvec_ndofn(inode)
          ievab = ievab+1
          temp_dual(ievab) = uelem_dual(idofn,inode)
        enddo
      enddo

      if (nnode.eq.3) then
        call gausstri(nGauss,ninte,sg,tg,wg)
      else
        call gauss2(nGauss,ninte,sg,tg,wg)
      endif

      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.false.,nnode,
     $      elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco

c       create prop_grad (contains prop gradient and value) 
        prop_grad(:,:)  = 0.0d0
        do inode = 1,nnode
          do idofn = 1,3! 3 = two gradients+the value
            do iset = 1,2
              prop_grad(iset,idofn) = prop_grad(iset,idofn)+
     $                    shap(idofn,inode)*elemdata_nodal(iset,inode)
            enddo
          enddo
        end do

        udiff(:) = 0.0d0

c       compute the deformation gradient at Gauss Point
        Fdef(1,1) = 1.0d0
        Fdef(2,2) = 1.0d0
        Fdef(1,2) = 0.0d0
        Fdef(2,1) = 0.0d0 
        do inode = 1,nnode
          do j = 1,ndim
            do i = 1,ndim
              Fdef(i,j)=Fdef(i,j)+uelem(i,inode)*shap(j,inode)
            end do
          end do
        end do
      
c       compute the Cauchy tensor and its inverse
        Ctens(1:2,1:2) = 0.0d0
        do j = 1,ndim
          do i = j,ndim! Ctens is symmetric
            do r = 1,ndim
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do
        Ctens(1,2)=Ctens(2,1)! Ctens is symmetric

        Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)
        CdetI = 1.0d0/Cdet

        Cinv(1,1) = CdetI*Ctens(2,2)
        Cinv(2,2) = CdetI*Ctens(1,1)
        Cinv(1,2) = -CdetI*Ctens(1,2)
        Cinv(2,1) = Cinv(1,2)

c       principal invariants Inv1 and Inv2
        Inv1 = Ctens(1,1)+Ctens(2,2)
        Inv2 = Cdet

c       compute second Piola Kirchhoff stress SecPK
        dK1(1:2,1:2) = ident(1:2,1:2)-(Cinv(1:2,1:2)*CdetI)
        dK2(1:2,1:2) = (ident(1:2,1:2)*CdetI) +
     $                 (Inv2-(Inv1*CdetI))*Cinv(1:2,1:2)

        K1 = Inv1+CdetI-3.0d0   
 
        do inode = 1, nnode
          bb(inode,1,1) = shap(1,inode)
          bb(inode,2,1) = 0.0d0
          bb(inode,3,1) = shap(2,inode)
          bb(inode,4,1) = 0.0d0
          bb(inode,1,2) = 0.0d0
          bb(inode,2,2) = shap(2,inode)
          bb(inode,3,2) = 0.0d0
          bb(inode,4,2) = shap(1,inode)
        end do

        do iset = 1,2! JFD waste of time to do the 2 indexes if we optimize only one parameter (mu or gamma)!
          do knode = 1,nnode
c           set variations in the shape function for the property JFD not clear, ask Sevan 
            var_shape_prop(:,:) = 0.0d0
            if (ielem_nod_grad(iset,knode).eq.1) then
              do idofn = 1,3! JFD idofn is not properly used ... idofn does not refer to a dof
                var_shape_prop(iset,idofn) = shap(idofn,knode)
              enddo
            endif

            if (iset.eq.1) then
              SecPK_grad(1:2,1:2) = prop_grad(2,3)*2.0d0*K1
     $              *var_shape_prop(iset,3)*exp(K1*prop_grad(iset,3))*
     $              dK1(1:2,1:2)
            else
              SecPK_grad(1:2,1:2) = var_shape_prop(iset,3)
     $              *(2.0d0*exp(K1*prop_grad(1,3))*dK1(1:2,1:2)-
     $              dK2(1:2,1:2))
            end if

            FS(:,:) = 0.0d0
            do j = 1, ndim
              do r = 1, ndim
                do i = 1, ndim
                  FS(i,j) = FS(i,j)+Fdef(i,r)*SecPK_grad(r,j)
                end do
              end do
            end do

            T(1) = FS(1,1) ! is FS a useful structure? use T from the start?
            T(2) = FS(2,2)
            T(3) = FS(1,2)
            T(4) = FS(2,1)

            ievab = 0
            eforc(:) = 0.0d0
            do inode = 1, nnode
              do i = 1,2! 2 = elemvec_ndofn(inode)
                ievab = ievab+1
                do j = 1,4
                  eforc(ievab)=eforc(ievab)+
     $                       bb(inode,j,i)*T(j)*wtjac
                end do
              end do
            end do! inode
            do ievab = 1,2*nnode
              egrad(iset,knode) = egrad(iset,knode) +
     $                            eforc(ievab)*temp_dual(ievab)
            end do

c           account for the regularization term (see the computation of regularization for more info)
c JFD            if (iset.eq.1) then
              if (ireg.eq.0) then
c               There is not supposed to be a regularization term: do nothing
              elseif (ireg.eq.1) then! H1 or Tikhonov reg.
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2))
              elseif ((ireg.eq.2).or.(ireg.eq.21)) then! TVD reg.
                deno = sqrt(beta(iset)*beta(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $              var_shape_prop(iset,2)*prop_grad(iset,2))/deno
              elseif (ireg.eq.3) then! power reg.
                deno = sqrt(beta(iset)*beta(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              ((1+deno-beta(iset))**(powe-1.0d0))*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2))/deno
              elseif (ireg.eq.31) then! power reg. (proposed by Paul Barbone)
                deno = beta(iset)*beta(iset)+prop_grad(iset,1)*
     $             prop_grad(iset,1)+prop_grad(iset,2)*prop_grad(iset,2)
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2))*
     $              (deno**(1-0.5d0*powe)-(deno-beta(iset)*beta(iset))*
     $            (1-0.5d0*powe)*(deno**(-0.5d0*powe)))/(deno**(2-powe))
              elseif (ireg.eq.4) then! logarithm reg.
                deno = sqrt(beta(iset)*beta(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2))/(deno*
     $              Treg(iset)*(1+(deno-beta(iset))/Treg(iset)))
              elseif (ireg.eq.41) then! logarithm reg. (second implementation)
                deno = Treg(iset)*Treg(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2)
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2))/deno
              elseif (ireg.eq.5) then! exponential reg. (contrast preserving)
                deno = sqrt(beta(iset)*beta(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2))*
     $              exp((beta(iset)-deno)/Treg(iset))/(Treg(iset)*deno)
              elseif (ireg.eq.6) then! fraction reg. (contrast preserving)
                deno = 1 + (prop_grad(iset,1)*prop_grad(iset,1)+
     $                      prop_grad(iset,2)*prop_grad(iset,2))/
     $                     (Treg(iset)*Treg(iset))
                deno = Treg(iset)*Treg(iset)*deno*deno
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              (var_shape_prop(iset,1)*prop_grad(iset,1)+
     $               var_shape_prop(iset,2)*prop_grad(iset,2))/deno
              elseif (ireg.eq.7) then! TVD reg. normalized by the field value with offset
                deno = sqrt(beta(iset)*beta(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              ((var_shape_prop(iset,1)*prop_grad(iset,1)+
     $                var_shape_prop(iset,2)*prop_grad(iset,2))/
     $                  (deno*prop_grad(iset,3))-
     $               sqrt(beta(iset)*beta(iset)+
     $                  prop_grad(iset,1)*prop_grad(iset,1)+
     $                  prop_grad(iset,2)*prop_grad(iset,2))
     $               *var_shape_prop(iset,3)/(prop_grad(iset,3)**2.0d0))
              elseif (ireg.eq.71) then! TVD reg. normalized by the field value without offset
                deno = sqrt(beta(iset)*beta(iset)+
     $              prop_grad(iset,1)*prop_grad(iset,1)+
     $              prop_grad(iset,2)*prop_grad(iset,2))
                egrad(iset,knode)=egrad(iset,knode) + wtjac*alpha(iset)*
     $              ((var_shape_prop(iset,1)*prop_grad(iset,1)+
     $                var_shape_prop(iset,2)*prop_grad(iset,2))/
     $                  (deno*prop_grad(iset,3))-
     $               (sqrt(beta(iset)*beta(iset)+
     $                  prop_grad(iset,1)*prop_grad(iset,1)+
     $                  prop_grad(iset,2)*prop_grad(iset,2))-beta(iset))
     $               *var_shape_prop(iset,3)/(prop_grad(iset,3)**2.0d0))
              else
                Print*,"elem606.f: ireg=",ireg,
     $                 " is not implemented: exiting"
                stop
              endif! ireg
c            elseif (iset.eq.2) then! JFD: the difference with iset.eq.1 is ...Mu instead of ...Gamma only!
c              if (ireg.eq.1) then
c                egrad(iset,knode) = egrad(iset,knode) + alpha(iset)*(
c     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
c     $              var_shape_prop(iset,2)*prop_grad(iset,2))*wtjac
c              elseif (ireg.eq.2) then
c                deno = sqrt(beta(iset)*beta(iset)+
c     $              prop_grad(iset,1)*prop_grad(iset,1)+
c     $              prop_grad(iset,2)*prop_grad(iset,2))
c                egrad(iset,knode) = egrad(iset,knode) + 
c     $              0.5d0*alpha(iset)*wtjac*(
c     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
c     $              var_shape_prop(iset,2)*prop_grad(iset,2))/deno
c              elseif (ireg.eq.3) then
c                deno = sqrt(beta(iset)*beta(iset)+
c     $              prop_grad(iset,1)*prop_grad(iset,1)+
c     $              prop_grad(iset,2)*prop_grad(iset,2))
c                egrad(2,inode) = egrad(2,inode) + 
c     $              0.5d0*alpha(iset)*Treg(iset)*wtjac/deno*(
c     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
c     $              var_shape_prop(iset,2)*prop_grad(iset,2))*
c     $              exp(-deno/Treg(iset))
c              elseif (ireg.eq.4) then
c                deno = sqrt(beta(iset)*beta(iset)+
c     $              prop_grad(iset,1)*prop_grad(iset,1)+
c     $              prop_grad(iset,2)*prop_grad(iset,2))
c                egrad(iset,knode) = egrad(iset,knode) + 
c     $              0.5d0*alpha(iset)*wtjac*((deno+1)**(powe-1.0d0))*(
c     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
c     $              var_shape_prop(iset,2)*prop_grad(iset,2))/deno
c              endif
c            endif! iset

            udiff(iset)=udiff(iset)+shap(3,knode)*uelem_diff(iset,knode)

          end do! knode
        end do! iset

        elemDataMatch(ielem) = elemDataMatch(ielem) + 0.5d0*wtjac*
     $        (udiff(1)*udiff(1)+udiff(2)*udiff(2))

        if (ireg.eq.1) then! H1 or Tikhonov reg.
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      (alpha(2)*(prop_grad(2,1)*prop_grad(2,1)+
     $                prop_grad(2,2)*prop_grad(2,2))+
     $      alpha(1)*(prop_grad(1,1)*prop_grad(1,1)+
     $                   prop_grad(1,2)*prop_grad(1,2)))
        elseif (ireg.eq.2) then! TVD reg. with offset
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)*sqrt(beta(2)*beta(2)+
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))+
     $       alpha(1)*sqrt(beta(1)*beta(1)+
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2)))
        elseif (ireg.eq.21) then! TVD reg. without offset
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)*(sqrt(beta(2)*beta(2)+
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))-beta(2))+
     $       alpha(1)*(sqrt(beta(1)*beta(1)+
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2))-beta(1)))
        elseif (ireg.eq.3) then! power reg.
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)/powe*((1+sqrt(beta(2)*beta(2)+
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))-beta(2))**(powe)-1.0d0)+
     $      alpha(1)/powe*((1+sqrt(beta(1)*beta(1)+
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2))-beta(1))**(powe)-1.0d0))
        elseif (ireg.eq.31) then! power regularization (proposed by Paul Barbone)
c JFD: for powe.eq.1, this is another implementation of TVD
          deno = prop_grad(2,1)*prop_grad(2,1)+prop_grad(2,2)*
     $           prop_grad(2,2)
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      (alpha(2)*deno/(beta(2)*beta(2)+deno)**(1-0.5d0*powe))
          deno = prop_grad(1,1)*prop_grad(1,1)+prop_grad(1,2)*
     $           prop_grad(1,2)
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      (alpha(1)*deno/(beta(1)*beta(1)+deno)**(1-0.5d0*powe))
c        elseif (ireg.eq.32) then! power with Tikhonov
c           1.0d0/powe * ((1.0d0+grad(u)^2)*(powe/2)-1.0d0) - not implemented 2011-05-17
        elseif (ireg.eq.4) then! logarithm reg.
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $       (alpha(2)*log(1+(sqrt(beta(2)*beta(2)+prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2))-
     $           beta(2))/Treg(2))+
     $       alpha(1)*log(1+(sqrt(beta(1)*beta(1)+prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2))-
     $           beta(1))/Treg(1)))
        elseif (ireg.eq.41) then! logarithm reg. (second implementation to avoid using a sqrt operator and the associated beta that has a significant influence on the convergence speed)
c [Herbert and Leahy, a generalized EM algorithm for 3D bayesian reconstructions from poisson data using Gibbs priors, IEEE trans Med imaging, MI-8, 1990, pp. 194-202]
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $       (alpha(2)*log(1+(prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))/(Treg(2)*Treg(2)))+
     $       alpha(1)*log(1+(prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2))/(Treg(1)*Treg(1))))
c        elseif (ireg.eq.42) then! logarithm with cosh [Green, Bayesian reconstructions from emission tomography data using a modified EM algorithm, IEEE Trans. Med. Imaging 9, 1990, pp. 84-93]
c           log(cosh(grad(u))) - not yet non-singular in 0 ... not implemented 2011-05-17 
        elseif (ireg.eq.5) then! exponential reg. (contrast preserving: no penalty for large jumps)
c JFD: proposed by Paul Barbone
c JFD: this regularization is very unstable
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)*(1.0d0 - exp((beta(2)-sqrt(beta(2)*beta(2)+
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2)))/Treg(2))) +
     $       alpha(1)*(1.0d0-exp((beta(1)-sqrt(beta(1)*beta(1)+
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2)))/Treg(1))))
        elseif (ireg.eq.6) then! fraction reg. (contrast preserving: no penalty for large jumps)
c JFD: this regularization is very unstable
          deno = (prop_grad(2,1)*prop_grad(2,1)+
     $            prop_grad(2,2)*prop_grad(2,2))/(Treg(2)*Treg(2))
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $                     (alpha(2)*deno/(1+deno))
          deno = (prop_grad(1,1)*prop_grad(1,1)+
     $            prop_grad(1,2)*prop_grad(1,2))/(Treg(1)*Treg(1))
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $                     (alpha(1)*deno/(1+deno))
        elseif (ireg.eq.7) then! TVD reg. normalized by the field with offset
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)*sqrt(beta(2)*beta(2)+
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))/prop_grad(2,3)+
     $       alpha(1)*sqrt(beta(1)*beta(1)+
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2))/prop_grad(1,3))
        elseif (ireg.eq.71) then! TVD reg. normalized by the field without offset
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      (alpha(2)*(sqrt(beta(2)*beta(2)+
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))-beta(2))/prop_grad(2,3)+
     $       alpha(1)*(sqrt(beta(1)*beta(1)+
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2))-beta(1))/prop_grad(1,3))
        endif
        l2grad1(ielem) = l2grad1(ielem) + sqrt(prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2))
        l2grad2(ielem) = l2grad2(ielem) + sqrt(prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2))
      end do! iinte
      dataMatch = dataMatch + elemDataMatch(ielem)
      regularization = regularization + elemRegul(ielem)
      l2grad1(ielem) = l2grad1(ielem)/dble(ninte)
      l2grad2(ielem) = l2grad2(ielem)/dble(ninte)
      suml2grad1 = suml2grad1 + l2grad1(ielem)*l2grad1(ielem)
      suml2grad2 = suml2grad2 + l2grad2(ielem)*l2grad2(ielem)

      return


 4    return
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

