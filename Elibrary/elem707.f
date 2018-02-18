c**********************************************************
      subroutine elem707 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     NOTES:
c     Material model by Sevan Goenezen, written by Sevan Goenezen
c     2D, plane strain, finite elasticity, incompressible
c     The strain energy density function has an exponential 
c     stress-strain with no deviatoric stress
c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem, s, ii
      parameter (mnode_elem = 9,mquad_elem=25)
      integer ielem, itask,ireg,iset
      integer nGauss, ninte, iinte, inode, jnode, knode
      integer ievab, jevab, jdofn, idofn, i, j, k, l, q, r
      double precision xjaco,wtjac,alpha,beta,h,h1,h2,h3,temp,temp_der
      double precision gamm, mu, Cdet, Inv1, Inv2, K1, K2(2,2), Fdet
      double precision shap(3,mnode_elem), pres(3),deno,tauMult
      double precision sg(mquad_elem),tg(mquad_elem),wg(mquad_elem)
      double precision DivQqGrad(2), DivQq(2), SecPK_grad(2,2)
      double precision Fdef(2,2), Ctens(2,2), dWdC(2,2), dJdC(2,2)
      double precision Cinv(2,2), SecPK(2,2), ident(2,2), Finv(2,2)
      double precision dK1(2,2),dK2(2,2),Ctang(2,2,2,2),d2JdC(2,2,2,2)
      double precision d2K1(2,2,2,2),d2K2(2,2,2,2),B(2,2,2,2)
      double precision Igeo(4,4),Lmat(2,2,2,2), Dgeomat(4,4)
      double precision T(4), Tp(4), FS(2,2), bb(mnode_elem,4,3)
      double precision Qq(mnode_elem,2,2), QqGrad(mnode_elem,2,2)
      double precision dCinvdC(2,2,2,2), dWdC_grad(2,2),udiff(2)
      double precision temp_dual(3*mnode_elem),prop_grad(2,3)
      double precision var_shape_prop(2,3)
      double precision Qlin(mnode_elem,2,2,2), DivQlin(2,2)
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

c---------------------------------------------------------------

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask  
      
c     --------------------------------------------------------
c     Initialize the elemental parameters
c     --------------------------------------------------------
 1    elemvec_ndofn(1:nnode) = 3
      buildSymmetricMatrices = .false.
      return
      
c     --------------------------------------------------------
c     Read and store the material properties
c     --------------------------------------------------------
c     nGauss     : no. of gauss pts/direction
c     ireg       : regularization type (0/1/2/3/4:none/H1/TVD/newTV/power)
c     alpha      : regularization weight for the non-linear parameter
c     beta       : extra parameter for TVD for the non-linear parameter - JFD ask Sevan for a better description
c     tauMult    : stabilization factor
c     s          : stabilization type (0/1:off/on) - JFD use stype instead of s?
 2    read(iread,*) nGauss,ireg,alpha,beta,tauMult,s
      write(iwrit,205) nGauss,ireg,alpha,beta,tauMult,s
 205  format(' finite elasticity (2D,2DOF) '/
     +       ' gauss pts/dir .........................',i12,/
     +       ' reg type(0/1/2:none/H1/TVD) ...........',i12,/
     +       ' regularization parameter ..............',1p,e16.4,/
     +       ' extra parameter (TVD) .................',1p,e16.4,/
     +       ' stabilization factor ..................',1p,e16.4,/
     +       ' stabilization terms 0/1 off/on ........',i12)
c      nprop = 6
      pelem(1) = dble(nGauss)
      pelem(2) = dble(ireg)
      pelem(3) = alpha
      pelem(4) = beta
      pelem(5) = tauMult
      pelem(6) = dble(s)

      return

c     --------------------------------------------------------
c     Build the elemental consistent tangent stiffness matrix (estif)
c       and the elemental RHS/residual (eforc)
c     --------------------------------------------------------
 3    nGauss = int(pelem(1))
      tauMult= pelem(5)
      s      = int(pelem(6))
      ndim = ndime


c     determine the characteristic element length h 
c        (for bilinear elements compute both diagonals and take the larger one)
      if (nnode.eq.4) then
        h1=0.0d0
        h=0.0d0
        do i=1, ndim
          h1= h1+ (xelem(i,3)-xelem(i,1))**2.0d0
          h = h + (xelem(i,4)-xelem(i,2))**2.0d0
        enddo
        if (h1 .gt. h) then
          h=sqrt(h1)
        else
          h=sqrt(h)
        endif
      elseif (nnode.eq.3) then
        h1=0.0d0
        h=0.0d0
        h3=0.0d0 
        do i=1, ndim
          h1= h1+ (xelem(i,1)-xelem(i,2))**2.0d0
          h = h + (xelem(i,1)-xelem(i,3))**2.0d0
          h3= h3+ (xelem(i,2)-xelem(i,3))**2.0d0
        enddo
        if (h1.gt.h) then
          h=h1    
        endif
        if (h3.gt.h) then
          h=h3
        endif
        h=sqrt(h)
      endif! if (nnode.eq.X)

c     set element stiffness to zero
      estif(:,:) = 0.0d0! elemental stiffness
      eforc(:) = 0.0d0! RHS
      Qq(:,:,:) = 0.0d0
      Qlin(:,:,:,:) = 0.0d0

c     identity matrix ident
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0  

      if (nnode.eq.3) then
         call gausstri(nGauss,ninte,sg,tg,wg)
      else
         call gauss2(nGauss,ninte,sg,tg,wg)
      endif
c
      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.false.,nnode,
     $      elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco
         
        mu    = 0.0d0
        gamm = 0.0d0
        do inode = 1,nnode
          mu    = mu   + shap(3,inode)*elemdata_nodal(2,inode)
          gamm  = gamm + shap(3,inode)*elemdata_nodal(1,inode)
        enddo

        if (s.eq.0) then
          temp  = 0.0d0! JFD better name than temp
        else if (s.eq.1) then
          temp = (0.5d0*tauMult*(h**2.0d0))/mu
        else
          write(iwrit,200) 
 200      format(4x,'Stabilization property must be 0 or 1')
          stop
        endif

c       Computing the deformation gradient at Gauss Point
        Fdef(1,1) = 1.0d0
        Fdef(2,2) = 1.0d0
        Fdef(1,2) = 0.0d0
        Fdef(2,1) = 0.0d0 
 
        do i = 1,ndim
          do j = 1,ndim
            do inode = 1,nnode
              Fdef(i,j)=Fdef(i,j)+uelem(i,inode)*shap(j,inode)
            end do
          end do
        end do
c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)

c       Finv is inverse of deformation gradient at Gauss Point
        Finv(1,1) = (1.0d0/Fdet)*Fdef(2,2)
        Finv(2,2) = (1.0d0/Fdet)*Fdef(1,1)
        Finv(1,2) = -(1.0d0/Fdet)*Fdef(1,2)
        Finv(2,1) = -(1.0d0/Fdet)*Fdef(2,1)

c       Computing pressure and spatial derivatives at Gauss Point
        pres(1:3) = 0.0d0

        do i = 1, 3
          do inode = 1, nnode
            pres(i) = pres(i) + uelem(3,inode)*shap(i,inode)
          end do
        end do

c       Computing the Cauchy tensor and its inverse
        Ctens(1:2,1:2) = 0.0d0
        do i = 1,ndim
          do j = 1,ndim
            do r = 1,ndim
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do

        Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)

        Cinv(1,1) = (1.0d0/Cdet)*Ctens(2,2)
        Cinv(2,2) = (1.0d0/Cdet)*Ctens(1,1)
        Cinv(1,2) = -(1.0d0/Cdet)*Ctens(1,2)
        Cinv(2,1) = -(1.0d0/Cdet)*Ctens(2,1)

c       Principal invariants Inv1 and Inv2
        Inv1 = Ctens(1,1)+Ctens(2,2)
        Inv2 = Cdet

c       dJdC is the derivative of the Jacobian with respect to 
c       the Cauchy Green tensor
        dJdC(1:2,1:2) = 0.5d0*Fdet*Cinv(1:2,1:2)

c       Computing second Piola Kirchhoff stress SecPK 
c       and the material tangent Ctang
        K1 = Inv1-2.0d0*log(Fdet)-2.0d0
        K2(1:2,1:2) = ident(1:2,1:2) - Cinv(1:2,1:2)   

        dWdC(1:2,1:2) = 0.5d0*mu*exp(K1*gamm)*K2(1:2,1:2)

        SecPK(1:2,1:2) = 2.0d0*(dWdC(1:2,1:2)-pres(3)*dJdC(1:2,1:2))

        do i = 1, 2
          do j = 1, 2
            do q = 1, 2
              do r = 1, 2
                dCinvdC(i,j,q,r) = -0.5d0*(Cinv(i,q)*Cinv(j,r)+
     $                                   Cinv(i,r)*Cinv(j,q)) 
              enddo
            enddo
          enddo
        enddo         

        do i = 1, 2
          do j = 1, 2
            do q = 1, 2
              do r = 1, 2
                Ctang(i,j,q,r) = 0.5d0*mu*exp(gamm*K1)*
     $                           (-dCinvdC(i,j,q,r)+gamm*
     $                           (ident(i,j)-Cinv(i,j))*
     $                           (ident(q,r)-Cinv(q,r)))
                d2JdC(i,j,q,r) = 0.25d0*Fdet*(Cinv(q,r)*Cinv(i,j)
     $                           +2.0d0*dCinvdC(i,j,q,r))
              enddo
            enddo
          enddo
        enddo         

c       construct the element stiffness matrix 
        Igeo(1,1) = SecPK(1,1)
        Igeo(1,2) = 0.0d0
        Igeo(1,3) = 0.5d0*SecPK(1,2)
        Igeo(1,4) = 0.5d0*SecPK(1,2)
        Igeo(2,2) = SecPK(2,2)
        Igeo(2,3) = 0.5d0*SecPK(2,1)
        Igeo(2,4) = -0.5d0*SecPK(2,1)
        Igeo(3,3) = 0.25d0*(SecPK(2,2)+SecPK(1,1))
        Igeo(3,4) = 0.25d0*(SecPK(2,2)-SecPK(1,1))
        Igeo(4,4) = 0.25d0*(SecPK(2,2)+SecPK(1,1))
        ! symmetric parts of Igeo:
        Igeo(2,1) = Igeo(1,2)
        Igeo(3,1) = Igeo(1,3)
        Igeo(4,1) = Igeo(1,4)
        Igeo(3,2) = Igeo(2,3)
        Igeo(4,2) = Igeo(2,4)
        Igeo(4,3) = Igeo(3,4)

        do i = 1, ndim
          do j = 1, ndim
            do q = 1, ndim
              do r = 1, ndim
                Lmat(i,j,q,r) = 4.0d0*(Fdef(i,1)*Fdef(q,1)*
     $                          (Ctang(j,1,r,1)-pres(3)*d2JdC(j,1,r,1))+
     $                           Fdef(i,1)*Fdef(q,2)*
     $                          (Ctang(j,1,r,2)-pres(3)*d2JdC(j,1,r,2))+
     $                           Fdef(i,2)*Fdef(q,1)*
     $                          (Ctang(j,2,r,1)-pres(3)*d2JdC(j,2,r,1))+
     $                           Fdef(i,2)*Fdef(q,2)*
     $                          (Ctang(j,2,r,2)-pres(3)*d2JdC(j,2,r,2)))
              end do
            end do
          end do
        end do

        Dgeomat(1,1) = Lmat(1,1,1,1)
        Dgeomat(1,2) = Lmat(1,1,2,2)
        Dgeomat(1,3) = 0.5d0*(Lmat(1,1,1,2)+Lmat(1,1,2,1))
        Dgeomat(1,4) = 0.5d0*(Lmat(1,1,1,2)-Lmat(1,1,2,1))
        Dgeomat(2,2) = Lmat(2,2,2,2)
        Dgeomat(2,3) = 0.5d0*(Lmat(2,2,1,2)+Lmat(2,2,2,1))
        Dgeomat(2,4) = 0.5d0*(Lmat(2,2,1,2)-Lmat(2,2,2,1))
        Dgeomat(3,3) = 0.25d0*Lmat(1,2,1,2)+0.5d0*Lmat(1,2,2,1)
     $                  +0.25d0*Lmat(2,1,2,1)
        Dgeomat(3,4) = 0.25d0*(Lmat(1,2,1,2)-Lmat(2,1,2,1))
        Dgeomat(4,4) = 0.25d0*Lmat(1,2,1,2)+0.25d0*Lmat(2,1,2,1)
     $                  -0.5d0*Lmat(1,2,2,1)
         ! symmetric part of the material tangent
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
          bb(inode,1,3) = shap(1,inode)
          bb(inode,2,3) = shap(2,inode)
          bb(inode,3,3) = 0.0d0
          bb(inode,4,3) = 0.0d0
        end do  

        ievab = 0
        do inode = 1, nnode
          do idofn = 1,3! 3 = elemvec_ndofn(inode)
            ievab = ievab+1
            jevab = 0
            do jnode = 1, nnode
              do jdofn = 1,3! 3 = elemvec_ndofn(jnode)
                jevab = jevab+1
                if (idofn.eq.1  .AND. jdofn.eq.3) then
                  estif(ievab,jevab)=estif(ievab,jevab)-Fdet*
     $                   (shap(1,inode)*Finv(1,1)+
     $                    shap(2,inode)*Finv(2,1))*shap(3,jnode)*wtjac
                endif

                if (idofn.eq.2  .AND. jdofn.eq.3) then
                  estif(ievab,jevab)=estif(ievab,jevab)-Fdet*
     $                   (shap(1,inode)*Finv(1,2)+
     $                    shap(2,inode)*Finv(2,2))*shap(3,jnode)*wtjac
                endif
                  
                if (idofn.eq.3 .AND. jdofn.eq.1) then
                  estif(ievab,jevab)=estif(ievab,jevab)+Fdet*
     $                    (shap(1,jnode)*Finv(1,1)+
     $                    shap(2,jnode)*Finv(2,1))*shap(3,inode)*wtjac
                  do j = 1,ndim
                   do i = 1,ndim
                    do k = 1,ndim
                     do l = 1,ndim
                       estif(ievab,jevab)=estif(ievab,jevab)+4.0d0*
     $                      temp *pres(i)*shap(j,inode)*d2JdC(i,j,k,l)
     $                      *Fdef(1,k)*shap(l,jnode)*wtjac
                     enddo
                    enddo
                   enddo
                  enddo
                endif
                  
                if (idofn.eq.3  .AND. jdofn.eq.2) then
                  estif(ievab,jevab)=estif(ievab,jevab)+Fdet*
     $                   (shap(1,jnode)*Finv(1,2)+
     $                    shap(2,jnode)*Finv(2,2))*shap(3,inode)*wtjac
                  do j = 1,ndim
                   do i = 1,ndim
                    do k = 1,ndim
                     do l = 1,ndim
                       estif(ievab,jevab)=estif(ievab,jevab)+4.0d0*
     $                      temp*pres(i)*shap(j,inode)*d2JdC(i,j,k,l)
     $                      *Fdef(2,k)*shap(l,jnode)*wtjac
                     enddo
                    enddo
                   enddo
                  enddo
                endif

                if (idofn.eq.3  .AND. jdofn.eq.3) then
                  do i = 1,ndim
                    do j = 1,ndim
                      estif(ievab,jevab)=estif(ievab,jevab)+2.0d0*
     $                   temp*dJdC(i,j)*shap(i,jnode)*shap(j,inode)
     $                   *wtjac
                    enddo
                  enddo
                endif

                if (idofn.lt.3  .AND. jdofn.lt.3) then
                  do i = 1, 4
                    do j = 1, 4
                      estif(ievab,jevab)=estif(ievab,jevab)+
     $                         bb(inode,i,idofn)*Dgeomat(i,j)*
     $                         bb(jnode,j,jdofn)*wtjac

                    enddo
                  enddo
                endif
              enddo! jdofn
            enddo! jnode
          enddo! idofn
        enddo! inode

c       Linearization of the stabilization, divergence term
        if (nnode.ne.3) then
          if (s.eq.1) then
            jevab = 1
            ievab = 0
            do jnode = 1, nnode
              ievab = 0
              call DivLin(ielem,jnode,Qlin,nnode,elnods,ndim,xelem,
     $                    uelem,elemdata_nodal)
              DivQlin(:,:) = 0.0d0
              do i = 1, ndim
                do j = 1, ndim
                  do inode = 1, nnode
          DivQlin(1,i) = DivQlin(1,i) + Qlin(inode,1,i,j)*shap(j,inode)
          DivQlin(2,i) = DivQlin(2,i) + Qlin(inode,2,i,j)*shap(j,inode)
                  enddo
                enddo
              enddo
     
              do knode = 1, nnode
                ievab = ievab + 3
                estif(ievab,jevab) = estif(ievab,jevab) -
     $                ((Finv(1,1)*DivQlin(1,1)+
     $                Finv(1,2)*DivQlin(1,2))*shap(1,knode)+
     $                (Finv(2,1)*DivQlin(1,1)+Finv(2,2)*DivQlin(1,2))*
     $                shap(2,knode))*wtjac*temp
     
                estif(ievab,jevab+1) = estif(ievab,jevab+1) -
     $                ((Finv(1,1)*DivQlin(2,1)+
     $                 Finv(1,2)*DivQlin(2,2))*shap(1,knode)+
     $                (Finv(2,1)*DivQlin(2,1)+Finv(2,2)*DivQlin(2,2))*
     $                 shap(2,knode))*wtjac*temp 
              
              end do
              jevab = jevab + 3
            end do! jnode
    
            call DivInterp(ielem,Qq,nnode,elnods,ndim,xelem,uelem,
     $                     elemdata_nodal)
    
            DivQq(1:2) = 0.0d0
            do i = 1, ndim
              do j = 1, ndim
                do inode = 1, nnode
                  DivQq(i) = DivQq(i) + Qq(inode,i,j)*shap(j,inode)
                enddo
              enddo
            enddo

            ievab = 0
            do jnode = 1, nnode
              jevab = 1
              ievab = ievab + 3
              do knode = 1, nnode
                do i = 1, ndim
                  do j = 1, ndim
                    do l = 1, ndim
                      estif(ievab,jevab) = estif(ievab,jevab) + temp*
     $                              Finv(i,1)*Finv(l,j)*shap(l,knode)*
     $                              DivQq(j)*shap(i,jnode)* wtjac
     
                      estif(ievab,jevab+1) = estif(ievab,jevab+1) +temp*
     $                              Finv(i,2)*Finv(l,j)*shap(l,knode)*
     $                              DivQq(j)*shap(i,jnode)* wtjac
                    enddo! l=
                  enddo! j=
                enddo! i=
                jevab = jevab + 3
              enddo! knode = 
            enddo! jnode =
          endif! if (s.eq.1)
        endif! if (nnode.ne.3)

c       create element residual for right hand side
        do inode = 1, nnode
          bb(inode,1,1) = shap(1,inode)
          bb(inode,2,1) = 0.0d0
          bb(inode,3,1) = shap(2,inode) 
          bb(inode,4,1) = 0.0d0
          bb(inode,1,2) = 0.0d0
          bb(inode,2,2) = shap(2,inode)
          bb(inode,3,2) = 0.0d0
          bb(inode,4,2) = shap(1,inode)
          bb(inode,1,3) = shap(1,inode)
          bb(inode,2,3) = shap(2,inode)
          bb(inode,3,3) = shap(3,inode)
          bb(inode,4,3) = 0.0d0
        end do

        Tp(1:4) = 0.0d0
        do i = 1, ndim
          do j = 1, ndim
            Tp(i) = Tp(i) + 2.0d0*temp*dJdC(j,i)*pres(j)
          enddo
        enddo

        Tp(3) = Fdet-1.0d0 

        if (nnode.ne.3) then
          if (iinte.eq.1 .AND. s.eq.1) then
            call DivInterp(ielem,Qq,nnode,elnods,ndim,xelem,uelem,
     $                     elemdata_nodal)
          end if

          if (s.eq.1) then
c           Compute the Divergence of Qq after interpolation at Gauss Points
            DivQq(1:2) = 0.0d0
            do i = 1, ndim
              do j = 1, ndim
                do inode = 1, nnode
                  DivQq(i) = DivQq(i) + Qq(inode,i,j)*shap(j,inode)
                enddo
              enddo
            enddo

            do i = 1, ndim
              do j = 1, ndim
                Tp(i) = Tp(i) - temp*Finv(i,j)*DivQq(j)
              enddo
            enddo
          endif! if()
        endif! if (nnode.ne.3)	 

        FS(:,:) = 0.0d0
        do i = 1, ndim
          do j = 1, ndim
            do r = 1, ndim
              FS(i,j) = FS(i,j)+Fdef(i,r)*SecPK(r,j)
            enddo
          enddo
        enddo

        T(1) = FS(1,1)
        T(2) = FS(2,2)
        T(3) = FS(1,2)
        T(4) = FS(2,1)
         
        ievab = 0
        do inode = 1, nnode
          do i = 1,3! 3 = elemvec_ndofn(inode)
            ievab = ievab+1
c            if (i.eq.elemvec_ndofn(inode)) then
            if (i.eq.3) then
              do j = 1, 4
                eforc(ievab)=eforc(ievab)+
     $                          bb(inode,j,i)*Tp(j)*wtjac
              enddo
            else
              do j = 1, 4
                eforc(ievab)=eforc(ievab)+
     $                          bb(inode,j,i)*T(j)*wtjac
               enddo
             endif
           enddo
         enddo! inode
c	 Following is for output purpose
c           do j=1, 12
c             Print*,"element force at",j, "  is", eforc(j)
c           end do
c           stop
      enddo! iinte
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
 6    nGauss  = int(pelem(1))
      eforc(:) = 0.0d0! elemental RHS/residual
      ndim = ndime! local copy

      if (nnode.eq.3) then
        call gausstri(nGauss,ninte,sg,tg,wg)
      else
        call gauss2(nGauss,ninte,sg,tg,wg)
      endif
c
      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.false.,nnode,
     $      elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco
c
        ievab = 0
        do inode = 1, nnode
          do idofn = 1,3! 3 = elemvec_ndofn(inode)
            ievab = ievab+1
            do jnode = 1, nnode
              do jdofn = 1,3! 3 = elemvec_ndofn(jnode)
                if (idofn.eq.jdofn) then
                  eforc(ievab) = eforc(ievab)-
     $                       shap(3,inode)*shap(3,jnode)
     $                       *uelem_diff(jdofn,jnode)*wtjac
                endif
              enddo! jdofn
            enddo! jnode
          enddo! idofn
        enddo! inode
      enddo! iinite

      return


c     --------------------------------------------------------
c     Compute the objective function (dataMatch+regularization)
c      and its gradient (egrad) on the element
c     --------------------------------------------------------
 7    nGauss  = int(pelem(1))
      ireg    = int(pelem(2))
      alpha   = pelem(3)
      beta    = pelem(4)
      tauMult = pelem(5)
      s       = int(pelem(6))
      ndim = ndime

c     identity matrix ident
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0 

      elemDataMatch(ielem) = 0.0d0
      elemRegul(ielem) = 0.0d0
      l2grad1(ielem) = 0.0d0
      l2grad2(ielem) = 0.0d0
      
      egrad(:,:) = 0.0d0! initialize egrad

      temp_dual(:) = 0.0d0
      ii=0

      do inode = 1,nnode
        do idofn = 1,3! 3 = elemvec_ndofn(inode)
          ii = ii+1
          temp_dual(ii)   = uelem_dual(idofn,inode)
        enddo
      enddo

c     determine the characteristic element legth h for bilinear elements 
c       compute both diagonals and take the bigger one

      if (nnode.eq.4) then
        h1=0.0d0
        h2=0.0d0
        do i=1, ndim
          h1 = h1 + (xelem(i,3)-xelem(i,1))**2.0d0
          h2 = h2 + (xelem(i,4)-xelem(i,2))**2.0d0
        enddo
        if (h1 .gt. h2) then
          h=sqrt(h1)
        else
          h=sqrt(h2)
        endif
      elseif (nnode.eq.3) then
        h1=0.0d0
        h2=0.0d0
        h3=0.0d0 
        do i=1, ndim
          h1 = h1 + (xelem(i,1)-xelem(i,2))**2.0d0
          h2 = h2 + (xelem(i,1)-xelem(i,3))**2.0d0
          h3 = h3 + (xelem(i,2)-xelem(i,3))**2.0d0
        enddo
        if (h1.gt.h2) then
          h=h1
        else
          h=h2
        endif
        if (h3.gt.h) then
          h=h3
        endif
        h=sqrt(h)
      endif!if (nnode.eq.X)
      
      if (nnode.eq.3) then
        call gausstri(nGauss,ninte,sg,tg,wg)
      else
        call gauss2(nGauss,ninte,sg,tg,wg)
      endif
c
      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.false.,nnode,
     $      elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco

c       create prop_grad (contains prop gradient and value) 
        prop_grad(:,:)  = 0.0d0
        do iset = 1,2
          do idofn = 1,3! two gradents+the value
            do inode = 1,nnode
               prop_grad(iset,idofn) = prop_grad(iset,idofn)+
     $                            shap(idofn,inode)*
     $                      elemdata_nodal(iset,inode)
            enddo
          enddo
        enddo

c       compute the gradient and the objective function
        udiff(:)      = 0.0d0
        QqGrad(:,:,:) = 0.0d0

c       Compute the deformation gradient at Gauss Point
        Fdef(1,1) = 1.0d0
        Fdef(2,2) = 1.0d0
        Fdef(1,2) = 0.0d0
        Fdef(2,1) = 0.0d0 

        do i = 1,ndim
          do j = 1,ndim
            do inode = 1,nnode
              Fdef(i,j)=Fdef(i,j)+uelem(i,inode)*shap(j,inode)
            enddo
          enddo
        enddo
c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)

c       Finv is inverse of deformation gradient at Gauss Point
        Finv(1,1) = (1.0d0/Fdet)*Fdef(2,2)
        Finv(2,2) = (1.0d0/Fdet)*Fdef(1,1)
        Finv(1,2) = -(1.0d0/Fdet)*Fdef(1,2)
        Finv(2,1) = -(1.0d0/Fdet)*Fdef(2,1)

c       Compute pressure and spatial derivatives at Gauss Point
        pres(1:3) = 0.0d0
        do i = 1, 3
          do inode = 1, nnode
            pres(i) = pres(i) + uelem(3,inode)*shap(i,inode)
          enddo
        enddo

c       Compute the Cauchy tensor and its inverse
        Ctens(1:2,1:2) = 0.0d0
        do i = 1,ndim
          do j = 1,ndim
            do r = 1,ndim
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            enddo
          enddo
        enddo

        Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)

        Inv2 = Cdet
        Inv1 = Ctens(1,1)+Ctens(2,2)

        Cinv(1,1) = (1.0d0/Cdet)*Ctens(2,2)
        Cinv(2,2) = (1.0d0/Cdet)*Ctens(1,1)
        Cinv(1,2) = -(1.0d0/Cdet)*Ctens(1,2)
        Cinv(2,1) = -(1.0d0/Cdet)*Ctens(2,1)

c       dJdC is the derivative of the Jacobian with respect to 
c         the Cauchy Green tensor

        dJdC(1:2,1:2) = 0.5d0*Fdet*Cinv(1:2,1:2)

        do inode = 1, nnode
          bb(inode,1,1) = shap(1,inode)
          bb(inode,2,1) = 0.0d0
          bb(inode,3,1) = shap(2,inode) 
          bb(inode,4,1) = 0.0d0
          bb(inode,1,2) = 0.0d0
          bb(inode,2,2) = shap(2,inode)
          bb(inode,3,2) = 0.0d0
          bb(inode,4,2) = shap(1,inode)
          bb(inode,1,3) = shap(1,inode)
          bb(inode,2,3) = shap(2,inode)
          bb(inode,3,3) = shap(3,inode)
          bb(inode,4,3) = 0.0d0
        end do

        do iset = 1,2
          do knode = 1,nnode
c           set variations in the shape function for the property 
            var_shape_prop(:,:) = 0.0d0
            if (ielem_nod_grad(iset,knode).eq.1) then
              do idofn = 1,3
                var_shape_prop(iset,idofn) = shap(idofn,knode)
              enddo
            endif

            if (s.eq.0) then
              temp  = 0.0d0
              temp_der  = 0.0d0
            endif 
    
            if (s.eq.1) then
              temp = 0.5d0*(tauMult*(h**2.0d0))/prop_grad(2,3)
              temp_der = 0.0d0
            endif
    
            if (s.eq.1 .AND. iset.eq.2) then  
              temp_der = -0.5d0*var_shape_prop(2,3)*
     $                (tauMult*(h**2.0d0))/(prop_grad(2,3)**2.0d0)
            endif

            Tp(1:4) = 0.0d0
            do i = 1, ndim
              do j = 1, ndim
                Tp(i) = Tp(i) + 2.0d0*temp_der*dJdC(j,i)*pres(j)
              enddo
            enddo

c            Tp(3) = 0.0d0 is zero here

            if (nnode.ne.3) then
              if (s.eq.1) then
                call DivInterpGrad(ielem,QqGrad,iset,nnode,elnods,
     $                             ndim,xelem,uelem,elemdata_nodal)
                call DivInterp(ielem,Qq,nnode,elnods,ndim,xelem,uelem,
     $                         elemdata_nodal)
c               Compute the Divergence of Qq after interpolation at Gauss Points
                DivQqGrad(1:2) = 0.0d0
                do i = 1, ndim
                  do j = 1, ndim
                    DivQqGrad(i) = DivQqGrad(i) + QqGrad(knode,i,j)
     $                        *shap(j,knode)
                  enddo
                enddo

                DivQq(1:2) = 0.0d0
                do i = 1, ndim
                  do j = 1, ndim
                    do inode = 1, nnode
                      DivQq(i) = DivQq(i)+Qq(inode,i,j)*shap(j,inode)
                    enddo
                  enddo
                enddo

                do i = 1, ndim
                  do j = 1, ndim
                    Tp(i) = Tp(i) - temp_der*Finv(i,j)*DivQq(j)
     $                     - temp*Finv(i,j)*DivQqGrad(j)
                  enddo
                enddo
              endif! (s.eq.1)
            endif! nnode.ne.3
 
c           Compute derivative of second Piola Kirchhoff stress SecPK 
c             with respect to the material parameters
         
            K1 = Inv1-2.0d0*log(Fdet)-2.0d0
            K2(1:2,1:2) = ident(1:2,1:2) - Cinv(1:2,1:2)

            if (iset.eq.1) then
              dWdC_grad(1:2,1:2) = 0.5d0*prop_grad(2,3)
     $                  *K1*var_shape_prop(iset,3)
     $                  *exp(K1*prop_grad(iset,3))*K2(1:2,1:2)
            elseif (iset.eq.2) then
              dWdC_grad(1:2,1:2) = 0.5d0*var_shape_prop(iset,3)
     $                    *exp(K1*prop_grad(1,3))*K2(1:2,1:2)
            endif

            SecPK_grad(1:2,1:2) = 2.0d0*dWdC_grad(1:2,1:2)

            FS(:,:) = 0.0d0
            do i = 1, ndim
              do j = 1, ndim
                do r = 1, ndim
                  FS(i,j) = FS(i,j)+Fdef(i,r)*SecPK_grad(r,j)
                enddo
              enddo
            enddo

            T(1) = FS(1,1)
            T(2) = FS(2,2)
            T(3) = FS(1,2)
            T(4) = FS(2,1)
c	     Print*,"Fdef",Fdef
c	     Print*,"dWdC",dWdC
c	     stop
         
            ievab = 0
            eforc(:) = 0.0d0
            do inode = 1, nnode
              do i = 1,3! 3 = elemvec_ndofn(inode)
                ievab = ievab+1
                if (i.eq.elemve c_ndofn(inode)) then
                  do j = 1, 4
                    eforc(ievab)=eforc(ievab)+
     $                          bb(inode,j,i)*Tp(j)*wtjac
                  enddo
                else
                  do j = 1, 4
                    eforc(ievab)=eforc(ievab)+
     $                          bb(inode,j,i)*T(j)*wtjac
                  enddo
                endif
              enddo
            enddo! inode
 
            do ii = 1, 3*nnode!for triangular element, change 12 to 3*nnode
              egrad(iset,knode) = egrad(iset,knode) +
     $                            eforc(ii)*temp_dual(ii)
            enddo

c     account for the regularization term
            if (iset.eq.2) then
              if (ireg.eq.1) then
                egrad(iset,knode) = egrad(iset,knode) + alpha*(
     $              var_shape_prop(iset,2)*prop_grad(iset,2))*
     $              wtjac
              elseif (ireg.eq.2) then
                deno = sqrt(
     $              beta*beta+
     $              prop_grad(2,1)*prop_grad(2,1)+
     $              prop_grad(2,2)*prop_grad(2,2))
                egrad(iset,knode) = egrad(iset,knode) + 
     $              0.5d0*alpha*(
     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
     $              var_shape_prop(iset,2)*prop_grad(iset,2))*
     $              wtjac/deno
              endif
            endif! iset.eq.2
           udiff(iset)=udiff(iset)+shap(3,knode)*uelem_diff(iset,knode)
          enddo! knode
        enddo! iset
      
        elemDataMatch(ielem) = elemDataMatch(ielem) + 
     $        0.5d0*(
     $        udiff(1)*udiff(1)+
     $        udiff(2)*udiff(2))*
     $        wtjac

        if (ireg.eq.1) then
          elemRegul(ielem) = elemRegul(ielem)+ 
     $           0.5d0*alpha*(
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))*
     $           wtjac
        elseif (ireg.eq.2) then
          elemRegul(ielem) = elemRegul(ielem)+ 
     $           0.5d0*alpha*sqrt(beta*beta+
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))*
     $           wtjac
        endif
        l2grad1(ielem) = l2grad1(ielem) + sqrt(prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2))
        l2grad2(ielem) = l2grad2(ielem) + sqrt(prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2))
      enddo! iinte

      dataMatch=dataMatch+elemDataMatch(ielem)
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


c*******************************************************************************
      subroutine DivInterp(ielem,Qq,nnode,elnods,ndim,xelem,uelem,
     $   elemdata_nodal)
c*******************************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none

      integer mnode_elem
      parameter (mnode_elem = 9)
      integer ielem, ndim, nnode
      double precision Qq(mnode_elem,2,2)
      integer elnods(mnode)
      double precision xelem(ndim,*)
      double precision uelem(mdofn,*)
      double precision elemdata_nodal(nset_nodal,*)

      integer ii,i,j,r,inode
      double precision K1,Inv1,Inv2,Cdet,xjaco,Fdet
      double precision Fdef(2,2),K2(2,2),Cinv(2,2),Ctens(2,2),ident(2,2)
      double precision ss(4), tt(4),dWdC(2,2)
      double precision shap(3,mnode_elem)
      
c     Last term in stabilization contains Divergence of a tensor which is spatially 
c       dependent. This will be linearly approximated to enable taking the Divergence

      data ss/-1.0d0,1.0d0,1.0d0,-1.0d0/,
     $     tt/-1.0d0,-1.0d0,1.0d0,1.0d0/

c     identity matrix ident

      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0 

      Qq(:,:,:) = 0.0d0
      do inode = 1, nnode
        call shape2(ielem,ss(inode),tt(inode),shap,xjaco,.false.,
     $     nnode,elnods,ndim,xelem)

c       Compute the deformation gradient at Gauss Point
         
         Fdef(1,1) = 1.0d0
         Fdef(2,2) = 1.0d0
         Fdef(1,2) = 0.0d0
         Fdef(2,1) = 0.0d0 
 
         do i = 1,ndim
           do j = 1,ndim
             do ii = 1,nnode
               Fdef(i,j)=Fdef(i,j)+uelem(i,ii)*shap(j,ii)
             end do
           end do
         end do
c        Fdet is the Jacobian, determinant of the deformation gradient
         Fdet = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)

c        Computing the Cauchy tensor and its inverse

         Ctens(1:2,1:2) = 0.0d0

         do i = 1,ndim
           do j = 1,ndim
             do r = 1,ndim
               Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
             end do
           end do
         end do

         Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)

         Cinv(1,1) = (1.0d0/Cdet)*Ctens(2,2)
         Cinv(2,2) = (1.0d0/Cdet)*Ctens(1,1)
         Cinv(1,2) = -(1.0d0/Cdet)*Ctens(1,2)
         Cinv(2,1) = -(1.0d0/Cdet)*Ctens(2,1)

c        Principal invariants Inv1 and Inv2

         Inv1 = Ctens(1,1)+Ctens(2,2)
         Inv2 = Cdet

     
         K1 = Inv1-2.0d0*log(Fdet)-2.0d0
         K2(1:2,1:2) = ident(1:2,1:2) - Cinv(1:2,1:2)   


         dWdC(1:2,1:2) = 0.5d0*elemdata_nodal(2,inode)
     $               *exp(K1*elemdata_nodal(1,inode))*K2(1:2,1:2)     

        do i = 1, ndim
          do j =  1, ndim
            do r = 1, ndim
              Qq(inode,i,j) = Qq(inode,i,j) + 2.0d0*Fdef(i,r)*dWdC(r,j)
            enddo
          enddo
        enddo
      enddo! inode

      end subroutine DivInterp



c*******************************************************************************
      subroutine DivLin(ielem,jnode,Qlin,nnode,elnods,ndim,xelem,
     $    uelem,elemdata_nodal)

c*******************************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      
      integer mnode_elem
      parameter (mnode_elem = 9)
      integer ielem, jnode, ndim, nnode
      double precision Qlin(mnode_elem,2,2,2)
      integer elnods(mnode)
      double precision xelem(ndim,*)
      double precision uelem(mdofn,*)
      double precision elemdata_nodal(nset_nodal,*)

      integer ii,i,j,r,q,l,m,inode
      double precision K1,Inv1,Inv2,Cdet,xjaco
      double precision Fdef(2,2),K2(2,2),Cinv(2,2),Ctens(2,2),ident(2,2)
      double precision ss(4),tt(4),dWdC(2,2),Fdet,Ctang(2,2,2,2)
c JFD      double precision elem_nodal_gamma(mnode_elem)! unused
      double precision shap(3,mnode_elem)
      double precision elem_nodal_mu(mnode_elem)
      double precision dCinvdC(2,2,2,2)
      

c     identity matrix
                                                        
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0

      data ss/-1.0d0,1.0d0,1.0d0,-1.0d0/,
     $        tt/-1.0d0,-1.0d0,1.0d0,1.0d0/


      Qlin(:,:,:,:) = 0.0d0
      do inode = 1, nnode
        call shape2(ielem,ss(inode),tt(inode),shap,xjaco,.false.,
     $      nnode,elnods,ndim,xelem)

c       Compute the deformation gradient at Gauss Point
        Fdef(1,1) = 1.0d0
        Fdef(2,2) = 1.0d0
        Fdef(1,2) = 0.0d0
        Fdef(2,1) = 0.0d0

        do i = 1,ndim
          do j = 1,ndim
            do ii = 1,nnode
              Fdef(i,j)=Fdef(i,j)+uelem(i,ii)*shap(j,ii)
            enddo
          enddo
        enddo
c        Fdet is the Jacobian, determinant of the deformation gradient
         Fdet = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)

c        Computing the Cauchy tensor and its inverse

         Ctens(1:2,1:2) = 0.0d0

         do i = 1,ndim
           do j = 1,ndim
             do r = 1,ndim
               Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
             end do
           end do
         end do

         Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)

         Cinv(1,1) = (1.0d0/Cdet)*Ctens(2,2)
         Cinv(2,2) = (1.0d0/Cdet)*Ctens(1,1)
         Cinv(1,2) = -(1.0d0/Cdet)*Ctens(1,2)
         Cinv(2,1) = -(1.0d0/Cdet)*Ctens(2,1)

c        Principal invariants Inv1 and Inv2

         Inv1 = Ctens(1,1)+Ctens(2,2)
         Inv2 = Cdet

         K1 = Inv1-2.0d0*log(Fdet)-2.0d0
         K2(1:2,1:2) = ident(1:2,1:2) - Cinv(1:2,1:2)   


         dWdC(1:2,1:2) = 0.5d0*elemdata_nodal(2,inode)*
     $          exp(K1*elemdata_nodal(1,inode))*K2(1:2,1:2)

     
         do i = 1, 2
           do j = 1, 2
             do q = 1, 2
               do r = 1, 2
                 dCinvdC(i,j,q,r) = -0.5*(Cinv(i,q)*Cinv(j,r)+
     $                                   Cinv(i,r)*Cinv(j,q))
               end do
             end do
           end do
         end do

         do i = 1, 2
           do j = 1, 2
             do q = 1, 2
               do r = 1, 2
                 Ctang(i,j,q,r) = 0.5*elemdata_nodal(2,inode)*
     $                            exp(elemdata_nodal(1,inode)*K1)*
     $                            (-dCinvdC(i,j,q,r)+
     $                            elemdata_nodal(1,inode)*
     $                            (ident(i,j)-Cinv(i,j))*
     $                            (ident(q,r)-Cinv(q,r)))
               end do
             end do
           end do
         end do

         Qlin(inode,1,1,1) = Qlin(inode,1,1,1) + 2.0d0*(shap(1,jnode)*
     $                              dWdC(1,1)+shap(2,jnode)*dWdC(2,1))
     
         Qlin(inode,1,1,2) = Qlin(inode,1,1,2) + 2.0d0*(shap(1,jnode)*
     $                              dWdC(1,2)+shap(2,jnode)*dWdC(2,2))
     
     
         Qlin(inode,2,2,1) = Qlin(inode,1,1,1) 
         Qlin(inode,2,2,2) = Qlin(inode,1,1,2)
                                         
  
         do i = 1, ndim
           do j =  1, ndim
             do m = 1, ndim
               do r = 1, ndim
                 do q = 1, ndim
                   do l = 1, ndim
                     Qlin(inode,m,i,j) = Qlin(inode,m,i,j) + 
     $                                  4.0d0*Fdef(i,r)*Ctang(r,j,q,l)*
     $                                  Fdef(m,l)*shap(q,jnode)
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo! inode
      end subroutine DivLin


c*******************************************************************************
      subroutine DivInterpGrad(ielem,QqGrad,iset,nnode,elnods,
     $  ndim,xelem,uelem,elemdata_nodal)
c*******************************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none

      integer mnode_elem
      parameter (mnode_elem = 9)
      integer ielem, iset, nnode, ndim
      double precision QqGrad(mnode_elem,2,2)
      integer elnods(mnode)
      double precision xelem(ndim,*)
      double precision uelem(mdofn,*)
      double precision elemdata_nodal(nset_nodal,*)

      integer ii,i,j,r,inode
      double precision K1,Inv1,Inv2,Cdet,xjaco,Fdet
      double precision Fdef(2,2),K2(2,2),Cinv(2,2),Ctens(2,2),ident(2,2)
      double precision ss(4), tt(4),dWdC_grad(2,2)
      double precision shap(3,mnode_elem)
      
c        Last term in stabilization contains Divergence of a tensor which is spatially 
c        dependent. This will be linearly approximated to enable taking the Divergence

         data ss/-1.0d0,1.0d0,1.0d0,-1.0d0/,
     $        tt/-1.0d0,-1.0d0,1.0d0,1.0d0/

c     identity matrix ident

      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0 

      QqGrad(:,:,:) = 0.0d0

      do inode = 1, nnode
        call shape2(ielem,ss(inode),tt(inode),shap,xjaco,.false.,
     $      nnode,elnods,ndim,xelem)

c       Compute the deformation gradient at Gauss Point
        Fdef(1,1) = 1.0d0
        Fdef(2,2) = 1.0d0
        Fdef(1,2) = 0.0d0
        Fdef(2,1) = 0.0d0 
        do i = 1,ndim
          do j = 1,ndim
            do ii = 1,nnode
              Fdef(i,j)=Fdef(i,j)+uelem(i,ii)*shap(j,ii)
            enddo
          enddo
        enddo
c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)

c       Compute the Cauchy tensor and its inverse
        Ctens(1:2,1:2) = 0.0d0
        do i = 1,ndim
          do j = 1,ndim
            do r = 1,ndim
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            enddo
          enddo
        enddo

        Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)

        Cinv(1,1) = (1.0d0/Cdet)*Ctens(2,2)
        Cinv(2,2) = (1.0d0/Cdet)*Ctens(1,1)
        Cinv(1,2) = -(1.0d0/Cdet)*Ctens(1,2)
        Cinv(2,1) = -(1.0d0/Cdet)*Ctens(2,1)

c       Principal invariants Inv1 and Inv2

        Inv1 = Ctens(1,1)+Ctens(2,2)
        Inv2 = Cdet
         
        K1 = Inv1-2.0d0*log(Fdet)-2.0d0
        K2(1:2,1:2) = ident(1:2,1:2) - Cinv(1:2,1:2)

        if (iset.eq.2) then
          dWdC_grad(1:2,1:2) = 0.5d0*exp(K1*elemdata_nodal(1,inode))
     $                          *K2(1:2,1:2)
        elseif (iset.eq.1) then
          dWdC_grad(1:2,1:2) = 0.5d0*elemdata_nodal(2,inode)
     $                   *K1*exp(K1*elemdata_nodal(1,inode))*K2(1:2,1:2)
        else
          Print*,"iset is wrong"
          stop
        end if

        do i = 1, ndim
          do j =  1, ndim
            do r = 1, ndim
              QqGrad(inode,i,j) = QqGrad(inode,i,j)+
     $                            2.0d0*Fdef(i,r)*dWdC_grad(r,j)
            end do
          end do
        end do

      end do! inode

      end subroutine DivInterpGrad
