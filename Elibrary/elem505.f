c**********************************************************
      subroutine elem505 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c      Material model by Sevan Goenezen, written by Sevan Goenezen   
c     2D (plane strain), triangles/quadrilaterals, finite elasticity, incompressible
c     The strain energy density (stress-strain) function is a modified
c      Veronda-Westman relation; it is used in the paper by Sevan Goenezen
c      and Assad Oberai (2010?).
c      W=\mu/(2.\gamma).(exp(\gamma.(J^(-2/3).I_1-3))-1)
c      S=-p.J.C^-1+\mu.J^(-2/3).(I_{dentity}-(I_1.C^-1)/3).exp(\gamma.J^(-2/3).I_1-3)
c      with this choice: trace(sigma) = -3 p
c     cite Maniatty's paper on stabilized FEM
c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem, s, ii
      parameter (mnode_elem = 9,mquad_elem=25)
      integer ielem, itask,k,ireg,iset
      integer l, ninte, iinte, inode, jnode,knode
      integer ievab, jevab, jdofn, idofn, i, j, q, r
      double precision xjaco,wtjac,temp,temp_der
      double precision h, h1, h2, tauMult! stabilization variables
      double precision powe, deno, alpha(2), beta(2), Treg(2)! regularization variables
      double precision gamm, mu, Cdet, Inv1, Inv2, K1, K2(2,2), Fdet
      double precision shap(3,mnode_elem), pres(3)
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
c     l          : no. of gauss pts/direction
c     ireg       : regularization type (0/1/2/3/4:none/H1/TVD/newTV/power)
c     alpha(1)   : regularization weight for the non-linear parameter
c     alpha(2)   : regularization weight for the shear modulus
c     beta(1)    : extra parameter for TVD for the non-linear parameter - JFD ask Sevan for a better description
c     beta(2)    : extra parameter for TVD for the shear modulus- JFD idem
c     Treg(1)    : second extra parameter for new TVD for the non-linear parameter
c     Treg(2)    : second extra parameter for new TVD for the shear modulus
c     tauMult    : stabilization factor
c     s          : stabilization type (0/1:off/on) - JFD use stype instead of s?
 2    read(iread,*) l,ireg,alpha(1),alpha(2),beta(1),beta(2),
     +     Treg(1),Treg(2),tauMult,s
      write(iwrit,205) l,ireg,alpha(1),alpha(2),beta(1),beta(2),
     +     Treg(1),Treg(2),tauMult,s
 205  format(/' finite elasticity (2D,plane strain approx.,3DOFs) '/
     +       ' gauss pts/dir .........................',i12,/
     +       ' reg type(0/1/2/4:none/H1/TVD/log, etc).',i12,/
     +       ' regularization parameter 1 ............',1p,e16.4,/
     +       ' regularization parameter 2 ............',1p,e16.4,/
     +       ' extra parameter (TVD) .................',1p,e16.4,/! modify description here JFD
     +       ' extra parameter (TVD) .................',1p,e16.4,/! idem JFD
     +       ' second extra parameter (TVD) ..........',1p,e16.4,/! modify description here JFD
     +       ' second extra parameter (TVD) ..........',1p,e16.4,/! idem JFD
     +       ' stabilization factor ..................',1p,e16.4,/
     +       ' stabilization type (0/1:off/on) .......',i12)
c      nprop = 10
c     check for error in input file
      if (alpha(1).lt.0.0d0) then
        print*,'elem505.f: alpha(1) is less than 0.0d0: exiting'
        stop
      elseif (alpha(2).lt.0.0d0) then
        print*,'elem505.f: alpha(2) is less than 0.0d0: exiting'
        stop
      elseif (beta(1).le.0.0d0) then
        print*,'elem505.f: beta(1) is less or equal to 0.0d0: exiting'
        stop
      elseif (beta(2).le.0.0d0) then
        print*,'elem505.f: beta(2) is less or equal to 0.0d0: exiting'
        stop
      elseif (Treg(1).le.0.0d0) then
        print*,'elem505.f: Treg(1) is less or equal to 0.0d0: exiting'
        stop
      elseif (Treg(2).le.0.0d0) then
        print*,'elem505.f: Treg(2) is less or equal to 0.0d0: exiting'
        stop
      endif
c     note: ireg is checked in itask.eq.7 for keeping things simple
c     note: s is checked in itask.eq.3 for keeping things simple
      pelem(1) = dble(l)
      pelem(2) = dble(ireg)
      pelem(3) = alpha(1)
      pelem(4) = alpha(2)
      pelem(5) = beta(1)
      pelem(6) = beta(2)
      pelem(7) = Treg(1)
      pelem(8) = Treg(2)
      pelem(9) = tauMult
      pelem(10) = dble(s)
      return

c     --------------------------------------------------------
c     Build the elemental consistent tangent stiffness matrix (estif)
c       and the elemental RHS/residual (eforc)
c     --------------------------------------------------------
 3    l      = int(pelem(1))
      tauMult= pelem(9)
      s      = int(pelem(10))
      ndim = ndime

c     determine the characteristic length h of the element
c      (for quadrilaterals: use the longer diagonal)
c      (for triangles: use the longest edge)
      if (nnode.eq.4) then
        h1=0.0d0
        h=0.0d0
        do i=1,ndim
          h1 = h1 + (xelem(i,3)-xelem(i,1))**2.0d0
          h = h + (xelem(i,4)-xelem(i,2))**2.0d0
        enddo
        if (h1 .gt. h) then
           h=h1
        endif
        h=sqrt(h)
      endif! nnode.eq.4
      
      if (nnode.eq.3) then
        h1=0.0d0
        h2=0.0d0
        h=0.0d0 
        do i=1, ndim
          h1 = h1 + (xelem(i,1)-xelem(i,2))**2.0d0
          h2 = h2 + (xelem(i,1)-xelem(i,3))**2.0d0
          h = h + (xelem(i,2)-xelem(i,3))**2.0d0
        enddo
        if (h1.gt.h) then
          h=h1
        endif
        if (h2.gt.h) then
          h=h2
        endif
        h=sqrt(h)
      endif! nnode.eq.3

c     initialize variables
      estif(:,:) = 0.0d0! elemental tangent stiffness
      eforc(:) = 0.0d0! elemental RHS
      Qq(:,:,:) = 0.0d0
      Qlin(:,:,:,:) = 0.0d0
      ident(1,1) = 1.0d0! identity matrix
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0  

      if (nnode.eq.3) then
         call gausstri(l,ninte,sg,tg,wg)
      else
         call gauss2(l,ninte,sg,tg,wg)
      endif
c
      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.false.,nnode,
     $      elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco
         
        mu = 0.0d0
        gamm = 0.0d0
        do inode = 1,nnode
          mu = mu + shap(3,inode)*elemdata_nodal(2,inode)
          gamm = gamm + shap(3,inode)*elemdata_nodal(1,inode)
        enddo

        if (s.eq.0) then
          temp  = 0.0d0
        else if (s.eq.1) then
          temp = (0.5d0*tauMult*(h**2.0d0))/mu
        else
          write(iwrit,200) 
 200      format(4x,'Stabilization property must be 0 or 1 ... exiting')
          stop
        endif

c       compute the deformation gradient at the current Gauss Point
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

c       compute pressure and spatial derivatives at Gauss Point
        pres(1:3) = 0.0d0
        do i = 1, 3
          do inode = 1, nnode
            pres(i) = pres(i) + uelem(3,inode)*shap(i,inode)
          end do
        end do

c       compute the Cauchy tensor and its inverse
        Ctens(1:2,1:2) = 0.0d0
        do j = 1,ndim
          do i = j,ndim! implement the symmetries
            do r = 1,ndim
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do
c       enforce the symmetries
        Ctens(1,2)=Ctens(2,1)

        Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)

        Cinv(1,1) = (1.0d0/Cdet)*Ctens(2,2)
        Cinv(2,2) = (1.0d0/Cdet)*Ctens(1,1)
        Cinv(1,2) = -(1.0d0/Cdet)*Ctens(1,2)
        Cinv(2,1) = Cinv(1,2)! = -(1.0d0/Cdet)*Ctens(2,1)

c       principal invariants Inv1 and Inv2
        Inv1 = Ctens(1,1)+Ctens(2,2)
        Inv2 = Cdet

c       dJdC is the derivative of the Jacobian with respect to 
c       the Cauchy Green tensor
        dJdC(1:2,1:2) = 0.5d0*Fdet*Cinv(1:2,1:2)

c       compute second Piola Kirchhoff stress SecPK 
c       and the material tangent Ctang
         
        K1 = (Fdet**(-2.0d0/3.0d0))*(Inv1+1.0d0)-3.0d0
        K2(1:2,1:2) =    (ident(1:2,1:2) - 
     $                (1.0d0/3.0d0)*Cinv(1:2,1:2)*(Inv1+1.0d0))
     $                 *(Fdet**(-2.0d0/3.0d0))


        dWdC(1:2,1:2) = 0.5d0*mu*exp(K1*gamm)*K2(1:2,1:2)

        SecPK(1:2,1:2) = 2.0d0*(dWdC(1:2,1:2)-pres(3)*dJdC(1:2,1:2))

        do i = 1, 2
          do j = 1, 2
            do q = 1, 2
              do r = 1, 2
                dCinvdC(i,j,q,r) = -0.5d0*(Cinv(i,q)*Cinv(j,r)+
     $                                   Cinv(i,r)*Cinv(j,q)) 
              end do
            end do
          end do
        end do         

        do i = 1, 2
          do j = 1, 2
            do q = 1, 2
              do r = 1, 2
                Ctang(i,j,q,r) = 0.5d0*mu*exp(K1*gamm)*
     $                           ((gamm*K2(i,j)*K2(q,r))-
     $                           (1.0d0/3.0d0)*K2(i,j)*Cinv(q,r)-
     $                           (1.0d0/3.0d0)*(Fdet**(-2.0d0/3.0d0))*
     $                           (dCinvdC(i,j,q,r)*(Inv1+1.0d0)
     $                             +Cinv(i,j)*ident(q,r)))
     
                d2JdC(i,j,q,r) = 0.25d0*Fdet*(Cinv(q,r)*Cinv(i,j)
     $                           +2.0d0*dCinvdC(i,j,q,r))
              end do
            end do
          end do
        end do         

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
     $                 +0.25d0*Lmat(2,1,2,1)
        Dgeomat(3,4) = 0.25d0*(Lmat(1,2,1,2)-Lmat(2,1,2,1))
        Dgeomat(4,4) = 0.25d0*Lmat(1,2,1,2)+0.25d0*Lmat(2,1,2,1)
     $                 -0.5d0*Lmat(1,2,2,1)
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
                end if

                if (idofn.eq.2  .AND. jdofn.eq.3) then
                   estif(ievab,jevab)=estif(ievab,jevab)-Fdet*
     $                   (shap(1,inode)*Finv(1,2)+
     $                    shap(2,inode)*Finv(2,2))*shap(3,jnode)*wtjac
                end if

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
                        end do
                      end do
                    end do
                  end do
                end if

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
                        end do
                      end do
                    end do
                  end do
                end if
                  
                if (idofn.eq.3  .AND. jdofn.eq.3) then
                  do i = 1,ndim
                    do j = 1,ndim
                      estif(ievab,jevab)=estif(ievab,jevab)+2.0d0*
     $                   temp*dJdC(i,j)*shap(i,jnode)*shap(j,inode)
     $                   *wtjac
                    end do
                  end do
                end if

                if (idofn.lt.3  .AND. jdofn.lt.3) then
                  do i = 1, 4
                    do j = 1, 4
                      estif(ievab,jevab)=estif(ievab,jevab)+
     $                         bb(inode,i,idofn)*Dgeomat(i,j)*
     $                         bb(jnode,j,jdofn)*wtjac
                    end do
                  end do
                end if
              end do! jdofn
            end do! jnode
          end do! idofn
        end do! inode

c       Linearization of the stabilization, divergence term
        if (nnode.ne.3) then
          if (s.eq.1) then
            jevab = 1
            ievab = 0
            do jnode = 1, nnode
              ievab = 0
              call DivLin2(ielem,jnode,Qlin,nnode,elnods,xelem,
     $               elemdata_nodal,uelem)
              DivQlin(:,:) = 0.0d0
              do i = 1, ndim
                do j = 1, ndim
                  do inode = 1, nnode
          DivQlin(1,i) = DivQlin(1,i) + Qlin(inode,1,i,j)*shap(j,inode)
          DivQlin(2,i) = DivQlin(2,i) + Qlin(inode,2,i,j)*shap(j,inode)
                  end do
                end do
              end do
     
              do knode = 1, nnode
                ievab = ievab + 3
                estif(ievab,jevab) = estif(ievab,jevab)-
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
    
            call DivInterp2 (ielem,Qq,nnode,elnods,xelem,elemdata_nodal,
     $             uelem)
            DivQq(1:2) = 0.0d0
            do i = 1, ndim
              do j = 1, ndim
                do inode = 1, nnode
                  DivQq(i) = DivQq(i) + Qq(inode,i,j)*shap(j,inode)
                end do
              end do
            end do

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
     
                      estif(ievab,jevab+1) = estif(ievab,jevab+1) +
     $                          temp*Finv(i,2)*Finv(l,j)*shap(l,knode)*
     $                              DivQq(j)*shap(i,jnode)* wtjac
                     end do
                   end do
                 end do
                 jevab = jevab + 3
              end do
            end do
          end if
        end if  ! nnode.ne.3

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
          end do
        end do

        Tp(3) = Fdet-1.0d0 

        if (nnode.ne.3) then
          if (iinte.eq.1 .AND. s.eq.1) then
            call DivInterp2 (ielem,Qq,nnode,elnods,xelem,elemdata_nodal,
     $             uelem)
          end if
          if (s.eq.1) then
c           compute the divergence of Qq after interpolation at Gauss Points
            DivQq(1:2) = 0.0d0
            do i = 1, ndim
              do j = 1, ndim
                do inode = 1, nnode
                  DivQq(i) = DivQq(i) + Qq(inode,i,j)*shap(j,inode)
                end do
              end do
            end do
            do i = 1, ndim
              do j = 1, ndim
                Tp(i) = Tp(i) - temp*Finv(i,j)*DivQq(j)
              end do
            end do
          end if
        end if   ! nnode.ne.3	 

        FS(:,:) = 0.0d0
        do i = 1, ndim
          do j = 1, ndim
            do r = 1, ndim
              FS(i,j) = FS(i,j)+Fdef(i,r)*SecPK(r,j)
            end do
          end do
        end do

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
              end do
            else
              do j = 1, 4
                eforc(ievab)=eforc(ievab)+
     $                          bb(inode,j,i)*T(j)*wtjac
              end do
            end if
          end do
        end do! inode

c	Following is for output purpose

c        do j=1, 12
c          Print*,"element force at",j, "  is", eforc(j)
c        end do
c        stop

      end do! iinte

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
 6    l  = int(pelem(1))
      ndim = ndime

      eforc(:) = 0.0d0! elemental RHS/residual

      if (nnode.eq.3) then
         call gausstri(l,ninte,sg,tg,wg)
      else
         call gauss2(l,ninte,sg,tg,wg)
      endif
c
      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,
     $               .false.,nnode,elnods,ndim,xelem)
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
     $                    shap(3,inode)*shap(3,jnode)
     $                    *uelem_diff(jdofn,jnode)*wtjac
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
 7    l       = int(pelem(1))
      ireg    = int(pelem(2))
      alpha(1)= pelem(3)
      alpha(2)= pelem(4)
      beta(1) = pelem(5)
      beta(2) = pelem(6)
      Treg(1) = pelem(7)
      Treg(2) = pelem(8)
      tauMult = pelem(9)
      s       = int(pelem(10))

      powe = 0.5d0! power of the power regularization (ireg.eq.4)

      egrad(:,:) = 0.0d0! initialize egrad

c     identity matrix ident
      ident(1,1) = 1.0d0
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0
      ndim = ndime

      elemDataMatch(ielem) = 0.0d0
      elemRegul(ielem) = 0.0d0
      l2grad1(ielem) = 0.0d0
      l2grad2(ielem) = 0.0d0

      temp_dual(:) = 0.0d0
      ii=0

      do inode = 1,nnode
        do idofn = 1,3! 3 = elemvec_ndofn(inode)
          ii = ii+1
          temp_dual(ii)   = uelem_dual(idofn,inode)
        enddo
      enddo

c     determine the characteristic length h of the element
c      (for quadrilaterals: use the longer diagonal)
c      (for triangles: use the longest edge)
      if (nnode.eq.4) then
        h1=0.0d0
        h=0.0d0
        do i=1, ndim
          h1 = h1 + (xelem(i,3)-xelem(i,1))**2.0d0
          h = h + (xelem(i,4)-xelem(i,2))**2.0d0
        enddo
        if (h1 .gt. h) then
          h=sqrt(h1)
        endif
      endif  !nnode.eq.4
      
      if (nnode.eq.3) then
        h1=0.0d0
        h2=0.0d0
        h=0.0d0 
        do i=1, ndim
          h1 = h1 + (xelem(i,1)-xelem(i,2))**2.0d0
          h2 = h2 + (xelem(i,1)-xelem(i,3))**2.0d0
          h = h + (xelem(i,2)-xelem(i,3))**2.0d0
        enddo
        if (h1.gt.h) then
          h=h1   
        endif
        if (h2.gt.h) then
          h=h2
        endif
        h=sqrt(h)
      endif        !nnode.eq.3      
      
      if (nnode.eq.3) then
        call gausstri(l,ninte,sg,tg,wg)
      else
        call gauss2(l,ninte,sg,tg,wg)
      endif
c
      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.false.,nnode,
     $        elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco

c       create prop_grad (contains prop gradient and value) 
        prop_grad(:,:)  = 0.0d0
        do iset = 1,2
          do idofn = 1,3!two gradents+the value
            do inode = 1,nnode
              prop_grad(iset,idofn) = prop_grad(iset,idofn)+
     $                            shap(idofn,inode)*
     $                      elemdata_nodal(iset,inode)
            end do
          end do
        end do

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

c       compute the spatial derivatives of the pressure and the pressure at the Gauss Point
c        (shap contain the derivative of the shape functions and its value)
        pres(1:3) = 0.0d0
        do i = 1, 3
          do inode = 1, nnode
            pres(i) = pres(i) + uelem(3,inode)*shap(i,inode)
          end do
        end do

c       compute the Cauchy tensor and its inverse
        Ctens(1:2,1:2) = 0.0d0
        do i = 1,ndim
          do j = 1,ndim
            do r = 1,ndim
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do

        Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)

        Inv2 = Cdet
        Inv1 = Ctens(1,1)+Ctens(2,2)

        Cinv(1,1) = (1.0d0/Cdet)*Ctens(2,2)
        Cinv(2,2) = (1.0d0/Cdet)*Ctens(1,1)
        Cinv(1,2) = -(1.0d0/Cdet)*Ctens(1,2)
        Cinv(2,1) = -(1.0d0/Cdet)*Ctens(2,1)

c       dJdC is the derivative of the Jacobian with respect to 
c       the Cauchy Green tensor
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
            end if 
            if (s.eq.1) then
              temp = 0.5d0*(tauMult*(h**2.0d0))/prop_grad(2,3)
              temp_der = 0.0d0
            end if
    
            if (s.eq.1 .AND. iset.eq.2) then  
              temp_der = -0.5d0*var_shape_prop(2,3)*
     $                (tauMult*(h**2.0d0))/(prop_grad(2,3)**2.0d0)
            end if

            Tp(1:4) = 0.0d0
            do i = 1, ndim
              do j = 1, ndim
                Tp(i) = Tp(i) + 2.0d0*temp_der*dJdC(j,i)*pres(j)
              end do
            end do

c           Tp(3) = 0.0d0 is zero here
            if (nnode.ne.3) then
              if (s.eq.1) then
                call DivInterpGrad2 (ielem,QqGrad,iset,nnode,elnods,
     $                 xelem,elemdata_nodal,uelem)
                call DivInterp2 (ielem,Qq,nnode,elnods,xelem,
     $                 elemdata_nodal,uelem)
              end if
              if (s.eq.1) then! JFD surprising .eq.1 cf above
c               Computing the Divergence of Qq after interpolation at Gauss Points

                DivQqGrad(1:2) = 0.0d0
                do i = 1, ndim
                  do j = 1, ndim
                    DivQqGrad(i) = DivQqGrad(i) + QqGrad(knode,i,j)
     $                        *shap(j,knode)
                  end do
                end do

                DivQq(1:2) = 0.0d0
                do i = 1, ndim
                  do j = 1, ndim
                    do inode = 1, nnode
                      DivQq(i) = DivQq(i) + Qq(inode,i,j)
     $                    *shap(j,inode)
                    end do
                  end do
                end do
                do i = 1, ndim
                  do j = 1, ndim
                    Tp(i) = Tp(i) - temp_der*Finv(i,j)*DivQq(j)
     $                     - temp*Finv(i,j)*DivQqGrad(j)
                  end do
                end do
              end if
            end if! nnode.ne.3
 
c           Compute derivative of second Piola Kirchhoff stress SecPK 
c           with respect to the material parameters
         
            K1 = (Fdet**(-2.0d0/3.0d0))*(Inv1+1.0d0)-3.0d0
            K2(1:2,1:2) = (ident(1:2,1:2) - 
     $                 (1.0d0/3.0d0)*Cinv(1:2,1:2)*(Inv1+1.0d0))
     $                  *(Fdet**(-2.0d0/3.0d0))

            if (iset.eq.1) then
              dWdC_grad(1:2,1:2) = 0.5d0*prop_grad(2,3)*
     $                    K1*var_shape_prop(iset,3)*
     $                    exp(K1*prop_grad(iset,3))*K2(1:2,1:2)
            elseif (iset.eq.2) then
              dWdC_grad(1:2,1:2) = 0.5d0*var_shape_prop(iset,3)
     $                      *exp(K1*prop_grad(1,3))*K2(1:2,1:2)
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
c                if (i.eq.elemvec_ndofn(inode)) then
                if (i.eq.3) then
                  do j = 1, 4
                    eforc(ievab)=eforc(ievab)+
     $                          bb(inode,j,i)*Tp(j)*wtjac
                  end do
                else
                  do j = 1, 4
                    eforc(ievab)=eforc(ievab)+
     $                          bb(inode,j,i)*T(j)*wtjac
                  end do
                end if
              end do
            end do! inode
            do ii = 1, 3*nnode
              egrad(iset,knode) = egrad(iset,knode) +
     $                            eforc(ii)*temp_dual(ii)
            end do

c           account for the regularization term
c JFD            if (iset.eq.1) then
              if (ireg.eq.0) then
c               there is not supposed to be a regularization term: do nothing
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
              else
                Print*,"elem505.f: ireg=",ireg,
     $                 " is not implemented: exiting"
                stop
              endif! ireg

c            elseif (iset.eq.2) then
c              if (ireg.eq.1) then
c                egrad(iset,knode) = egrad(iset,knode) + alpha(iset)*(
c     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
c     $              var_shape_prop(iset,2)*prop_grad(iset,2))*
c     $              wtjac
c              else if (ireg.eq.2) then
c                deno = sqrt(
c     $              beta(iset)*beta(iset)+
c     $              prop_grad(iset,1)*prop_grad(iset,1)+
c     $              prop_grad(iset,2)*prop_grad(iset,2))
c                egrad(iset,knode) = egrad(iset,knode) + 
c     $              0.5d0*alpha(iset)*(
c     $              var_shape_prop(iset,1)*prop_grad(iset,1)+
c     $              var_shape_prop(iset,2)*prop_grad(iset,2))*
c     $              wtjac/deno
c              endif
c            endif! iset
            udiff(iset) = udiff(iset)+shap(3,knode)*
     $                      uelem_diff(iset,knode)
          enddo! knode
        enddo! iset
      
        elemDataMatch(ielem) = elemDataMatch(ielem) + 0.5d0*wtjac*
     $       (udiff(1)*udiff(1)+udiff(2)*udiff(2))

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
        elseif (ireg.eq.4) then! logarithm reg.
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $       (alpha(2)*log(1+(sqrt(beta(2)*beta(2)+prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2))-
     $           beta(2))/Treg(2))+
     $       alpha(1)*log(1+(sqrt(beta(1)*beta(1)+prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2))-
     $           beta(1))/Treg(1)))
        elseif (ireg.eq.41) then! logarithm reg. (second implementation to avoid using a sqrt operator and the associated beta that has a significant influence on the convergence speed)
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $       (alpha(2)*log(1+(prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))/(Treg(2)*Treg(2)))+
     $       alpha(1)*log(1+(prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2))/(Treg(1)*Treg(1))))
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


c*******************************************************************************
      subroutine DivInterp2 (ielem,Qq,nno,elnods,xelem,elemd,uel)
c     Not updated for this material model - JFD ask sevan 
c     This subroutine is used when dealing with quadrilaterals and not triangles
c*******************************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none

      integer ii,i,j,r,inode,mnode_elem,ielem
      parameter (mnode_elem = 9)
      double precision K1,Inv1,Inv2,Cdet,xjaco,Fdet
      double precision Fdef(2,2),K2(2,2),Cinv(2,2),Ctens(2,2),ident(2,2)
      double precision ss(4), tt(4),dWdC(2,2)
      double precision shap(3,mnode_elem)
      double precision Qq(mnode_elem,2,2)
      integer nno
      integer elnods(*)
      double precision xelem(ndime,*)
      double precision elemd(nset_nodal,*)
      double precision uel(mdofn,*)
      
c     Last term in stabilization contains Divergence of a tensor which is spatially 
c      dependent. This will be linearly approximated to enable taking the Divergence

      data ss/-1.0d0,1.0d0,1.0d0,-1.0d0/,
     $     tt/-1.0d0,-1.0d0,1.0d0,1.0d0/

      ident(1,1) = 1.0d0! identity matrix
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0 

      Qq(:,:,:) = 0.0d0
      do inode = 1, nno
        call shape2(ielem,ss(inode),tt(inode),shap,xjaco,.false.,nno,
     $         elnods,ndime,xelem)

c       compute the deformation gradient at Gauss Point
        Fdef(1,1) = 1.0d0
        Fdef(2,2) = 1.0d0
        Fdef(1,2) = 0.0d0
        Fdef(2,1) = 0.0d0 
        do i = 1,ndime
          do j = 1,ndime
            do ii = 1,nno
              Fdef(i,j)=Fdef(i,j)+uel(i,ii)*shap(j,ii)
            end do
          end do
        end do
c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)

c       compute the Cauchy tensor and its inverse
        Ctens(1:2,1:2) = 0.0d0
        do i = 1,ndime
          do j = 1,ndime
            do r = 1,ndime
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do

        Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)

        Cinv(1,1) = (1.0d0/Cdet)*Ctens(2,2)
        Cinv(2,2) = (1.0d0/Cdet)*Ctens(1,1)
        Cinv(1,2) = -(1.0d0/Cdet)*Ctens(1,2)
        Cinv(2,1) = -(1.0d0/Cdet)*Ctens(2,1)

c       principal invariants Inv1 and Inv2
        Inv1 = Ctens(1,1)+Ctens(2,2)
        Inv2 = Cdet

        K1 = Inv1-2.0d0*log(Fdet)-2.0d0
        K2(1:2,1:2) = ident(1:2,1:2) - Cinv(1:2,1:2)   

        dWdC(1:2,1:2) = 0.5d0*elemd(2,inode)
     $                 *exp(K1*elemd(1,inode))*K2(1:2,1:2)

        do i = 1, ndime
          do j =  1, ndime
            do r = 1, ndime
              Qq(inode,i,j) = Qq(inode,i,j) + 2.0d0*Fdef(i,r)*dWdC(r,j)
            end do
          end do
        end do

      end do! inode
      end subroutine DivInterp2

c*******************************************************************************
      subroutine DivLin2(ielem,jnode,Qlin,nno,elnods,xelem,elemd,uel)
c     Not updated for this material model - JFD ask sevan
c     This subroutine is used when dealing with quadrilaterals and not triangles
c*******************************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none

      integer ii,i,j,r,q,l,m,inode,jnode,mnode_elem,ielem
      parameter (mnode_elem = 9)
      double precision K1,Inv1,Inv2,Cdet,xjaco
      double precision Fdef(2,2),K2(2,2),Cinv(2,2),Ctens(2,2),ident(2,2)
      double precision ss(4),tt(4),dWdC(2,2),Fdet,Ctang(2,2,2,2)
      double precision elem_nodal_gamma(mnode_elem),shap(3,mnode_elem)
      double precision Qlin(mnode_elem,2,2,2),elem_nodal_mu(mnode_elem)
      double precision dCinvdC(2,2,2,2)
      integer nno
      integer elnods(*)
      double precision xelem(ndime,*)
      double precision elemd(nset_nodal,*)
      double precision uel(mdofn,*)
      
      ident(1,1) = 1.0d0! identity matrix
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0

      data ss/-1.0d0,1.0d0,1.0d0,-1.0d0/,
     $     tt/-1.0d0,-1.0d0,1.0d0,1.0d0/

      Qlin(:,:,:,:) = 0.0d0
      do inode = 1, nno
        call shape2(ielem,ss(inode),tt(inode),shap,xjaco,.false.,nno,
     $       elnods,ndime,xelem)

c       compute the deformation gradient at Gauss Point
        Fdef(1,1) = 1.0d0
        Fdef(2,2) = 1.0d0
        Fdef(1,2) = 0.0d0
        Fdef(2,1) = 0.0d0
        do i = 1,ndime
          do j = 1,ndime
            do ii = 1,nno
              Fdef(i,j)=Fdef(i,j)+uel(i,ii)*shap(j,ii)
            end do
          end do
        end do
c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)

c       compute the Cauchy tensor and its inverse
        Ctens(1:2,1:2) = 0.0d0
        do i = 1,ndime
          do j = 1,ndime
            do r = 1,ndime
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do

        Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)

        Cinv(1,1) = (1.0d0/Cdet)*Ctens(2,2)
        Cinv(2,2) = (1.0d0/Cdet)*Ctens(1,1)
        Cinv(1,2) = -(1.0d0/Cdet)*Ctens(1,2)
        Cinv(2,1) = -(1.0d0/Cdet)*Ctens(2,1)

c       principal invariants Inv1 and Inv2
        Inv1 = Ctens(1,1)+Ctens(2,2)
        Inv2 = Cdet

        K1 = Inv1-2.0d0*log(Fdet)-2.0d0
        K2(1:2,1:2) = ident(1:2,1:2) - Cinv(1:2,1:2)   


        dWdC(1:2,1:2) = 0.5d0*elemd(2,inode)*
     $         exp(K1*elemd(1,inode))*K2(1:2,1:2)

     
        do i = 1, 2
          do j = 1, 2
            do q = 1, 2
              do r = 1, 2
                dCinvdC(i,j,q,r) = -0.5d0*(Cinv(i,q)*Cinv(j,r)+
     $                                   Cinv(i,r)*Cinv(j,q))
              end do
            end do
          end do
        end do

        do i = 1, 2
          do j = 1, 2
            do q = 1, 2
              do r = 1, 2
                Ctang(i,j,q,r) = 0.5d0*elemd(2,inode)*
     $                           exp(elemd(1,inode)*K1)*
     $                           (-dCinvdC(i,j,q,r)+
     $                           elemd(1,inode)*
     $                           (ident(i,j)-Cinv(i,j))*
     $                           (ident(q,r)-Cinv(q,r)))
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
  
        do i = 1, ndime
          do j =  1, ndime
            do m = 1, ndime
              do r = 1, ndime
                do q = 1, ndime
                  do l = 1, ndime
                    Qlin(inode,m,i,j) = Qlin(inode,m,i,j) + 
     $                                  4.0d0*Fdef(i,r)*Ctang(r,j,q,l)*
     $                                  Fdef(m,l)*shap(q,jnode)
                  end do
                end do
              end do
            end do
          end do
        end do
      end do! inode
      end subroutine DivLin2


c*******************************************************************************
      subroutine DivInterpGrad2 (ielem,QqGrad,iset,nno,elnods,xelem,
     $     elemd,uel)
c     Not updated for this material model  - JFD ask sevan
c     This subroutine is used when dealing with quadrilaterals and not triangles
c*******************************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none

      integer ii,i,j,r,inode,mnode_elem,ielem,iset
      parameter (mnode_elem = 9)
      double precision K1,Inv1,Inv2,Cdet,xjaco,Fdet
      double precision Fdef(2,2),K2(2,2),Cinv(2,2),Ctens(2,2),ident(2,2)
      double precision ss(4), tt(4),dWdC_grad(2,2)
      double precision shap(3,mnode_elem)
      double precision QqGrad(mnode_elem,2,2)
      integer nno
      integer elnods(*)
      double precision xelem(ndime,*)
      double precision elemd(nset_nodal,*)
      double precision uel(mdofn,*)
      
c     Last term in stabilization contains Divergence of a tensor which is spatially 
c      dependent. This will be linearly approximated to enable taking the Divergence

      data ss/-1.0d0,1.0d0,1.0d0,-1.0d0/,
     $     tt/-1.0d0,-1.0d0,1.0d0,1.0d0/

      ident(1,1) = 1.0d0! identity matrix
      ident(2,2) = 1.0d0
      ident(1,2) = 0.0d0
      ident(2,1) = 0.0d0 

      QqGrad(:,:,:) = 0.0d0

      do inode = 1, nno
        call shape2(ielem,ss(inode),tt(inode),shap,xjaco,.false.,nno,
     $         elnods,ndime,xelem)

c       compute the deformation gradient at Gauss Point
        Fdef(1,1) = 1.0d0
        Fdef(2,2) = 1.0d0
        Fdef(1,2) = 0.0d0
        Fdef(2,1) = 0.0d0 
        do i = 1,ndime
          do j = 1,ndime
            do ii = 1,nno
              Fdef(i,j)=Fdef(i,j)+uel(i,ii)*shap(j,ii)
            end do
          end do
        end do
c       Fdet is the Jacobian, determinant of the deformation gradient
        Fdet = Fdef(1,1)*Fdef(2,2)-Fdef(1,2)*Fdef(2,1)

c       compute the Cauchy tensor and its inverse
        Ctens(1:2,1:2) = 0.0d0
        do i = 1,ndime
          do j = 1,ndime
            do r = 1,ndime
              Ctens(i,j)= Ctens(i,j)+Fdef(r,i)*Fdef(r,j)
            end do
          end do
        end do

        Cdet = Ctens(1,1)*Ctens(2,2)-Ctens(1,2)*Ctens(2,1)

        Cinv(1,1) = (1.0d0/Cdet)*Ctens(2,2)
        Cinv(2,2) = (1.0d0/Cdet)*Ctens(1,1)
        Cinv(1,2) = -(1.0d0/Cdet)*Ctens(1,2)
        Cinv(2,1) = -(1.0d0/Cdet)*Ctens(2,1)

c       principal invariants Inv1 and Inv2
        Inv1 = Ctens(1,1)+Ctens(2,2)
        Inv2 = Cdet
         
        K1 = Inv1-2.0d0*log(Fdet)-2.0d0
        K2(1:2,1:2) = ident(1:2,1:2) - Cinv(1:2,1:2)

        if (iset.eq.2) then
          dWdC_grad(1:2,1:2) = 0.5d0*exp(K1*elemd(1,inode))
     $                          *K2(1:2,1:2)
        else if (iset.eq.1) then
          dWdC_grad(1:2,1:2) = 0.5d0*elemd(2,inode)
     $                *K1*exp(K1*elemd(1,inode))*K2(1:2,1:2)
        else
          Print*,"iset is wrong"
          stop
        end if

        do i = 1, ndime
          do j =  1, ndime
            do r = 1, ndime
              QqGrad(inode,i,j) = QqGrad(inode,i,j)+
     $                            2.0d0*Fdef(i,r)*dWdC_grad(r,j)
            end do
          end do
        end do

      end do! inode
      end subroutine DivInterpGrad2

