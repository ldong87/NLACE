c**********************************************************
      subroutine elem66 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     NOTES:
c     1) inverse element
c     2) 2D, plane strain, compressible, linear elasticity
c     3) accepts as properies (lambda and mu) or (gamma = lambda/mu and mu)
c     note: by fixing gamma, the element can become incompressible plane stress
c       JFD gamma is not the non-linear parameter
c     Modified by Olalekan Babaniyi
c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem
      parameter (mnode_elem = 9,mquad_elem=26)
      integer ielem, itask
      integer nGauss, ninte, iinte, inode, jnode, iset
      integer ievab, jevab, jdofn, idofn, ii, jj
      integer ireg, propset, logflag
      double precision xjaco, wtjac
      double precision xforc 
      double precision lambda, mu, lambda_ref, mu_ref
      double precision psi, phi ! in case of logarithimic approach, psi=log(mu/mu_ref), phi=log(lambda/lambda_ref)
      double precision alpha, beta, gamma, deno! JFD do no like gamma cf non-linear code (variable gamm)
      double precision shap(3,mnode_elem)
      double precision sg(mquad_elem),tg(mquad_elem),wg(mquad_elem)
      double precision dmatrix(3,3),bb(mnode_elem,3,2)
      double precision dmatrix_prop1(3,3),dmatrix_prop2(3,3)
      double precision log_prop1,log_prop2
      double precision temp_primal(3), temp_dual(3)
      double precision udiff(2), prop_grad(2,3), prop_grad1(2,3)
      double precision var_shape_prop(2,3)
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
c----------------------------------------------------------

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask
      
c     -----------------------------
c     initialize the elemental parameters
c     -----------------------------
 1    elemvec_ndofn(1:nnode) = 2
c      buildSymmetricMatrices = .true.
      return
     

c     --------------------------------------------------------
c     read and store the material properties
c     --------------------------------------------------------
c     nGauss     : no. of gauss pts/direction
c     ireg  : regularization type
c     alpha : regularization parameter
c     beta: : another regularization parameter JFD modify names here
c     propset: property set. 1-> lambda and mu; 
c                            2-> gamma and mu, where gamma = \lambda/mu
c     logflag: flag for log element (0/1:off/on)
c     lambda_ref: dimensional factor to scale lambda
c     mu_ref: dimensional factor to scale mu
 2    read(iread,*) nGauss,ireg,alpha,beta,propset,logflag,
     $     lambda_ref,mu_ref
      write(iwrit,205) nGauss,ireg,alpha,beta,propset,logflag,
     $     lambda_ref,mu_ref
 205  format(/,4x,'LINEAR ELASTICITY(2D,2DOF)      '/
     +       4x,'gauss pts/dir                     ',i10,/
     +       4x,'reg type(0/1/2:none/H1/TVD)       ',i10,/
     +       4x,'regularization parameter          ',e15.6,/
     +       4x,'extra parameter (TVD)             ',e15.6,/
     +       4x,'property set(1/2)                 ',i10,/
     +       4x,'flag for log element(0/1:off/on)  ',i10,/
     +       4x,'dimensional factor to scale lambda',e15.6,/
     +       4x,'dimensional factor to scale mu    ',e15.6)
c      nprop = 8
c     check for error in input file - JFD to do
      pelem(1) = dble(nGauss)
      pelem(2) = dble(ireg)
      pelem(3) = alpha
      pelem(4) = beta
      pelem(5) = dble(propset)
      pelem(6) = dble(logflag)
      pelem(7) = lambda_ref
      pelem(8) = mu_ref
      return


c     --------------------------------------------------------
c     Build the elemental consistent stiffness matrix (estif)
c       and the elemental RHS/internal force (eforc)
c     --------------------------------------------------------
 3    nGauss  = int(pelem(1))
      propset = int(pelem(5))
      logflag = int(pelem(6))
      lambda_ref = pelem(7)
      mu_ref = pelem(8)
      ndim = ndime

      estif(:,:) = 0.0d0! elemental tangent stiffness
      eforc(:) = 0.0d0! elemental RHS/residual

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
         
        lambda = 0.0d0
        mu     = 0.0d0
        gamma  = 0.0d0
        psi    = 0.0d0
        phi    = 0.0d0

        if (propset.eq.1) then
           if (logflag.eq.0) then
              do inode = 1,nnode
                 lambda = lambda + shap(3,inode)*elemdata_nodal(1,inode)
                 mu     = mu     + shap(3,inode)*elemdata_nodal(2,inode)
              enddo
           else if (logflag.eq.1) then
              do inode = 1,nnode
                 phi    = phi + shap(3,inode)*elemdata_nodal(1,inode)
                 psi    = psi + shap(3,inode)*elemdata_nodal(2,inode)   
              enddo
              lambda = lambda_ref*exp(phi)
              mu = mu_ref*exp(psi)
           endif                 
        else if (propset.eq.2) then
           if (logflag.eq.0) then
              do inode = 1,nnode
                 gamma  = gamma + shap(3,inode)*elemdata_nodal(1,inode)
                 mu     = mu    + shap(3,inode)*elemdata_nodal(2,inode)
              enddo
              lambda = gamma*mu
           else if (logflag.eq.1) then
              do inode = 1,nnode
                 phi    = phi + shap(3,inode)*elemdata_nodal(1,inode) 
                 psi    = psi + shap(3,inode)*elemdata_nodal(2,inode)  
              enddo
              gamma = lambda_ref*exp(phi)
              mu = mu_ref*exp(psi)
              lambda = gamma*mu
           endif               
        endif

c       define D-matrix
        dmatrix(:,:) = 0.0d0
        dmatrix(1,1) = lambda+2.0d0*mu
        dmatrix(2,2) = lambda+2.0d0*mu
        dmatrix(3,3) = mu
        dmatrix(1,2) = lambda
        dmatrix(2,1) = lambda
         
c       create  the b-matrix
        do inode = 1,nnode
          bb(inode,1,1) = shap(1,inode)
          bb(inode,1,2) = 0.0d0
          bb(inode,2,1) = 0.0d0
          bb(inode,2,2) = shap(2,inode)
          bb(inode,3,1) = shap(2,inode)
          bb(inode,3,2) = shap(1,inode)
        enddo

c       construct the stiffness matrix 
        ievab = 0
        do inode = 1, nnode
          do idofn = 1,2!2 = elemvec_ndofn(inode)
            ievab = ievab +1 
            jevab = 0
            do  jnode = 1, nnode
              do jdofn = 1,2!2 = elemvec_ndofn(jnode)
                jevab = jevab+1
                do ii = 1,3
                  do jj = 1,3
                    estif(ievab,jevab) = estif(ievab,jevab)+
     $                          bb(inode,ii,idofn)*dmatrix(ii,jj)*
     $                          bb(jnode,jj,jdofn)*wtjac
                  enddo
                enddo
              enddo! jdofn
            enddo! jnode
          enddo! idofn
        enddo! inode
      enddo! iinte

c     compute element internal force vector
      do ievab =  1,nnode*2!nevab! nevab does not exist as a global variable for openmp JFD
        jevab = 0
        do jnode = 1,nnode
          do jdofn = 1,2!2 = elemvec_ndofn(jnode)
            jevab = jevab + 1
            xforc = estif(ievab,jevab)*uelem(jdofn,jnode)
            eforc(ievab) = eforc(ievab) + xforc
          enddo
        enddo
      enddo      
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
          do idofn = 1,2!2 = elemvec_ndofn(inode)
            ievab = ievab+1
            do jnode = 1, nnode
              do jdofn = 1,2!2 = elemvec_ndofn(jnode)
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
c     and its gradient (egrad) on the element
c     --------------------------------------------------------
 7    nGauss  = int(pelem(1))
      ireg    = int(pelem(2))
      alpha   = pelem(3)
      beta    = pelem(4)
      propset = int(pelem(5))
      logflag = int(pelem(6))
      lambda_ref = pelem(7)
      mu_ref = pelem(8)
      ndim = ndime

      elemDataMatch(ielem) = 0.0d0
      elemRegul(ielem) = 0.0d0
      l2grad1(ielem) = 0.0d0
      l2grad2(ielem) = 0.0d0
      egrad(:,:) = 0.0d0! initialize egrad

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
         
c       create  the b-matrix
        do inode = 1,nnode
          bb(inode,1,1) = shap(1,inode)
          bb(inode,1,2) = 0.0d0
          bb(inode,2,1) = 0.0d0
          bb(inode,2,2) = shap(2,inode)
          bb(inode,3,1) = shap(2,inode)
          bb(inode,3,2) = shap(1,inode)
        enddo

c       create Bu and Bw
        do ii = 1,3
          temp_primal(ii) = 0.0d0
          temp_dual(ii) = 0.0d0

          do inode = 1,nnode
            do idofn = 1,2!2 = elemvec_ndofn(inode)
              temp_primal(ii) = temp_primal(ii)+
     $                 bb(inode,ii,idofn)*uelem(idofn,inode)
              temp_dual(ii) = temp_dual(ii)+
     $                 bb(inode,ii,idofn)*uelem_dual(idofn,inode)
            enddo
          enddo
        enddo

c       create prop_grad and prop_grad1 (contains prop gradient and value) 
        prop_grad(1:2,1:3)  = 0.0d0
        prop_grad1(1:2,1:3) = 0.0d0
        do iset = 1,2
          do idofn = 1,3      !two gradents+the value
            do inode = 1,nnode
              prop_grad(iset,idofn) = prop_grad(iset,idofn)+
     $                 shap(idofn,inode)*
     $                 elemdata_nodal(iset,inode)
              prop_grad1(iset,idofn) = prop_grad1(iset,idofn)+
     $                 shap(idofn,inode)*
     $                 elemdata_nodal(iset,inode)
            enddo
c        compute lambda and mu if user is working with phi and psi as optimization variables
            if ((logflag.eq.1) .and. (iset.eq.1)) then
               prop_grad1(iset,idofn) = lambda_ref
     $              *exp(prop_grad1(iset,idofn))
            else if ((logflag.eq.1) .and. (iset.eq.2)) then
               prop_grad1(iset,idofn) = mu_ref
     $              *exp(prop_grad1(iset,idofn))
            end if
          enddo
        enddo


c         write(*,*)'propgrad',ielem,prop_grad(1,1),prop_grad(1,2),
c     $        prop_grad(2,1),prop_grad(2,2)
c
c         write(*,*)'elemdata',ielem,
c     $        elemdata_nodal(1,1),
c     $        elemdata_nodal(1,2),
c     $        elemdata_nodal(2,1),
c     $        elemdata_nodal(2,2)
c
c         write(*,*)'nod_grad',ielem,
c     $        ielem_nod_grad(1,1),
c     $        ielem_nod_grad(1,2),
c     $        ielem_nod_grad(2,1),
c     $        ielem_nod_grad(2,2)
c
c         write(*,*)'propgrad',ielem,prop_grad(1,1),prop_grad(1,2),
c     $        prop_grad(2,1),prop_grad(2,2)

c     create log_prop to help calculate derivaive of D-matrices
        phi = 0.0d0
        psi = 0.0d0

        if (logflag.eq.0) then
           log_prop1 = 1
           log_prop2 = 1
        else if (logflag.eq.1) then
           do inode = 1,nnode
              phi = phi + shap(3,inode)*elemdata_nodal(1,inode)
              psi = psi + shap(3,inode)*elemdata_nodal(2,inode)
           enddo
           log_prop1 = lambda_ref*exp(phi)
           log_prop2 = mu_ref*exp(psi)   
        end if

c     compute the gradient and the functional 
        udiff(:) = 0.0d0

        do inode = 1,nnode
c         set variations in the shape function for the property 
          var_shape_prop(:,:) = 0.0d0
          do iset = 1,2
            if (ielem_nod_grad(iset,inode).eq.1) then
              do idofn = 1,3
                var_shape_prop(iset,idofn) = 
     $                    shap(idofn,inode)
              enddo
            endif
          enddo

c         create the derivative of  D-matrices
          dmatrix_prop1(:,:) = 0.0d0
          dmatrix_prop2(:,:) = 0.0d0
          if (propset.eq.1) then
            dmatrix_prop1(1,1) = var_shape_prop(1,3)*log_prop1
            dmatrix_prop1(2,2) = var_shape_prop(1,3)*log_prop1
            dmatrix_prop1(1,2) = var_shape_prop(1,3)*log_prop1
            dmatrix_prop1(2,1) = var_shape_prop(1,3)*log_prop1

            dmatrix_prop2(1,1) = 2.0d0*var_shape_prop(2,3)*log_prop2
            dmatrix_prop2(2,2) = 2.0d0*var_shape_prop(2,3)*log_prop2
            dmatrix_prop2(3,3) = var_shape_prop(2,3)*log_prop2
          elseif (propset.eq.2) then
            dmatrix_prop1(1,1) = prop_grad1(2,3)*var_shape_prop(1,3)
            dmatrix_prop1(2,2) = prop_grad1(2,3)*var_shape_prop(1,3)
            dmatrix_prop1(1,2) = prop_grad1(2,3)*var_shape_prop(1,3)
            dmatrix_prop1(2,1) = prop_grad1(2,3)*var_shape_prop(1,3)

            dmatrix_prop2(1,1) = (2.0d0*var_shape_prop(2,3)+
     $            prop_grad1(1,3)*var_shape_prop(2,3))*log_prop2
            dmatrix_prop2(2,2) = (2.0d0*var_shape_prop(2,3)+
     $            prop_grad1(1,3)*var_shape_prop(2,3))*log_prop2
            dmatrix_prop2(3,3) = var_shape_prop(2,3)*log_prop2
            dmatrix_prop2(1,2) = prop_grad1(1,3)*var_shape_prop(2,3)
     $           *log_prop2
            dmatrix_prop2(2,1) = prop_grad1(1,3)*var_shape_prop(2,3)
     $           *log_prop2
          endif

          do ii = 1,3
            do jj = 1,3
              egrad(1,inode) = egrad(1,inode)+
     $                 temp_dual(ii)*dmatrix_prop1(ii,jj)
     $                 *temp_primal(jj)*wtjac
              egrad(2,inode) = egrad(2,inode)+
     $                 temp_dual(ii)*dmatrix_prop2(ii,jj)
     $                 *temp_primal(jj)*wtjac
            enddo
          enddo
            
c         account for the regularization term
          if (ireg.eq.0) then
c           there is not supposed to be a regularization: do nothing
          elseif (ireg.eq.1) then
            egrad(1,inode) = egrad(1,inode) + alpha*(
     $              var_shape_prop(1,1)*prop_grad(1,1)+
     $              var_shape_prop(1,2)*prop_grad(1,2))*
     $              wtjac
               
            egrad(2,inode) = egrad(2,inode) + alpha*(
     $              var_shape_prop(2,1)*prop_grad(2,1)+
     $              var_shape_prop(2,2)*prop_grad(2,2))*
     $              wtjac
          elseif ((ireg.eq.2).or.(ireg.eq.21)) then
            deno = dsqrt(
     $              beta*beta+
     $              prop_grad(1,1)*prop_grad(1,1)+
     $              prop_grad(1,2)*prop_grad(1,2)+
     $              prop_grad(2,1)*prop_grad(2,1)+
     $              prop_grad(2,2)*prop_grad(2,2))
            egrad(1,inode) = egrad(1,inode) + 
     $              0.5d0*alpha*(
     $              var_shape_prop(1,1)*prop_grad(1,1)+
     $              var_shape_prop(1,2)*prop_grad(1,2))*
     $              wtjac/deno
     $              
            egrad(2,inode) = egrad(2,inode) + 
     $              0.5d0*alpha*(
     $              var_shape_prop(2,1)*prop_grad(2,1)+
     $              var_shape_prop(2,2)*prop_grad(2,2))*
     $              wtjac/deno
          else
            print*,"ireg has a wrong value ... exiting"
            stop
          endif
            
          udiff(1) = udiff(1)+ shap(3,inode)*uelem_diff(1,inode)
          udiff(2) = udiff(2)+ shap(3,inode)*uelem_diff(2,inode)
        enddo! inode

        elemDataMatch(ielem) = elemDataMatch(ielem) + 0.5d0*wtjac*
     $        (udiff(1)*udiff(1)+udiff(2)*udiff(2))

        if (ireg.eq.1) then
          elemRegul(ielem) = elemRegul(ielem)+ 
     $           0.5d0*alpha*(
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))*
     $           wtjac
        elseif (ireg.eq.2) then
          elemRegul(ielem) = elemRegul(ielem)+ 
     $           0.5d0*alpha*dsqrt(beta*beta+
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))*
     $           wtjac
        elseif (ireg.eq.21) then
          elemRegul(ielem) = elemRegul(ielem)+ 
     $           0.5d0*alpha*(dsqrt(beta*beta+
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(2,1)*prop_grad(2,1)+
     $           prop_grad(2,2)*prop_grad(2,2))-beta)*
     $           wtjac
        endif
        l2grad1(ielem) = l2grad1(ielem) + sqrt(prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2))
        l2grad2(ielem) = l2grad2(ielem) + sqrt(prop_grad(2,1)*
     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2))
      enddo! iinte
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
