c**********************************************************
      subroutine elem77 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     Written by Sevan Goenezen
c     NOTES:
c     1) inverse element
c     2) 2D, plane strain, incompressible, linear elasticity
c     3) Added stabilization

c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem
      parameter (mnode_elem = 9,mquad_elem=26)
      integer ielem, itask
      integer nGauss, ninte, iinte, inode, jnode,knode, iset
c ij is added new
      integer ievab, jevab, jdofn, idofn, ii, jj, ij, i, s 
      integer ireg, propset
      double precision xjaco, wtjac
      double precision xforc 
      double precision mu, mux, muy
      double precision alpha, beta, gamma, deno! JFD do no like gamma cf non-linear element (gamm variable)
      double precision shap(3,mnode_elem), shapehod(3,mnode_elem)
      double precision sg(mquad_elem),tg(mquad_elem),wg(mquad_elem)
      double precision temp_primal(3*mnode_elem),temp_dual(3*mnode_elem)
      double precision udiff(2), prop_grad(1,3), var_shape_prop(1,3)
      double precision alpha2, h, h1, temp, temp1
      double precision temp_der, temp1_der
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
c----------------------------------------------------------

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask
      
c     -----------------------------
c     initialize the elemental parameters
c     -----------------------------
 1    elemvec_ndofn(1:nnode) = 3
      buildSymmetricMatrices = .false.
      return
      

c     --------------------------------------------------------
c     read and store the material properties
c     --------------------------------------------------------
c     nGauss: no. of gauss pts/direction
c     ireg  : regularization type
c     alpha : regularization parameter
c     beta  : another regularization parameter JFD modify names here
c     alpha2: stabilization factor
c     s     : type of stabilization
 2    read(iread,*) nGauss,ireg,alpha,beta, alpha2, s
      write(iwrit,205) nGauss,ireg,alpha,beta, alpha2, s
 205  format(/,4x,'LINEAR ELASTICITY(2D,2DOF) '/
     +       4x,'gauss pts/dir              ',i10,/
     +       4x,'reg type(0/1/2:none/H1/TVD)',i10,/
     +       4x,'regularization parameter   ',e15.6,/
     +       4x,'extra parameter (TVD)      ',e15.6,/
     +       4x,'stabilization factor    ',e15.6,/
     +       4x,'stabilization terms 0/1/2 off/laplacian stabilization
     $ /stabilization with higher order terms   ',i10)

c      nprop = 6
      pelem(1) = dble(nGauss)
      pelem(2) = dble(ireg)
      pelem(3) = alpha
      pelem(4) = beta
      pelem(5) = alpha2
      pelem(6) = dble(s)
      return


c     --------------------------------------------------------
c     Build the elemental consistent stiffness matrix (estif)
c       and the elemental RHS/internal force (eforc)
c     --------------------------------------------------------
 3    nGauss  = int(pelem(1))
      alpha2  = pelem(5)
      s       = int(pelem(6))
      ndim = ndime

      estif(:,:) = 0.0d0! elemental tangent stiffness
      eforc(:) = 0.0d0! elemental RHS/residual

c     Determine the characteristic element length h for bilinear elements 
c      compute both diagonals and take the bigger one
      h1=0.0d0
      h=0.0d0
      do i=1, ndim
        h1= h1+ (xelem(i,3)-xelem(i,1))**2.0d0
        h = h + (xelem(i,4)-xelem(i,2))**2.0d0
      enddo
      if (h1 .gt. h) then
        h=dsqrt(h1)
      else
        h=dsqrt(h)
      endif

      if (nnode.eq.3) then
        print*,"elem77.f: nnode.eq.3 is not implemented yet: exiting"
        stop
        call gausstri(nGauss,ninte,sg,tg,wg)
      else
        call gauss2(nGauss,ninte,sg,tg,wg)
      endif

      do iinte = 1, ninte
        call shapehigher_deriv(ielem,sg(iinte),tg(iinte),shap,
     $           shapehod,xjaco,.false.,nnode,elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco
         
        mu  = 0.0d0
        mux = 0.0d0
        muy = 0.0d0
        do inode = 1,nnode
          mu     = mu     + shap(3,inode)*elemdata_nodal(1,inode)
          mux     = mux   + shap(1,inode)*elemdata_nodal(1,inode)
          muy     = muy   + shap(2,inode)*elemdata_nodal(1,inode)
        enddo

        if (s.eq.0) then
          temp  = 0.0d0 ! JFD choose other name than temp...
          temp1 = 0.0d0 ! JFD choose other name than temp1...
        else if (s.eq.1) then
          temp  = 0.0d0
          temp1 = (alpha2*(h**2.0d0))/mu
        else if (s.eq.2) then
          temp  = (alpha2*(h**2.0d0))/mu
          temp1 = (alpha2*(h**2.0d0))/mu
        else
          write(iwrit,200) 
 200      format(4x,'Stabilization property must be 0, 1, or 2')
          stop
        endif

c       construct the stiffness matrix 
        do inode = 1, nnode
          do  jnode = 1, nnode
            estif((inode*3)-2,(jnode*3)-2) = 
     $        estif((inode*3)-2,(jnode*3)-2)+
     $        2*mu*(shap(1,inode)*shap(1,jnode)+
     $        0.5*shap(2,inode)*shap(2,jnode))*wtjac

            estif((inode*3)-2,(jnode*3)-1) =
     $               estif((inode*3)-2,(jnode*3)-1)+ 
     $               mu*shap(2,inode)*shap(1,jnode)*wtjac
    
            estif((inode*3)-2,jnode*3) =
     $                            estif((inode*3)-2,jnode*3)
     $                                    -shap(3,jnode)*
     $                                      shap(1,inode)*wtjac

            estif((inode*3)-1,(jnode*3)-2) = 
     $            estif((inode*3)-1,(jnode*3)-2)+
     $            mu*shap(1,inode)*shap(2,jnode)*wtjac

            estif((inode*3)-1,(jnode*3)-1) =
     $           estif((inode*3)-1,(jnode*3)-1)+
     $           2*mu*(shap(2,inode)*
     $           shap(2,jnode)+0.5*shap(1,inode)*shap(1,jnode))*wtjac
                 
            estif((inode*3)-1,jnode*3) =
     $           estif((inode*3)-1,jnode*3)
     $                                     -shap(3,jnode)*
     $                                      shap(2,inode)*wtjac

            estif((inode*3),(jnode*3)-2) =
     $            estif((inode*3),(jnode*3)-2)+
     $            (shap(3,inode)*
     $            shap(1,jnode)-temp*mu*(0.5*shap(2,inode)*
     $            shapehod(3,jnode)+shap(1,inode)*(shapehod(1,jnode)+
     $            0.5*shapehod(2,jnode)))-
     $            temp*(shap(1,inode)*(shap(1,jnode)*mux+
     $            0.5*shap(2,jnode)*muy)+0.5*shap(2,inode)*
     $            shap(2,jnode)*mux))*wtjac

            estif((inode*3),(jnode*3)-1) = 
     $            estif((inode*3),(jnode*3)-1)+
     $            (shap(3,inode)*
     $            shap(2,jnode)-temp*mu*(0.5*shap(1,inode)*
     $            shapehod(3,jnode)+shap(2,inode)*(shapehod(2,jnode)+
     $            0.5*shapehod(1,jnode)))-
     $            temp*(shap(2,inode)*(shap(2,jnode)*muy+
     $            0.5*shap(1,jnode)*mux)+0.5*shap(1,inode)*
     $            shap(1,jnode)*muy))*wtjac

            estif((inode*3),(jnode*3)) = 
     $           estif((inode*3),(jnode*3))+
     $           0.5d0*temp1*(shap(1,inode)*
     $           shap(1,jnode)+shap(2,inode)*shap(2,jnode))*wtjac
               
          enddo! jnode
        enddo! inode
      enddo! iinte

c     compute element internal force vector
      do ievab =  1,nnode*3!nevab ! nevab does not exist as it was a global variable
        jevab = 0
        do jnode = 1,nnode
          do jdofn = 1,3! 3 = elemvec_ndofn(jnode)
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
 6    nGauss = int(pelem(1))
      ndim = ndime
      eforc(:) = 0.0d0! elemental RHS/residual

      if (nnode.eq.3) then
        call gausstri(nGauss,ninte,sg,tg,wg)
      else
        call gauss2(nGauss,ninte,sg,tg,wg)
      endif

      do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,
     $               .false.,nnode,elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco

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
      enddo! iinte
      return

c     --------------------------------------------------------
c     Compute the objective function (dataMatch+regularization)
c      and its gradient (egrad) on the element
c     --------------------------------------------------------
 7    nGauss  = int(pelem(1))
      ireg    = int(pelem(2))
      alpha   = pelem(3)
      beta    = pelem(4)
      alpha2  = pelem(5)
      s       = int(pelem(6))
      ndim = ndime
      
      elemDataMatch(ielem) = 0.0d0
      elemRegul(ielem) = 0.0d0
      l2grad1(ielem) = 0.0d0
      l2grad2(ielem) = 0.0d0
      egrad(:,:) = 0.0d0! initialize egrad

      temp_primal(:) = 0.0d0
      temp_dual(:) = 0.0d0
      ii=0

      do inode = 1,nnode
        do idofn = 1,3! 3 = elemvec_ndofn(inode)
          ii = ii+1
          temp_primal(ii) = uelem(idofn,inode)
          temp_dual(ii)   = uelem_dual(idofn,inode)
        enddo
      enddo

c     determine the characteristic element length h for bilinear elements 
c      compute both diagonals and take the bigger one
      h1=0.0d0
      h=0.0d0
      do i=1, ndim
        h1= h1+ (xelem(i,3)-xelem(i,1))**2.0d0
        h = h + (xelem(i,4)-xelem(i,2))**2.0d0
      enddo
      if (h1 .gt. h) then
        h=dsqrt(h1)
      else
        h=dsqrt(h)
      endif

c     compute the gradient and the objective function value
      if (nnode.eq.3) then
        call gausstri(nGauss,ninte,sg,tg,wg)
      else
        call gauss2(nGauss,ninte,sg,tg,wg)
      endif
c
      do iinte = 1, ninte
        call shapehigher_deriv(ielem,sg(iinte),tg(iinte),shap,
     $           shapehod,xjaco,.false.,nnode,elnods,ndim,xelem)
        wtjac = wg(iinte)*xjaco
 
c       create prop_grad (contains prop gradient and value) 
        prop_grad(1,1:3)  = 0.0d0
        do idofn = 1,3! two gradents+the value
          do inode = 1,nnode
            prop_grad(1,idofn) = prop_grad(1,idofn)+
     $                 shap(idofn,inode)*
     $                 elemdata_nodal(1,inode)
          enddo
        enddo

c       compute the gradient and the objective function value
        udiff(:) = 0.0d0

        do knode = 1,nnode
c         set variations in the shape function for the property 
          var_shape_prop(:,:) = 0.0d0
          estif(:,:) = 0.0d0
          if (ielem_nod_grad(1,knode).eq.1) then
            do idofn = 1,3
              var_shape_prop(1,idofn) = shap(idofn,knode)
            enddo
          endif
c         take out pf knode loop the following! JFD ask sevan
          mu  = 0.0d0
          mux = 0.0d0
          muy = 0.0d0
          do inode = 1,nnode
            mu     = mu     + shap(3,inode)*elemdata_nodal(1,inode)
            mux     = mux   + shap(1,inode)*elemdata_nodal(1,inode)
            muy     = muy   + shap(2,inode)*elemdata_nodal(1,inode)
          enddo

          if (s.eq.0) then
            temp  = 0.0d0
            temp1 = 0.0d0
            temp_der = 0.0d0
            temp1_der = 0.0d0
          elseif (s.eq.1) then
            temp  = 0.0d0
            temp1 = (alpha2*(h**2.0d0))/mu
            temp_der = 0.0d0
            temp1_der = -var_shape_prop(1,3)*(alpha2*(h**2.0d0))/(mu**2)
          elseif (s.eq.2) then
            temp = (alpha2*(h**2.0d0))/mu
            temp1 = (alpha2*(h**2.0d0))/mu
            temp_der  = -var_shape_prop(1,3)*(alpha2*(h**2.0d0))/(mu**2)
            temp1_der = -var_shape_prop(1,3)*(alpha2*(h**2.0d0))/(mu**2)
          else
            write(iwrit,300) 
 300        format(4x,'Stabilization property must be 0, 1, or 2')
            stop
          endif

          do inode = 1, nnode     
            do  jnode = 1, nnode
              estif((inode*3)-2,(jnode*3)-2) = 
     $          estif((inode*3)-2,(jnode*3)-2)+
     $          2*var_shape_prop(1,3)*(shap(1,inode)*shap(1,jnode)+
     $          0.5*shap(2,inode)*shap(2,jnode))*wtjac

              estif((inode*3)-2,(jnode*3)-1) =
     $           estif((inode*3)-2,(jnode*3)-1)+ 
     $           var_shape_prop(1,3)*shap(2,inode)*shap(1,jnode)*wtjac
    
              estif((inode*3)-2,jnode*3) = 0.0d0

              estif((inode*3)-1,(jnode*3)-2) = 
     $           estif((inode*3)-1,(jnode*3)-2)+
     $           var_shape_prop(1,3)*shap(1,inode)*shap(2,jnode)*wtjac

              estif((inode*3)-1,(jnode*3)-1) =
     $           estif((inode*3)-1,(jnode*3)-1)+
     $           2*var_shape_prop(1,3)*(shap(2,inode)*
     $           shap(2,jnode)+0.5*shap(1,inode)*shap(1,jnode))*wtjac
                 
              estif((inode*3)-1,jnode*3) = 0.0d0

              estif((inode*3),(jnode*3)-2) =
     $          estif((inode*3),(jnode*3)-2)-
     $          temp_der*(shap(1,inode)*(shap(1,jnode)*mux+
     $          0.5*shap(2,jnode)*muy)+0.5*shap(2,inode)*
     $          shap(2,jnode)*mux)*wtjac-
     $          temp*(shap(1,inode)*(shap(1,jnode)*var_shape_prop(1,1)+
     $          0.5*shap(2,jnode)*var_shape_prop(1,2))+0.5*
     $          shap(2,inode)*shap(2,jnode)*var_shape_prop(1,1))*wtjac

              estif((inode*3),(jnode*3)-1) = 
     $          estif((inode*3),(jnode*3)-1)-
     $          temp_der*(shap(2,inode)*(shap(2,jnode)*muy+
     $          0.5*shap(1,jnode)*mux)+0.5*shap(1,inode)*
     $          shap(1,jnode)*muy)*wtjac-
     $          temp*(shap(2,inode)*(shap(2,jnode)*var_shape_prop(1,2)+
     $          0.5*shap(1,jnode)*var_shape_prop(1,1))+0.5*
     $          shap(1,inode)*shap(1,jnode)*var_shape_prop(1,2))
     $          *wtjac

              estif((inode*3),(jnode*3)) = 
     $           estif((inode*3),(jnode*3))+
     $           0.5d0*temp1_der*(shap(1,inode)*
     $           shap(1,jnode)+shap(2,inode)*shap(2,jnode))*wtjac
               
            enddo! jnode
          enddo! inode

          do ii = 1, nnode*3!nevab! nevab was a global variable
            do jj = 1, nnode*3!nevab! nevab was a global variable
              egrad(1,knode) = egrad(1,knode)+
     $                 temp_dual(ii)*estif(ii,jj)*temp_primal(jj)
            enddo
          enddo
            
c         account for the regularization term
          if (ireg.eq.1) then
            egrad(1,knode) = egrad(1,knode) + alpha*(
     $              var_shape_prop(1,1)*prop_grad(1,1)+
     $              var_shape_prop(1,2)*prop_grad(1,2))*
     $              wtjac
          elseif (ireg.eq.2) then
            deno = dsqrt(
     $              beta*beta+
     $              prop_grad(1,1)*prop_grad(1,1)+
     $              prop_grad(1,2)*prop_grad(1,2))
            egrad(1,knode) = egrad(1,knode) + 
     $              0.5d0*alpha*(
     $              var_shape_prop(1,1)*prop_grad(1,1)+
     $              var_shape_prop(1,2)*prop_grad(1,2))*
     $              wtjac/deno
          endif
          udiff(1) = udiff(1)+ shap(3,knode)*uelem_diff(1,knode)
          udiff(2) = udiff(2)+ shap(3,knode)*uelem_diff(2,knode)
        enddo! knode

        elemDataMatch(ielem) = elemDataMatch(ielem) + 0.5d0*wtjac*
     $        (udiff(1)*udiff(1)+udiff(2)*udiff(2))

        if (ireg.eq.1) then
          elemRegul(ielem) = elemRegul(ielem)+ 
     $           0.5d0*alpha*(
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2))*
     $           wtjac
        elseif (ireg.eq.2) then
          elemRegul(ielem) = elemRegul(ielem)+ 
     $           0.5d0*alpha*dsqrt(beta*beta+
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2))*
     $           wtjac
        endif
        l2grad1(ielem) = l2grad1(ielem) + sqrt(prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2))
        l2grad2(ielem) = 0.0d0
c        l2grad2(ielem) = l2grad2(ielem) + sqrt(prop_grad(2,1)*
c     $           prop_grad(2,1)+prop_grad(2,2)*prop_grad(2,2))

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
