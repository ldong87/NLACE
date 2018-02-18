c**********************************************************
      subroutine elem37 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     Written by Sevan Goenezen
c     NOTES:
c     1) inverse element, use only linear tetrahedrals
c     2) 3D, incompressible material, linear elasticity
c     3) Added stabilization

c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem
      parameter (mnode_elem = 4,mquad_elem=30)
      integer ielem, itask, ireg
      integer nGauss, ninte, iinte, inode, jnode,knode, iset
      integer ievab, jevab, jdofn, idofn, ii, jj, i, s 
      double precision xjaco, wtjac, xforc, h, tmp1, tmp2
      double precision tmp3, tmp4, tmp5, mu
      double precision alpha, beta, deno
      double precision shape(4,mnode_elem),wg(mquad_elem)
      double precision sg(mquad_elem),tg(mquad_elem), zg(mquad_elem)
      double precision temp_primal(4*mnode_elem),temp_dual(4*mnode_elem)
      double precision udiff(3), prop_grad(1,4), var_shape_prop(1,4)
      double precision alpha2, temp
      double precision temp_der
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
 1    elemvec_ndofn(1:nnode) = 4
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
 2    read(iread,*) nGauss,ireg,alpha,beta,alpha2,s
      write(iwrit,205) nGauss,ireg,alpha,beta,alpha2,s
 205  format(/,4x,'LINEAR ELASTICITY(2D,2DOF) '/
     +       4x,'gauss pts/dir              ',i10,/
     +       4x,'reg type(0/1/2:none/H1/TVD)',i10,/
     +       4x,'regularization parameter   ',e15.6,/
     +       4x,'extra parameter (TVD)      ',e15.6,/
     +       4x,'stabilization factor    ',e15.6,/
     +       4x,'stabilization terms 0/1 off/
     $           laplacian stabilization  ',i10)

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

c     determine the characteristic length h of the element
c       (for tetrahedra: h is the length of the longest edge)
      tmp1=(xelem(1,1)-xelem(1,2))**2.0d0
      tmp2=(xelem(1,1)-xelem(1,3))**2.0d0
      tmp3=(xelem(1,1)-xelem(1,4))**2.0d0
      tmp4=(xelem(1,2)-xelem(1,3))**2.0d0
      tmp5=(xelem(1,2)-xelem(1,4))**2.0d0
      h=(xelem(1,3)-xelem(1,4))**2.0d0
      do i=2,ndim! 1, ndim
        tmp1 = tmp1 + (xelem(i,1)-xelem(i,2))**2.0d0
        tmp2 = tmp2 + (xelem(i,1)-xelem(i,3))**2.0d0
        tmp3 = tmp3 + (xelem(i,1)-xelem(i,4))**2.0d0
        tmp4 = tmp4 + (xelem(i,2)-xelem(i,3))**2.0d0
        tmp5 = tmp5 + (xelem(i,2)-xelem(i,4))**2.0d0
        h = h + (xelem(i,3)-xelem(i,4))**2.0d0
      enddo
      if (tmp1.gt.h) then
        h=tmp1
      endif
      if (tmp2.gt.h)then
        h=tmp2
      endif
      if (tmp3.gt.h) then
        h=tmp3
      endif
      if (tmp4.gt.h) then
        h=tmp4
      endif
      if (tmp5.gt.h) then
        h=tmp5
      endif
      h=sqrt(h)
      

      call gausstet(nGauss,ninte,sg,tg,zg,wg)

      do iinte = 1, ninte
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shape,
     $              xjaco,.false.,nnode,ndim,elnods,xelem)
        wtjac = wg(iinte)*xjaco
         
        mu  = 0.0d0
        do inode = 1,nnode
          mu = mu + shape(4,inode)*elemdata_nodal(1,inode)
        enddo

        if (s.eq.0) then
          temp  = 0.0d0 
        else if (s.eq.1) then
          temp = 0.5d0*(alpha2*(h**2.0d0))/mu
        else
          write(iwrit,200) 
 200      format(4x,'Stabilization property must be 0 or 1')
          stop
        endif

c     construct the stiffness matrix 	  

          do inode = 1, nnode
	    do jnode = 1, nnode
	    
c            first block is upper left 3x3 matrix 
	      estif((inode*4)-3,(jnode*4)-3) = 
     $         estif((inode*4)-3,(jnode*4)-3) + 
     $          mu*(2.0d0*shape(1,inode)*shape(1,jnode) +
     $              shape(2,inode)*shape(2,jnode) +
     $              shape(3,inode)*shape(3,jnode))*wtjac
     
              estif((inode*4)-3,(jnode*4)-2) = 
     $         estif((inode*4)-3,(jnode*4)-2) +
     $          mu*shape(2,inode)*shape(1,jnode)*wtjac
     
              estif((inode*4)-3,(jnode*4)-1) = 
     $         estif((inode*4)-3,(jnode*4)-1) +
     $          mu*shape(3,inode)*shape(1,jnode)*wtjac
     
              estif((inode*4)-2,(jnode*4)-3) =
     $         estif((inode*4)-2,(jnode*4)-3) +
     $          mu*shape(1,inode)*shape(2,jnode)*wtjac
     
              estif((inode*4)-2,(jnode*4)-2) = 
     $         estif((inode*4)-2,(jnode*4)-2) + 
     $          mu*(2.0d0*shape(2,inode)*shape(2,jnode) +
     $              shape(1,inode)*shape(1,jnode) +
     $              shape(3,inode)*shape(3,jnode))*wtjac
     
              estif((inode*4)-2,(jnode*4)-1) =
     $         estif((inode*4)-2,(jnode*4)-1) +
     $          mu*shape(3,inode)*shape(2,jnode)*wtjac
     
              estif((inode*4)-1,(jnode*4)-3) =
     $         estif((inode*4)-1,(jnode*4)-3) +
     $          mu*shape(1,inode)*shape(3,jnode)*wtjac
     
              estif((inode*4)-1,(jnode*4)-2) =
     $         estif((inode*4)-1,(jnode*4)-2) +
     $          mu*shape(2,inode)*shape(3,jnode)*wtjac
     
              estif((inode*4)-1,(jnode*4)-1) =
     $         estif((inode*4)-1,(jnode*4)-1) +
     $          mu*(2.0d0*shape(3,inode)*shape(3,jnode) +
     $              shape(1,inode)*shape(1,jnode) +
     $              shape(2,inode)*shape(2,jnode))*wtjac
     
c            Contribution from incompressibility constraint

              estif(inode*4,(jnode*4)-3) = 
     $         estif(inode*4,(jnode*4)-3) +
     $          shape(4,inode)*shape(1,jnode)*wtjac
     
              estif(inode*4,(jnode*4)-2) =
     $         estif(inode*4,(jnode*4)-2) +
     $          shape(4,inode)*shape(2,jnode)*wtjac
     
              estif(inode*4,(jnode*4)-1) =
     $         estif(inode*4,(jnode*4)-1) +
     $          shape(4,inode)*shape(3,jnode)*wtjac
     
c            Contribution from pressure term

              estif((inode*4)-3,jnode*4) = 
     $         estif((inode*4)-3,jnode*4) -
     $          shape(1,inode)*shape(4,jnode)*wtjac
     
              estif((inode*4)-2,jnode*4) =
     $         estif((inode*4)-2,jnode*4) -
     $          shape(2,inode)*shape(4,jnode)*wtjac
     
              estif((inode*4)-1,jnode*4) =
     $         estif((inode*4)-1,jnode*4) -
     $          shape(3,inode)*shape(4,jnode)*wtjac
     
c           Stabilization contribution     
     
              estif(inode*4,jnode*4) =
     $         estif(inode*4,jnode*4) + 
     $          temp*(shape(1,inode)*shape(1,jnode) +
     $               shape(2,inode)*shape(2,jnode) +
     $               shape(3,inode)*shape(3,jnode))*wtjac
     
            end do    ! jnode
          end do       ! inode

      enddo! iinte
      
c     compute element internal force vector
      do ievab =  1,nnode*4
        jevab = 0
        do jnode = 1,nnode
          do jdofn = 1,4! 4 = elemvec_ndofn(jnode)
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

      call gausstet(nGauss,ninte,sg,tg,zg,wg)

      do iinte = 1, ninte
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shape,xjaco,
     $              .false.,nnode,ndim,elnods,xelem)
        wtjac = wg(iinte)*xjaco

        ievab = 0
        do inode = 1, nnode
          do idofn = 1,4! 4 = elemvec_ndofn(inode)
            ievab = ievab+1
            do jnode = 1, nnode
              do jdofn = 1,4! 4 = elemvec_ndofn(jnode)
                if (idofn.eq.jdofn) then
                  eforc(ievab) = eforc(ievab)-
     $                       shape(4,inode)*shape(4,jnode)
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
        do idofn = 1,4! 4 = elemvec_ndofn(inode)
          ii = ii+1
          temp_primal(ii) = uelem(idofn,inode)
          temp_dual(ii)   = uelem_dual(idofn,inode)
        enddo
      enddo

c     determine the characteristic length h of the element
c       (for tetrahedra: h is the length of the longest edge)
      tmp1=(xelem(1,1)-xelem(1,2))**2.0d0
      tmp2=(xelem(1,1)-xelem(1,3))**2.0d0
      tmp3=(xelem(1,1)-xelem(1,4))**2.0d0
      tmp4=(xelem(1,2)-xelem(1,3))**2.0d0
      tmp5=(xelem(1,2)-xelem(1,4))**2.0d0
      h=(xelem(1,3)-xelem(1,4))**2.0d0
      do i=2,ndim! 1, ndim
        tmp1 = tmp1 + (xelem(i,1)-xelem(i,2))**2.0d0
        tmp2 = tmp2 + (xelem(i,1)-xelem(i,3))**2.0d0
        tmp3 = tmp3 + (xelem(i,1)-xelem(i,4))**2.0d0
        tmp4 = tmp4 + (xelem(i,2)-xelem(i,3))**2.0d0
        tmp5 = tmp5 + (xelem(i,2)-xelem(i,4))**2.0d0
        h = h + (xelem(i,3)-xelem(i,4))**2.0d0
      enddo
      if (tmp1.gt.h) then
        h=tmp1
      endif
      if (tmp2.gt.h)then
        h=tmp2
      endif
      if (tmp3.gt.h) then
        h=tmp3
      endif
      if (tmp4.gt.h) then
        h=tmp4
      endif
      if (tmp5.gt.h) then
        h=tmp5
      endif
      h=sqrt(h)

c     compute the gradient and the objective function value

      call gausstet(nGauss,ninte,sg,tg,zg,wg)
c
      do iinte = 1, ninte
        call shape3(ielem,sg(iinte),tg(iinte),zg(iinte),shape,xjaco,
     $              .false.,nnode,ndim,elnods,xelem)
        wtjac = wg(iinte)*xjaco
 
c       create prop_grad (contains prop gradient and value) 
        prop_grad(1,1:4)  = 0.0d0
        do idofn = 1,4! three gradents+the value
          do inode = 1,nnode
            prop_grad(1,idofn) = prop_grad(1,idofn)+
     $                 shape(idofn,inode)*
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
            do idofn = 1,4
              var_shape_prop(1,idofn) = shape(idofn,knode)
            enddo
          endif
c         take out pf knode loop the following! JFD ask sevan
          mu  = 0.0d0
          do inode = 1,nnode
            mu     = mu     + shape(4,inode)*elemdata_nodal(1,inode)
          enddo

          if (s.eq.0) then
            temp  = 0.0d0
            temp_der = 0.0d0
          elseif (s.eq.1) then
            temp = (alpha2*(h**2.0d0))/mu
            temp_der = -var_shape_prop(1,4)*(alpha2*(h**2.0d0))/(mu**2)
          else
            write(iwrit,300) 
 300        format(4x,'Stabilization property must be 0, 1, or 2')
            stop
          endif

          do inode = 1, nnode
	    do jnode = 1, nnode
	    
c            first block is upper left 3x3 matrix 
	      estif((inode*4)-3,(jnode*4)-3) = 
     $         estif((inode*4)-3,(jnode*4)-3) + 
     $               var_shape_prop(1,4)*
     $             (2.0d0*shape(1,inode)*shape(1,jnode) +
     $              shape(2,inode)*shape(2,jnode) +
     $              shape(3,inode)*shape(3,jnode))*wtjac
     
              estif((inode*4)-3,(jnode*4)-2) = 
     $         estif((inode*4)-3,(jnode*4)-2) +
     $              var_shape_prop(1,4)*
     $             shape(2,inode)*shape(1,jnode)*wtjac
     
              estif((inode*4)-3,(jnode*4)-1) = 
     $         estif((inode*4)-3,(jnode*4)-1) +
     $               var_shape_prop(1,4)*
     $             shape(3,inode)*shape(1,jnode)*wtjac
     
              estif((inode*4)-2,(jnode*4)-3) =
     $         estif((inode*4)-2,(jnode*4)-3) +
     $              var_shape_prop(1,4)*
     $             shape(1,inode)*shape(2,jnode)*wtjac
     
              estif((inode*4)-2,(jnode*4)-2) = 
     $         estif((inode*4)-2,(jnode*4)-2) + 
     $              var_shape_prop(1,4)*
     $             (2.0d0*shape(2,inode)*shape(2,jnode) +
     $              shape(1,inode)*shape(1,jnode) +
     $              shape(3,inode)*shape(3,jnode))*wtjac
     
              estif((inode*4)-2,(jnode*4)-1) =
     $         estif((inode*4)-2,(jnode*4)-1) +
     $             var_shape_prop(1,4)*
     $             shape(3,inode)*shape(2,jnode)*wtjac
     
              estif((inode*4)-1,(jnode*4)-3) =
     $         estif((inode*4)-1,(jnode*4)-3) +
     $             var_shape_prop(1,4)*
     $             shape(1,inode)*shape(3,jnode)*wtjac
     
              estif((inode*4)-1,(jnode*4)-2) =
     $         estif((inode*4)-1,(jnode*4)-2) +
     $              var_shape_prop(1,4)*
     $             shape(2,inode)*shape(3,jnode)*wtjac
     
              estif((inode*4)-1,(jnode*4)-1) =
     $         estif((inode*4)-1,(jnode*4)-1) +
     $              var_shape_prop(1,4)*
     $             (2.0d0*shape(3,inode)*shape(3,jnode) +
     $              shape(1,inode)*shape(1,jnode) +
     $              shape(2,inode)*shape(2,jnode))*wtjac
     
c           Stabilization contribution     
     
              estif(inode*4,jnode*4) =
     $         estif(inode*4,jnode*4) + 
     $          temp_der*(shape(1,inode)*shape(1,jnode) +
     $               shape(2,inode)*shape(2,jnode) +
     $               shape(3,inode)*shape(3,jnode))*wtjac
     
            end do    ! jnode
          end do       ! inode


          do ii = 1, nnode*4!nevab! nevab was a global variable
            do jj = 1, nnode*4!nevab! nevab was a global variable
              egrad(1,knode) = egrad(1,knode)+
     $                 temp_dual(ii)*estif(ii,jj)*temp_primal(jj)
            enddo
          enddo
            
c         account for the regularization term
          if (ireg.eq.1) then
            egrad(1,knode) = egrad(1,knode) + alpha*(
     $              var_shape_prop(1,1)*prop_grad(1,1)+
     $              var_shape_prop(1,2)*prop_grad(1,2)+
     $              var_shape_prop(1,3)*prop_grad(1,3) )*
     $              wtjac
          elseif (ireg.eq.2) then
            deno = dsqrt(
     $              beta*beta+
     $              prop_grad(1,1)*prop_grad(1,1)+
     $              prop_grad(1,2)*prop_grad(1,2)+
     $              prop_grad(1,3)*prop_grad(1,3))
            egrad(1,knode) = egrad(1,knode) + 
     $              0.5d0*alpha*(
     $              var_shape_prop(1,1)*prop_grad(1,1)+
     $              var_shape_prop(1,2)*prop_grad(1,2)+
     $              var_shape_prop(1,3)*prop_grad(1,3) )*
     $              wtjac/deno
          endif
          udiff(1) = udiff(1)+ shape(4,knode)*uelem_diff(1,knode)
          udiff(2) = udiff(2)+ shape(4,knode)*uelem_diff(2,knode)
	  udiff(3) = udiff(3)+ shape(4,knode)*uelem_diff(3,knode)
        enddo! knode

        elemDataMatch(ielem) = elemDataMatch(ielem) + 0.5d0*wtjac*
     $        (udiff(1)*udiff(1)+udiff(2)*udiff(2)+udiff(3)*udiff(3))

        if (ireg.eq.1) then
          elemRegul(ielem) = elemRegul(ielem)+ 
     $           0.5d0*alpha*(
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(1,3)*prop_grad(1,3) )*
     $           wtjac
        elseif (ireg.eq.2) then
          elemRegul(ielem) = elemRegul(ielem)+ 
     $           0.5d0*alpha*dsqrt(beta*beta+
     $           prop_grad(1,1)*prop_grad(1,1)+
     $           prop_grad(1,2)*prop_grad(1,2)+
     $           prop_grad(1,3)*prop_grad(1,3))*
     $           wtjac
        endif
        l2grad1(ielem) = l2grad1(ielem) + sqrt(prop_grad(1,1)*
     $           prop_grad(1,1)+prop_grad(1,2)*prop_grad(1,2)
     $           +prop_grad(1,3)*prop_grad(1,3))
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
