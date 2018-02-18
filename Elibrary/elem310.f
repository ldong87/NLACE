c**********************************************************************
      subroutine elem310 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c*********************************************************************
c treat traction as an optimization variable 
c tetrahedral element where interior node defines the triangle
c Assad Oberai and Li Dong
c*********************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer l
      integer ielem, itask, nnode, nno, ndim, iinte, ninte
      integer ii, inode, jnode, inter, iset, knode,linter
      integer ievab, idofn, jevab,jdofn,mnode_elem, mquad_elem
      integer ireg, inode2d
      parameter (mnode_elem = 4,mquad_elem=20)
      integer elnods(mnode) 
      integer elnods2d(3), node_num_2d(4)
      double precision mag_pos1, mag_pos2, pos2dotuni1
      double precision pelem(*), eforc(mevab), xelem(ndime,*)
      double precision pos1(3), pos2(3), uni1(3)
      double precision shap(4,mnode_elem),shap2d(3,3),xel(2,3)
      double precision elemdata_nodal(nset_nodal,*), uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
      double precision tmpCoor(3,3)
      double precision sg(mquad_elem),tg(mquad_elem), wg(mquad_elem)
      double precision estif(mevab,mevab)
      double precision xjaco,wtjac
      double precision deno, alpha, beta
      double precision trac(3,3),dual_loc(3),var_shape_prop(3,3)
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
c     alpha : reg. param
c     ireg: flag for regularizatio(1-L2)
c           flag for regularizatio(2-H1)
c           flag for regularizatio(3-TVD)
c     beta : extra param. for TVD
 2    read(iread,*) l, inter, alpha, ireg, beta
      write(iwrit,205) l, inter, alpha, ireg, beta
 205  format(/' Tractions as opt. vars.(3D,4DOF) '/
     +       ' gauss pts/dir .........................',i12,/
     +       ' interior node .........................',i12,/
     +       ' reg. param.  .......................',1p,e16.4,/
     +       ' reg. type...........................',i12,/
     +       ' extra param. for TVD..................',1p,e16.4,/)
      pelem(1) = dble(l)
      pelem(2) = dble(inter)
      pelem(3) = alpha 
      pelem(4) = dble(ireg)
      pelem(5) = beta
      return
      
c     -----------------------------------------------------------------
c     Build the elemental lhs and rhs arising from pressure term 
c     -----------------------------------------------------------------
 3    l     = int(pelem(1))
      inter = int(pelem(2))
      alpha = pelem(3)
      ireg = int(pelem(4))
      beta = pelem(5)
      estif(:,:)=0.0d0
      eforc(:) = 0.0d0
      elnods2d(:)=0 
      node_num_2d(:)=0 
      tmpCoor(:,:) = 0.0d0
      pos1(:)=0.0d0
      uni1(:)=0.0d0
      pos2(:)=0.0d0
      xel(:,:) =0.0d0
c----------------------------------------------------------     
c associating flag to each node of the tet , flag is 0 , if it is an interior node
c and 1,2,3 value if it is a triangular face node
      
      inode2d=0
      do inode =1,nnode
       if (elnods(inode).eq.inter) then
        node_num_2d(inode)=0
                           ! interior node is flagged 0
       else
        inode2d=inode2d+1
        node_num_2d(inode)=inode2d
        elnods2d(inode2d)=elnods(inode) 
        tmpCoor(1:ndime,inode2d)=xelem(1:ndime,inode)
       endif
      enddo

       nno = 3 
               ! reduced number of nodes for triangular elements
       ndim = 2
               ! reduce dimension for triangular elements

c Compute coordinates in a surface coordinate system
       pos1(1:3)=tmpCoor(1:3,2)-tmpCoor(1:3,1)
       pos2(1:3)=tmpCoor(1:3,3)-tmpCoor(1:3,1)
       mag_pos1=dsqrt(pos1(1)**2.0d0+pos1(2)**2.0d0+pos1(3)**2.0d0)
       mag_pos2=dsqrt(pos2(1)**2.0d0+pos2(2)**2.0d0+pos2(3)**2.0d0)
       uni1(:) = pos1(:)/mag_pos1
       pos2dotuni1 = pos2(1)*uni1(1)+pos2(2)*uni1(2)+pos2(3)*uni1(3)
       pos2(1:3) = pos2(1:3) - pos2dotuni1*uni1(1:3)
       mag_pos2=dsqrt(pos2(1)**2.0d0+pos2(2)**2.0d0+pos2(3)**2.0d0)
       xel(1,1) = 0.0d0
       xel(2,1) = 0.0d0
       xel(1,2) = mag_pos1
       xel(2,2) = 0.0d0
       xel(1,3) = pos2dotuni1
       xel(2,3) = mag_pos2 
c done calculating coordinates

      call gausstri(l,ninte,sg,tg,wg)

       do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap2d,xjaco,
     $              .false.,nno,elnods2d,ndim,xel)

        trac(1:3,1:3) = 0.0d0
        do inode = 1,nnode
           inode2d = node_num_2d(inode)
           if (inode2d.ne.0) then
              do ii = 1,3
                 trac(ii,3) = trac(ii,3)+
     $            shap2d(3,inode2d)*elemdata_nodal(2+ii,inode)
	      enddo
	    endif
	enddo

        wtjac = wg(iinte)*xjaco
c construct element force vector( rhs ) for primal problem
      ievab =0

       do inode =1,nnode
         do idofn = 1,4
            ievab = ievab +1
            inode2d = node_num_2d(inode)
               if (inode2d.ne.0.and.idofn.ne.4) then
                  eforc(ievab) = eforc(ievab)
     $            -shap2d(3,inode2d)*trac(idofn,3)
     $            *wtjac
               endif
            enddo
           ! idofn
         enddo
         ! inod
      enddo
      ! iinte
      return
 4    return
 5    continue
      eforc(:) = 0.0d0
      return
 6    continue
      eforc(:) = 0.0d0
      return
c     -----------------------------------------------------------------
c     Build the elemental gradient matrix
c     -----------------------------------------------------------------
 7    l     = int(pelem(1))
      inter = int(pelem(2))
      alpha = pelem(3)
      ireg = int(pelem(4))
      beta = pelem(5)

      egrad(:,:) = 0.0d0
      elemDataMatch(ielem) = 0.0d0
      elemRegul(ielem) = 0.0d0
      l2grad1(ielem) = 0.0d0
      l2grad2(ielem) = 0.0d0

      elnods2d(:)=0 
      node_num_2d(:)=0 
      tmpCoor(:,:) = 0.0d0
      pos1(:)=0.0d0
      uni1(:)=0.0d0
      pos2(:)=0.0d0
      xel(:,:) =0.0d0
c----------------------------------------------------------     
c associating flag to each node of the tet , flag is 0 , if it is an interior node
c and 1,2,3 value if it is a triangular face node
      
      inode2d=0
      do inode =1,nnode
       if (elnods(inode).eq.inter) then
        node_num_2d(inode)=0
                           ! interior node is flagged 0
       else
        inode2d=inode2d+1
        node_num_2d(inode)=inode2d
        elnods2d(inode2d)=elnods(inode) 
        tmpCoor(1:ndime,inode2d)=xelem(1:ndime,inode)
       endif
      enddo

       nno = 3 
               ! reduced number of nodes for triangular elements
       ndim = 2
               ! reduce dimension for triangular elements

c Compute coordinates in a surface coordinate system
       pos1(1:3)=tmpCoor(1:3,2)-tmpCoor(1:3,1)
       pos2(1:3)=tmpCoor(1:3,3)-tmpCoor(1:3,1)
       mag_pos1=dsqrt(pos1(1)**2.0d0+pos1(2)**2.0d0+pos1(3)**2.0d0)
       mag_pos2=dsqrt(pos2(1)**2.0d0+pos2(2)**2.0d0+pos2(3)**2.0d0)
       uni1(:) = pos1(:)/mag_pos1
       pos2dotuni1 = pos2(1)*uni1(1)+pos2(2)*uni1(2)+pos2(3)*uni1(3)
       pos2(1:3) = pos2(1:3) - pos2dotuni1*uni1(1:3)
       mag_pos2=dsqrt(pos2(1)**2.0d0+pos2(2)**2.0d0+pos2(3)**2.0d0)
       xel(1,1) = 0.0d0
       xel(2,1) = 0.0d0
       xel(1,2) = mag_pos1
       xel(2,2) = 0.0d0
       xel(1,3) = pos2dotuni1
       xel(2,3) = mag_pos2 
c done calculating coordinates
c--------------------------------------------

      call gausstri(l,ninte,sg,tg,wg)

       do iinte = 1, ninte
        call shape2(ielem,sg(iinte),tg(iinte),shap2d,xjaco,
     $              .false.,nno,elnods2d,ndim,xel)

        trac(1:3,1:3) = 0.0d0
        dual_loc(1:3) = 0.0d0
        do inode = 1,nnode
          inode2d = node_num_2d(inode)
          if (inode2d.ne.0) then
            do iset = 1,3
             do idofn = 1,3
              trac(iset,idofn) = trac(iset,idofn)+
     $            shap2d(idofn,inode2d)*elemdata_nodal(2+iset,inode)
             enddo
            enddo
            do idofn = 1,3
             dual_loc(idofn) = dual_loc(idofn)+
     $            shap2d(3,inode2d)*uelem_dual(idofn,inode)
	    enddo
	  endif
	enddo

        wtjac = wg(iinte)*xjaco
c construct element gradient

       do inode =1,nnode
         do iset = 1,3
            inode2d = node_num_2d(inode)
            if (inode2d.ne.0) then

            var_shape_prop(:,:) = 0.0d0
            if (ielem_nod_grad(iset+2,inode).eq.1) then
              do idofn = 1,3
                var_shape_prop(iset,idofn) = shap2d(idofn,inode2d)
              enddo
            endif

            egrad(2+iset,inode) = egrad(2+iset,inode)
     $              -var_shape_prop(iset,3)*dual_loc(iset)
     $              *wtjac

             if (ireg.eq.1) then ! l2-regularization
               egrad(2+iset,inode) = egrad(2+iset,inode)
     $              +var_shape_prop(iset,3)*trac(iset,3)
     $              *wtjac*alpha
              else if (ireg.eq.2) then ! h1-regularization
               egrad(2+iset,inode) = egrad(2+iset,inode)
     $              +(var_shape_prop(iset,1)*trac(iset,1)
     $              +var_shape_prop(iset,2)*trac(iset,2)
     $              )*wtjac*alpha
              else if ((ireg.eq.3).or.(ireg.eq.31)) then! TVD reg.
c                deno = sqrt(beta*beta+trac(iset,1)*trac(iset,1)+
c     $                                trac(iset,2)*trac(iset,2))
                deno = sqrt(beta*beta+trac(1,1)*trac(1,1)+
     $                                trac(1,2)*trac(1,2)+
     $                                trac(2,1)*trac(2,1)+
     $                                trac(2,2)*trac(2,2)+
     $                                trac(3,1)*trac(3,1)+
     $                                trac(3,2)*trac(3,2))
                egrad(2+iset,inode)=egrad(2+iset,inode) + wtjac*alpha*
     $              (var_shape_prop(iset,1)*trac(iset,1)+
     $               var_shape_prop(iset,2)*trac(iset,2))/deno
              endif !ireg
            endif
           ! (inode_2d ...)
         enddo
           ! iset
       enddo
         ! inod
c	evaluate the regularization contribution to the
c	objective function
       if (ireg.eq.1) then   !L2 
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      alpha*(
     $		trac(1,3)*trac(1,3)+trac(2,3)*trac(2,3)+
     $		trac(3,3)*trac(3,3))
       else if (ireg.eq.2) then! H1 
          elemRegul(ielem) = elemRegul(ielem) + wtjac*0.5d0*
     $      alpha*(
     $		trac(1,1)*trac(1,1)+trac(1,2)*trac(1,2)+
     $		trac(2,1)*trac(2,1)+trac(2,2)*trac(2,2)+
     $      trac(3,1)*trac(3,1)+trac(3,2)*trac(3,2))
        elseif (ireg.eq.3) then! TVD reg. with offset
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      alpha*sqrt(beta*beta+
     $		trac(1,1)*trac(1,1)+trac(1,2)*trac(1,2)+
     $		trac(2,1)*trac(2,1)+trac(2,2)*trac(2,2)+
     $		trac(3,1)*trac(3,1)+trac(3,2)*trac(3,2))
        elseif (ireg.eq.31) then! TVD reg. without offset
          elemRegul(ielem) = elemRegul(ielem) + wtjac*
     $      alpha*(sqrt(beta*beta+
     $		trac(1,1)*trac(1,1)+trac(1,2)*trac(1,2)+
     $		trac(2,1)*trac(2,1)+trac(2,2)*trac(2,2)+
     $		trac(3,1)*trac(3,1)+trac(3,2)*trac(3,2))-beta)
       endif
      enddo
      ! iinte
      regularization = regularization + elemRegul(ielem)
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
