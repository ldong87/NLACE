c*************************************************************
      subroutine elem611 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     Written by Tom Seidl     
c     traction / spring boundary element (2D/3D)
c     Notes: boundary_springs.pdf  
c     ** Elemental Tractions not yet present in 3D 
c     ** 2D uses 1 deformation, 3D can use 1 or 2 

c*************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem
      parameter (mnode_elem = 3,mquad_elem=16)
      integer ielem, itask, iset 
      integer nGauss, ninte, iinte, inode, jnode, deg_p 
      integer ievab, jevab, jdofn, idofn, rdofn, i, j  
      double precision xjaco, wtjac
      double precision xk11, xk12, xk13, xk21, xk22 ! spring tensor entries
      double precision xk23, xk31, xk32, xk33
      double precision xk11b, xk12b, xk13b, xk21b, xk22b ! b spring tensor entries
      double precision xk23b, xk31b, xk32b, xk33b
      double precision xkelem(ndime,ndime) ! elemental spring tensor
      double precision shap(ndime,nnode)
      double precision sg(mquad_elem),wg(mquad_elem),tg(mquad_elem)
      double precision t1,t2 ! elemental traction data
      integer nnode
      double precision pelem(*)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,*)
      double precision elemdata_nodal(nset_nodal,*)
      double precision uelem(mdofn,*)
      double precision felem(mdofn) ! elemental traction vector
      double precision uelem_meas(mdofn,*)
      logical deriv
c      integer TID,OMP_GET_THREAD_NUM! uncomment for openmp
c-------------------------------------------------------------

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask
      
c     --------------------------------------------------------
c     Initialize the elemental parameters
c     --------------------------------------------------------
      
 1     if (ndime.eq.2) then 
        elemvec_ndofn(1:nnode) = 2
       else
        elemvec_ndofn(1:nnode) = 4
       endif
c      buildSymmetricMatrices = .true.
      return


c     --------------------------------------------------------
c     Read and store the material properties
c     --------------------------------------------------------
c     nGaussPtsPerDir xk11 xk12 xk21 xk22
c     nGauss     : no. of gauss pts/direction
c     xk11        : spring tensor component
c     xk12        : spring tensor component
c     xk21        : spring tensor component
c     xk22        : spring tensor component
c     t1          : surface normal x component
c     t2          : surface normal y component

c     -Add 2 deformations capability to 2D -TS

     
 2    if (ndime.eq.2) then
        read(iread,*) nGauss,xk11,xk12,xk21,xk22,t1,t2
        write(iwrit,205) nGauss,xk11,xk12,xk21,xk22,t1,t2

 205  format(/' boundary springs (2D, 2DOFs)'/
     +       ' gauss pts/dir .........................',i12,/
     +       ' K_11....................................',1p,e16.4,/
     +       ' K_12....................................',1p,e16.4,/
     +       ' K_21....................................',1p,e16.4,/
     +       ' K_22....................................',1p,e16.4,/
     +       ' traction_1..............................',1p,e16.4,/
     +       ' traction_2..............................',1p,e16.4)

        pelem(1) = dble(nGauss)
        pelem(2) = xk11
        pelem(3) = xk12
        pelem(4) = xk21
        pelem(5) = xk22
        pelem(6) = t1
        pelem(7) = t2
      
      else ! 2D / 3D 

       if (nmeas.eq.1) then 
        read(iread,*) deg_p,xk11,xk12,xk13,xk21,xk22,xk23,xk31,xk32,xk33
        write(iwrit,206) deg_p,xk11,xk12,xk13,xk21,xk22,xk23,xk31,xk32,
     +xk33
206   format(/' boundary springs (3D, 4DOFs)'/
     +       ' degree of precision..........................',i12,/
     +       ' K_11....................................',1p,e16.4,/
     +       ' K_12....................................',1p,e16.4,/
     +       ' K_13....................................',1p,e16.4,/
     +       ' K_21....................................',1p,e16.4,/
     +       ' K_22....................................',1p,e16.4,/
     +       ' K_23....................................',1p,e16.4,/
     +       ' K_31....................................',1p,e16.4,/
     +       ' K_32....................................',1p,e16.4,/
     +       ' K_33....................................',1p,e16.4)

        pelem(1) = dble(deg_p)
        pelem(2) = xk11
        pelem(3) = xk12
        pelem(4) = xk13
        pelem(5) = xk21
        pelem(6) = xk22
        pelem(7) = xk23
        pelem(8) = xk31
        pelem(9) = xk32
        pelem(10) = xk33
       
      else  ! 2 deformations

       read(iread,*) deg_p,xk11,xk12,xk13,xk21,xk22,xk23,xk31,xk32,xk33,
     + xk11b,xk12b,xk13b,xk21b,xk22b,xk23b,xk31b,xk32b,xk33b 
        write(iwrit,207) deg_p,xk11,xk12,xk13,xk21,xk22,xk23,xk31,xk32,
     +xk33,xk11b,xk12b,xk13b,xk21b,xk22b,xk23b,xk31b,xk32b,xk33b
207   format(/' boundary springs (3D,finite elasticity,3DOFs)'/
     +       ' degree of precision..........................',i12,/
     +       ' K_11....................................',1p,e16.4,/
     +       ' K_12....................................',1p,e16.4,/
     +       ' K_13....................................',1p,e16.4,/
     +       ' K_21....................................',1p,e16.4,/
     +       ' K_22....................................',1p,e16.4,/
     +       ' K_23....................................',1p,e16.4,/
     +       ' K_31....................................',1p,e16.4,/
     +       ' K_32....................................',1p,e16.4,/
     +       ' K_33....................................',1p,e16.4,/
     +       ' K_11b...................................',1p,e16.4,/
     +       ' K_12b...................................',1p,e16.4,/
     +       ' K_13b...................................',1p,e16.4,/
     +       ' K_21b...................................',1p,e16.4,/
     +       ' K_22b...................................',1p,e16.4,/
     +       ' K_23b...................................',1p,e16.4,/
     +       ' K_31b...................................',1p,e16.4,/
     +       ' K_32b...................................',1p,e16.4,/
     +       ' K_33b...................................',1p,e16.4)

        pelem(1) = dble(deg_p)
        pelem(2) = xk11
        pelem(3) = xk12
        pelem(4) = xk13
        pelem(5) = xk21
        pelem(6) = xk22
        pelem(7) = xk23
        pelem(8) = xk31
        pelem(9) = xk32
        pelem(10) = xk33
        pelem(11) = xk11b
        pelem(12) = xk12b
        pelem(13) = xk13b
        pelem(14) = xk21b
        pelem(15) = xk22b
        pelem(16) = xk23b
        pelem(17) = xk31b
        pelem(18) = xk32b
        pelem(19) = xk33b

       endif ! nmeas

      endif ! 2D/3D

      return
      

c     --------------------------------------------------------
c     Build the elemental consistent tangent stiffness matrix (estif)
c       and the elemental RHS (eforc)
c     --------------------------------------------------------
      
 3    if (ndime.eq.2) then 
       nGauss = int(pelem(1))
       xkelem(1,1) = pelem(2)
       xkelem(1,2) = pelem(3)
       xkelem(2,1) = pelem(4)
       xkelem(2,2) = pelem(5)
       t1          = pelem(6)
       t2          = pelem(7)
c     initialize some variables
       estif(:,:) = 0.0d0! elemental tangent stiffness
       eforc(:) = 0.0d0! elemental RHS/residual 
       felem(1) = t1
       felem(2) = t2

       call gauss1(nGauss,sg,wg)

       do iinte = 1,nGauss! for all Gauss integration points
         call shape1(sg(iinte),shap,xjaco,.false.,nnode,ndime,xelem)
         wtjac = wg(iinte)*xjaco

               ! edit the parts below to create estif and eforc   -TS

c       build the tangent stiffness matrix: estif is symmetric
         ievab = 0
         do inode = 1, nnode
           do idofn = 1,2 ! 2 = elemvec_ndofn(inode)
             ievab = ievab+1
             jevab = 0! jevab=0 when not using the symmetries
             do jnode = 1, 2
                do jdofn = 1,2 ! 2 = elemvec_ndofn(inode)
                 jevab = jevab+1
                     estif(ievab,jevab)= estif(ievab,jevab) +
     +               xkelem(idofn,jdofn)*shap(2,inode)
     +               *shap(2,jnode)*wtjac ! spring stiffness
                end do! jdofn
             end do! jnode
             eforc(ievab) = eforc(ievab) - 
     +       felem(idofn)*shap(2,inode)*wtjac  ! traction forc
           end do! idofn
         end do! inode
      end do! iinte

c       eforc  (springs) loop     
         ievab = 0
         do inode = 1, 2 ! 2 nodes
           do idofn = 1,2 ! 2 dofs
             ievab =  ievab + 1
             jevab = 0
             do jnode = 1, 2 ! 2 nodes
                do jdofn = 1, 2 ! 2 dofs
                    if (nset_nodal.eq.4) then
                       rdofn = jdofn + 2 + 2*( ifield - 1) 
                    else
                       rdofn = jdofn + 1 + 2*(ifield - 1)
                    endif
                    jevab = jevab + 1

          eforc(ievab) = eforc(ievab) - estif(ievab,jevab)
     +    *(bc_step*elemdata_nodal(rdofn,jnode) - uelem(jdofn,jnode))

                end do! jdofn
             end do! jnode
           end do! idofn
         end do! inode
    
      
      else  ! ndim = 3; triangle surface element. 

       deg_p = int(pelem(1))

        if (ifield.eq.1)  then
          xkelem(1,1) = pelem(2)
          xkelem(1,2) = pelem(3)
          xkelem(1,3) = pelem(4)
          xkelem(2,1) = pelem(5)
          xkelem(2,2) = pelem(6)
          xkelem(2,3) = pelem(7)
          xkelem(3,1) = pelem(8)
          xkelem(3,2) = pelem(9)
          xkelem(3,3) = pelem(10)
        else
          xkelem(1,1) = pelem(11)
          xkelem(1,2) = pelem(12)
          xkelem(1,3) = pelem(13)
          xkelem(2,1) = pelem(14)
          xkelem(2,2) = pelem(15)
          xkelem(2,3) = pelem(16)
          xkelem(3,1) = pelem(17)
          xkelem(3,2) = pelem(18)
          xkelem(3,3) = pelem(19)
        endif 

c     initialize some variables
       estif(:,:) = 0.0d0! elemental tangent stiffness
       eforc(:) = 0.0d0! elemental RHS/residual 
       call gausstri(deg_p,ninte,sg,tg,wg)
       deriv = .true.
       do iinte = 1,ninte ! for all Gauss integration points
         call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.true.
     +         ,nnode,elnods,ndime,xelem)
         wtjac = wg(iinte)*xjaco
         ievab = 0
         do inode = 1, 3 ! 3 nodes
           do idofn = 1,3 ! 4 dofs (skip pressure)
             ievab = (inode-1)*4 + idofn
             jevab = 0 
             do jnode = 1, 3 ! 3 nodes
                do jdofn = 1, 3 ! 4 dofs (skip pressure)
                     jevab = (jnode-1)*4 + jdofn 
                     estif(ievab,jevab)= estif(ievab,jevab) +
     +               xkelem(idofn,jdofn)*shap(3,inode)
     +               *shap(3,jnode)*wtjac

                end do! jdofn
             end do! jnode
           end do! idofn
         end do! inode


       end do! iinte

c       eforc loop     
         ievab = 0
         do inode = 1, 3 ! 3 nodes
           do idofn = 1,3 ! 4 dofs (skip pressure)
             ievab = (inode-1)*4 + idofn
             jevab = 0
             do jnode = 1, 3 ! 3 nodes
                do jdofn = 1, 3 ! 4 dofs (skip pressure)
                    if (nset_nodal.eq.5) then
                       rdofn = jdofn + 2 + 3 * (ifield - 1)  !for nodal data fields
                    else
                       rdofn = jdofn + 1 + 3 * (ifield - 1)
                    endif
                    jevab = (jnode-1)*4 + jdofn 

          eforc(ievab) = eforc(ievab) - estif(ievab,jevab)
     +    *(bc_step*elemdata_nodal(rdofn,jnode) - uelem(jdofn,jnode))

                end do! jdofn
             end do! jnode
           end do! idofn
         end do! inode
    
     
      endif ! 2D/3D if 
      
      return
cc     --------------------------------------------------------
cc     Evaluate the elemental RHS/residual (eforc) from a source term
cc     --------------------------------------------------------
 5    continue
      eforc(:) = 0.0d0! elemental RHS/residual
cc     turned off for now
      return
c     -------------------------------------------------------
c     Evaluate the elemental RHS (eforc) for the dual problem
c     -------------------------------------------------------
 6    continue ! no contribution from the boundary springs
      return    
 7    return
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

