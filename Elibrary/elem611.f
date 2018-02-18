c*************************************************************
      subroutine elem611 (ielem,itask,pelem,nnode,estif,eforc,elnods, 
     $    xelem,elemdata_nodal,uelem,uelem_meas,TT_nodal)
c     spring boundary element (2D/3D)
c*************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer  mnode_elem, mquad_elem
      parameter (mnode_elem = 3,mquad_elem=16)
      integer ielem, itask, iset, Tflag
      integer nGauss, ninte, iinte, inode, jnode, deg_p
      integer ievab, jevab, jdofn, idofn, i, j
      double precision xjaco, wtjac, kmag ! kmag -> magnitude of K tensor
      double precision xk11, xk12, xk13, xk21, xk22 ! spring tensor entries
      double precision xk23, xk31, xk32, xk33
      double precision xkelem(mdofn,mdofn) ! elemental spring tensor
      double precision shap(ndime,nnode)
      double precision sg(mquad_elem),wg(mquad_elem),tg(mquad_elem)
c      !real(kind=8) :: t1,t2 ! elemental traction data
      integer nnode
      double precision pelem(*)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,*)
      double precision elemdata_nodal(nset_nodal,*)
      double precision uelem(mdofn,*)
c     !real(kind=8) :: felem(mdofn) ! elemental traction vector
      double precision uelem_meas(mdofn,*)
      double precision TT_nodal(mnode,mdofn,mdofn)
      logical deriv
c      integer TID,OMP_GET_THREAD_NUM! uncomment for openmp
c-------------------------------------------------------------

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask
          
c     --------------------------------------------------------
c     Initialize the elemental parameters
c     --------------------------------------------------------
          
 1     elemvec_ndofn(1:nnode) = mdofn
c        buildSymmetricMatrices = .true.
      return


c     --------------------------------------------------------
c     Read and store the material properties
c     --------------------------------------------------------
c     nGaussPtsPerDir xk11 xk12 xk21 xk22
c     nGauss      : no. of gauss pts/direction
c     kmag        : magnitude of the spring tensor
c     xk11        : spring tensor component
c     xk12        : spring tensor component
c     xk21        : spring tensor component
c     xk22        : spring tensor component
c     Tflag       : Use T from objective function

         
 2     if (ndime.eq.2) then   
          read(iread,*) nGauss,kmag,xk11,xk12,xk22,Tflag
          write(iwrit,205) nGauss,kmag,xk11,xk12,xk22,Tflag

 205      format(/' boundary springs (2D, 2DOFs)'/ 
     +     ' gauss pts/dir .........................',i12,/ 
     +     ' kmag....................................',1p,e16.4,/ 
     +     ' K_11....................................',1p,e16.4,/ 
     +     ' K_12....................................',1p,e16.4,/ 
     +     ' K_22....................................',1p,e16.4,/ 
     +     ' Tflag...................................',i12)

          pelem(1) = dble(nGauss)
          pelem(2) = kmag
          pelem(3) = xk11
          pelem(4) = xk12
          pelem(5) = xk22
          pelem(6) = dble(Tflag)
               
      else ! 2D / 3D

            read(iread,*) deg_p,kmag,xk11,xk12,xk13,xk22,xk23,xk33,Tflag
            write(iwrit,206) deg_p,kmag,xk11,xk12,xk13,xk22,xk23,xk33,
     $                     Tflag
 206      format(/' boundary springs (3D, 4DOFs)'/ 
     +       ' degree of precision..........................',i12,/ 
     +       ' magnitude of spring tensor..............',1p,e16.4,/ 
     +       ' K_11....................................',1p,e16.4,/ 
     +       ' K_12....................................',1p,e16.4,/ 
     +       ' K_13....................................',1p,e16.4,/ 
     +       ' K_22....................................',1p,e16.4,/ 
     +       ' K_23....................................',1p,e16.4,/ 
     +       ' K_33....................................',1p,e16.4,/ 
     +       ' Tflag...................................',i12)

            pelem(1) = dble(deg_p)
            pelem(2) = kmag
            pelem(3) = xk11
            pelem(4) = xk12
            pelem(5) = xk13
            pelem(6) = xk22
            pelem(7) = xk23
            pelem(8) = xk33
            pelem(9) = dble(Tflag)
                   

      endif ! 2D/3D
 
      return
          

c     --------------------------------------------------------
c     Build the elemental consistent tangent stiffness matrix (estif)
c       and the elemental RHS (eforc)
c     --------------------------------------------------------
          
 3     if (ndime.eq.2) then
          nGauss = int(pelem(1))
          kmag = pelem(2)
          Tflag = int(pelem(6))
          if (Tflag /= 1) then
            xkelem(:,:) = 0.0d0
            xkelem(1,1) = kmag*pelem(3)
            xkelem(1,2) = kmag*pelem(4)
            xkelem(2,1) = kmag*pelem(4)
            xkelem(2,2) = kmag*pelem(5)
        endif
c     initialize some variables
          estif(:,:) = 0.0d0! elemental tangent stiffness
          eforc(:) = 0.0d0! elemental RHS/residual

          call gauss1(nGauss,sg,wg)
 
          do iinte = 1,nGauss! for all Gauss integration points
              call shape1(sg(iinte),shap,xjaco,.false.,nnode,
     $                          ndime,xelem)
              wtjac = wg(iinte)*xjaco

              if (Tflag == 1) then
                xkelem(:,:) = 0.0d0
                do inode = 1,nnode
                  do idofn = 1,ndime
                    do jdofn = 1,ndime
                      xkelem(idofn,jdofn) = xkelem(idofn,jdofn) 
     $                       + TT_nodal(inode,idofn,jdofn)       
     $                                *shap(2,inode)
                    enddo
                  enddo
                enddo
                xkelem = kmag*matmul(xkelem,xkelem)
              endif

              ievab = 0
              do inode = 1, nnode
                  do idofn = 1, 2
                      ievab = (inode-1)*mdofn + idofn
                      do jnode = 1, nnode
                          do jdofn = 1, 2
                              jevab = (jnode-1)*mdofn + jdofn
                              estif(ievab,jevab)= estif(ievab,jevab) + 
     $                      xkelem(idofn,jdofn)*shap(2,inode) 
     $                         *shap(2,jnode)*wtjac ! spring stiffness
                          end do! jdofn
                      end do! jnode
                  end do! idofn
              end do! inode
          end do! iinte

          ievab = 0
          do inode = 1, nnode
              do idofn = 1, 2
                  ievab =  (inode-1)*mdofn + idofn
                  do jnode = 1,nnode
                      do jdofn = 1, 2
                          jevab =  (jnode-1)*mdofn + jdofn
                          eforc(ievab) = eforc(ievab) 
     $                             - estif(ievab,jevab) 
     $                     *(bc_step*uelem_meas(jdofn,jnode)
     $                     - uelem(jdofn,jnode))
                      end do! jdofn
                  end do! jnode
              end do! idofn
          end do! inode
            
              
      else  ! ndim = 3; triangle surface element.

          deg_p = int(pelem(1))
          kmag  = pelem(2)
          Tflag = int(pelem(9))

         if (Tflag /= 1) then
           xkelem(:,:) = 0.0d0
           xkelem(1,1) = kmag*pelem(3)
           xkelem(1,2) = kmag*pelem(4)
           xkelem(1,3) = kmag*pelem(5)
           xkelem(2,1) = kmag*pelem(4)
           xkelem(2,2) = kmag*pelem(6)
           xkelem(2,3) = kmag*pelem(7)
           xkelem(3,1) = kmag*pelem(5)
           xkelem(3,2) = kmag*pelem(7)
           xkelem(3,3) = kmag*pelem(8)
         endif
c     initialize some variables
          estif(:,:) = 0.0d0! elemental tangent stiffness
          eforc(:) = 0.0d0! elemental RHS/residual
          call gausstri(deg_p,ninte,sg,tg,wg)
          deriv = .TRUE.
          do iinte = 1,ninte ! for all Gauss integration points
              call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.true. 
     +         ,nnode,elnods,ndime,xelem)
              wtjac = wg(iinte)*xjaco


              if (Tflag == 1) then
                xkelem(:,:) = 0.0d0
                do inode = 1,nnode
                  do idofn = 1,ndime
                    do jdofn = 1,ndime
                      xkelem(idofn,jdofn) = xkelem(idofn,jdofn) 
     $                            + TT_nodal(inode,idofn,jdofn)
     $                                      *shap(3,inode)
                    enddo
                  enddo
                enddo
                xkelem = kmag*matmul(xkelem,xkelem)
              endif

              ievab = 0
              do inode = 1, 3 ! 3 nodes
                  do idofn = 1,3 ! 4 dofs (skip pressure)
                      ievab = (inode-1)*4 + idofn
                      jevab = 0
                      do jnode = 1, 3 ! 3 nodes
                          do jdofn = 1, 3 ! 4 dofs (skip pressure)
                              jevab = (jnode-1)*4 + jdofn
                              estif(ievab,jevab)= estif(ievab,jevab) + 
     $                         xkelem(idofn,jdofn)*shap(3,inode) 
     $                         *shap(3,jnode)*wtjac
                          end do! jdofn
                      end do! jnode
                  end do! idofn
              end do! inode


          end do! iinte

          ievab = 0
          do inode = 1, 3 ! 3 nodes
              do idofn = 1,3 ! 4 dofs (skip pressure)
                  ievab = (inode-1)*4 + idofn
                  jevab = 0
                  do jnode = 1, 3 ! 3 nodes
                      do jdofn = 1, 3 ! 4 dofs (skip pressure)
                          jevab = (jnode-1)*4 + jdofn
                          eforc(ievab) = eforc(ievab) 
     $                          - estif(ievab,jevab) 
     $                     *(bc_step*uelem_meas(jdofn,jnode) 
     $                     - uelem(jdofn,jnode))
                      end do! jdofn
                  end do! jnode
              end do! idofn
          end do! inode

             
      endif ! 2D/3D if
          
      return
c!c     --------------------------------------------------------
c!c     Evaluate the elemental RHS/residual (eforc) from a source term
c!c     --------------------------------------------------------
 5     continue
       eforc(:) = 0.0d0! elemental RHS/residual
c!c     turned off for now
       return
c!     -------------------------------------------------------
c!     Evaluate the elemental RHS (eforc) for the dual problem
c!     -------------------------------------------------------
 6      continue ! no contribution from the boundary springs
        eforc(:) = 0.0d0! 
        return
 7      continue
        egrad(:,:) = 0.0d0
        return
 4      return
 8      return
 9      return
 10     return
 11     return
 12     return
 13     return
 14     return
 15     return
 16     return
 17     return
        end  !subroutine elem611

