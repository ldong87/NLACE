c**********************************************************
      subroutine elem618 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c     Projec measured displacement to calculate stress components
c     2D, Linear, Incompressible, Plane-Stress
c     input: displacement ux(node)= elemdata_nodal(1,node)
c                         uy(node)= elemdata_nodal(2,node) 
c                         mu(node)= elemdata_nodal(3,node) 
c     output: at each nodes dofn = s11 s22 s12 s21
c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem
      parameter (mnode_elem = 9,mquad_elem=26)
      integer ielem, itask, iset
      integer l, ninte, iinte, inode, jnode, nodei, nodej, nGauss
      integer ievab, jevab, jdofn, idofn
      double precision xjaco, wtjac, mu
      double precision xforc, e11, e22, e12
      double precision shap(3,mnode_elem)
      double precision sg(mquad_elem),tg(mquad_elem),wg(mquad_elem)
      double precision temp_primal(3), temp_dual(3)
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

      go to (1,2,3,3,5,6,6,8,9,10,11,12,13,14,14,16,17), itask
      
c     -----------------------------
c     initialize element parameters
c     -----------------------------
 1    elemvec_ndofn(1:nnode) = 4

      return
      
c     --------------------------------------------------------
c     input material properties
c     --------------------------------------------------------
c     l     : no. of gauss pts/direction

 2    read(iread,*) l
      write(iwrit,205) l
 205  format(4x,'Projected displacement'//
     +       4x,'gauss pts/dir        ',i10)
      pelem(1) = dble(l)
      return

c     ------------------------
c     from element stiffness
c     ------------------------
 3    nGauss  = int(pelem(1))
      estif(:,:) = 0.0d0 !
      eforc(:) = 0.0d0 !
       if (nnode.eq.3) then
       call gausstri(nGauss,ninte,sg,tg,wg)! Gausstri for cells
       else
       call gauss2(nGauss,ninte,sg,tg,wg)
       endif
   
c
      do iinte = 1, ninte
         call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.false.,nnode,
     $      elnods,ndime,xelem)
         wtjac = wg(iinte)*xjaco
c
         ievab = 0
         do inode = 1, nnode
            do idofn = 1,elemvec_ndofn(inode)
               ievab = ievab +1 
               jevab = 0
               do  jnode = 1, nnode
                  do jdofn = 1,elemvec_ndofn(jnode)
                     jevab = jevab+1
                     if (idofn.eq.jdofn) then
                        estif(ievab,jevab) = 
     $                       estif(ievab,jevab) + 
     $                       shap(3,inode)*shap(3,jnode)
     $                       *wtjac
                     endif
                  enddo         !jdofn
               enddo            !jnode
            enddo               !idofn
         enddo                  !inode
c     
      enddo                     !iinte
c     -------------------------------------
c     compute element internal force vector
c     -------------------------------------

      do  ievab =  1,4*nnode
         jevab = 0
         do  jnode = 1,nnode
            do  jdofn = 1,4
               jevab = jevab + 1
               xforc = estif(ievab,jevab)*uelem(jdofn,jnode)
               eforc(ievab) = eforc(ievab) + xforc
            enddo
         enddo
      enddo 
        

      return
cc     --------------------------------------------------------
cc     Evaluate the elemental RHS/residual (eforc) from a source term
cc     --------------------------------------------------------
      
 5    nGauss  = int(pelem(1))
      eforc(:) = 0.0d0
c     the right hand side is entered as a body force
        if (nnode.eq.3) then
        call gausstri(nGauss,ninte,sg,tg,wg)
        else
        call gauss2(nGauss,ninte,sg,tg,wg)
        endif
       do iinte = 1, ninte
         call shape2(ielem,sg(iinte),tg(iinte),shap,xjaco,.false.,nnode,
     $      elnods,ndime,xelem)
           wtjac = wg(iinte)*xjaco
           e11=0.0d0
           e22=0.0d0
           e12=0.0d0
            mu=0.0d0
          do inode = 1, nnode
             e11=e11 + shap(1,inode)*elemdata_nodal(1,inode)
             e22=e22 + shap(2,inode)*elemdata_nodal(2,inode)
             e12=e12 + 0.50d0*(shap(2,inode)*elemdata_nodal(1,inode)+ 
     $                         shap(1,inode)*elemdata_nodal(2,inode))
              mu= mu + shap(3,inode)*elemdata_nodal(3,inode)
          enddo
          do jnode=1,nnode
              ievab= 4*(jnode-1) 
             eforc(ievab+1) = eforc(ievab+1) + shap(3,jnode)*wtjac
     $           * mu *( 4*e11 + 2*e22 ) 
             eforc(ievab+2) = eforc(ievab+2) + shap(3,jnode)*wtjac
     $           * mu *( 2*e11 + 4*e22 ) 
             eforc(ievab+3) = eforc(ievab+3) + shap(3,jnode)*wtjac
     $           * mu *(  2*e12 ) 
             eforc(ievab+4) = eforc(ievab+4) + shap(3,jnode)*wtjac
     $           * mu *(  2*e12 ) 
          enddo !jnode
         

        enddo ! iinte

      return
 6    return
 8    return
 9    return
 10   return
 11   return
 12   return
 13   return
 14   return
 16   return
 17   return
      end
