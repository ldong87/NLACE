c**********************************************************************
      subroutine elem309 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c*********************************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer ielem, itask, nnode
      integer l, inter
      integer logflag
      integer elnods(mnode) 
      double precision pelem(*), eforc(*), xelem(ndime,*)
      double precision elemdata_nodal(nset_nodal,*), uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
      double precision estif(mevab,*)
      double precision fpress
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
c     fpress : fluid pressure on each patch
c     logflag: flag for log approach(0/1:off,on)
 2    read(iread,*) l, inter, fpress, logflag
      write(iwrit,205) l, inter, fpress, logflag
 205  format(/' finite elasticity (3D,4DOF) '/
     +       ' gauss pts/dir .........................',i12,/
     +       ' interior node .........................',i12,/
     +       ' fluid pressure  .......................',1p,e16.4,/
     +       ' logflag ...............................',i12,/)
      pelem(1) = dble(l)
      pelem(2) = dble(inter)
      pelem(3) = fpress
      pelem(4) = dble(logflag)
      return
      
c     -----------------------------------------------------------------
c     Build the elemental lhs and rhs arising from pressure term 
c     -----------------------------------------------------------------
 3    l     = int(pelem(1))
      inter = int(pelem(2))
      fpress = pelem(3)
      logflag = int(pelem(4))
      eforc(1:mevab) = 1.0d0;
      return
 4    return
 5    continue
      eforc(1:mevab) = 0.0d0
      return
 6    return
 7    return 
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
