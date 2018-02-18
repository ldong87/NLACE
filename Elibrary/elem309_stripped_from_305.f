c**********************************************************
      subroutine elem309 (ielem,itask,pelem,nnode,estif,eforc,elnods,
     $     xelem,elemdata_nodal,uelem,uelem_meas)
c**********************************************************
      USE IOUNIT
      USE MAINMEM
      implicit none
      integer mnode_elem, mquad_elem
      parameter (mnode_elem = 4,mquad_elem=30)
      integer ielem,itask,nnode
      integer l,logflag, int_node 
      double precision fpress
c     Localized information
      double precision pelem(*)
      double precision estif(mevab,mevab)
      double precision eforc(mevab)
      integer elnods(mnode)
      double precision xelem(ndime,*)
      double precision elemdata_nodal(nset_nodal,*)
      double precision uelem(mdofn,*)
      double precision uelem_meas(mdofn,*)
c-------------------------------------------------------------

      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17), itask  
      
c     --------------------------------------------------------
c     Initialize the elemental parameters
c     --------------------------------------------------------
 1    elemvec_ndofn(1:nnode) = 4
      buildSymmetricMatrices=.false.
      return


c     --------------------------------------------------------
c     Read and store the material properties
c     --------------------------------------------------------
c     l          : no. of gauss pts/direction
c     ireg       : regularization type (0/1/2/3/4:none/H1/TVD/newTV/power)
c     alpha(1)   : regularization weight for the non-linear parameter
c     alpha(2)   : regularization weight for the shear modulus
c     beta(1)    : extra parameter for TVD - JFD ask Sevan for a better name
c     beta(2)    : extra parameter for TVD - JFD ask Sevan for a better name
c     Treg(1)    : second extra parameter
c     Treg(2)    : second extra parameter
c     tauMult    : stabilization factor
c     s          : stabilization type (0/1:off/on) - JFD use stype instead of s?
 2    read(iread,*) l,int_node,fpress,
     +     logflag
      write(iwrit,205) l,int_node,fpress,logflag 
 205  format(/' finite elasticity (3D,4DOF) '/
     +       ' gauss pts/dir .........................',i12,/
     +       ' interior node ',i12,/
     +       ' fluid pressure ',1p,e16.4,/
     +       ' log-flag(0/1:off/on) .......',i12)
c     check for error in input file

      pelem(1) = dble(l)
      pelem(2) = dble(int_node)
      pelem(3) = fpress
      pelem(4) = dble(logflag)
      return


c     --------------------------------------------------------
c     Build the elemental consistent tangent stiffness matrix (estif)
c      and the elemental RHS/residual (eforc)
c     --------------------------------------------------------
 3    l      = int(pelem(1))

      eforc(:) = 1.0d0

      return


c     --------------------------------------------------------
c     Evaluate the elemantal RHS/residual from a source term
c     --------------------------------------------------------
 5    continue
      eforc(:) = 0.0d0! elemental RHS/residual
c     turned off for now
      return


c     --------------------------------------------------------
c     Build the elemental RHS/residual (eforc) for the dual problem
c     --------------------------------------------------------
 6    continue 
      eforc(:) = 0.0d0! elemental RHS/residual

      return


c     --------------------------------------------------------
c     Compute the objective function (dataMatch+regularization)
c      and its gradient (egrad) on the element
c     --------------------------------------------------------
 7    continue
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

