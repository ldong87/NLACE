c**********************************************************
c List of the elements:
c**********************************************************
      subroutine elem305
c     3D, linear elasticity
c     modified Veronda-Westman
c----------------------------------------------------------
      subroutine elem505
c     2D, plane strain, compressible, linear elasticity
c     modified Veronda-Westman
c----------------------------------------------------------
      subroutine elem606
c     2D, plane stress, incompressible, linear elasticity
c     Veronda-Westman
c----------------------------------------------------------
      subroutine elem607
c     2D, plane stress, incompressible, linear elasticity
c     Veronda-Westman with H1 norm for the objective function
c----------------------------------------------------------
      subroutine elem707
c     old element used by Sevan Goenezen

c**********************************************************
c List of the regularizations implemented:
c**********************************************************
ireg.eq.1 : H1 or Tikhonov                                 0.5*grad(p)^2
ireg.eq.2 : Total Variation (Diminuishing) with offset     sqrt(beta^2+grad(p)^2)
ireg.eq.21: Total Variation (Diminuishing) without offset  sqrt(beta^2+grad(p)^2)-beta
ireg.eq.3 : power reg.                                     ((1+sqrt(beta^2+grad(p)^2)-beta)^powe-1)/powe
ireg.eq.31: power reg. (by Paul Barbone)                   grad(p)^2/(beta^2+grad(p)^2)^(1-powe/2)
ireg.eq.4 : logarithmic reg.                               ln(1+sqrt(beta^2+grad(p)^2)-beta)
ireg.eq.41: logarithmic reg.                               ln(1+(grad(p)/Treg)^2)
ireg.eq.5 : exponential reg. (by Paul Barbone)             1-exp(-(sqrt(beta^2+grad(p)^2)-beta)/Treg)
ireg.eq.6 : fraction reg.                                  0.5*(grad(p)/Treg)^2/(1+(grad(p)/Treg)^2)


