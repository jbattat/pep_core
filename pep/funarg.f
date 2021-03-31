      subroutine FUNARG(tjd)
 
      implicit none
 
c
c m.e.ash  sept 1967   subroutine funarg
c fundamental arguments of brown lunar theory are evaluated at
c julian ephemeris date tjd
c
      real*10 tjd
      include 'argfun.inc'
      real*10 FUNCOF
c
c times from fundamental epochs
      T1 = tjd - 2415020.0_10
c
c mean longitude of moon ascending node (big omega)
      Ascm = FUNCOF(T1, 0.71995354167_10, -1, -4.7094228332E-5_10,
     .       +0.432630E-14_10, +1.266E-22_10)
c
c longm-perm (little l)
      L = FUNCOF(T1, 0.82251280093_10, +362, +9.1645684716E-5_10,
     .    +1.913865E-14_10, +8.203E-22_10)
c
c longs-pers (little l prime)
      Lp = FUNCOF(T1, 0.99576620370_10, +27, +3.7778519279E-5_10,
     .     -0.031233E-14_10, -1.900E-22_10)
c
c longm-ascm (big f)
      F = FUNCOF(T1, 0.03125246914_10, +367, +4.8195691688E-5_10,
     .    -0.668609E-14_10, -0.190E-22_10)
c
c longm-longs (big d)
      D = FUNCOF(T1, 0.97427079475_10, +338, +6.3192198393E-5_10,
     .    -0.299023E-14_10, +1.077E-22_10)
      return
      end
