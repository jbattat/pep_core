      subroutine NUTAT(tjd, dpsi, deps)
 
      implicit none
 
c
c m.e.ash   oct 1967    subroutine nutat
c evaluate earth nutation in seconds of arc
c (explanatory supplement to the ephemeris, pp.44-45)
c
      real*10 tjd, dpsi, deps
c tjd =julian ephemeris date
c dpsi=nutation in longitude
c deps=nutation in obliquity
      real*10 t, LDPSI, LDEPS, SDPSI, SDEPS
 
      t    = (tjd - 2415020.0_10)/36525.0_10
      dpsi = (LDPSI(t) + SDPSI(t))*1.0E-4_10
      deps = (LDEPS(t) + SDEPS(t))*1.0E-4_10
      return
      end
