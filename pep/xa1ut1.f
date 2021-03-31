      real*10 function XA1UT1(jd, fract)

      implicit none

c
c ash/amuchastegui - august 1969 - real*8 function xa1ut1
c determines a.1-ut1 by extrapolation beyond table
c
      integer   jd
      real*10 fract, ETUTF, UT2UT1

c USNO preliminary times and coordinates of pole stopped
c giving weekly values of ut1-utc
c

c extrapolation
      XA1UT1 = (ETUTF(jd,fract) - 32.15_10) + UT2UT1(jd, fract)

      return
      end
