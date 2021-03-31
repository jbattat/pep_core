      subroutine WOBILS(jd, fract, x, y)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   jd
      real      x, y
 
c*** end of declarations inserted by spag
 
 
c
c     m.ash   nov 1970    dummy subroutine wobils
c     compute earth wobble in seconds of arc from international polar
c     motion service (formerly known as international latitude service)
c     36 day tabular interval table second difference interpolation
c     1891 to 1956
c
      real*10 fract
 
      x = 0.
      y = 0.
      return
      end
