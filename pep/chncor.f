      subroutine CHNCOR(x, y)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash  nov 1966   subroutine chncor
c change of coordinates from true equinox and equator of date
c to the mean equinox and equator of the reference epoch
c
      real*10 x(3), y(3)
      include 'nutprc.inc'
 
      do i = 1, 3
         x(i) = Nutpr(1, i)*y(1) + Nutpr(2, i)*y(2) + Nutpr(3, i)*y(3)
      end do
      return
      end
