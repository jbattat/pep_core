      subroutine CORCHN(x, y)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash  jan 1967   subroutine corchn
c change of coordinates from the mean equinox and equator of the
c reference epoch to the true equinox and equator of date
c
      real*10 x(3), y(3)
      include 'nutprc.inc'
 
      do i = 1, 3
         x(i) = Nutpr(i, 1)*y(1) + Nutpr(i, 2)*y(2) + Nutpr(i, 3)*y(3)
      end do
      return
      end
