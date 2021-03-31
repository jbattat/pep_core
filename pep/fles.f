      subroutine FLES(azim, elev, rlam, i8or9, f)
 
      implicit none
c
c c.moulton   aug 1975   dummy subroutine fles

c arguments
      real*10 azim, elev, f(3), rlam
      integer i8or9
c
      f(1) = 0.
      f(2) = 0.
      f(3) = 0.
      return
      end
