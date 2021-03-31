      subroutine LESROT(c, s, e1, e2)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 c, s, w
      integer   i
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   aug 1975    subroutine lesrot
c perform les-8/9 attitude error and offset rotations
c
      real*10 e1(3), e2(3)
 
      do i = 1, 3
         w     = e1(i)
         e1(i) = s*e2(i) + c*e1(i)
         e2(i) = c*e2(i) - s*w
      end do
 
      return
      end
