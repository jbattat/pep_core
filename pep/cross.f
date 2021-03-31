      subroutine CROSS(v1, v2, v3)
 
      implicit none
 
c
c
c---purpose:     to find the cross product of vector v1 with vector v2
c                and store the result in vector v3. (v3 = v1 x v2)
c
c---references:  none
c
c---input
c   arguments:   v1     = first vector in the cross product
c                         (double precision)
c                v2     = second vector in the cross product
c                         (double precision)
c
c---output
c   arguments:   v3     = vector cross product (v3 = v1 x v2)
c                         (double precision)
c
c---common
c   blocks:      none
c
c---routines
c   called:      none
c
c---ver./date/
c   programmer:  v1.0/08-90/jlh (usno/na)
c
c---notes:       none
c
c-----------------------------------------------------------------------
 
      real*10 v1(3), v2(3), v3(3)
 
c perform cross product
 
      v3(1) = v1(2)*v2(3) - v1(3)*v2(2)
      v3(2) = v1(3)*v2(1) - v1(1)*v2(3)
      v3(3) = v1(1)*v2(2) - v1(2)*v2(1)
 
      return
      end
