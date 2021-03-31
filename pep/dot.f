      real*10 function DOT(v1, v2)
 
      implicit none
 
c
c
c---purpose:     to find the dot product of vector v1 with vector v2.
c
c---references:  none
c
c---input
c   arguments:   v1     = first vector in the dot product
c                         (double precision)
c                v2     = second vector in the dot product
c                         (double precision)
c
c---output
c   arguments:   dot    = dot product
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
c---notes:       vectors must be of length 3
c
c-----------------------------------------------------------------------
 
      real*10 v1(3), v2(3)
 
      DOT = v1(1)*v2(1) + v1(2)*v2(2) + v1(3)*v2(3)
 
      end
