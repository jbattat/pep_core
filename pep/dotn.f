      real*10 function DOTN(v1, v2, n)
 
      implicit none
 
c
c
c---purpose:     to find the dot product of the first n elements of
c                vector v1 with the first n elements of vector v2.
c
c---references:  none
c
c---input
c   arguments:   v1     = first vector in the dot product
c                         (double precision)
c                v2     = second vector in the dot product
c                         (double precision)
c                n      = vector dimension 
c                         (integer)
c
c---output
c   arguments:   dotn   = n-dimensional dot product
c                         (double precision)
c
c---common
c   blocks:      none
c
c---routines
c   called:      none
c
c---ver./date/
c   programmer:  v1.0/10-91/jlh (usno/omd)
c
c---notes:       11-25-91 mam (sao)  restriction to n <= 9 removed
c
c-----------------------------------------------------------------------
 
      real*10 v1(*), v2(*)
      integer   n, i
 
      DOTN = 0.0
      if( n .gt. 0 ) then
         do i = 1, n
            DOTN = DOTN + v1(i)*v2(i)
         enddo
      endif
 
      end
