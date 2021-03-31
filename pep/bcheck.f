      subroutine BCHECK(l, nice, n)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, m, n, nice
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   oct 1967   subroutine bcheck
c increment nparam for coordinates,biases,etc plus body initial
c conditions, call acheck for adjustable body parameters
c
      integer*2 l(100)
      include 'fcntrl.inc'
 
      m = 6
      if( n .le. 0 ) m = iabs(n)
      do i = 1, m
         if( l(i) .gt. 0 ) Nparam = Nparam + 1
      end do
      m = n - m
      if( m .gt. 0 ) call ACHECK(l(7), nice, m)
      return
      end
