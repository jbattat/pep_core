      subroutine ACHECK(l, nice, n)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, n, n1, nice
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   aug 1966   subroutine acheck
c check consistency of l-vectors controling adjustment to parameters
c (solar system or body)
c
      integer*2 l(n)
      include 'fcntrl.inc'
      nice = 0
      n1   = 1
      if( l(1) .gt. 0 ) then
         Nparam = Nparam + 1
         if( n .ge. 2 ) then
            do i = 2, n
               if( l(i) .gt. 0 ) then
                  Nparam = Nparam + 1
                  if( l(i) .le. l(i-1) ) nice = nice + 1
               else
                  n1 = i
                  go to 100
               endif
            end do
         endif
         return
      endif
  100 do i = n1, n
         if( l(i) .gt. 0 ) nice = nice + 1
      end do
      return
      end
