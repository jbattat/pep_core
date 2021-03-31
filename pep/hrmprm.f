      subroutine HRMPRM(lvec, cycle, prefix, k, names, iskale, type,
     .                  order, xnom, harm)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ibase, j, l, m, n
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1973  subroutine hrmprm
c
c generates names for harmonics
c
c parameters
      integer*2 lvec(1), order
      integer*4 cycle, type, k
      real*10 xnom(1), harm(1)
      real*10 iskale(1)
      character*8 prefix, names(2, 1)
c
c external
      real*10 SKALE
c
c locals
      character*8 zharm1 /'J( )'/, zharm2 /'J(  )'/,
     .     tharm1 /'C( , )'/, tharm2 /'C(  , )'/, tharm3 /'C(  ,  )'/,
     .     temp
      character*1 tem(8)
      equivalence(tem, temp)
      character*1 nums(10)/'1','2','3','4','5','6','7','8','9','0'/
      character*1 ess/'S'/
c
c zonals
      if( type .gt. 1 ) then
c
c tesserals
         m     = 0
         ibase = 0
         do l = 2, order
            ibase = ibase + m
            m     = cycle*l
            n     = 0
            do i = 1, m, cycle
               n = n + 1
               if( lvec(ibase+i) .ne. 0 ) then
                  k = k + 1
                  names(1, k) = prefix
                  iskale(k)   = SKALE(-(l*100+n))
                  xnom(k)     = harm(ibase + i)/iskale(k)
                  if( l .ge. 10 ) then
                     if( n .ge. 10 ) then
                        temp = tharm3
                        j    = mod(n, 10)
                        if( j .eq. 0 ) j = 10
                        tem(7) = nums(j)
                        tem(6) = nums(n/10)
                     else
                        temp   = tharm2
                        tem(6) = nums(n)
                     endif
                     j = mod(l, 10)
                     if( j .eq. 0 ) j = 10
                     tem(4) = nums(j)
                     tem(3) = nums(l/10)
                  else
                     temp   = tharm1
                     tem(3) = nums(l)
                     tem(5) = nums(n)
                  endif
                  if( type .eq. 3 ) tem(1) = ess
                  names(2, k) = temp
               endif
            end do
         end do
      else
         m = cycle*(order - 1)
         l = 1
         do i = 1, m, cycle
            l = l + 1
            if( lvec(i) .ne. 0 ) then
               k = k + 1
               names(1, k) = prefix
               xnom(k)     = harm(i)
               if( l .ge. 10 ) then
                  temp = zharm2
                  j    = mod(l, 10)
                  if( j .eq. 0 ) j = 10
                  tem(4) = nums(j)
                  tem(3) = nums(l/10)
               else
                  temp   = zharm1
                  tem(3) = nums(l)
               endif
               names(2, k) = temp
            endif
         end do
      endif
      return
      end
