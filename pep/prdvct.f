      subroutine PRDVCT(n0, n1, n2, l, m, l1, l2, n)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, n, nn
 
c*** end of declarations inserted by spag
 
 
c
c ash/friedman    nov 1968    subroutine prdvct
c increment subset of normal equations
c
      integer*2 n0, n1, n2, l(100), m(100), l1(100), l2(100)
 
c common
      include 'wrkcompm.inc'
 
      if( n2 .ge. n1 ) then
         if( n0 .ge. n1 ) then
c
c these parameters exist and are in normal equations
            Jfk = l1(n0) - 1
            call PRDBDY(l, m, n)
            return
         endif
      endif
c
c these parameters do not exist or are not in normal equations
      nn = iabs(n)
      do i = 1, nn
         if( m(i) .gt. 0 ) Lbj = Lbj + 1
      end do
 
      return
      end
