      subroutine PRDBDY(lbd, mbd, length)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, len, length
 
c*** end of declarations inserted by spag
 
 
c
c ash/friedman    nov 1968    subroutine prdbdy
c
      integer*2 lbd(7), mbd(7)
c
c common
      include 'wrkcompm.inc'
 
      len = 6
      if( length .le. 0 ) len = iabs(length)
 
      do i = 1, len
         if( lbd(i) .gt. 0 ) then
            Jfk = Jfk + 1
            if( mbd(i) .gt. 0 ) then
               Lbj = Lbj + 1
               Ipts(Jfk) = Lbj
            endif
         else if( mbd(i) .gt. 0 ) then
            Lbj = Lbj + 1
         endif
      end do
 
      len = length - len
      if( len .gt. 0 ) call PRDPRM(lbd(7), mbd(7), len)
      return
      end
