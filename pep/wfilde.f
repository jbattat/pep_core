      subroutine WFILDE(b, side, ipoch, nparam, temp, lhsflg, rhsflg)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ipoch, ist, j, k, nparam
 
c*** end of declarations inserted by spag
 
 
c
c subroutine wfilde
c modified for real*16 b from form paul macneil, sept., 1978
c d. white  april 1974  subroutine wfilda
c
c write out a rhs and b to proper place
c start of ipoch epoch is at ipoch*(lhsflg*nparam + rhsflg)
c
c parameters
      real*10   b(1)
      real*10 side(nparam), temp(nparam)
      integer*4 lhsflg, rhsflg
c
c common
      include 'filtda.inc'
c
c zero buffer
      do i = 1, nparam
         temp(i) = 0.0_10
      end do
c
c get starting point
      ist = ipoch*(lhsflg*nparam + rhsflg) + 1
      if( rhsflg .eq. 1 ) write(Mfile, rec = ist) side
      if( lhsflg .eq. 1 ) then
         do i = 1, nparam
            k = i*(i - 1)/2
            do j = 1, i
               temp(j) = b(k + j)
            end do
            write(Mfile, rec = ist + i) temp
         end do
      endif
      return
      end
