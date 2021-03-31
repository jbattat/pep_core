      subroutine WFILDA(b, side, ipoch, nparam, temp, lhsflg, rhsflg)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ipoch, ist, j, k, nparam
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1974  subroutine wfilda
c
c write out a rhs and b to proper place
c start of ipoch epoch is at ipoch*(lhsflg*nparam + rhsflg)
c
c parameters
      real*10 side(nparam), temp(nparam), b(1)
      logical*4 lhsflg, rhsflg
c
c common
      include 'filtda.inc'
c
c zero buffer
      do i = 1, nparam
         temp(i) = 0._10
         end do
c
c get starting point
c     ist = ipoch*(lhsflg*nparam + rhsflg) + 1
      ist = 0
      if(rhsflg) ist = 1
      if(lhsflg) ist = ist + nparam
      ist = ipoch*ist + 1
      if(rhsflg) write(Mfile, rec = ist) side
      if(lhsflg) then
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
