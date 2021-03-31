      subroutine RFILDA(b, side, ipoch, nparam, temp, lhsflg, rhsflg)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ipoch, ist, j, k, nparam
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1974  subroutine rfilda
c
c read in rhs and b from proper place
c start of ipoch epoch is at ipoch*(lhsflg*nparam + rhsflg)
c
c parameters
      real*10 side(nparam), temp(nparam), b(1)
      logical*4 lhsflg, rhsflg
c
c common
      include 'filtda.inc'
c
c get starting point
c     ist = ipoch*(lhsflg*nparam + rhsflg) + 1
      ist = 0
      if(rhsflg) ist = 1
      if(lhsflg) ist = ist + nparam
      ist = ipoch*ist + 1
      if(rhsflg) read(Mfile, rec = ist) side
      if(lhsflg) then
         do i = 1, nparam
            k = i*(i - 1)/2
            read(Mfile, rec = ist + i) temp
            do j = 1, i
               b(k + j) = temp(j)
               end do
            end do
      endif
      return
      end
