      subroutine PLPTRP(kt)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, k, kk, l
      real*10 s1, s2
 
c*** end of declarations inserted by spag
 
 
c everett interpolation for central body position and partials
c j.f. chandler - 1976 sep, revised 1977 may
c input is array yplp, output is array ytp (ptrs in kq/kqt)
c * * * implicit real*8 * * * * *
c
c  nq = number of quantities to be interpolated
c  kq(1-nq) = pointers to output array for each quantity
c  lmvlt = number of coordinates for each quantity (3 or 6)
c
c argument
      integer kt
c kt - code for needed body: 0=>central body, other=>target(kt)

c array dimensions
      include 'globdefs.inc'
c        common
      include 'prtpin.inc'
      include 'tapdtplp.inc'
      include 'yvectplp.inc'
 
      do j = 1,Nqt(kt)
         if(Krt(j,kt).gt.0) then
            l = Kqt(j,kt)
            do i = 1,Lmvlt(kt)
               s1 = 0._10
               s2 = 0._10
               do k = 1, 5
                  kk = 6 - k
                  s1 = s1*P(2) + Yplp(kk,i,2,j,kt)
                  s2 = s2*P(4) + Yplp(kk,i,1,j,kt)
               end do
               Ytp(i,l,kt) = s1*P(1) + s2*P(3)
            end do
         endif
      end do
      return
      end
