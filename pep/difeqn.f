      subroutine DIFEQN(ndif,cdf)
 
      implicit none

c subr. difeqn - j.f.chandler - 1980 may
c compare input eq.-eq. corrections with saved (see difnom)
c
c arguments
      integer*4 ndif(2)
      logical*4 cdf(5)
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'anctrl.inc'
      include 'eqenox.inc'
      include 'eqenxm.inc'
      include 'restor.inc'
      include 'wrkcomrs.inc'

c local
      integer*4 i,j,k
      real*4 eq(3)
      character*8 sitser
 
      do k = 1,Mumeqn
         do j = 1,Numeqn
            if(Eqnsit(j).eq.Eqnst1(k) .and. Eqnser(j).eq.Eqnsr1(k)) then
               Nsav = Leqn1(j) - 1
               do i = 1,3
                  eq(i) = deqnx(j,i)
                  end do
               sitser(1:4) = Eqnsit(j)
               sitser(5:8) = Eqnser(j)
               call DIFPM4(eq,Deqnx1(1,k),Leqn(1,j),3,ndif,
     .                     cdf,sitser,'EQEQ  ')
               goto 50
            endif
         end do
 
   50 end do
 
      return
      end
