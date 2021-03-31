      subroutine DIFSIT(ndif,cdf)
 
      implicit none

c subr. difsit - j.f.chandler - 1980 may
c compare input site coordinates with saved (see difnom)
c (see difnom)
c
c arguments
      integer*4 ndif(2)
      logical*4 cdf(5)

c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'anctrl.inc'
      include 'restor.inc'
      include 'stcord.inc'
      include 'stcrdm.inc'
      include 'wrkcomrs.inc'

c local
      integer   j,k

      do k = 1, Mumsit
         do j = 1, Numsit
            if(Site1(1,k).eq.Site(1,j)) then
               Nsav = Lscrd1(j) - 1
               call DIFBDY(Scord(1,j),Scord1(1,k),Lscrd(1,j),-6,
     .          ndif,cdf,sitd(j),'COORD ')
               goto 50
            endif
         end do
   50 end do
 
      return
      end
