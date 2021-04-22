      subroutine DIFSPT(mspot, mplnt, ndif, cdf)
 
      implicit none

c subr. difspt - j.f.chandler - 1980 may
c compare input spot coordinates with saved (see difnom)

c arguments
      integer*2 mplnt
      integer*4 mspot,ndif(2)
      logical*4 cdf(5)

c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'anctrl.inc'
      include 'restor.inc'
      include 'sptcrd.inc'
      include 'sptcdm.inc'
      include 'wrkcomrs.inc'

c local
      integer   i
      character*8 sptn/'SPT ****'/
 
      do while( mspot.lt.Mumspt )
         if(mplnt.ne.Msplnt(mspot+1)) return
         mspot = mspot + 1
 
         do i = 1, Numspt
            if(Nsplnt(i).eq.mplnt .and. Spot(i).eq.Spot1(mspot)) then
               sptn(5:8)= Spot(i)
               Nsav     = Lspt1(i) - 1
               call DIFBDY(Spcord(1,i), Spcrd1(1,mspot),
     .                     Lspcrd(1,i), -6, ndif, cdf, sptn, 'COORD ')
               goto 100
            endif
            end do
  100    end do
 
      return
      end
