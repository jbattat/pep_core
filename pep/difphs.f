      subroutine DIFPHS(mphase, mplnt, ndif, cdf)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, m
 
c*** end of declarations inserted by spag
 
 
c subr. difphs - j.f.chandler - 1980 may
c compare input phase corrections with saved (see difnom)

c arguments
      integer mphase
      integer*2 mplnt
      integer*4 ndif(2)
      logical*4 cdf(5)
 
c array dimensions
      include 'globdefs.inc'
c commons
      include 'anctrl.inc'
      include 'phase.inc'
      include 'phasem.inc'
      include 'restor.inc'
      include 'wrkcomrs.inc'
 
      character*8 sitser
      character*6 phsn/'...PHS'/
c
 
      do while( mphase.lt.Mumphs )
         if(mplnt.ne.Mplphs(mphase+1)) return
         mphase = mphase + 1
 
         do i = 1, Numphs
            if(Nplphs(i).eq.mplnt .and. Phsit(i).eq.Phsit1(mphase)
     .               .and. Phser(i).eq.Phser1(mphase)) then
               sitser(1:4) = Phsit(i)
               sitser(5:8) = Phser(i)
               m = mplnt
               call EBCDI(m, phsn, 3)
               Nsav = Lphs1(i) - 1
               call DIFPM4(Aphase(1,i), Aphs1(1,mphase),
     .                     Lphs(1,i), 9, ndif, cdf, sitser, phsn)
               go to 100
            endif
 
            end do
  100    end do
 
      return
      end
