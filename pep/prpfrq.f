      subroutine PRPFRQ(kick, kobj, f)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 freq2
      integer   i
 
c*** end of declarations inserted by spag
 
 
      integer*4 kick, kobj
      real*10 f(2)
c
c     r.w.king and r.b. goldstein may 1978
c     routine to determine the frequency to be used in ionsphere correct
c     correction calculations
c
c        common
      include 'obscrd.inc'
      real*10 freqtr(2)
      equivalence(Dstf(5), freq2), (Dstf(7), freqtr)
      include 'sitcrd.inc'
      include 'spqind.inc'
c
c if received frequency not continually monitored, use
c nominal for series
      do i = 1, 2
         f(i) = Freq
         if( kobj .eq. 2 ) f(i) = freq2
      end do
      if( Idumob .eq. -2 ) then
c
c radar-link observations
         if( kick .le. 1 ) then
            f(1) = Save(29)*Obscon(1)
            f(2) = Save(29)*Obscon(2)
c
c fermtr-link observations
         else if( Numsav .ge. 32 ) then
            if( kobj .ne. 2 ) then
               f(1) = Save(31)
               f(2) = Save(32)
            else if( Numsav .ge. 39 ) then
               f(1) = Save(38)
               f(2) = Save(39)
            endif
         endif
      endif
c
c
      return
      end
