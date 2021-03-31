      subroutine FRMPHS(rst,mphase,mplnt,ntop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i
 
c*** end of declarations inserted by spag
 
 
c
c m.e. ash - february 1970 - subroutine frmphs
c planetary optical observation phase correction routine for
c forming normal equations from saved normal equations
c
c arguments
      real*10 rst(1000)
      integer*4 mphase,ntop
      integer*2 mplnt

c array dimensions
      include 'globdefs.inc'
c commons
      include 'anctrl.inc'
      include 'phase.inc'
      include 'phasem.inc'
      include 'restor.inc'
 
      do while( mphase.lt.Mumphs )
 
         if(mplnt.ge.0) then
 
            if(mplnt.ne.Mplphs(mphase+1)) return
            mphase = mphase + 1
         else
            mphase = mphase + 1
            if(Mplphs(mphase).ge.0) goto 100
         endif
 
         if(Numphs.gt.0) then
            do i = 1, Numphs
               if(mplnt.ge.0 .or. Nplphs(i).ge.0) then
                  if(Nplphs(i).ne.mplnt) goto 20
               endif
               if(Phsit(i).eq.Phsit1(mphase)) then
                  if(Phser(i).eq.Phser1(mphase)) then
                     Nrst = Lphs1(i) - 1
                     call FRMBDY(rst,Lphs(1,i),Mphs(1,mphase),-9,ntop)
                     goto 100
                  endif
               endif
   20       end do
         endif
 
         do i = 1, 9
            if(Mphs(i,mphase).gt.0) Nsav = Nsav + 1
         end do
  100 end do
 
      return
      end
