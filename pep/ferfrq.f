      subroutine FERFRQ(nrvfrq)
 
      implicit none

c
c           r. king    june 1977  subroutine ferfrq
c           calculate transmitter frequencies for vlbi n-count
c           observable, either for alsep observations using the mit
c           differential doppler reciever (ddr), or pioneer-venus
c           multiprobe observations using the jpl digital recording
c           assembly (dra) and bandwidth reduction assembly (bra).
c
c parameter
      integer*4 nrvfrq
c           nrvfrq= 0  instantaneous received frequencies not available
c                      at the epoch of observation, calculate the
c                       average over the count interval from phase
c           nrvfrq= 1  instantaneous received frequencies available
c                      at the epoch of observation. use these to calcula
c                      the transmitter frequencies.  however, if
c                      con(16) for a spacecraft is.gt.0, the
c                      xmtr. freq. calculated using received freq.
c                      is overridden by the input power series.

c array dimensions
      include 'globdefs.inc'
c
c        commons
      include 'difnct.inc'
      real*10 tc
      equivalence (tc,Cnttim)
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'number.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct,freq2,freqtr(2)
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (Dstf(5),freq2),(Dstf(7),freqtr)
      integer*2 nddiff
      equivalence (Kob(13),nddiff)
      real*10 freql(2),xfactr,dbias
      equivalence (Save(29),freql),(Save(33),xfactr),
     .            (Save(34),dbias)
      include 'sitcrd.inc'
c
c quantities internal to this routine
      real*10 dtf,epsf,epsf0
      integer   ifrq2, imf
      real*10 freqt(2),freqr(2),count(2),ft1dot,ft2dot,freqt0
      real*10 frqcr1,frqcr2
      character*2 mfrq(2)/'  ','**'/
 
      ifrq2 = 0
c if undifferenced n-count, object must be spacecraft with
c power series model for transmitted frequency
      if(nddiff.ne.-1) then
c
c setup for calculating xmtr drift rate
         if(Nk1.le.0) then
 
c set xmtr freq drift rate = 0. for first point
            ft1dot = 0._10
            ft2dot = 0._10
 
c set difnct correction from drift = 0. for first pt (for print out)
            frqcr1 = 0._10
            frqcr2 = 0._10
c jct(25)=-2 do not calculate drift rate
c =-1 calculate and print out but do not apply drift rate
c = 0 apply drift rate correction to difnct
            if(Jct(25).eq.-1) write(Iout,20) Jct(25)
   20       format(1x,
     .    ' CORRECTIONS TO DIFNCT FOR XMTR DRIFT NOT APPLIED, JCT(25)='
     .    , i3)
            Line = Line + 1
         endif
 
c set initial xmtr freq equal value at first included observation
         if(Nobs(1).eq.1) freqt0 = Freqts(1)
c
c------------calculate transmitter frequency for the first object-------
c
c
c     save(31-32) are the recorded or calculated received
c     frequencies from the first object at sites 1 and 2
         freqr(1) = Save(31)
         freqr(2) = Save(32)
 
c save the current freqr--freqrs set = freqrl in fermtr if pt kept
         Freqrl(1) = freqr(1)
         Freqrl(2) = freqr(2)
c
c calculate average received frequency over count interval
c (beta was calculated from change in delay over interval)
         if(nrvfrq.le.0) then
            if(Nobs(1).ge.1) then
 
c refer the received freq to the midpt of the ddr count interval
               freqr(1) = (freqr(1) + Freqrs(1))/2._10
               freqr(2) = (freqr(2) + Freqrs(2))/2._10
            endif
         endif
c
c calculate a transmitted frequency from the received
c frequency at each station
         freqt(1) = freqr(1)/(1._10 - Beta(1,1))
         freqt(2) = freqr(2)/(1._10 - Beta(2,1))
c
c decide which transmitted frequency to use
         if(Klanb.le.0 .or. Pcond(22,Klanb).le.0._10) then
            imf  = 1
            Ifrq = 1
            if(Nplnt0.ne.10) then
               freqtr(1) = freqt(Ifrq)
               goto 100
 
c checks on hand-copied freq.s for alsep ddr observations
            else if(Nobs(1).ge.1) then
               if(ABS(freqt(1)-freqt0).lt.epsf) then
                  freqtr(1) = freqt(Ifrq)
                  goto 100
               else
                  if(ABS(freqt(2)-freqt0).ge.epsf) goto 300
                  Ifrq = 2
                  freqtr(1) = freqt(Ifrq)
                  goto 100
               endif
            else
c xmtr freq at initial point should be within epsf0 hz of nominal
c and drift by no more than epsf ht during a track
               epsf0 = 50.E3_10
               epsf  = 3.E3_10
               if(ABS(freqt(1)-Freq).lt.epsf0) then
                  freqtr(1) = freqt(Ifrq)
                  goto 100
               else
                  if(ABS(freqt(2)-Freq).ge.epsf0) goto 300
                  Ifrq = 2
                  freqtr(1) = freqt(Ifrq)
                  goto 100
               endif
            endif
         endif
      else
         freqt(1) = 0._10
         freqt(2) = 0._10
      endif
 
c ct - utc ignored for now in freq. calculations
      freqtr(1) = Pcond(22,Klanb) + Pcond(23,Klanb)
     .            *Tfrqsb + 0.5_10*Pcond(24,Klanb)*Tfrqsb**2 +
     .            Pcond(25,Klanb)*Tfrqsb**3/6._10
      Ifrq = 3
      if(nddiff.eq.-1) goto 400
c     set ifrq=3 for printout to signal that input constants, not
c     printed freqt(i) are used for freqtr(1)
c
c           correct xmtr freq for drift, either over the count interval
c           or over the light time difference between stations
c     save xmtr freq for ft1dot calculation at next observation
c     ---freqts(1) set eq. freqtl(1) in fermtr if pt not deleted
  100 Freqtl(1) = freqtr(1)
      if(Nobs(1).ge.1) then
         if(Jct(25).le.-2) goto 200
         dtf = tc
         if(nrvfrq.gt.0) dtf = Difdly(2,1)
c
c correct xmtr freq for drift since midpt of interval
         ft1dot = (freqtr(1) - Freqts(1))/((tc+Tcs)*.5_10)
c calculate the difnct correction due to drift (used only for
c test and printout).
         frqcr1 = .5_10*ft1dot*tc*Difdly(2,1)*360._10
         if(Error(1).gt.0.) then
            if(abs(frqcr1).gt.Error(1)) then
               imf = 2
               goto 200
            endif
         endif
      endif
      if(Jct(25).ge.0) freqtr(1) = freqtr(1) + .5_10*ft1dot*dtf
c flag pt if xmtr freq calculations from the two stations differ
c zero out second xmtr frequency for print out
  200 if(nddiff.le.0) freqtr(2) = 0._10
      if(Nplnt2.gt.0) then
c
c
c-----------calculate transmitter frequency for second object-----------
c
         ifrq2 = Ifrq
         if(Klans1.gt.0 .and. Pcond(22,Klans1).gt.0._10) then
            freqtr(2) = Pcond(22,Klans1) + Pcond(23,Klans1)
     .                  *Tfrqsc + 0.5_10*Pcond(24,Klans1)*Tfrqsc**2 +
     .                  Pcond(25,Klans1)*Tfrqsc**3/6._10
            ifrq2     = 3
         else if(nrvfrq.gt.0) then
 
c freqtr(2)= freq2 in ferctl for dummy mode
            if(Idumob.ne.1) freqtr(2) = Save(37 + Ifrq)
     .          /(1._10 - Beta(Ifrq,2))
         else
            count(1)  = Result(1) - Reslts(1)
            count(2)  = Result(2) - Reslts(2)
            freqtr(2) = (freqr(Ifrq) + (count(Ifrq)/tc-dbias)/xfactr)
     .                  /(1._10 - Beta(Ifrq,2))
         endif
c
c
         Freqtl(2) = freqtr(2)
c
c correct freqtr(2) for drift
         if(Nobs(1).ge.1) then
            if(Jct(25).gt.-2) then
               dtf = tc
               if(nrvfrq.gt.0) dtf = Difdly(2,2)
               ft2dot = (freqtr(2) - Freqts(2))/((tc+Tcs)*.5_10)
               frqcr2 = -.5_10*ft2dot*dtf*Difdly(2,2)*360._10
               if(Error(1).gt.0) then
                  if(abs(frqcr1+frqcr2).ge.Error(1)) imf = 2
                  if(abs(frqcr2)/Error(1).gt.Eps(7)) then
                     imf = 2
                     goto 400
                  endif
               endif
               if(Jct(25).ge.0) freqtr(2) = freqtr(2)
     .             + .5_10*ft2dot*dtf
            endif
         endif
      endif
      goto 400
c     transmitted frequencies are measured in terms of clock time at
c     site no. ifrq.  therefore they should logically be corrected for
c     fdev and rate of ct-at.  this amounts to a small error in the
c     overall scale factor and hence is ignored.
c
  300 Ifrq = 0
      imf  = 2
 
  400 if(Line.gt.57) call OBSPAG
      write(Iout,500) freqt(1),freqt(2),Ifrq,ifrq2,freqtr(1),
     .                 freqtr(2)
      Line = Line + 1
  500 format(' FREQT=', 2F12.0, '  IFRQ=', 2I2, '  FREQTR=', 2F14.2)
      if(Jct(25).ge.-1) write(Iout,600) ft1dot,ft2dot,frqcr1,
     .                            frqcr2, mfrq(imf)
  600 format('+', t90, 'FTDOT=', f6.3, 1x, f6.3, ' FRQCR=', f6.1, 1x,
     .       f6.1, 2x, a2)
      return
c
c
      end
