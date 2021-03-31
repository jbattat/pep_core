      subroutine FERCNT(nrvfrq)
 
      implicit none

c
c           r. king   march 1978  subroutine fercnt
c           calculate quantities used in forming counted-cycle
c           vlbi observables.  replaces code formerly in sub-
c           routine mnterf
c
c parameter
      integer*4 nrvfrq
 
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'coord.inc'
      real*10 dutrec
      equivalence (Angdum(10),dutrec)
      include 'difnct.inc'
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'ltrapobs.inc'
      real*10 difnct
      equivalence (Deriv(2,1),difnct)
      include 'number.inc'
      include 'obscrd.inc'
      real*10 freql(2),xfactr,dbias
      equivalence (Save(29),freql),(Save(33),xfactr),
     .            (Save(34),dbias)
      real*10 ctrecf,freq2,freqtr(2),fdev
      equivalence (Dstf,ctrecf),(Dstf(5),freq2),
     .            (Dstf(7),freqtr),(Dstf(10),fdev)
      integer*2 nintrf, nddiff
      equivalence (Kob(12),nintrf),(Kob(13),nddiff)
c
c quantitites internal to this routine
      real*10 count(2),delt,dfnct0,dfnctf,tx,tz
      integer   i, nst, numsit
c
c----------------------------------------------------------------------
c
c           reslts(i) is n-count from site i at beginning of count
c              interval (last included observation)
c           reslt1(i) is n-count saved for next observation
c           if observation kept, reslts set = reslt1 in fermtr
c           reslt1 always written on iabs2
c
c-----------initialization-------------------------------------------
      if(Nk1.ge.0) then
c
c
c-----------setup for dummy observations
c
         if(Idumob.eq.1) then
            Ifrq = 1
            do i = 1, numsit
 
c store received freqs in save(31-32)
               Save(i + 30) = (1._10 - Beta(i,1))*freqtr(1)
               count(i)     = Save(i + 30)
 
c store received freqs from obj. 2 in save(38-39)
               if(nddiff.eq.1) Save(i + 37) = (1._10 - Beta(i,2))
     .             *freqtr(2)
               if(nddiff.eq.-1) count(i)    = count(i) - freqtr(1)
               if(nddiff.eq.1) count(i)     = freqtr(2)
     .                 *(1._10 - Beta(i,2)) - count(i)
               count(i)  = Cnttim*(dbias + xfactr*count(i))
               Result(i) = Reslts(i) + count(i)
            end do
            nst = 30 + numsit
            if(nddiff.eq.1) nst = 37 + numsit
            if(Numsav.lt.nst) Numsav = nst
         endif
c
c save n-counts from both sites to write on iabs2
         do i = 1, numsit
            Reslt1(i) = Result(i)
         end do
      else
         numsit = 2
         if(nddiff.lt.0) numsit = 1
 
c n-counts at beginning of series are zero by construction
         Reslts(1) = 0._10
         Reslts(2) = 0._10
      endif
c
c setup receive time offsets
      delt = 0._10
      if(Numsav.ge.48 .and. Save(48).ne.1._10) delt = Save(48)
c delt= utrec2 - utrec
c
      if(nrvfrq.le.0) then
         if(Nk1.lt.0) return
         if(Nk1.ne.0) goto 100
      else if(Nk1.ge.0) then
         goto 100
      endif
c
c calculate theoretical value of initial phase difference
      if(nddiff.eq.-1) Difdly(1,1) = Rsave(1,1,1)
      if(Idumob.ne.1) call FERFRQ(nrvfrq)
      dfnct0 = freqtr(1)*(delt + Difdly(1,1))*xfactr
      if(nddiff.eq.1) dfnct0 = dfnct0 - freqtr(2)
     .                             *(delt + Difdly(1,2))*xfactr
c rbsx(1) or eqnx(1), representing initial sampling error, are
c now added to difnct at each observation in ferctl rather than
c subtracted once at the initial epoch in fercnt
c if observable is phase, subtract integer value at initial point
      if(nintrf.eq.1) then
         dfnctf = MOD(dfnct0,1._10)
         dfnct0 = dfnct0 - dfnctf
      endif
      if(Nk1.lt.0) return
c
c calculate accumulated n-count
  100 if(nddiff.eq.-1) Difdly(2,1) = Rsave(2,1,1)
      if(nrvfrq.ne.0 .or. Nk1.ne.0) then
         if(Idumob.ne.1) call FERFRQ(nrvfrq)
      endif
      difnct = freqtr(1)*(delt + Difdly(2,1))*xfactr
      if(nddiff.eq.1) difnct = difnct - freqtr(2)
     .                             *(delt + Difdly(2,2))*xfactr
      difnct = difnct - dfnct0
 
c additional coding for undifferenced n-count (one-way doppler)
      if(Idumob.ne.1) then
         if(nddiff.lt.0) then
            if(Klanb.le.0) call SUICID(
     .' CANNOT CALCULATE UNDIFFERENCED N-COUNT TERMS FOR NON-SPACECRAFT,
     . STOP IN FERCNT', 20)
            tx = Tfrqsb
 
c tz= (t-sbcom(3)) - (t-utrec) = utrec0-sbcom(3)
            tz     = Tfrqsb - dutrec
            difnct = difnct +
     .       ((Save(29)-Pcond(22,Klanb))*dutrec
     .       - Pcond(23,Klanb)*(tx**2-tz**2)/2._10
     .       - Pcond(24,Klanb)*(tx**3-tz**3)/6._10
     .       - Pcond(25,Klanb)*(tx**4-tz**4)/24._10)*xfactr
         endif
         if(Ifrq.eq.0) difnct = 0._10
c this causes deletion of point if xmtr freq no good
c
c
c form differenced observable from n-counts at each station
         if(Nice.gt.0) call SUICID(
     .' CANNOT CALCULATE RATES FOR N-COUNT OBSERVABLE, STOP IN FERCNT 61
     .0  ', 17)
         if(numsit.eq.2) Result(1) = Result(1) - Result(2)
         Result(2) = 0._10
      endif
      return
      end
