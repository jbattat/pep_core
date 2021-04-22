      subroutine FERMTR(mocpar)
 
      implicit none
 
c
c           m. ash and b. preston, feb. 1970, subroutine fermtr
c     processes interferometer observations - modified dec 1971
c     by r.king to process differential n-count (quasi-vlbi) observable
c
c           fermtr is analogous to subroutine radar (see radar for
c              further comment cards)
c
c arguments
      integer*4 mocpar
c     mocpar = 0 compar called in midst of least squares iteration
c     read all three data sets iobcon, iobs, iabs1 which become in
c     subroutines of compar iiobcn, iiobs, iiabs1
c     mocpar = 1 compar called at end of least squares iteration to
c                calculate dummy observations
c     read only iobcon which become in subroutines of compar iiobcn
c
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdat.inc'
      include 'difnct.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'ltrapx.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 fdev,freq2,tmdly0
      equivalence (Result(1),tmdly0)
      equivalence (Dstf(5),freq2),(Dstf(10),fdev)
c dstf(5) =freq2=reference frequency for second object
c dstf(10)= fdev =1+fractional frequency offset from atomic
c           time for unit of time delay measurement
      include 'kobequiv.inc'
      include 'obstap.inc'
      include 'plndta.inc'
      include 'prpgat.inc'
      include 'redobs.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'statsrad.inc'
      include 'tidal.inc'
      real*10 tidstf(72)
      equivalence (tidstf(1),Dxdhe(1,1))
c
c  nprec  = 0 not used
c  nprec  = 1 nutation-precession determined at receiving time at
c             first site
c
c  ndprec = 0 derivative of nutation-precession matrix not used in
c             determining velocities
c  ndprec = 1 derivative of nutation-precession matrix used in
c             determining velocities
c
c  kob(3) = ntmdly in radar
c  kob(4) =   ndop in radar
c
c  nmedia =-1 effect of interplanetary media on signal is ignored
c  nmedia = 0 effect of interplanetary media on delay included
c  nmedia = 1 effect of interplanetary media on delay and delay
c             rate included
c
c  neatm  =-1 effect of earth atmosphere on signal is ignored
c  neatm  = 0 effect of earth atmosphere on delay only is
c             included
c  neatm  = 1 effect of earth atmosphere on delay and delay rate
c             is included
c
c  npatm  =-1 effect of planet atmosphere on signal is ignored
c  npatm  = 0 effect of planet atmosphere on delay only is
c             included
c  npatm  = 1 effect of planet atmosphere on delay and delay rate
c             is included
c
c  npshp  =-1 effect of planet shape on signal is ignored
c  npshp  = 0 effect of planet shape on signal is included
c
c  neion  =-1 effect of earth ionosphere on probe signal ignored
c  neion  = 0 effect of earth ionosphere on probe delay is calc.
c  neion  = 1 effect of earth ionosphere on probe delay rate calc.
c
c  iwob   = 0 wobble effect not included
c  iwob   = 1 wobble effect is included
c
c  nintrf =-1 interferometry observables are delay and/or delay
c             rate (conventional vlbi)
c  nintrf = 0 interferometry observable is accumulated cycle
c             count
c  nintrf = 1 interferometry observable is phase
c
c  nddiff =-1 undifferenced n-count (one site only)
c  nddiff = 0 regular observable
c  nddiff = 1 doubly-differenced observable
c
c  nlibpr =-2,-1,0,1,2 control integers for determination of lunar
c             libration and partials w.r.t. libration parameters
c             (see mnspt for documentation)
c
c  kob(j) = possible additional control integers, j=15,16)
c
c external functions
      integer*4 JULDAY
      real*4 STORNE

c internal to subroutine fermtr only
      character*24 obsnam(3)/'    DIFF. DELAY (SEC)   ',
     .  '        N-COUNT ( DEG.)', '         PHASE  (CYCLES)'/
      character*16 delrat/'DIFF. DELAY RATE'/
      character*2 astrik(2)/'  ','* '/
      real*10 secbx
      integer*4 i, idopob, iprlog, jfernm, ns
      integer*2 iyr19
c-----------initialization----------------------------------------------
c
c
c           initialize quantities in subroutine dguess (entry dginit)
c           and setup initial guess for delay iterations
      call DGINIT
c           nk1 used to indicate first observation of a series for
c           partial calculations - set = 0 after first call to partl
c           set = 1 after second call to partl
c
c
c-----------determine data about observation series---------------------
c
      nprec  = 1
      ndprec = 1
      if(Ict(29).eq.-1) ndprec = 0
c kob(1-16) initialized =-1 in cmpar3
c zero tide partials since etide & mtide are not called by fersb
      do i = 1, 72
         tidstf(i) = 0._10
         end do
 
c set nmedia, neatm, neion for printout only from jcal
      if(Jclone(9).ge.2) nmedia  = 0
      if(Jclone(10).ge.2) nmedia = 1
      if(Jclone(1).ge.2 .or. Jclone(3).ge.2) neatm = 0
      if(Jclone(2).ge.2 .or. Jclone(4).ge.2) neatm = 1
      if(Jclone(5).ge.2 .or. Jclone(7).ge.2) neion = 0
      if(Jclone(6).ge.2 .or. Jclone(8).ge.2) neion = 1
      if(Ncodf.eq.12 .or. Ncodf.eq.32) nintrf = 0
      if(Ncodf.eq.10 .or. Ncodf.eq.30) nintrf = 1
      jfernm = nintrf + 2
      nddiff = 0
      if(Ncodf.gt.20) nddiff = 1
      if(Nsite2.eq.0) nddiff = -1
      Iwob = 1
      if(Lmrx(7).gt.0) nlibpr = -2
      if(Ilib.gt.0) nlibpr = -nlibpr
      if(nlibpr.lt.0 .and. Jct(26).eq.0) nlibpr = 0
c
c
c-----------determine site quantities for both sites--------------------
c
      ns = 2
      if(Nsite2.eq.0) ns = 1
      call ESHAPE(ns)
c
c
c-----------write first page heading------------------------------------
c
      call OBSBEG(' INTERFEROMETER ', Sitf(1),obsnam(jfernm),4,
     .            delrat, 4)
c
c write out setup time
      write(Intern,200) obsnam(jfernm)
  200 format('  GRNWCH  JULIAN  REC UTC (WWV)  UT1-UTC AT-UTC  CT-AT',
     .       8x, a24, 12x, '    DELAY RATE               '/
     .       '   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC',
     .       7x, 'OBSERVED     ERROR     OBS-TH', 9x,
     .       'OBSERVED     ERROR     OBS-TH')
      call OBSBGT
c
c
c*  start=300
c-----------read observation data card----------------------------------
c
      goto 500
c
c if err on iiobs, dummy read necessary
  300 write(Iout,400) Iiobs
  400 format(' **** ERROR ON DATA SET IIOBS=', i2,
     .', OBSERVATION CARD RECORD SKIPPED IN FERMTR LABEL=310 *********')
      Line = Line + 1
      read(Iiobs,700)
      Nerrra = Nerrra + 1
      goto 600
c
c zero propco
  500 call PRPZRO
      if(mocpar.gt.0 .or. Niobs.ne.0 .or. Iiobs.le.0) goto 800
  600 if(nintrf.ge.0) then
c
c
c observable is differential n-count
         read(Iiobs,650,err=300) Ncodeb(1,1),Ihrb(1,1),Iminb(1,1),
     .        secbx, Resltb(1,1),Errorb(1,1),Resltb(2,1),Errorb(2,1),
     .        Imnthb(1,1),Idayb(1,1),Iyearb(1,1)
  650    format(i1,2I2,f5.1,d24.16,e8.0,d24.16,e8.0,3I2)
         Atutsb(1,1) = 0._10
         Ututsb(1,1) = 0._10
      else
c
c observable is differential delay and/or delay rate
         read(Iiobs,700,err=300) Ncodeb(1,1),Ihrb(1,1),Iminb(1,1),
     .    secbx,Resltb(1,1),Errorb(1,1),Resltb(2,1),Errorb(2,1),
     .    Atutsb(1,1),Ututsb(1,1),Imnthb(1,1),Idayb(1,1),Iyearb(1,1)
 
  700    format(i1,2I3,f8.4,d15.8,e7.0,d15.8,e7.0,f8.4,f7.4,3I2)
      endif
      Niobs = -1
      Secb(1,1) = STORNE(secbx)
      if(Ncodeb(1,1).gt.0) then
         Fdsb(1,1) = Ihrb(1,1)*3600._10 + Iminb(1,1)*60._10 + Secb(1,1)
         iyr19=Iyearb(1,1)+ctime*100
         Jdsb(1,1) = JULDAY(Imnthb(1,1),Idayb(1,1),iyr19)
         Jdb(1,1)  = Jdsb(1,1)
         Niobs = 1
      endif
c
c*  start=600
c           read observation library tape and process iobcon
c           data set to alter error weightings or setup
c           dummy observations
  800 call OBSRED(mocpar)
      if(Nice.lt.-1) then
c
c-----------------------------------------------------------------------
c*  start=9000
c
c           printout at end of observation series
         call OBSEND
 
         return
      endif
c
c move second reference freq from /obstap/ to /obscrd/
c (would not be needed if dstf(5) weren't clobbered each
c observation)
      freq2 = Freqa2(Ntape)
 
c set up counting interval for dummy observations
      if(Idumob.eq.1) Save(28) = Intscc
c
c setup controls for propagation corrections
      Izctl = 0
      if(Jct(2).ge.0) then
         iprlog = 1
         if(Ncode.eq.3) iprlog = 2
         if(nintrf.ge.0) iprlog = -2
         call PRPLOG(iprlog,4)
      endif
c
c*  start=800
c set reflct used by dguess if not saved from obslib tape
      if(Dstf(4).eq.0._10 .or. Dstf(4).eq.1._10) Dstf(4) = Tguess
c dstf(4)=tmdly1 is overwwritten by initial epoch partials for pvm
c observables, so don't use to get initial guess
      if(Idumob.eq.1) tmdly0 = 2._10*Tguess
c  dguess expects tmdly0 to be two-way time delay
c  tguess set for observed body in cmpar3
c
c     logic used in calling ferctl and partl
c
c  (1) conventional vlbi observations (nintrf=-1)
c  observable is instantaneous delay and/or delay rate.  idopob
c  remains = 0.  ferctl and partl called only once.
c
c  (2) counted-cycle vlbi (or one-way doppler) observations
c     (nintrf=0, ncode=1)
c  observable is accumulated cycle count over the interval from utrec
c  - tc at the initial observation to the current utrec. if the cal-
c  culation of transmit frequency for each observation requires
c  use of the delay difference over a count interval to get
c  beta (=v/c), then ferctl is called twice, once at the start
c  and once at the end of the current count interval.  with idopob
c  = 0, fermtr calls ferctl to setup for beginning of count interval.
c  ferctl calls fermn or fersb for calculations, sets idopob = 1, and
c  returns to fermtr.  fermtr then calls ferctl a second time
c  for calculations at the end of the count interval; ferctl calls
c  fermn or fersb and this time returns idobob= -1.  the first
c  pass through ferctl is eliminated if the begin time of the current
c  observation interval is the same as the final time of the previous
c  interval.  in this case, ferctl merely updates the old information
c  and returns from the first call  with idopob =-1.  even if beta
c  is calculated using instantaneous velocities, this logic is still
c  used for the first observation in the series.
c  in either case, partl is called twice for
c  the first observation (nk1=-1,0) but only once thereafter (nk1=1),
c  the partials being calculated using the current delay minus the
c  initial delay with the transmit frequency assumed constant.
c
c  (3) observations are phase (nintrf=1, ncode=1)
c  observable is absolute but ambiguous (modulo 1 cycle) phase at
c  utrec.  since the received frequency may need to be calculated
c  from differenced phases, this observable is calculated in the
c  same way as the counted-cycle observable except that the
c  quantity subtracted at the initial epoch is zero.
c
      idopob = 0
  850 call FERCTL(idopob)
c
c        idopob= 0 when ferctl first called for observation
c        idopob= -1 set in ferctl on second call (interf called
c                  second time only if idopob was set = 1 on
c                  first call
      if(Jd.le.-10) goto 500
      if(Jd.gt.0) then
c
c see if this is differential n-count observable
         if(idopob.gt.0) then
            if(Ict(1).gt.0) then
               if(Idumob.ne.1 .or. Ict(3).ge.0) call PARTL(4)
            endif
            goto 850
 
c idumob=1  dummy mode            idumob=-1  not dummy mode
         else if(Idumob.eq.1) then
 
c for dummy mode only
            Result(1) = Deriv(2,1)
            Result(2) = Deriv(2,2)
         endif
      endif
c*  start=1000
c
c calculation of error quantities and page logic
      call OBSCNT
      if(Jd.le.0) goto 500
c
c save result from last included observation for freqtr calcl.
      if(nintrf.ge.0 .and. Nice.le.0 .and. Nast11.ne.2) then
         do i = 1, 2
            Freqrs(i) = Freqrl(i)
            Freqts(i) = Freqtl(i)
            Reslts(i) = Reslt1(i)
            end do
         Tcs = Cnttim
      endif
c
c
c-----------printout observed minus theory------------------------------
c
      if(mocpar.gt.0 .or. Ict(2).le.0) then
         if(Nice.lt.0) then
            write(Iout,810) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,
     .             Ututs,Atuts,Ctat,Result(1),(Deriv(i,1),i=1,2),
     .             astrik(Nast11)
         else if(Nice.eq.0) then
            write(Iout,810) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,
     .       Ututs,Atuts,Ctat,Result(1),(Deriv(i,1),i=1,2),
     .       astrik(Nast11),Result(2),(Deriv(i,2),i=1,2),astrik(Nast12)
  810       format(i3,'/',i2,'/',i2,i8,2I3,4F8.4,
     .             1pd16.9,1pe8.1,1pe13.5,1A1,1pd17.9,
     .             1pe8.1,1pe13.5,1A1)
         else
            write(Iout,820) Imonth,Iday,Iyear,Jds,Ihr,
     .            Imin,Sec,Ututs,Atuts,Ctat,Result(2),
     .            (Deriv(i,2),i = 1,2),astrik(Nast12)
  820       format(i3,'/',i2,'/',i2,i8,2I3,4F8.4,
     .             38x,1pd17.9,1pe8.1,1pe13.5,1A1)
         endif
 
c*  start=1500
         Line = Line + 1
      endif
c
c
c-----------write(obs-th), partials buffer-----------------------------
c*  start=2000
c
      if(Ict(1).lt.0) goto 500
      if(Ict(1).gt.0) then
         if(Idumob.eq.1 .and. Ict(3).lt.0) goto 500
 
c calculate partial derivatives
         call PARTL(4)
      endif
c
c*  start=2200
      if(Idumob.ne.1 .and. nintrf.ge.0) then
         do i = 1, 2
            Result(i) = Reslt1(i)
            end do
      endif
      call COMRIT(0)
      goto 500
      end
