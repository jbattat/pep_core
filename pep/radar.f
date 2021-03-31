      subroutine RADAR(mocpar)
 
      implicit none

c m.e.ash   jan 1967    subroutine radar
c radar observations of the sun,moon,planets are processed

c arguments
      integer*4 mocpar
c     mocpar = 0 compar called in midst of least squares iteration
c    read all three data sets iobcon, iobs, iabs1 which become in
c     subroutines of compar iiobcn, iiobs, iiabs1
c     mocpar = 1 compar called at end of least squares iteration to
c                calculate dummy observations
c     read only iobcon which become in subroutines of compar iiobcn

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdateq.inc'
c
c           ncodf=1   radar observation
c           ncodf=2   transponder
c           ncodf=3   differential radar observation (decampli aug 1973)
c                     (feature relative to subradar point if nspot.gt.0)
c           ncodf=19  les-8/9 one way doppler count observable
c                     (first point of observation series is reference,
c                     no effect on orbit fit) (ncode must be 1)
c
c           ncode=1   time delay only
c           ncode=2   time delay and doppler
c           ncode=3   doppler only
c
      include 'difnct.inc'
      include 'empcnd.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'lfix.inc'
      include 'ltrapx.inc'
      real*10 deriv1(296)
      include 'mtrapx.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev
      equivalence (ctrecf,Dstf),(fdev,Dstf(10))
c        dstf(1) =ctrecf=coordinate time at reception (fraction of day)
c        dstf(2) = longitude of subradar point (deg)
c        dstf(3) = latitude of subradar point (deg)
c        dstf(4) = time delay from receive to reflection (sec)
c        dstf(5) = longitude of observed spot away from subradar point
c                  if spot name='&&&&' and nspot=0 (deg)
c        dstf(6) = latitude of observed spot away from subradar point
c                  if spot name='&&&&' and nspot=0 (deg)
c        dstf(10)= fdev =1+fractional frequency offset from atomic
c                        time for unit of time delay measurement
      real*10 tmdly0,dop0
      equivalence (Result(1),tmdly0),(Result(2),dop0)
      include 'kobequiv.inc'
      include 'obstap.inc'
      include 'obstuf.inc'
      include 'param.inc'
      include 'plndta.inc'
      include 'prpgat.inc'
      include 'redobs.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'statsrad.inc'
      include 'tidal.inc'
      real*10 tidstf(72)
      equivalence (tidstf,Dxdhe)
c
c  nprec  = 0 nutation-precession determined midway between
c             receiving and sending times
c  nprec  = 1 nutation-precession determined at both receiving and
c             sending times
c
c  ndprec = 0 derivative of nutation-precession matrix not used in
c             determining velocities
c  ndprec = 1 derivative of nutation-precession matrix used in
c             determining velocities
c
c  ntmdly =-1 classical determination of time delay
c  ntmdly = 0 special relativistic determination of time delay
c             (same as classical because deviation of earth
c             velocity from its average is small)
c  ntmdly = 1 general relativistic determination of time delay
c  ntmdly = 2 general rel., including effect of earth
c
c  ndop   =-1 classical determination of doppler shift
c  ndop   = 0 special relativistic determination of doppler shift
c             (same as classical because deviation of earth
c             velocity from its average is small)
c  ndop   = 1 general relativistic determination of doppler shift
c
c  nmedia =-1 effect of interplanetary media on radar signal is
c             ignored
c  nmedia = 0 effect of interplanetary media on time delay only is
c             included
c  nmedia = 1 effect of interplanetary media on time delay and
c             doppler shift is included
c
c  neatm  =-1 effect of earth atmosphere on radar signal is
c             ignored
c  neatm  = 0 effect of earth atmosphere on time delay only is
c             included
c  neatm  = 1 effect of earth atmosphere on time delay and doppler
c             shift is included
c
c  npatm  =-1 effect of planet atmosphere on radar signal is
c             ignored
c  npatm  = 0 effect of planet atmosphere on time delay only is
c             included
c  npatm  = 1 effect of planet atmosphere on time delay and
c             doppler shift is included
c
c  npshp  =-1 effect of planet shape on radar signal is ignored
c  npshp  = 0 effect of planet shape on time delay is included
c
c  nplng  = -1,0 no long. of subradar pt. calculated
c  nplng  = 1 long of subradar pt. is calc. but not printed out
c  nplng  = 2  long of subradar point is calc. and printed out
c
c  neion  =-1 effect of earth ionosphere on probe signal ignored
c  neion  = 0 effect of earth ionosphere on probe delay is calc.
c  neion  = 1 effect of earth ionosphere on probe dopler is calc.
c
c  iwob   = 0 wobble effect not included
c  iwob   = 1 wobble effect is included
c
c  nlibpr =-2,-1,0,1,2 control integers for determination of lunar
c             libration and partials w.r.t. libration parameters
c             (see mnspt for documentation)
c
c  lopler =-1 phase delay doppler
c  lopler = 1 instantaneous doppler
c
c  calcvl = -1,0 don't calculate transponder and relative
c                velocities unless needed for doppler
c  calcvl = 1    calculate velocities (needed for either ctvary
c                partial or light time iteration correction)
c
c  kob(j) = possible additional control integer, (j=11)
c
c savb(i,1), i=1,numsav,  saved precession-nutation,etc. from input
c                         observation library tape for radar, optical
c                         and first part of transit observation
c savb(i,2) i=1,numsav   these saved quantities from iiobs
c
c
c in following, j=1 is for quantity from observation card data set iiobs
c               j=2 is for quantity from input observation library
c                      tape iiabs1
c
c resltb(1,j) = time delay for radar observation (seconds)
c               right ascension for optical observation (sec.of time)
c               time that first limb of planet cuts first limb of sun
c                 for transit observation (first part) (seconds)
c resltb(2,j) = doppler shift for radar observation (cycles/sec)
c               declination for optical observatoon (sec.of arc)
c               time that second limb of planet cuts first limb of sun
c                 for transit observation (first part) (seconds)
c
c errorb(i,j) = error in resltb(i,j) in same units (i=1,4)
c
c
c           set up observation data set quantities
c           niobc,niobs,niabs1  indicate status of data sets  iiobc,
c           iiobs, iiabs1
c              -1  not to be read
c              0  ready to be read
c              1  already read
c           if  iiobc, iiobs, or iiabs1  are initially zero,  data set
c           is not read in any case.
c
c external functions
      real*10 ANCORS,DOT
      real*4 STORNE
      integer*4 JULDAY

c local variables
      character*4 pound9/'####'/,amper9/'&&&&'/,blnk4/'    '/
      character*2 astrik(2)/'  ', '* '/
      character*24 obsnam(3)/'     TIME DELAY (SEC)   ',
     .     '  DOPPLER COUNT (CYCLES)', '   PULSAR PHASE (CYCLES)'/
      character*16 dopshf/' DOPPLER SHIFT  '/
      character*16 typobs(4)/' DEL/DOP RADAR  ', ' RADIO TRACKING ',
     .                       ' DP.COUNT RADAR ', '   PULSAR PHASE '/
      logical*4 cislun
      real*10 dummy,secbx
      integer*4 i,ictst,idopob,iprlog,iradnm,j,jradnm,
     . klpdl,lfirst,ns
      integer*2 iyr19
c
c*  start=100
c initialization
      Jdsav  = 0
      cislun = Ncp0.eq.3 .or. Ncp0.eq.10
      iradnm = 1
      if(Ncodf.eq.2) iradnm  = 2
      if(Ncodf.eq.19) iradnm = 3
      if(Ncodf.eq.18) iradnm = 4
      jradnm = 1
      if(Ncodf.eq.19) jradnm = 2
      if(Ncodf.eq.18) jradnm = 3
c
c determine data about observation series
      nddiff = 0
      if(Nsite2.eq.0) nddiff = -1
      nprec  = 0
      ndprec = 0
      ictst  = 1
      if(Ncodf.eq.2) ictst = 0
      if(Ict(29).ge.ictst) ndprec = 1
      if(Ict(28).ge.ictst .and. .not.cislun) nprec = 1
      if(nddiff.lt.0) nprec = 1
      Ntrmct = 2
      if(Ncodf.eq.18) Ntrmct  = 4
      if(Jct(71).gt.0) Ntrmct = Jct(71)
 
c calculate factor in relativistic corrections
      Reltrm(1) = 0._10
      Reltrm(2) = 0._10
      if(Reldel.ne.0._10) then
         if(Ncodf.ne.3) then
            ntmdly = 1
            if(cislun) ntmdly = 2
            Reltrm(1) = Reldel*Gmc2*(1._10 + prmter(42))/2._10
            if(Klanb.gt.0 .or. Reldop.gt.0._10) then
               ndop = 1
               Reltrm(2) = Reltrm(1)
            endif
         endif
      endif
 
c zero tide partials since etide & mtide called only by radmn
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
      if(Klanb.le.0) then
         if(Klan.gt.0 .and. Klan.le.u_mxpl) then
            do i = 8, 10
               if(Pcond(i,Klan).ne.0._10) then
                  npshp = 1
                  if(Ncode.eq.3) nplng = 1
               endif
            end do
            if(Ncodf.eq.1 .or. Ncodf.eq.3) then
               if(Spotf.eq.pound9 .or.
     .          ((Spotf.eq.blnk4 .or. Ncodf.eq.3) .and.
     .          (Iabs2.gt.0 .or. npshp.gt.0)) .or.
     .          ABS(1._10-DOT(Pcom(2),Pcom(2))).lt.1.E-3_10) nplng = 1
               if(Ict(30).gt.0) nplng = 2
               if(Ict(2).gt.0.and.mocpar.le.0.and.nplng.gt.1) nplng = 1
 
c make sure required parameters were specified
               if(nplng.gt.0 .and.
     .            Pcom(1)*Pcond(13,Klan)*Pcond(14,Klan)*Pcond(15,Klan)
     .           .eq.0._10)
     .            call SUICID(
     .'MISSING PARAMETERS REQUIRED FOR COMPUTING SUB-RADAR POINT POSITIO
     .N, WARNING IN RADAR', -21)
            endif
         endif
      endif
      Iwob = 1
      do i = 7, 30
         if(Lmrx(i).le.0) goto 100
         if(Lmrx(i).eq.3 .or. Lmrx(i).eq.4) nlibpr = -2
      end do
  100 if(Ilib.gt.0) nlibpr = -nlibpr
      if(nlibpr.le.-1 .and. Jct(26).eq.0) nlibpr = 0
c
c force calculation of velocities if ctvary partial is to
c be computed
      if(Ict(1).gt.0) then
         if(npshp.ge.0 .and. Npshp1.le.0)
     .        call SUICID('NO SHAPE MODEL INPUT, STOP IN RADAR ', 9)
         do i = 1, 100
            if(Lprmx(i).le.0) goto 200
            if(Lprmx(i).eq.72 .or. Lprmx(i).eq.81) then
               if(Ict(4).le.0) then
                  do j = 1, 100
                     if(Mprmx(j).le.0) goto 110
                     if(Mprmx(j).eq.Lprmx(i)) goto 150
                  end do
               endif
  110          if(Lprmx(i).ne.72) then
                  if(Atctsc.eq.0._10) call SUICID(
     .           'CANNOT COMPUTE PARTIAL WHEN ATCTSC=0, STOP IN RADAR ',
     .           13)
                  if(Jct(67).eq.0 .and.
     .                (Ncodf.eq.1 .or. Ncodf.eq.2)) call SUICID(
     .           'SIGNIFICANT ERRORS IN ATCTSC PARTIAL SINCE JCT(67)=0',
     .           13)
                  goto 200
               else
                  calcvl = 1
               endif
            endif
  150    end do
      endif
  200 if(Jct(67).gt.0) calcvl = 1
      lopler = Ict(21)
      if(Ncodf.ne.2) lopler = 1
      if(Ncodf.eq.2 .and. lopler.ne.-1) lopler = 1
 
c set up solar system barycenter computations
      Nswcns = 0
      Nvlcns = 0
      if(Ict(27).gt.0) Nswcns = 1
      if(Ncodf.ne.2) Nswcns   = 0
      if(Ncodf.eq.18) Nswcns  = 1
      if(Ncodf.eq.18 .and. Ict(27).gt.1) Nvlcns = 1
      if(cislun .and. Ict(27).ne.0) Nswcns = 1
c
c set up for binary pulsar
      if(Ncodf.eq.18 .and. Psrprm(10).gt.0._10) call PSRBIN
c
c       calculate obscon(1)
c          obscon(1)=(uplink freq/exciter freq. ) typically 96.
c          obscon(2)=obscon(1)*turnaround ratio
c            turnaround ratio=240/221  s band
c                             880/221  x band
c
      if(Obscon(1).eq.0._10) then
         if(Ncodf.eq.1) Obscon(1) = Obscon(2)
         if(Ncodf.eq.2) Obscon(1) = (221._10/240._10)*Obscon(2)
      endif
 
c set xpdly to planet radius or transponder delay
      Xpdly = 0._10
      if(Klap.gt.0 .and. Klap.le.u_mxpl+1) then
         if(Ncodf.ne.1 .or. Nspot.le.0) then
            klpdl = Klap
            if(klpdl.eq.u_mxpl+1) klpdl = -2
            Xpdly = Pcond(7,klpdl)/Ltvel
            if(Ncodf.eq.2) Xpdly = -.5_10*Pcond(29,klpdl)
         endif
      endif
      if(Nplnt0.eq.10) Pradls = Mcond(7)/Ltvel
c
c determine receiving and sending site quantities
      ns = 2
      if(Nsite2.eq.0) ns = 1
      call ESHAPE(ns)
      dummy = ANCORS(ns)
 
      call DGINIT
      Vswrte(1) = 7.2921164E-5_10*Rc(1)
      Vswrte(2) = 7.2921164E-5_10*Rc(2)
c
c switch for first point of series for les-8/9 doppler
c count observable
      lfirst = -1
c
c write first page heading
      call OBSBEG(typobs(iradnm),Sitf(1),obsnam(jradnm),4,
     .            dopshf, 4)
c
c write efg header record
      if(Jct(30).ne.0) call EFG1
c
c write out setup time
      write(Intern,300) obsnam(jradnm)
  300 format('  GRNWCH  JULIAN  REC UTC (WWV)  UT1-UTC AT-UTC  CT-AT',
     .       8x, a24, 12x, 'DOPPLER SHIFT (CYCLES/SECOND)'/
     .       '   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC',
     .       7x, 'OBSERVED     ERROR     OBS-TH', 9x,
     .       'OBSERVED     ERROR     OBS-TH')
 
c add headers for rngmod if needed
      if(Ncodf.eq.2) then
         if(Ict(67).ne.0) write(Intern,350)
  350    format(87x,'M   BITS   XOR       SUBCODE  PREFIX O-C'/97x,
     .          '(HEX)        USEC      USEC')
         if(Ict(25).ne.0) write(Intern,400)
  400    format(106x,'NRNG  CDLMAX    RAWRNG'/113x, 'USEC      USEC')
      endif
      if(Ncodf.eq.18) write(Intern,500)
  500 format(97x,'DLY FROM SSBC PRLX DLY I-S FREQ'/102x, 2('SEC',8x),
     .       'HZ')
      call OBSBGT
c
c consistency check for observing away from subradar point
      if(Spotf.ne.amper9) goto 800
      if(itime.ne.2) goto 800
      call SUICID(
     .' OBSERVATIONS AWAY FROM SUBRADAR POINT WITH ITIME BAD, STOP IN RA
     .DAR', 17)
c     itime.eq.2 means a1-utc,ut1-utc are read in from cells that
c     longitude and latitude of observed spot come from
c
c     if observed spot name is '&&&&', spots away from subradar point
c     are observed with input longitude and latitude on observation card
c     in fields 58-65,66-72 instead of at-uts,ut-uts. this longitude and
c     latitude away from the subradar point  are in dstf(5),dstf(6)
c     if observation is from observation library tape. if from
c     observation card data in fields 58-65,66-72 are moved into
c     dstf(5),dstf(6).
c
c*  start=400
c           read observation data card
c
c           if err on iiobs, dummy read necessary
  600 write(Iout,700) Iiobs
  700 format(' **** ERROR ON DATA SET IIOBS=', i2,
     .', OBSERVATION CARD RECORD SKIPPED IN RADAR LABEL=310 **********'
     .)
      Line = Line + 1
      read(Iiobs,1000)
      Nerrra = Nerrra + 1
      goto 900
 
c zero propco here
  800 call PRPZRO
      if(mocpar.gt.0 .or. Niobs.ne.0 .or. Iiobs.le.0) goto 1100
  900 read(Iiobs,1000,err=600) Ncodeb(1,1),Ihrb(1,1),
     . Iminb(1,1),secbx,Resltb(1,1),Errorb(1,1),Resltb(2,1),
     . Errorb(2,1),Atutsb(1,1),Ututsb(1,1),Imnthb(1,1),
     . Idayb(1,1),Iyearb(1,1)
 1000 format(i1,2I3,f8.4,f13.7,e7.0,d15.8,e7.0,f8.4,f7.4,3I2)
      Niobs = -1
      Secb(1,1) = STORNE(secbx)
      if(Ncodeb(1,1).gt.0) then
         Fdsb(1,1) = Ihrb(1,1)*3600._10 + Iminb(1,1)*60._10 + Secb(1,1)
         iyr19 = Iyearb(1,1)+ctime*100
         Jdsb(1,1) = JULDAY(Imnthb(1,1),Idayb(1,1),iyr19)
         Jdb(1,1)  = Jdsb(1,1)
         Niobs = 1
      endif
c
c*  start=600
 1100 call OBSRED(mocpar)
      if(Nice.lt.-1) then
c
c printout at end of observation series
         if(Jct(30).ne.0) call EFG2
         call OBSEND
         return
      endif
 
      if(Ncodf.eq.19 .and. Ncode.ne.1) call SUICID(
     .   ' NCODF=19 AND NCODE.NE.1, STOP IN RADAR ', 10)
c
c set up controls for propagation corrections
      Izctl = 0
      if(Jct(2).ge.0) then
         iprlog = 2
         if(lopler.eq.-1) iprlog = -2
         call PRPLOG(iprlog,1)
      endif
 
c check if extra headers needed
      if(Nk1.lt.0 .and. Nlhd.ne.2) then
         if(Nice.ge.0) then
 
c not needed - doppler
            Nlhd = 2
         else if(Nlhd.gt.2) then
 
c no doppler, merge headers
            Nlhd = Nlhd - 2
            do i = 1, 2
               do j = 13, 16
                  Hdr(j,i) = Hdr(j,i + Nlhd)
                  end do
               end do
         endif
      endif
c
c        setup longitude and latitude of observed spot for radar
c        observations away from the subradar point in case that
c        observations are from cards (ib=1). if observations from
c        observation library tape, dstf(5),dstf(6) setup in
c        subroutine obsred
      if(Spotf.eq.amper9 .and. Idumob.eq.-1) then
         Dstf(5) = Atutsb(1,1) + Pcom(6)
         Dstf(6) = Ututsb(1,1)
      endif
c
c*  start=800
c set tmdly0  used by radctl
      if(Idumob.eq.1 .or. Nice.gt.0) tmdly0 = Tguess*2.01_10
c
c calculation of theoretical value of time delay and/or
c doppler shift
      idopob = 0
 1200 call RADCTL(idopob)
      if(Jd.le.-10) goto 800
      if(Jd.gt.0) then
         if(Ict(55).eq.1 .and. Klanb.gt.0) call SUBSCP(idopob)
c
c is this les-8/9 doppler count observable
         if(Ncodf.eq.19) then
            if(lfirst.eq.0) lfirst = 1
            if(lfirst.lt.0) then
               lfirst    = 0
               deriv1(2) = Deriv(2,1)
            endif
            Deriv(2,1) = Freq*(deriv1(2) - Deriv(2,1))
c
c see if this is phase delay doppler observable
         else if(idopob.gt.0) then
c     idopob is set to 0 in radar before radctl called
c     if this is phase delay doppler observable, radctl sets up
c     for start of counting interval and calls mdeldp or sbdldp for
c     calculations. upon return from these routines radctl sets idopob=1
c     and returns to radar. if idopob=1 radar calls partl to calculate
c     partials at start of counting interval. then radar calls radctl
c     again, which if idopob=1 sets idopob=-1 and does calculations at
c     end of counting interval. upon return to radar usual things are
c     done and partl called a second time (actually first time if not
c     phase delay doppler observable) to evaluate partials at the end
c     of the counting interval. the partials of the phase delay doppler
c     observable are the difference of the two calculations times the
c     frequency divided by the counting interval.
c     this first pass through partial is eliminated in radctl if the
c     begin time of the current observation is the same as the final
c     time of the previous observation. in this case old information
c     is merely updated
c
            if(Ict(1).gt.0) then
               if(Idumob.ne.1 .or. Ict(3).ge.0) then
                  Nice = -1
                  call PARTL(1)
                  Nice = 1
                  do i = 1, Numpar
                     deriv1(i) = Deriv(i,1)
                     end do
               endif
            endif
            goto 1200
         endif
c
c idumob=1  dummy mode            idumob=-1  not dummy mode
         if(Idumob.ne.1) then
c
c demod and 'fix' range observable
            if(Nice.lt.0 .and. Ncodf.eq.2) call RNGMOD
         else
 
c for dummy mode only
            Result(1) = Deriv(2,1)
            Result(2) = Deriv(2,2)
         endif
      endif
c
c*  start=1000
c page logic
      call OBSCNT
      if(Jd.le.0) goto 800
      if(Ict(2).le.0 .or. mocpar.gt.0) then
c
c printout observed minus theory
         if(Nice.lt.0) then
            if(Ncodf.ne.19) then
               write(Iout,1110) Imonth,Iday,Iyear,Jds,Ihr,Imin,
     .               Sec, Ututs, Atuts, Ctat, tmdly0,
     .               (Deriv(i,1),i=1,2),astrik(Nast11)
            else
               write(Iout,1120) Imonth,Iday,Iyear,Jds,Ihr,Imin,
     .               Sec, Ututs, Atuts, Ctat, tmdly0,
     .               (Deriv(i,1),i=1,2),astrik(Nast11)
               if(lfirst.le.0 .and.
     .            ABS(tmdly0).gt.1.E-5_10) call SUICID(
     .   ' FIRST DOPPLER COUNT OBSERVABLE NOT ZERO, STOP IN RADAR ',14)
            endif
         else if(Nice.eq.0) then
            write(Iout,1110) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,
     .            Ututs, Atuts, Ctat, tmdly0, (Deriv(i,1),i=1,2),
     .            astrik(Nast11),dop0,(Deriv(i,2),i=1,2),
     .            astrik(Nast12)
 1110       format(i3,'/', i2, '/', i2, i8, 2I3, 4F8.4,
     .             f16.10, 1pe8.1, 1pe13.5, 1A1, 1pd17.9,
     .             1pe8.1, 1pe13.5, 1A1)
 1120       format(i3,'/', i2, '/', i2, i8, 2I3, 4F8.4,
     .             f16.3, 1pe8.1, 1pe13.5, 1A1)
         else
            if(Ncodf.ge.3 .and. Nspot.le.0) then
               write(Iout,1130) Imonth,Iday,Iyear,Jds,Ihr,Imin,
     .               Sec, Ututs, Atuts, Ctat, dop0, (Deriv(i,2),i=1,2),
     .               astrik(Nast12)
 
 1130          format(i3,'/', i2, '/', i2, i8, 2I3,
     .                   4F8.4, 20x, 'BANDWIDTH', 9x,
     .                   1pd17.9, 1pe8.1, 1pe13.5, 1A1)
            else
 
               write(Iout,1140) Imonth,Iday,Iyear,Jds,Ihr,Imin,
     .               Sec, Ututs, Atuts, Ctat, dop0, (Deriv(i,2),i=1,2),
     .               astrik(Nast12)
 
 1140          format(i3,'/', i2, '/', i2, i8, 2I3, 4F8.4,
     .             38x, 1pd17.9, 1pe8.1, 1pe13.5, 1A1)
            endif
         endif
c
c*  start=1500
         Line = Line + 1
      endif
c
c*  start=2000
c calculate partial derivatives w.r.t. quantities to be adj.
      if(Ict(1).lt.0) goto 800
      if(Ict(1).ne.0) then
         if(Idumob.eq.1 .and. Ict(3).lt.0) goto 800
         if(idopob.lt.0) Nice = -1
         call PARTL(1)
c
c is this les-8/9 doppler count observable
         if(Ncodf.eq.19) then
            if(lfirst.le.0) then
               lfirst = 1
               do i = 3, Numpar
                  deriv1(i) = Deriv(i,1)
                  end do
            endif
            do i = 3, Numpar
               if(Lold(i).ne.1) Deriv(i,1)= Freq*(deriv1(i)-Deriv(i,1))
               end do
c
c see if this is phase delay doppler observable
         else if(idopob.lt.0) then
            Nice = 1
            if(Numpar.ge.3) then
               do i = 3, Numpar
                  if(Lold(i).ne.1) Deriv(i,2)
     .                = -(Deriv(i,1) - deriv1(i))*Freq/Tcs
                  end do
 
c save results for next observation
               do i = 3, Numpar
                  deriv1(i) = Deriv(i,1)
                  end do
            endif
         endif
      endif
c
c*  start=2200
c write(obs-th),partials buffer
      if(idopob.lt.0) then
         call TIMINC(Jd1s,Ctrcs,Jd,ctrecf,Cnttm2)
         Dstf(4) = (Rsave(1,1,1) + Rsave(2,1,1))/2._10
      endif
      call COMRIT(0)
      goto 800
      end
