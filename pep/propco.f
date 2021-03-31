      subroutine PROPCO(ipct, kick)

      implicit none


c*** start of declarations inserted by spag
      integer   i, ic, icall, ipa, ips, istrt, jct3, kobj, lps
      real      rtrn

c*** end of declarations inserted by spag



      integer*4 ipct, kick

c     subroutine to calculate and sum propogation corrections to
c     radio observations as needed
c     t. marshall eubanks october 1976
c     r.b. goldstein  and  r.w. king   may 1978
c     propco(-2) = first call of propco for phase delay doppler or
c                  counted-cycle vlbi corrections, results stored
c                  in tcal2 or fcal0
c     propco(-1) = second call for phase delay doppler or counted-cycle
c                  vlbi correction, stored in tcal1 or fcal
c     propco(0) = phase delay or vlbi phase correction applied
c     propco(1) = propagation correction to delay
c                 tmdly or dfdly equivalenced to deriv(2,1) in /ltrap/
c     propco(2) = propagation correction to doppler or diff. delay rate
c                 dop equivalenced to deriv(2,2) in /ltrap/
c
c     a note on phase delay  and ipct
c     npath = phadop or fercnt control
c     npath = 1  : ipct = -2
c     npath = 2  : ipct = -1
c
c     kick = 1  propco called by radar link
c     kick = 4  propco called by fermtr link
c
c
c
c     /prpgat/ stores the active propagation corrections
c     cal(i) is the ith correction to delay or rate
c            in the appropriate units
c     scal(i) is the scatter of the ith correction
c     ical(i) is an integer flag to control use of cal(i) and scal(i)
c     cal,scal,ical are all on the obslib tape type four record
c     jcal(i) is an integer flag to be set in the pep input stream
c     jcal is not on the obslib tape
c
c     ------------standard logic for correction application -----------
c     unless overridden in input stream :
c        standard logic for calculation
c     1) if a correction exists on tape,do not recalculate it
c     2) calculate a correction if  you are going to use it
c        standard logic for use:
c     1) use only the highest ranked corrections
c     2) use every compatible class of correction
c
c        the ranking of corrections is actually done in prplog and
c        is inherent in the order of the pointer vectors contained there
c
c
c     use of jcal
c     jcal is an override vector set in the pep input stream
c     jcal (i) is a two digit number
c     each digit of jcal is independent of the other
c     the first (ones) digit of jcal(i) refers
c     to the calculation of corrections
c     the second (tens) digit of jcal(i) refers to the use of correction
c     ones digit =
c                  1 : do not calculate the ith correction
c                  2 : use standard logic  for correction calculation
c                  3 : calculate the correction
c
c     second digit =
c                   1 : do not use the ith correction
c                   2 : use standard logic for correction use
c                   3 : use the correction regardless of rank
c
c     phase delay correction in propco :
c     phase delay dop is a method of calculating
c     the dop theoritical through differenced range measurements
c    similar logic is used in the interferometry link to calculate
c     the accumulated cycle-count observable, except that the
c     observable (and its corrections) are stored in the delay
c     slot, not the rate slot
c     when :
c     ipct = -2 then tcal2 (i) = cal(i)
c     ipct = -1 then tcal1 (i) = cal(i)
c     ipct = 0  if kick=1, then dop= dop + (tcal2-tcal1)/count time
c               if kick=4, then difnct= difnct + fcal - fcal0
c
c
c      - - - - important note- - - - - - - - - - -
c      by convention tmdly units are seconds (of range)
c      dop units are not hz but hz/frequency
c
c
c     the assignment of  the cal vectors
c     cal(i),scal(i),ical(i),jcal(i)
c     all refer to the same correction
c     i =    Correction Type            Rank   Units
c     1    Static neutral atmosphere      2     delay  (sec)
c     2    Static neutral atmosphere      2     rate  (sec/sec)
c
c     3    Active  terrestrial            2     delay
c     4    neutral atmosphere             2     rate
c
c     5    Passive terrestrial            3     delay
c     6    ionosphere                     3     rate
c
c     7    Active terrestrial             3     delay
c     8    ionosphere                     3     rate
c
c     9    Solar plasma                   1     delay
c    10    calculated in MEDIA                  rate
c
c    11    Static SX from rngns          20     delay
c    12    Calibrations from dopns       20     rate
c
c    13    Extrapolated                  17     delay
c    14    plasma corrections (SX)       17     delay
c
c    15    Range cal. via integrated doppler 18  delay
c    16    not used
c
c    17    RANCAL                         5     delay
c
c    19    Static neutral                 2     delay
c    20    planetary atmosphere                 rate (not coded)
c
c    21    Solar plasma                  25     delay
c    22    M. Eubanks model              25     rate
c
c    23    Waveform distortion            5     delay
c
c    25    Interstellar plasma            1     delay
c
c
c
c     notes on extrapolated corrections :
c     starting with viking 76 corrections derived
c     from an experiment on one spacecraft may be
c     applied to other spacecraft,and corrections
c     may be extrapolated in time (due to a
c     temporary loss of the ranging code,for example)
c
c
c     use of sumcor :
c     sumcor stores the total correction applied to the pep observables
c     on the obslib tape type 4 data record
c     ifmany runs are being made with the same corrections,
c     these corrections can be stored on the obslib tape and do not
c     have to be recalculated
c     sumcor (1) = total delay correction
c     sumcor (2) = total rate correction
c
c        the correction applied by propco is always put
c        in the sumcor vector:
c             jct(4)=1  use sumcor that exists
c             jct(4)=0  recalculate sumcor
c
c     notes on calibrations:
c     calibrations are defined as those corrections having a rank
c     ical(i) = 5
c     calibrations are not included in the search for the
c     highest ranked corrections
c     calibrations,like all other corrections,are added
c     to the computed value of the observable
c
c     commons prpgat contains the propogation
c     corrections stored in the type 4 records,
c     the jcal vectors entered in the pep input stream ,
c     and communicates between various pep sub-routines
c
c     note : prmred gives details on use of jct(3)
c
c        commons
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obscrd.inc'
      include 'kobequiv.inc'
      include 'prpgat.inc'
c
c local variables
      real*10 f(2)
      real*4    patcor(2), eatcor(2), eiocor(2), cor(2)
      real*4    dumcor(2)
      logical*1 init/.false./
      logical*1 awarn/.false./
      logical*1 iwarn/.false./
      integer*2 naus, iaus(20)
      logical   active
      character*1 blank/' '/
c
c initial logic
c set up of p/o controls
      if( .not. init) then
         jct3 = Jct(3)
         ips  = mod(jct3, 2)
         ipa  = mod(jct3/16, 2)
         init = .true.
      endif
c
c jct(2) and jct(4) logic:
c       propco off:         => zero sumcor
c       propco on: jct(4)=0 => recalculate sumcor
c                  jct(4)=1 => use existing sumcor
c

c decide whether to make corrections
      if(Jct(2).le.0 .or. Ncal.le.0) then

c propco turned off
         Sumcor(1) = 0.
         Sumcor(2) = 0.
         return
      endif
      if(Jct(4).eq.1) then
         if(ipct.lt.0) return
         goto 200
      endif

      icall = 1
      if(ipct.eq.2) icall = 2

c begin loop kobj= 1,2
      kobj = 1

  100 if(ipct.eq.0) then
c
c third call for phase-delay difference, form cals
         do i = icall, 10, 2
            if(Calcor(i)) call PRPSTR(ipct, i, kick, kobj, dumcor)
         end do
         i = 18 + icall
         if(Calcor(i)) call PRPSTR(ipct, i, kick, kobj, dumcor)
      else
c
c calculate various quantities needed
c   a. zenith angle and rate
c   b. frequencies for dispersive media corrections
c
         if(Izctl.gt.0) call PRPZEN(kobj)
         call PRPFRQ(kick, kobj, f)
c
c loop over calcor for all internally calculated correction
c for tmdly and doppler
c
c first neutral atmosphere
c
         do i = icall, 4, 2
            if(Calcor(i)) then
               if(itime.ne.5) then
                  active = (i.ge.3)
                  call EATCTL(icall, active, kick, kobj, eatcor)
                  call PRPSTR(ipct, i, kick, kobj, eatcor)
               else
                  if( .not. awarn) call SUICID(
     . ' WARNING: PROPCO CAN''T CALC. ATMOSPHERE/IONOSPHERE CORR.''S FOR
     . NORMAL POINTS', -19)
                  awarn = .true.
               endif
            endif
         end do
c
c now ionosphere
         istrt = 4 + icall
         do i = istrt, 8, 2
            if(Calcor(i)) then
               if(itime.ne.5) then

                  active = (i.ge.7)
                  call EIOCTL(icall, active, kick, kobj, f, eiocor)
                  call PRPSTR(ipct, i, kick, kobj, eiocor)
               else
                  if( .not. awarn) call SUICID(
     . ' WARNING: PROPCO CAN''T CALC. ATMOSPHERE/IONOSPHERE CORR.''S FOR
     . NORMAL POINTS', -19)
                  awarn = .true.
               endif
            endif
         end do
c
c now interplanetary media
         i = 8 + icall
         if(Calcor(i)) then
            if(kick.ne.4) then
               call MEDIA(icall, rtrn)
               cor(1) = rtrn
               cor(2) = 0.
               call PRPSTR(ipct, i, kick, kobj, cor)
            else
               if( .not. iwarn) call SUICID(
     .' WARNING:  PROPCO CANNOT CALCULATE INTERPLANETARY MEDIA CORRECTIO
     .NS FOR INTERFEROMETRIC OBSERVABLES ', -25)
               iwarn = .true.
            endif
         endif
c
c now planetary atmospheres
         i = 18 + icall
         if(Calcor(i)) then
            call PATCTL(icall, kobj, patcor)
            call PRPSTR(ipct, i, kick, kobj, patcor)
         endif
      endif

      if(nddiff.le.0 .or. kobj.eq.2) then
c
c use corrections to form sumcor
         if(ipct.lt.0) return
         naus = 0
         ic   = 2
         if(ipct.eq.1) ic = 1
         if(ipct.eq.0 .and. kick.eq.4) ic = 1
         Sumcor(ic) = 0.
         do i = ic, Ncal, 2
            if(Uscor(i)) then
               Sumcor(ic) = Sumcor(ic) + Cal(i)
               if(ipa.ne.0) then
                  naus = naus + 1
                  iaus(naus) = i
               endif
            endif
         end do

         if(ipa.ne.0 .and. naus.gt.0) then
            lps = (naus - 1)/6 + 1
            if(Line + lps.gt.57) call OBSPAG
            write(Iout, 110) (blank, iaus(i), Cal(iaus(i)), i=1,naus)
  110       format(' SUMCOR FORMATION:',
     .             (t20,6(a1,'CAL(',i2,')=',1pe9.2)))
            Line = Line + lps
         endif
      else
         kobj = 2
         goto 100
      endif

  200 if(ips.ne.0) then
         if(Line.gt.58) call OBSPAG
         icall = 2
         if(ipct.eq.1 .or. (ipct.eq.0 .and. kick.eq.4)) icall = 1
         write(Iout, 250) icall, Sumcor(icall)
  250    format(' SUMCOR(', i1, ')=', 1pe20.10)
         Line = Line + 1
      endif

      return
      end
