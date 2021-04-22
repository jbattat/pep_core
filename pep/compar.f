      subroutine COMPAR(mocpar)
 
      implicit none
c
c m.e.ash july 1969  subroutine compar  (main program for comparison
c of theory and observation and computation of partial derivatives)
c
 
c arguments
      integer*4 mocpar
c     mocpar = 0 compar called in midst of least squares iteration
c     read all three data sets iobcon, iobs, iabs1 which become in
c     subroutines of compar iiobcn, iiobs, iiabs1
c     mocpar = 1 compar called at end of least squares iteration to
c                calculate dummy observations
c     read only iobcon which become in subroutines of compar iiobcn

c array dimensions
      include 'globdefs.inc'

c commons
      include 'bdctrl.inc'
      include 'comdat.inc'
c        nplnt0 = planet number of observed planet. allowed numbers are
c                   0 sun
c                   1 mercury
c                   2 venus
c                   4 mars
c                   5 jupiter
c                   6 saturn
c                   7 uranus
c                   8 neptune
c                   9 pluto
c                  10 moon
c                  11,...  possible additional bodies
c                  -4 stars
c        ncode = external observation code. allowed numbers are
c               0 end of series
c                  1 radar observation of time delay
c                       or optical obs. of 1st planet limb contact time
c                       or right ascension or azimuth
c                  2 radar observation of time delay and doppler shift
c                       or optical obs. of 1st & 2nd limb times
c                       or r.a. & dec. or az. & el.
c                  3 radar observation of doppler shift
c                       or optical obs. of 2nd planet limb contact time
c                       or declination or elevation
c            note: 4,5,6 correspond to 1,2,3 but may be used only
c                  for optical observations
c        mnplnt= 0 observed body is moon
c        mnplnt= positive,negative observed body is not moon
c        klanb = 1,...,u_mxpl observed body is probe or satellite
c                  nplnt(klanb) (central body is sun if klan=0)
c        klan  = 1,...,u_mxpl observed body is planet nplnt(klan)
c                  unless klanb.gt.0 in which case central body of
c                  observed body nplnt(klanb) is nplnt(klan)
c        klan  =u_mxpl+1 observed body is moon
c                  unless klanb.gt.0 in which case central body of
c                  observed body nplnt(klanb) is moon
c        klan  =u_mxpl+2 observed body is sun
c                  unless klanb.gt.0 in which case central body of
c                  observed body nplnt(klanb) is earth
c        klans1= 0 observing site is on earth
c        klans1= 1,...,u_mxpl observing site is on body nplnt(klans1)
c        klans1=u_mxpl+1 observing site is on moon
c        klanr = 1,...,u_mxpl  input rotation data for observed body
c                  nplnt(klan) are in storage for nplnt(klanr).  in
c                  particular, observed planet rotation tape is
c                  iplnt(klanr).  for earth or moon, klanr=0 with inut
c                  or ilib nonzero.
c     nrewnd controls non-normal rewinding of planet,embary,moon tapes.
c           nrewnd =0,-1 only those tapes which are not in a rewound
c                        state and are not needed in processing
c                        observations of body nplnt0 are rewound.
c                        (this is normal situation)
c           nrewnd = 1   all tapes which are not in a rewound state are
c                        rewound.
c
c-----------------------------------------------------------------------
c
c            observation library tapes
c
c list of observations to be processed are on up to 10 observation
c library magnetic tapes plus input cards, with new library tapes
c written with input card data inserted at appropriate points. besides
c information on cards, library tapes have observed minus theoretical
c values of observations computed in run which produced tapes, data
c concerned with observations which stay essentially constant from run
c to run (e.g., precession-nutation matrix), and partial derivatives of
c observations. partial derivatives on newly written observation library
c tapes include all that were on input library tapes plus others needed
c in forming normal equations for given run. partial derivatives on
c observation library tapes whose values do not change signifigantly
c between iterations are not recalculated, whereas those whose values
c might change signifigantly (e.g., those with respect to planetary
c shape parameters) are recalculated before incrementing normal
c equations. there is option (ict(80)=1 with iterat=1) of using observed
c minus theoretical values on library tapes in forming normal equations
c rather than computing the values. data set iobcon contains data
c conrolling (1) alterations in error weightings for use in the normal
c equations of whole observation series or of individual observations
c within an observation series, (2) deletion of specific observations
c or whole observation series from the newly written library tapes,
c (3) overriding constants concerned with observation series in addition
c to error weighting constants, and (4) calculation of theoretical value
c of dummy observations and calculation of covariance matrix for dummy
c observations with or without additional real observations.
c
c iobs = observation card data set. first card contains information
c        about observation series, perhaps followed by additional such
c        cards. then come cards giving observation data for the
c        observation series, one card per observation. end of observa-
c        tion series is signaled by a blank card. then comes next obser-
c        vation series, etc. end of all the observation series (and end
c        of data set iobs) is signaled by two blank cards, one termina-
c        ting last observation series and an extra one terminating all
c        the observation series.
c iobs = 0, no observation card data set. in addition iobs is not read
c        if ict(80).gt.0. if ict(80).eq.0 and iobs.gt.0, then iobs is
c        read on first least squares iteration (iterat=1), but iobs
c        will be set equal to zero when
c        end is reached on iteration iterat=1 if ict(80).eq.0 so that
c        iobs will not be read on subsequent iterations, since the data
c        on iobs will have been written on the observation library tapes
c        to be read on subsequent iterations. ict(80).lt.0 is signal
c        that there are no input observation library tapes, program is
c        to read observation data only from iobs, rewinding iobs between
c        iterations, in which case partial derivatives have to be
c        recalculated each iteration.
c
c iabs1= input observation library data set. if 0, no such input.
c iabs2= output observation library data set. if 0, no such output.
c there is a loop on jtape from 1 to 10 to define values of iabs1,iabs2.
c if iterat=1, then      iabs1=iobs0(jtape)
c                        iabs2=iobs1(jtape)
c if iterat.gt.1 is even,iabs1=iobs1(jtape)
c                        iabs2=iobs2(jtape)
c if iterat.gt.1 is odd, iabs1=iobs2(jtape)
c                        iabs2=iobs1(jtape)
c jtape is incremented either (1) if end has been reached on iabs1 with
c iabs1.gt.0 and either iobs=0 or ntapa(2).gt.ntapa(3), or (2) if
c iabs1=0 and ntapa(2) makes a jump to a larger number. if ict(80).ge.0,
c for the first least squares analysis iteration, input library tapes
c iobs0 are read(along with cards iobs if ict(80).eq.0) and iobs1 is
c written. for the second
c iteration, iobs1 (written in first iteration) is read and iobs2 is
c written. in the third iteration iobs2 is read and iobs1 is written.
c and so it continues. the original input library tapes iobs0 are
c unchanged, but the input-output library tapes iobs1,iobs2 are
c alternately read and written, the program eating its tail from
c iteration to iteration.
c if ict(80).lt.0, then only observation cards iobs are read and
c iabs2=iobs1(1) is written from iteration to iteration.
c if ict(80).ge.0 saved partial derivatives are used from input
c observation library tapes where possible.
c if ict(80).lt.0 partial derivatives are recalculated each iteration.
c
      include 'crdbuf.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      integer*4 i2bod
      equivalence (i2bod, Jpert)
      include 'namtim.inc'
      integer*2 np3,np10
      equivalence (np3,Nplnt(-3)),(np10,Nplnt(-2))
      include 'number.inc'
      include 'obsdta.inc'
      include 'obstap.inc'
c in the following, i=1 denotes quantity belonging to iobcon
c                   i=2 denotes quantity belonging to iobs
c                   i=3 denotes quantity belonging to iabs1
c
c ntapa(i) = positive tape sequence number for observation series,
c            execpt that for iobcon it can be zero or negative to denote
c            dummy observations. for i = 1 or 2, ntapa(i) is greater
c            than or equal to tape sequence number for observation
c            series on iobcon or iobs, respectively, immediately
c            preceeding the given one if such exists. ntapa(3) for iabs1
c            corresponding to value of jtape greater than jtape
c            belonging to some previous iabs1 is greater than ntapa(3)
c            corresponding to this latter iabs1 if neither of the iabs1
c            in question are zero.
c
c nseqa(i) = positive observation series sequence number, greater than
c            previous observation series sequence numbers on the same
c            data set.
c
c ncoda(i) = code indicating type of observation series
c              1 radar observations of time delay and doppler shift
c              4 meridian circle obs of right ascension and declination
c                  referred to true equinox and equator of date reduced
c                  to center of earth, unless overriden.
c              5 photographic obs of right ascension and declination,
c                  topocentric, referred to mean equinox and equator of
c                  ref. epoch, circular aberration removed (unless
c                  overidden)
c
c nplnta(i)= planet number of observed body
c             0 sun       1 mercury   2 venus
c             4 mars      5 jupiter   6 saturn
c             7 uranus    8 neptune   9 pluto
c            10 moon
c            11,...,30 other natural planets,asteroids,satellites
c            31,... artificial space probes
c
c sita1(i) = first four characters of 8 character receiving site name
c sita1(4) = last four characters of this name for iabs1 (i=4)
c
c sera(i)  = four characters giving observation series name
c
c sita2(i) = first four characters of 8 character sending site name
c sita2(4) = last four characters of this name for iabs1 (i=3)
c
c spota(i) = four characters giving name of spot actually observed on
c            given body, if not center of body or sub-radar point
c
c erwgta(j,i)=factor by which quoted error of measurement j is
c             multiplied for use in the normal equations, j=1,2
c                j=1 radar measurement of time delay or optical
c                     measurement of right ascension
c                j=2  radar measurement of doppler shift or
c                     optical measurement of declination
c
c acctma(i)= accuracy constant for observation series
c            radar, accuracy constant for delay iteration
c            meridian circle, acc.constant for meridian crossing iterat.
c
c itima(i) = time code for observation series
c          = 10*c + f  (for c .ge. 0)
c          = c         (for c .lt. 0)
c         For radar observations, c is forced to be non-negative, so
c         f may be -1 in this case (means same as f=2).
c         Only last two digits of year are given for observations, so
c         an indication must be given at start of series what first
c         two digits are.  A series must not cross the start of a
c         century year, except for dummy observations.
c                 c is the century code, i.e., YEAR/100 - 19
c              1 21st century observations
c              0 20th century observations
c             -1 19th century observations
c             -2 18th century observations,etc.
c                 f is a flag for the type of observable
c                 (treated as 0 for c .lt. 0) 
c              0 optical: before atomic time
c              0 radar: observation time tag is UT2
c              1 atomic time exists, offset to be computed
c                observation time tag is UTC
c              2 AT - UT offset suppied with data
c              3 radar: observation time tag is UTC send time
c              5 normal point pseudo-obs at coordinate time
c
c fdeva(i) = fractional offset in units of 10**10 from a.1 time of unit
c            of time used in time delay measurement for radar obs series
c
c freqa(i) = frequency of sending radar for radar observation series
c
c nrewna(i)= 0,negative, at start of observation series only those tapes
c            not used in the series are rewound
c nrewna(i)= positive,  at start of observation series all tapes are
c            rewound
c nrewna(i) less than 0 or greater than 1 indicates in addition to the
c            above options that additional control data are to be read
c            from iobs or iobcon or exists on first record of series on
c            iabs1.
c
c ctlga(i) = 8-character name of reference star catalog used in
c            reducing the observation series
c
c niobc,niobs,niabs1 indicate the state of data sets iobcon,iobs,iabs1
c   -1 data set has reached end
c    0 data set is positioned at first record of observation series
c    1 data set has read first record of observation series and is
c      positioned to read subsequent records of series
c life=0 nothing yet written on iabs2
c life=1 first 2 records of iabs2 have been written
c
c ntaps(i)= old value of observation tape sequence number, i=1,4
c nseqs(i)= old value of observation series sequence number, i=1,4
c    i=1 from iobcon data set
c    i=2 from iobs data set
c    i=3 from iabs1 data set
c    i=4 for output library tape iabs2
c
c-----------------------------------------------------------------------
c
      include 'plndta.inc'
      include 'rotdta.inc'
      include 'sbdta.inc'
      include 'scdta.inc'
      include 'sitcrd.inc'
      include 'stats.inc'
      include 'tapdta.inc'
      include 'zeroes.inc'

c external functions
      integer*4 ITYPOB,JBDTST
c
c local
      integer   i,ict31,ictsav,j,jtape,
     .          jtape1,jtypob,jut1,jwob
      integer*2 ione/1/
c
c setup computations
      if(Iterat.gt.1 .or. mocpar.eq.1) Ict(41) = 0
      ictsav = Ict(1)
      if(mocpar.le.0) then
         if(Ict(80).gt.0) Iobs = 0
      else
         Ict(1) = 0
         Iobs   = 0
      endif
      call COMSET
      Iabs2 = Iobs1(1)
c
c*  start=100
c read first records of moon, embary, n-body data sets
      if(Nbody.gt.0) then
         call BDRD1(Libdy)
         Libdy = 1
      endif
      if(i2bod.gt.0) then
         call B2RD1(Lib2y)
         Lib2y = 1
      endif
      if(Imn.gt.0 .or. JBDTST(np10).ne.0) then
         call MNRD1(Limn)
         Limn = 1
      endif
      if(Iem.gt.0 .or. JBDTST(np3).ne.0) then
         call EMRD1(Liem)
         Liem = 1
      endif
      if(Inut.gt.0) then
         call RTRD1(Linut, 3)
         Linut = 1
      endif
      if(Ilib.gt.0) then
         call RTRD1(Lilib, 10)
         Lilib = 1
      endif
 
c ut1, wobble data sets
      if(Jct(33).gt.0) then
         call UT1RD1
         call WOBRD1
      endif
 
c ct-at data set
      if(Ictat.gt.0) call CTRD1(0)
c
c initialize control integers for observations
c*  start=1000
      if(Iobs.gt.0 .and. Ieof.eq.1) Iobs = 0
      Life  = -1
      Niobc = 0
      Niobs = 0
      if(Iobs.le.0) Niobs = -1
      Iabs2 = 0
      jtape = 0
c
c see if there is to be checkpoint restart
      if(mocpar.gt.0) then
c
c*  start=1100
c setup reading,writing of observation library tape
         jtape = jtape + 1
      else if(Ict(31).eq.0) then
         jtape = jtape + 1
      else
         jtape = Ict(31)
         jtape = iabs(jtape)
         if(Ict(31).ge.0) Life = -2
      endif
      do while( jtape.le.Numobt )
         jtape1 = jtape
         Iabs1s = 0
         Jiabs1 = 9999
         Iabs2  = Iobs1(jtape)
         Jiabs2 = 1
         if(Ict(80).ge.0 .or. mocpar.gt.0) then
            Iabs1s = Iobs0(jtape)
            Jiabs1 = 0
            if(Iterat.gt.1) then
               Iabs1s = Iobs2(jtape)
               Jiabs1 = 2
               if(Iterat .eq. (2*(Iterat/2))) then
                  Iabs1s = Iobs1(jtape)
                  Jiabs1 = 1
                  Iabs2  = Iobs2(jtape)
                  Jiabs2 = 2
                  if(Ict(80).gt.0 .and. Iterat.eq.2) then
                     Iabs1s = Iobs0(jtape)
                     Jiabs1 = 0
                  endif
               endif
            endif
c
c read and print first two records of observation library tape
            if(mocpar.gt.0) then
               Iabs2  = 0
               Iabs1s = 0
            endif
            Iiabs1 = Iabs1s
            call OBSRD1(Ntaps(3), Jiabs1, jtape, mocpar)
c
c position tapes for checkpoint restart
            if(Ict(31).ne.0) then
               ict31 = Ict(31)
               ict31 = iabs(ict31)
               if(jtape.eq.ict31) then
                  Iabs1 = Iabs1s
                  call CMPRST
               endif
            endif
         endif
         Niabs1 = 0
         do while( .true. )
c
c*  start=1200
            call CMPAR1(jtape1, mocpar)
            if(jtape1.le.0) then
c
c rewind output observation library tape
               if(Iabs2.gt.0) then
                  write(Iabs2) (Izero(1),i=1,110),(Zero(1),i=1,81),
     .             izero,izero,izero,izero,izero,izero,izero,izero
c end file iabs2
                  rewind Iabs2
                  Life = -1
               endif
c
c
               Ncodf = 0
               if(Ict(80).lt.0 .and. mocpar.le.0) go to 200
               jtape = jtape + 1
               go to 100
            else
               if(Ncodf.le.0) go to 200
 
               call CMPAR2(mocpar)
               call CMPAR3(mocpar, jtape)
c
c dispatch to observation routines
               jtypob = ITYPOB(Ncodf)
 
               if(jtypob.eq.1) then
c
c radar-type observables
                  if(Ksite(1).gt.0 .and. Ksite(1).ne.3) then
                     call STRADR(mocpar)
                  else
                     call RADAR(mocpar)
                  endif
c
c optical-type observables
               else if(jtypob.eq.2) then
                  if(Ksite(1).gt.0 .and. Ksite(1).ne.3) then
                     call STOPTC(mocpar)
                  else
                     call OPTIC(mocpar)
                  endif
c
c occultation-type observables
               else if(jtypob.eq.3) then
                  call TRNSIT(mocpar)
c
c interferometry-type observables
               else if(jtypob.eq.4) then
                  call FERMTR(mocpar)
c
c two-spacecraft delay observables
               else if(jtypob.eq.5) then
                  call STMRDR(mocpar)
 
               else
                  call SUICID(
     .           'UNKNOWN OBSERVATION TYPE, STOP IN COMPAR',10)
               endif
            endif
         end do
  100 end do
c
c end of comparison of theory and observation
c restore tapes
  200 call ERRTOT(Ncodg)
      if(Iobs.gt.0 .and. Iobs.ne.5) then
         rewind Iobs
         Itrwnd(Iobs) = 0
      endif
      if(i2bod.gt.0) then
         rewind i2bod
         Itrwnd(i2bod) = 0
      endif
      if(Nbody.gt.0) then
         rewind Ibody
         Itrwnd(Ibody) = 0
         if(Kpert.gt.0) then
            rewind Kpert
            Itrwnd(Kpert) = 0
         endif
      endif
      if(Inut.gt.0) then
         rewind Inut
         Itrwnd(Inut) = 0
      endif
      if(Imn.gt.0 .and. Itrwnd(Imn).gt.0) then
         rewind Imn
         Itrwnd(Imn) = 0
      endif
      if(Ilib.gt.0 .and. Itrwnd(Ilib).gt.0) then
         rewind Ilib
         Itrwnd(Ilib) = 0
      endif
      if(Iem.gt.0 .and. Itrwnd(Iem).gt.0) then
         rewind Iem
         Itrwnd(Iem) = 0
      endif
      if(Klan.gt.0 .and. Klan.le.u_mxpl) then
         if(Jplnt.gt.0) then
            rewind Jplnt
            Itrwnd(Jplnt) = 0
         endif
         if(Klanr.gt.0 .and. Jpr.gt.0) then
            rewind Jpr
            Itrwnd(Jpr) = 0
         endif
      endif
      if(Klanb.gt.0) then
         if(JBDTST(Nplnt(Klanb)).eq.0) then
            rewind Jsb
            Itrwnd(Jsb) = 0
         endif
      endif
      if(Klans1.gt.0) then
         if(JBDTST(Nplnt(Klans1)).eq.0) then
            rewind Jsc
            Itrwnd(Jsc) = 0
         endif
      endif
      if(Klssb.gt.0) then
         do i=1,9
            if(Ssbkl(i).gt.0) then
               j=Iplss(i)
               rewind j
               Itrwnd(j)=0
            endif
         end do
      endif
      if(Jct(33).gt.0) then
         jut1 = Jct(33)
         jwob = Jct(33) + 1
         rewind jut1
         Itrwnd(jut1) = 0
         rewind jwob
         Itrwnd(jwob) = 0
      endif
      if(Libhoc.gt.0) then
         rewind Libhoc
         Itrwnd(Libhoc)=0
      endif
      if(Jct(49).gt.0) call FLURD9
      if(Ictat.gt.0) then
         rewind Ictat
         Itrwnd(Ictat) = 0
      endif
      if(Jct(45).ne.0) call ANGOT(2)
 
c*  start=3000
      call ERRTOT(izr2)
      Ict(1) = ictsav
 
c test for unused iobs input
      if(Iobs.gt.0 .and. Niobs.ge.0) write(Iout, 300) Iobs
  300 format('0* * * WARNING. UNUSED OBSERVATION CARDS REMAIN ON IOBS=',
     .       i3, ' * * *')
c
      return
      end
