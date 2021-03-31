      subroutine PRMRED(in0,sepeat,moprnt,iseq1,nstop,jdpad)

      implicit none

c arguments
      integer*4 in0,moprnt,iseq1,nstop,jdpad
      integer*2 sepeat

c array dimensions
      include 'globdefs.inc'

c
c-----------data sets (see perdta for default values)------
c
      include 'inodta.inc'
c
c IN    =input data and control constants
c IOUT  =print output
c JOUT  =additional print output for analiz link to get extra printout
c        of correlation matrix and parameter adjustments (if zero, no
c        such printout)
c KOUT  =additional print output for extended print output in subroutine
c        SBOUT for SBFN integration of probe or planet (if zero, no
c        such printout)
c LOUT  =enabling switch for extra input link print of adjustable
c        parameter summary and a priori information (if zero, no such)
c MOUT  =extra "print" data set for summary of program flow, essentially
c        a digest of important data on iout, designed for speedy
c        on-line examination via time-sharing terminals. (if zero, no
c        such output)
c NOUT  =extra "print" data set for list of observations (with errors
c        and residuals) processed in compar and prdict links.  designed
c        for on-line processing and/or plotting via "smart" video
c        terminals. (if zero, no such output)
c IPUNCH=punch output
c JPUNCH=extra punched output -- new PEP input deck after each
c        iteration other than the last iteration
c KPUNCH
c IGRAPH=graphical output
c INTERN=internal buffer
c IBUF  =general purpose buffer
c IMAT  =saved inverse of coefficient matrix of normal equations and
c        solution vector calculated in SNORML (output)
c        (IMAT=0 don't save, IMAT.GT.0 do save)
c        used for series-by-series adjust
c JMAT  =same as IMAT but input to SOLFRM if ICT(5).gt.1
c IPERT =pertubing planets and Moon positions
c JPERT =additional perturbing body positions (asteroids or satellites)
c KPERT =positions of asteroid system c-of-m
c IPERT2=moon tape used as source of moon for nbody integration,
c        superseding IPERT
c IMAT0 = input saved normal equations and right hand side
c       = JMAT0+ 100 * JMATA  with JMAT0 used as data set number
c         and with JMATA a flag for type of restoration of
c         normal equations - 0 normal restore, 1 skip rhs, 2 skip lhs
c         data sets IMAT0(j),j=1,NUMMT0 are used if NUMMT0 .gt. 0
c         .and. NUMMT0 .le. maxmt0
c IMAT1 = output saved normal equations and right hand side
c         (series by series)
c        used for series-by-series adjust   (temp)
c IMAT2 = output total saved normal equations and right hand side
c         (must be assigned if restoring normal equations with imbedded
c          a priori constraint matrices, unless JCT(60)>0)
c IMAT3 = output total saved partially pre-reduced normal eqns
c         also used as input for de-correlation of partials. see ICT(11)
c IMAT4   unassigned
c LIBHOC= input data set of lunar libration corrections
c ICTAT = input integrated CT-AT data set
c IENG  =data set for gas leak and other man-made-force data
c JENG  =data set for radiation pressure data
c KENG
c
      include 'plndta.inc'
c
c IEM  =Earth-Moon barycenter coordinates,partials
c IMN   =Moon coordinates,partials
c INUT  =Earth nutation
c ILIB  =Moon physical libration, partials
c IPLNT(i)=planet NPLNT(i) coordinates,partials (i=1,numpln)
c IPLCON=Earth-Moon barcenter,Moon,Earth rotation,Moon rotation,planet
c        input data not in labeled commons
c IVCNT =velocity of center of mass of solar system relative to sun
c       (INUT,IVCNT not used.  Nutation & libration come from Moon or
c        n-body tape normally.)
c
      include 'obsdta.inc'
c
c
c IOBS  =observation cards (in sysin stream following input data and
c        control constants if IOBS=5, which is default value. if
c        ICT(1)=-1 or IOBS input 0, no observation cards read by
c        program. if IOBS input as 8 or some other number not 0 or 5,
c        then observation cards on tape)
c
c Format for first card of observation series for observation card data
c set IOBS.
c
c cols  1-2  Integer (i2 format) giving type of observation series.
c
c                 1 radar observations of time delay and/or Doppler
c                   shift
c                 2 radio tracking observations of time delay and/or
c                   Doppler shift
c                 3 differential radar observation
c                   (feature relative to subradar point if spot name
c                   not blank, bandwidth if spot name is blank)
c            (If observing sites on Earth, observables 1,2,3 processed
c            in radar link.  If observing sites on non-Earth body,
c            observables 1,2,3 processed in STRADR link.)
c
c            For observations made from the Earth
c            (processed in OPTIC link)
c                 4 meridian circle observations of right ascension and
c                   declination referred to the true equinox and equator
c                   of date for observed body not Earth,Moon satellite
c                   (geocentric)
c                 4 azimuth and elevation observations of Earth
c                   or Moon satellites (topocentric)
c                 5 photographic observations of right ascension and
c                   declination referred to the mean equinox and equator
c                   of reference epoch (elliptic aberration not removed:
c                   astrometric) (topocentric)
c                 6 same as 5 but with all aberration removed
c                   (astrographic)
c            For observations made from a non-Earth body
c            (processed in STOPTC link)
c                 4 look angles in spacecraft coordinate system
c               5,6 angles relative to star backgound (astometric,
c                   astrographic -- but abberation not quite right?)
c
c                 7 occultation
c                 8 transit of planet across face of sun (limbs)
c                 9 transit of planet across face of sun (black dot)
c
c                10 interferometer phase observation
c                11 conventional interferometer observation
c                12 counted cycle observation
c
c                13 time delay measurement for signal sent from ground
c                   to a spacecraft to another spacecraft and then
c                   back to the ground
c                14 time delay measurement for signal sent from ground
c                   to a spacecraft to another spacecraft back to the
c                   first spacecraft and then back to the ground
c
c                18 pulsar phase observation
c
c                19 les-8/9 one way Doppler count
c
c                20 azimuth-elevation (for sun,Moon,planets for which
c                   ncodf=4 means meridian circle)
c
c                21,22,23,... two observed objects difference observable
c                             of type ncodf-20
c
c cols  3-5  Integer (i3 format) giving planet number of observed body
c                 0 Sun        6 Saturn
c                 1 Mercury    7 Uranus
c                 2 Venus      8 Neptune
c                 4 Mars       9 Pluto
c                 5 Jupiter   10 Moon
c                11,...,30 natural asteroids and satellites
c                31,...    artificial space probes
c                negative  stars
c
c col   6    Flag for extra header card(s) if using 72 column format
c
c cols  7-10 First four characters (a4 format) of eight character
c            receiving site name (ignored for transit observations)
c
c cols 12-15 Four characters (a4 format) giving observation series name
c
c cols 17-20 First four characters (a4 format) of eight character
c            sending site name (ignored for optical observations)
c            for satellite mutual events, distinguishes between eclipse
c            and occultation.
c
c cols 22-25 Four characters (a4 format) giving name spot observed on
c            given body (if blank, subradar point for rader observations
c            or center of body for optical observations). However, if
c            planet number of observed body is negative, then this is
c            four character star name.
c            For radar observations of a planet, spot name = '&&&&'
c            denotes observations away from subradar point with spot
c            longitude and latitude on observation card for each obs.
c            For NCODF=13,14 spot name is first four characters of
c            second spacecraft name.
c            Also, for ncodf=7-9 spot name is first four characters
c            of second body name.
c
c cols 26-31 single precision number (e6.0 format) by which error of
c            first measurement (time delay or right ascension) of every
c            observation in series is to be multiplied for use in
c            forming normal equations.
c
c cols 32-37 Single-precision number (e6.0 format) by which error of
c            second measurement (Doppler shift or declination) of every
c            observation in series is to be multiplied for use in
c            forming normal equations.
c
c cols 38-43 Single-precision number (e6.0 format) giving accuracy
c            criterion for iterations in processing observations in
c            series (iterations to determine reflection and send times
c            given receive time for radar observations, iteration to get
c            body over meridian for meridian cicle observations, etc.).
c            The units and meaning depend on the observable type.  For
c            most types, the value represents a change of distance in
c            light-seconds, but it gives a change in the time of transit
c            or occultation (in seconds) or a change in Earth rotation
c            phase (in radians) for meridian-circle observations.  If
c            the value is omitted, 0.001 is used as a default.
c
c cols 44-45 Integer (i2 format) giving type of time epoch for
c            observations in series.
c            Radar observations:    1 WWV (UTC) epoch
c                                   0 UT2 epoch
c                                  -1 (A.1-epoch) and (UT1-epoch) are
c                                     on observation card
c                                   2 same as -1
c                                   3 epoch is send time
c                                   5 epoch is CT, normal point
c                                     pseudo-observations.
c            Optical observations:  1 epoch in 20th century for which
c                                     atomic time exists
c                                   0 epoch in 20th century for which
c                                     atomic time does not exist
c                                  -1 epoch in 1800's, etc.
c                                     (needed because only last two
c                                      digits of year are on card)
c            All observations: 10*i+j First two digits of the year are
c                                     19+i (for i.ge.0), j interpreted
c                                     as above
c
c From this point on, 2 formats are possible, one using full 80 columns,
c the other stopping at column 72.  For IOBS cards, the format selection
c is based on JCT(69); for dummy obs and error overrides, the format
c selection can be locally changed by insertion of '#72' or '#80'
c cards before any series.
c
c cols 46-52 Single-precision number (f7.1) format) giving offset*1e10
c            from A.1 atomic time of unit of time delay measurement in
c            radar observation series.
c            In the 72 column format this field appears on an optional
c            2nd header card (present if column 6 contains ':').
c
c cols 53-70 Double-precision number (d18.11 format) giving frequency of
c  or  46-63 sending site for radar observation series.
c
c cols 71-72 Integer (i2 format) indicating if ephemeris tapes (usually
c  or  64-65 backward in time) are to be rewound before processing
c            observations in series.
c                                        -1,0 tapes not to be rewound
c                                         1,2 tapes to be rewound
c             if  0,1 next card after header card(s) of series is
c                     observation data
c             if -1,2 next card after header card(s) of series is
c                     site card(s) then come observation data
c
c cols 73-75 Integer (i3 format) giving tape sequence number for output
c  or  66-68 observed minus theory and partial derivative tape for
c            observation series. (Must be positive and non-decreasing
c            from one observation series to the next.)
c
c cols 76-80 Integer (i5 format) giving series sequence number (must be
c  or  69-72 positive and increasing from one observation series to the
c            next for given value of tape sequence number).
c            In the 72 column format this field is i4.
c
c Format for optional second header card ('colon' card).
c
c col   1    Flag for optional further header cards.
c cols  2- 8 Same as cols 46-52 of 80 column format above.
c cols  9-16 Eight-character name (a8 format) of star catalog used in
c            reducing the observations in this series.
c
c
c Format for radar observation cards for observation data set IOBS.
c
c col   1    Integer (i1 format) giving type of observation
c                 1 radar observation of time delay
c                 2 radar observation of time delay and Doppler shift
c                 3 radar observation of Doppler shift
c                 3 bandwidth if ncodf=3, spot name blank
c cols  2-4  Integer (i3 format) giving hour of receive time epoch
c cols  5-7  Integer (i3 format) giving minute of receive time epoch
c cols  5-15 Single-precision number (f8.4 format) giving second of
c            receive time epoch.
c cols 16-28 Double-precision number (f13.7 format) giving observed time
c            delay in seconds.
c cols 29-35 Single-precision number (e7.0 format) giving error of time
c            delay measurement in seconds.
c cols 36-50 Double-precision number (d15.8 format) giving observed
c            Doppler shift in cycles per second.
c cols 51-57 Single-precision number (e7.0 format) giving error of
c            Doppler shift measurement in cycles per second.
c cols 58-65 Single-precision number (f8.4 format) giving A.1 time minus
c            epoch (if needed).
c cols 66-72 Single-precision number (f7.4 format) giving UT1 time minus
c            epoch (if needed).
c cols 73-74 Integer (i2 format) giving month observation (from 1 to 12)
c cols 75-76 Integer (i2 format) giving day of observation (from 1 to
c            31)
c cols 77-78 Integer (i2 format) giving last two digits of year of
c            observation (see also cols 44-45 on the series header).
c col  79    Character (a1 format) indicating observing site (not read,
c            for labeling purposes only).
c col  80    Character (a1 format) indicating observed body (not read,
c            for labeling purposes only).
c
c
c Format of optical observation cards for observation data set IOBS.
c
c col   1    integer (i1 format) giving type of observation
c              1 or 4 observation of right ascension
c              2 or 5 observation of right ascension and declination
c              3 or 6 observation of declination
c cols  2-15 give the approximate universal time of transit of the
c            observed body over the meridian of the observatory (within
c            twelve hours of the actual time of transit).  the center is
c            observed unless first card of series indicates the limb is
c            observed.
c cols  2-4  integer (i3 format) giving approximate hour of transit
c cols  5-7  integer (i3 format) giving approximate minute of transit
c cols  8-15 Single-precision number (f8.4 format) giving approximate
c            second of transit
c col  16    blank (not used)
c cols 17-30 give the right ascension of center from observation unless
c            first card of series indicates the limb is observed
c cols 17-19 integer (i3 format) giving hour of right ascension
c cols 20-22 integer (i3 format) giving minute of right ascension
c cols 23-30 Single-precision number (f8.4 format) giving second of
c            right ascension
c cols 31-37 Single-precision number (f7.4 format) giving error of
c            right ascension in seconds of time
c cols 38-50 give the geocentric declination of the center from
c            observation unless first card of series indicates the limb
c            is observed
c cols 38-40 integer (i3 format) giving degree of declination
c cols 41-43 integer (i3 format) giving minute of declination
c cols 44-50 Single-precision number (f7.3 format) giving second of
c            declination
c cols 51-57 Single-precision number (f7.3 format) giving error of
c            declination in seconds of arc
c cols 58-66 blank (not used)
c col  67    character (a1 format) indicating clamp
c cols 68-69 characters (2a1 format) indicating limb of right ascension
c            and declination (used in computing o-c if published
c            observation not reduced to center)
c cols 70-71 characters (2a1 format) indicating observer
c col  72    blank (not used)
c cols 73-74 integer (i2 format) giving month of observation
c            (from 1 to 12)
c cols 75-76 integer (i2 format) giving day of observation
c            (from 1 to 31)
c cols 77-78 integer (i2 format) giving last two digits of year of
c            observation (the century of observation is given on first
c            card of observation series)
c col  79    character (a1 format) indicating observing site (not read,
c            for labeling purposes only)
c col  80    character (a1 format) indicating observed body (not read,
c            for labeling purposes only)
c
c
c transit and occultation observation data cards have their own format
c interferometer and pulsar observation data cards have thier own format
c
c
c IOBCON=data for generation of dummy observations and deletion and
c        error weight alteration of real observations (written in
c        subroutine dltred of input link)
c numobt=number of observation library tapes (must be between 1 and 10
c        inclusive). then for i=1,numobt, we have
c IOBS0(i)   = input observation library tapes for first least squares
c              iteration
c IOBS1(i)   = output observation library tapes on odd least squares
c              iterations, input on even iterations
c IOBS2(1-10) =output observation library tapes on even least squares
c              iterations, input on odd iterations greater than first
c
c Observations to be processed are on up to 10 observation library files
c plus input cards.  New libraries are written with input card data
c inserted at appropriate points.  Besides information on cards, library
c files have observed minus theoretical values of observations computed
c in run which produced tapes, data concerned with observations which
c stay essentially constant from run to run (e.g., precession-nutation
c matrix), and partial derivatives of observations.  Partial derivatives
c on newly written observation library tapes include all that were on
c input library tapes plus others needed in forming normal equations for
c given run.  Partial derivatives on observation library tapes whose
c values do not change signifigantly between iterations are not
c recalculated, whereas those whose values might change signifigantly
c (e.g., those with respect to planetary shape parameters) are
c recalculated before incrementing normal equations.  There is option
c (ICT(80)=1 with iterat=1) of using observed minus theoretical values
c on library tapes in forming normal equations rather than computing the
c values.  Data set IOBCON contains data conrolling (1) alterations in
c error weightings for use in the normal equations of whole observation
c series or of individual observations within an observation series, (2)
c deletion of specific observations or whole observation series from the
c newly written library tapes, (3) overriding constants concerned with
c observation series in addition to error weighting constants, and (4)
c calculation of theoretical value of dummy observations and calculation
c of covariance matrix for dummy observations with or without additional
c real observations.
c
c IOBS = observation card data set. first card contains information
c        about observation series, perhaps followed by additional such
c        cards. then come cards giving observation data for the
c        observation series, one card per observation. end of observa-
c        tion series is signaled by a blank card. then comes next obser-
c        vation series, etc. end of all the observation series (and end
c        of data set IOBS) is signaled by two blank cards, one termina-
c        ting last observation series and an extra one terminating all
c        the observation series.
c IOBS = 0, no observation card data set. in addition IOBS is not read
c        if ICT(80).gt.0. if ICT(80).eq.0 and IOBS.gt.0, then IOBS is
c        read on first least squares iteration (iterat=1), but IOBS
c        will be set equal to zero when
c        end is reached on iteration iterat=1 if ICT(80).eq.0 so that
c        IOBS will not be read on subsequent iterations, since the data
c        on IOBS will have been written on the observation library tapes
c        to be read on subsequent iterations. ICT(80).lt.0 is signal
c        that there are no input observation library tapes, program is
c        to read observation data only from IOBS, rewinding IOBS between
c        iterations, in which case partial derivatives have to be
c        recalculated each iteration.
c
c in subroutine compar
c iabs1= input observation library data set. if 0, no such input.
c iabs2= output observation library data set. if 0, no such output.
c there is a loop on jtape from 1 to 10 to define values of iabs1,iabs2.
c if iterat=1, then      iabs1=IOBS0(jtape)
c                        iabs2=IOBS1(jtape)
c if iterat.gt.1 is even,iabs1=IOBS1(jtape)
c                        iabs2=IOBS2(jtape)
c if iterat.gt.1 is odd, iabs1=IOBS2(jtape)
c                        iabs2=IOBS1(jtape)
c jtape is incremented either (1) if end has been reached on iabs1 with
c iabs1.gt.0 and either IOBS=0 or ntapa(2).gt.ntapa(3), or (2) if
c iabs1=0 and ntapa(2) makes a jump to a larger number. if ICT(80).ge.0,
c for the first least squares analysis iteration, input library tapes
c IOBS0 are read(along with cards IOBS if ICT(80).eq.0) and IOBS1 is
c written. for the second
c iteration, IOBS1 (written in first iteration) is read and IOBS2 is
c written. in the third iteration IOBS2 is read and IOBS1 is written.
c and so it continues. the original input library tapes IOBS0 are
c unchanged, but the input-output library tapes IOBS1,IOBS2 are
c alternately read and written, the program eating its tail from
c iteration to iteration.
c if ICT(80).lt.0, then only observation cards IOBS are read and
c iabs2=IOBS1(1) is written from iteration to iteration.
c if ICT(80).ge.0 saved partial derivatives are used from input
c observation library tapes where possible.
c if ICT(80).lt.0 partial derivatives are recalculated each iteration.
c
c Global deletions are performed when Tdlt0, Tdlton, and Tdltof are
c specified.  Beginning at epoch Tdlt0 (julian ephemeris date) input
c data are alternately processed normally for Tdlton days and skipped
c for Tdltof days regardless of the other data editing controls in PEP.
c This periodic windowing extends both forward and backward in time
c from Tdlt0.  For input observation cards, the skipped data are given
c an overriding weight of zero, but such data from observation libraries
c are skipped altogether.
c
      include 'aprtbf.inc'
c
c epsa    no longer used
c ibuf1 = data set for parameter names.  if > 0, then names will appear
c         next to correlations in analiz link printout.
c ibuf2 = apriori b matrix
c ibuf3  saved right hand sides of normal equations for
c        series-by-series adjust
c ibuf4  partial solutions of series-by-series adjust
c ibuf5  multiple parameter set requests, pseudo-iptr's (iptrx's)
c ibuf6  parameter(and associated) values (by run) for
c        multiple parameter set (or multiple epoch) runs
c        also used for pprne generation (see arrwrt)
c ibuf7  temporary data set for partial prereduction of normal eqns
c ibuf8  temporary data set for partial prereduction of normal eqns
c ibuf9  temporary data set for partial prereduction of normal eqns
c ibuf10 unassigned

      include 'bdctrl.inc'
      integer*4 zbdctr/304/   !r8=304,r10=304

      include 'dtparm.inc'
      integer*4 zdtprm/6010/   !r8=6010,r10=6010

      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'lcntrl.inc'
      integer*4 zlcntr/1400/   !r8=1400,r10=1400

      include 'nrmwgt.inc'
      include 'mascon.inc'
      integer*4 zmascn/6404/   !r8=3204,r10=6404
      include 'maxmt0dt.inc'
c WGTOBS - Array of weighting factors applied to obslib tapes (IOBS0,
c          IOBS1, or IOBS2) in forming normal equations (on top of any
c          error weight information).  The sense of weighting is such
c          that wgtobs>1 increases the weight of the data (as if all
c          measurement uncertainties were divided by 'wgtobs').  A
c          value of zero causes the corresponding obslib to be skipped.
c WGTMT0 - Array of weighting factors applied to input saved normal
c          equations (on top of any error weight information).
c          The sense of weighting is the same as for wgtobs above.
c          A value of 0 causes the corresponding imat0 to be skipped.
c          A negative value causes the data on the corresponding imat0
c          to be deducted from the normal equations and error
c          statistics, i.e., reduces the count of measurements and
c          sum-squared residuals and subtracts information from the
c          normal equations.  This is meaningful only if some superset
c          of the same data is being included with a positive weight.
      include 'param.inc'
      include 'prpgat.inc'
c     use of jcal :
c     jcal is an override vector set in the PEP input stream
c  jcal(i) = 00 is the default
c any digit of jcal that is = 0 is replaced by the corresponding
c digit of JCT(2)
c
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
c     tens digit =
c                 1 : do not use the ith correction
c                 2 : use standard logic for correction use
c                 3 : use the correction regardless of rank
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           read following variables from above labeled commons
c from /INODTA/ and /NRMWGT/
c from /PLNDTA/ and /OBSDTA/
c from /FCNTRL/, /LCNTRL/ and /APRTBF/
c from /PARAM/
c from /BDCTRL/
c from /DTPARM/
c from /MASCON/
c maximum number of planets allowed
c integer control varible for propagation corrections
      namelist /NMLST1/Iout, Jout, Kout, Ipunch, Igraph, Intern, Ibuf,
     .         Imat, Ipert, Jpert, Kpert, Imat1, Imat2, Imat0, Nummt0,
     .         Jmat, repeat, noprnt, iseq, Lout, Mout, Nout, Jpunch,
     .         Kpunch, Ieng, Jeng, Keng, Imat3, Imat4, Ictat, Wgtobs,
     .         Wgtmt0, Iplcon, Ivcnt, Iobs, Iobs0, Iobs1, Iobs2, Iobcon,
     .         Numobt, Eps, Ict, Npage, Lprm, Epsa, Ipert1, Ipert2,
     .         Jpert1, Jpert2, Kpert1, Kpert2, Ibuf6, Ibuf7, Ibuf8,
     .         Ibuf9, Ibuf10, Ibuf1, Ibuf2, Ibuf3, Ibuf4, Ibuf5, Extprc,
     .         Typout, Jct, prmter, Mass, Relfct, Gmvary, Sunhar,
     .         Sunpoe, Aultsc, Ltvary, Reldel, Reldop, Ctvary, Ecinc,
     .         Seqinc, Seqasc, Sunrad, Mdstau, Mdstsc, Ltvel, Ibody,
     .         Jdbdy1, Jdbdy0, Jdbdy2, Jvlbdy, Epsbdy, Intbdy, Kbdy,
     .         Nplbdy, Nbody, Dt, Jddt, Ldt, Numdt, Jddt0, Kkbdy,
     .         maxpln, Jcal, Tdlt0,Tdlton,Tdltof,Masmsc,Rkmmsc,
     .         Lngmsc,Latmsc,Nummsc,Nplmsc,Libhoc,jdpad,Gauss
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      include 'loadid.inc'
      include 'timstf.inc'

c local variables
      integer*4 i,in1,iout2,iseq,maxpln,noprnt,t1,t2
      integer*2 repeat
      character*8 specs

      Dat0   = Dat1
      Ireal0 = Ireal1
      Itotsk = 0
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c initialize fundamental data set numbers
      In   = 5
      Iout = 6
      in0  = 98
c in0 is for PEP spooling of input stream
c
c read heading describing program run
      read(In,100) Heding,specs
  100 format(18A4,a8)
      if(specs(1:4).eq.'    ') then
         t1=ICHAR(specs(5:5))-ICHAR('0')
         t2=ICHAR(specs(7:7))-ICHAR('0')
         if((specs(5:5).eq.' '.or.(t1.ge.0.and.t1.le.9)) .and.
     .    (specs(7:7).eq.' '.or.(t2.ge.0.and.t2.le.9))) then
            read(specs,110) iout2,in1
  110       format(4x,2i2)
            if(iout2.gt.0) Iout = iout2
            if(in1.gt.0) in0 = in1
         endif
      endif
      call PLOG(Heding)
      write(Iout, 200) Npage, Lnkdat, Lnkdsn, Heding
  200 format('1INPUT STREAM TO THE PLANETARY EPHEMERIS PROGRAM (PEP)',
     .       t122, 'PAGE', i5/'0PEP VERSION= ', A8, 1x, A44/'0',
     .       18A4, 9x, 'TITLE CARD')

c npage initialized in main
      Npage = Npage + 1
c
c initialize peripheral data set numbers
      call PERDTA(in0)
c
c           initialize remainder of fcntrl labeled common (except for
c           iterat which is set equal to zero in main program, nparam
c           which is calculated in subroutine check, and date which was
c           determined in main program)
      Itimst = Ireal1

c npage =1       (initialized in main program)
      do i = 1, 80
         Ict(i) = 0
      end do
      Ict(1)  = -1
      Ict(5)  = -1
      Ict(10) = -2
      Ict(26) = 1
      Ict(34) = 2
      Ict(46) = -1
      do i = 1, 100
         Jct(i) = 0
      end do
      Jct(2)  = 11
      Jct(61) = -1
      Jct(62) = -1
      Eps(1)  = 1.0E3
      Eps(2)  = 1.0E2
      Eps(3)  = 1.0E1
      Eps(4)  = 1.0E1
      Eps(5)  = 1.0E1
      Eps(6)  = 1.0E1
      Eps(7)  = 1.0E2
      Eps(8)  = 1.0E1
      Eps(9)  = 0.1E0
      Eps(10) = 0.01E0
      Eps(11) = 1.0E-3
      Eps(12) = 1.0E-6
      Eps(13) = 1.0E-32
      do i = 14, 30
         Eps(i) = 0.0E0
      end do
      Eps(14) = 1._10
      Ntptop  = 99
      do i = 1, Ntptop
         Itrwnd(i) = 0
      end do
      Nsoltn = 0
c
c HEDING  =Input page heading describing program run.
c DATE    =Calendar date determined from operating system.
c ITIMST  =Real time at start of program run (binary integer number of
c          hundreths of seconds, determined from operating system).
c NPARAM  =Number of parameters to be adjusted in least squares analysis
c          (calculated in subroutine CHECK and again in ANSET).
c ITERAT  =Iteration number (initialized and incremented in main prog.)
c NPAGE   =Page number (initialized in subroutine input, incremented
c          after each new page heading is written).
c
c N.B. EPS in /FCNTRL/ read in &NMLST1 is not the same as EPS in &NMLST2
c
c Observed minus theory for measurement j deleted in normal equations
c if abs(obs-th).ge.EPS(j)*error
c j=1    time delay
c j=2    Doppler
c j=3    right ascension
c j=4    declination
c j=5    first limb of body crossing for transit or occultation
c j=6   second limb of body crossing for transit or occultation
c j=7    interferometric differential delay
c j=8    interferometric differential delay rate
c
c EPS(9) = Partial convergence criterion for adjustments to parameters.
c EPS(10)= Full convergence criterion for adjustments to parameters
c          least squares iteration terminates if either ITERAT.ge.
c          ICT(1)  or if all adjustments to parameters are less than
c          EPS(10) times corresponding standard deviations.
c
c EPS(11)= Partial convergence criterion for cleanup of matrix inverse
c          in subroutine SNORML.
c EPS(12)= Full convergence criterion for cleanup of matrix inverse
c          in subroutine SNORML.
c
c EPS(13)= Sine of elevation angle below which observed body is below
c          horizon at transmitting site or 2nd observing site, if any,
c          for dummy observations in subroutine HORIZN.
c
c EPS(14)= Fractional parameter adjustment to be used in least squares
c          iteration.
c
c EPS(15)= Convergence criterion for parameter adjustments when using
c          norm**2 metric (if JCT(1) = 1)
c
c EPS(17)= Sine of elevation angle below which observed body is below
c          horizon at observing site for dummy observations in
c          subroutine HORIZN.
c
c EPS(18)= Impact distance in kilometers below which observed body is
c          occulted by its central body for dummy observations in
c          subroutine HORIZN (entry point CNTOCC).
c
c EPS(19)= Sine of elevation angle below which the sun is below
c          horizon at observing site for dummy observations in
c          subroutine SUNHOR.
c
c EPS(20)= Sine of elevation angle below which the sun is below
c          horizon at transmitting site or 2nd observing site, if any,
c          for dummy observations in subroutine SUNHOR.
c
c EPS(21)= Sine of elevation angle below which the sun is below
c          horizon at observed spot on observed body for dummy
c          observations in subroutine SUNHRSP.
c
c ICT(1)=-99 Just print input parameters and stop.
c ICT(1) =-1 No comparison of theory and observation, no least squares
c            analysis.
c ICT(1) = 0 Comparison of theory and observation, no parameters
c            adjusted.
c ICT(1) > 0 Maximum number of least squares iterations involving
c            integration,comparison of theory and observation,solution
c            of normal equations,adjustment of parameters.
c
c ICT(2) =-1 Printout supressed for optical observed minus theory.
c ICT(2) = 0 No printout supression for observed minus theory.
c ICT(2) = 1 Printout supressed for radar and optical observed minus
c            theory.
c
c ICT(3) =-1 Dummy observations are not written on output observation
c            library tapes and partials are not calculated.
c ICT(3) = 0 Dummy observations are not written on output observation
c            library tapes but partials are calculated.
c ICT(3) = 1 Dummy observations are so written with partials calculated
c ICT(3) = 2 Same as 1, but also input observation library tapes with
c            ntape<0 (dummy data) are treated as real data, i.e.,
c            ntape -> |ntape|.
c
c ICT(4)     Governs the computing/copying of partials on obslib tapes.
c            Much of the logic is consolidated in routines pcops ('slow'
c            copying, called from everywhere) and lvt... (merging of L
c            vectors, called from cmpar1,2,3, shpcnt, trgcnt).
c            A 'quick' copy, if done, occurs at the beginning of PARTL
c            and preempts all the rest of the PARTL routines.
c            Note:  some partials (such as satellite/probe parameters
c            and certain non-motion parameters) are always (re)computed
c            unless there is a quick copy.
c            In the following, 'M vector' refers to an L vector read in
c            from tape.
c ICT(4) =-3 Use only saved partial derivatives and do quick copy.
c            Bug:  output tape has the merged L vectors, not the copied
c            M vectors.
c ICT(4) =-2 Use only saved partial derivatives if, for every L vector
c            set, corresponding M vector is also set.  If true, quick
c            copy made, otherwise behaves like ICT(4) = 0.
c ICT(4) =-1 Recompute only those partials that are specified in the
c            input L vector. Saved partials that are not in L vector
c            should be copied and not recomputed.
c            (Not implemented because PARTL has only the merged L vector
c            and not the original input L vector.)
c ICT(4) = 0 Must used saved partial from tape even if the input L
c            vector for that partial is set. If the L vector is set but
c            the partial is not on tape, then compute the partial.
c ICT(4) = 1 Partial derivatives of observations are computed, no saved
c            values are used, and the M vectors are ignored.
c
c ICT(5) =-1 No initialization or forming of normal equations from
c            subset of saved normal equations.
c ICT(5) = 0 Initialize normal equations from subset of saved normal
c            equations, process additional observations and solve
c            normal equations for ICT(1).gt.0 and iterat=1,2,3,...
c ICT(5) = 1 For iterat=1, form normal equations from subset of saved
c            normal equations and solve them if ICT(1).gt.0 (for
c            iterat.gt.1, behaves like ICT(5)=-1 if ICT(1).gt.1)
c ICT(5) = 2 For iterat=1, restore statistics, scale factors, solution
c            and inverse of normal equations from saved quantities on
c            data set jmat if ICT(1).gt.0 (for
c            iterat.gt.1, behaves like ICT(5)=-1 if ICT(1).gt.1)
c ICT(5) = 3 Same as ICT(5)=2, iterat=1 except go to compar link
c            first to generate an observation library tape
c ICT(5) = 4 Restore solution from dataset jmat plus zbar and
c            fbar-adjoint from partially prereduced normal equations
c            on imat0.  Enables computation of adjusts but not sigmas
c            for "reduced away" parameters; uncertainty of predict
c            cannot be calculated.
c
c ICT(6) =   Number of least squares analysis iterations during which
c            perturbing planet data set is to be updated with just
c            completed integrations before starting new integrations.
c
c ICT(7) = 0 No overriding error weighting for dummy observation tape
c            (ntape=negative) in forming normal equations in NRMICT or
c            NRMFRM or predicting residuals in PRDICT.
c ICT(7) = 1 Such overriding is allowed.
c
c ICT( 8)= 0 Series sequencing checked in subroutine compar if
c            ICT(80).ge.0 (if ICT(80).lt.0, no sequence checking)
c ICT( 8)= 1 Series sequencing not checked in subroutine COMPAR
c            for obs.library tape, but it is for cards IOBS and over-
c            riding observation series iobcon unless ICT(80).lt.0
c
c ICT( 9)= 0 Perturbing planet data set for all integrations is ipert
c ICT( 9)= 1 Perturbing planet data set for individual body integrations
c            is n-body data set ibody if nbody.gt.0 rather than ipert
c            (perturbing planet data set for n-body integration is still
c            ipert)
c
c ICT(10)=-2 No prediction of observed minus theory or de-correlation
c            of partial derivatives
c ICT(10)=-1 Such prediction for radio data only if ICT(11).gt.-2
c ICT(10)= 0 Such prediction for radar and optical data if ICT(11).gt.-2
c ICT(10)= 1 Such prediction for optical data only if ICT(11).gt.-2
c
c ICT(11)=-3 No prediction, but de-correlation of partial derivatives
c            using f-bar-adjoint from imat3 as projection matrix (with
c            output tapes determined by the same rules as prediction)
c            if ICT(10).gt.-2
c ICT(11)=-2 No prediction of observed minus theory or de-correlation,
c            even if ICT(10).gt.-2
c ICT(11)=-1 Prediction of observed minus theory with printout,no tape
c            if ICT(10).gt.-2
c ICT(11)= 0 Prediction of observed minus theory with printout and tape
c            if ICT(10).gt.-2
c ICT(11)= 1 Prediction of observed minus theory with tape,no printout
c            if ICT(10).gt.-2
c
c ICT(13)= 0 Normal sequencing of observation library data sets in
c            PRDICT link.
c ICT(13).gt.0  & iterat=1 in PRDICT link, input observation library
c            data set iabs1 is IOBS0(j) for j=1,...,ICT(13) and
c            IOBS1(j) for j=ICT(13)+1,...,numobt
c
c ICT(14)= 0 No uncertainty of prediction calculated in PRDICT link
c ICT(14)= 1 Calculate uncertainty of prediction in PRDICT link
c            need ICT(10).ge.-1, imat.gt.0 if ICT(5).lt.2 or
c            jmat.gt.0 if ICT(5)=2
c
c ICT(15)=-1 Nothing done in ANALIZ link even if ICT(1).gt.0.
c            Program is stopped after call to ANALIZ.
c ICT(15)= 0 Everything done in ANALIZ link which other control integers
c            say to do.
c ICT(15)= 1 Normal equations restored in ANALIZ link (if ICT(1).gt.0),
c            but there is no matrix inversion and solution.
c            Program is stopped after call to ANALIZ.
c
c ICT(16)= 0 Number of tapes read in PRDICT link is numobt
c ICT(16)> 0 Number of tapes read in PRDICT link is ICT(16) which could
c            be less than or greater than numobt (but not greater than
c            ten) there must be appropriate data set numbers in
c            IOBS1(i) or IOBS2(i) as the case may be.
c
c ICT(17)= 0 Form both the left hand side (lhs) and the right hand
c            side (rhs) of the normal equations in subroutine NRMICT.
c ICT(17)= 1 Form the rhs only in subroutine NRMICT.
c ICT(17)=-1 Form the lhs only in subroutine NRMICT.
c
c ICT(18)= 0 No series-by-series adjust.
c        = 1 Series-by-series adjust on. needs: imat, imat1, ibuf3,
c        ibuf4.  imat1 would not be needed if NRMICT were changed.
c
c ICT(19) = 0  Do not use correlated info.
c ICT(19) = 1  Make use of correlated observation information, if
c              available, in NRMICT.
c
c ICT(20) = 0 No printout of partials of observations in COMRIT.
c ICT(20).gt.0 Printout of partials in COMRIT for first ICT(20) obs.
c ICT(20).lt.0 Printout of partials in COMRIT for first iabs(ICT(20))
c              observations in each observation series.
c
c ICT(21)=-1 Doppler observable for probe is phase delay.
c ICT(21)= 0 Doppler observable for probe is averaged frequency.
c ICT(21)= 1 Doppler observable for probe is instantaneous freq.
c
c ICT(23)  binary coded control for atmospheric delay model selection
c     1 bit=1:  use new mapping function code (Mendes & Pavlis, 2002)
c          =0:  use old mapping function code
c     2 bit=1:  use new zenith delay code (Mendes & Pavlis, 2004)
c          =0:  use old zenith delay code
c
c ICT(24)=-1 Print out all elevation angles for radio observations if
c            the effect of the atmosphere or ionosphere is included.
c ICT(24)= 0 Do not print out elevation angles for radio observations.
c ICT(24)= n Print out all elevation angles .le. abs(n) degrees.
c ICT(24)<-1 Print out temp,pressure and humidity at both sites.
c ICT(24)=-3 Print out value for utrec in ATMION.
c
c ICT(25)= 0 No adjustment of range mod codelength to range
c ICT(25).gt.0  2.**(-ICT(25)) is the smallest codelength for which
c               range mod codelength is adjusted.
c
c ICT(26)= 0 Mean sun subroutine used in subroutine STANGL for
c            satellite based look angle observations (obsolete).
c ICT(26)= 1 Normal interpolation for such observations.
c
c ICT(27)= 0 No center of mass of solar system used in radar and
c            interferometer calculations.
c ICT(27)=-1 Position of Earth relative to center of mass of solar
c            system used in radar and interferometer calculations in
c            cislunar space.
c ICT(27)=-2 Same as ICT(27)=-1 with velocity also used.
c ICT(27)= 1 Same as ICT(27)=-1 plus position of center of mass of solar
c            system relative to sun used in planetary and mariner-type
c            space probe radar and interferometer calculations.
c ICT(27)= 2 Same as ICT(27)=1 with velocity also used.
c
c ICT(28)=-1 Nutation precesion evaluated midway between receive
c            and send time for radar and radio tracking data.
c ICT(28)= 0 Nutation precesion evaluated midway for radar data,
c            both ends for radio data.
c ICT(28)= 1 Nutation precesion evaluated at both send and receive
c            times for radar and radio tracking data.
c
c ICT(29)=-1 Derivative of nutation precesion not used for
c            radar or radio Doppler.
c ICT(29)= 0 Derivative of nutation precesion used for radio
c            Doppler but not radar Doppler.
c ICT(29)= 1 Derivative of nutation precesion used for radar
c            and radio Doppler.
c
c ICT(30)= 0 No printout of longitude,latitude of subradar point.
c ICT(30)= 1 There is such printout in sub.deldop for planets for which
c            the longitude,latitude of the subradar pt.are calculated.
c
c ICT(31)= 0 No checkpoint restart in COMPAR link
c ICT(31)= 1,...,10 Checkpoint restart in COMPAR link with jtape=ICT(31)
c ICT(31)=-1,...,-10 Checkpoint restart in COMPAR link with
c                    jtape=iabs(ICT(31)).
c                    Restart at start of tape.
c            (can only have checkpoint restart if ICT(80)=0)
c
c ICT(32)=   Value of ntape (tape sequence number) for COMPAR
c            checkpoint restart.
c
c ICT(33)=   Value of nseq (series sequence number) at which compar
c            checkpoint restart occurs if ICT(33).gt.0.
c            If ICT(33).lt.0, compar checkpoint restart occurs at series
c            following observation series with nseq=-ICT(33).
c
c ICT(34)   Seasonal terms in Earth rotation.
c              ut2-ut1 before 1956, uses erotat con(14-21)
c              a1-ut1 after 1956, uses erotat con(14-21)
c                                 see JCT(29)
c         = 0 no terms included
c         > 0 number of terms included
c             initial epoch is ercon1(1), also used for polynomial
c             epoch.  see JCT(34)
c
c ICT(35)=-1 stop program if error on n-body,planet,embary,Moon data
c            set in compar link
c ICT(35)= 0 skip optical observation, otherwise stop program if error
c            on one of above body data sets
c ICT(35)= 1 skip all execpt probe observations, which stop program,
c            if error on body data set
c ICT(35)= 2 skip all observations for which error on body data set
c
c ICT(36)=   similar logic for s-body data set
c
c ICT(37)=-1 do not skip dummy observation which is below horizon
c            of observing site
c ICT(37)= 0 skip such a dummy observation
c ICT(37).gt.0 same as ICT(37)=0 plus increment time of observation by
c            ICT(37) minutes for first occurance
c
c ICT(38)=-1 do not skip dummy observation which is is occulted by
c            central body of observed body
c ICT(38)= 0 skip such a dummy observation
c ICT(38).gt.0 same as ICT(38)=0 plus increment time of observation by
c            ICT(38) minutes for first occurance
c
c ICT(39)=-2 take nothing from n-body tape, just there for center of
c            mass of solar system
c ICT(39)=-1 take only Moon from n-body tape if there is one
c ICT(39)= 0 take everything from n-body tape that is there (compar l.)
c ICT(39)= 1 use any individual body tapes which are there
c            (note: must set iem=0 to override default if no embary tape
c
c ICT(40)= 0 Earth-Moon barycenter integration uses subroutine FN.
c ICT(40)= 1 Earth-Moon barycenter integration uses subroutine SBFN.
c
c ICT(41) = 0 do not skip observations in compar link
c ICT(41) > 0 skip observations in compar link on first least-squares
c             iteration, namely process every ICT(41) observations
c
c ICT(42) = 1 filter input in input stream
c ICT(42) = 0 no Kalman-Bucy filter (just maximum liklihood,
c             least-squares paramter estimation)
c ICT(42) > 1 filter input on data set no. ICT(42)
c   filter input described in subroutine FILTIN
c
c note: old code was removed in 1983 march
c ICT(43) = 0 use the old code in the radar and interferometry links
c             (DELDOP, SBDLDP, MDELDP; INTERF, MNTERF, etc.)
c ICT(43) = 1 use the new code in the radar and interferometry links
c             (RADCTL, RADSB, RADLND, RADMN; FERCTL, FERMN, etc.)
c
c ICT(44) = 1  apriori covariance matrices are in input stream
c ICT(44) =-1  apriori covariance matrices are in input stream, with no
c              right-hand-side contribution, and the right-hand side is
c              not to be adjusted after each iteration
c ICT(44) = 0  no apriori parameter covariance matrix
c ICT(44) > 1  apriori covar. matrices are on data set no. ICT(44)
c
c ICT(45) = 0 no fractional statistics printout in subroutine frstat
c             in midst of adjustments to parameters
c ICT(45) = 1 such printout on IOUT and JOUT
c
c ICT(46)= 0 normal equations solved and coefficient matrix inverted
c            with double precision arithmetic using fact that it is
c            symmetric in subroutine SNORML
c ICT(46)= 1 same as ICT(46)=0 but extended precision used  in
c            intermediate calculations
c ICT(46)= 2 same as ICT(46)=0 but extended precision used for
c            matrix inversion, and the results are then returned
c            as double precision
c ICT(46)>10 solve subset of the normal equations by eliminating the
c            ICT(46)-10 smallest eigenvalues and inverting the reduced
c            matrix that results
c
c ICT(47)=-1 no printout of scaled or unscaled normal equations
c ICT(47)= 0 printout of unscaled but not scaled normal equations
c ICT(47)= 1 printout of unscaled and scaled normal equations
c
c ICT(48)=-2 no printout of covariance or correlation matrices
c ICT(48)=-1 no printout of inverse of coefficient matrix of normal
c            equations (covariance matrix)
c ICT(48)= 0 printout of inverse of coefficient matrix of normal eqs.
c
c ICT(48)= maximum number of iterations for Gauss-Seidel cleanup of
c            solution in subroutine SNORML (if 0, no cleanup)
c
c ICT(49)= maximum number of iterations for cleanup of matrix inverse
c            in subroutine SNORML (if 0, no cleanup)
c ICT(49)=-1 no iterative cleanup nor calculation of error matrix
c            (product of calcuated inverse and orignal matrix minus the
c            identity matrix)
c
c ICT(50)=0 Moon orbit to be read from nbody or moon tape (+ nutation)
c ICT(50)=1 Moon orbit to be computed from Brown theory if no nbody
c           or moon tape
c
c ICT(55): controls computation of subspacecraft point in subscp
c        0: nothing done
c        1: do computations
c
c
c ICT(66)  binary coded control for mars rotation model and code
c     1 bit=0:  turns off new code. use old code for mars rot.
c          =1:  use new code
c     2 bit=1:  print conversion to proper time
c     4 bit=1:  print out i, psi, phi as calculated in spotpl
c     8 bit=1:  print out spot position vector in mars fixed and
c               1978.0 frame (mars centered)
c    16 bit=1:  print out spot logic vector after call to rotlog.
c               controls drop through sequence that computes partials
c               of rv w.r.t. various parameters
c    32 bit=1:  print out spot position vector in 1978.0 and reference
c               frame (mars centered)
c    64 bit=1:  print out partials of rv w.r.t. various parameters
c   128 bit=1:  print out final partial of position vector in frame
c               of reference epoch w.r.t. parameters
c   256 bit=1:  print out full rotation matrix
c
c ICT(67): control for "fixing" range observable based on the
c          assumption that a bit error occurred in the high order
c          bits that make up the range observable
c     blim=2**(-abs(ICT(67)))
c        blim is the maximum allowed range residual (seconds)
c        blim is also the minimum correction to range allowed
c  if(ICT(67) is less than zero, then the predicted residual
c  is used in the computations
c  this is all done in subroutine rngmod
c
c ICT(76)    control for disparate parameter sne combination
c        = 0 no DPSNEC operations
c        = 2 flag disparities for adjustable nominals
c        = 4 flag disparities for all nominals
c        = 6 copy nominals from imat0(1)
c      =10+n same as n, but also correct rhs for differences
c
c ICT(77)= 0 no punch output of input data at end of program run
c ICT(77)= 1 punch output of input data at end of program run so that
c            job can restarted at termination point with last adjusted
c            values of parameters if ICT(1).gt.0
c ICT(77).gt.1 write punch output to data set no. ICT(77)
c            after every iteration plus do same as ICT(77)=1
c
c ICT(78)= 0 no tapes are saved at end of least squares iteration
c ICT(78)=positive   at end of least squares iterations 1 to ICT(78),
c            inclusive, there is a pause with typewriter message
c            'remove and mount tapes as specified' for operator to
c            follow special instructions (except for last iteration)
c
c ICT(79)=-2 dummy observations not calculated each iteration nor
c            after convergence
c ICT(79)=-1 dummy observations not calculated each iteration, but are
c            after convergence
c ICT(79)= 0 dummy observations calculated each iteration, but not
c            after convergence
c ICT(79)= 1 dummy observations calculated each iteration and after
c            convergence
c
c ICT(80)=-1 input observation library tapes not used
c ICT(80)= 0 input observation library tapes are used, program eats its
c            tail from iteration to iteration
c ICT(80)= 1 input observation library tapes only are used to form the
c            normal equations for first iteration with no computation
c            of observed minus theory.
c            subsequent iterations same as ICT(80)=0
c
c itrwnd(i) = 0 or 1 according to whether data set i is in a rewound
c             state or not, i=1,ntptop
c             actually, any positive number indicates non-rewound.
c             in some cases, values other than 0 or 1 may have
c             significance as to the contents of the data set.
c
c note that if itrwnd in /fcntrl/ labeled common is given a length
c longer than 99, then subroutines xload (referring to extprc) and
c oprmsg (referring to typout) must have dimension of itrwnd changed
c besides all fortran routines that use the JCT vector.
c ntptop should be increased to the new length of itrwnd
c
c
c JCT(1) = 0  use eps(10) as convergence criterion
c JCT(1) = 1  use eps(15) rather than eps(10) as convergence criterion
c
c JCT(2) = jcal meta default-
c any digit of jcal that is = 0 is replaced by the corresponding
c digit of JCT(2)
c JCT(2) = 11 is the default
c if JCT(2).lt.0    propco is shut off
c
c JCT(3): binary coded printout control for propco
c        1 bit=1: print sumcor
c        2 bit=1: print uscor vector
c        4 bit=1: print calcor vector
c        8 bit=1: not used (formerly, print elevation angle and/or rate
c                           see ICT(24) )
c       16 bit=1: print each cal as it is added to sumcor
c       32 bit=1: print cal formation for phase delay dop.
c
c JCT(4): controls use of sumcor vector in propco
c        0: recalculate sumcor
c        1: use sumcor on obslib tape
c
c JCT(6): binary coded flags for compar link debug
c           1 bit: print delays in routine tsav
c           2 bit: print vectors & magnitudes in uvectr
c           4 bit: print positions in interpolators
c           8 bit: print velocities in interpolators
c          64 bit: print initial distance guess in obsred
c         128 bit: print nutation, precession and wobble angles
c         256 bit: print vector partials in cpartl
c         512 bit: print nutation,precession matrices in prcnut,preces
c        1024 bit: print time argument passed to integration reader
c        2048 bit: print relativistic time delay correction
c        4096 bit: print eccentric telescope mount correction
c        8192 bit: print CTAT and A1-UT1 corrections
c
c JCT(10):  binary coded control for solid Earth tides
c           1 bit: calculate effect of lunar tides on radio and laser
c                  ranging observables (coded only in radmn, fermn,
c                  and ferstr)
c           2 bit: calculate the effect of solar tides on these
c                  observables
c           4 bit: print out love numbers and lag angle used
c           8 bit: print out radial,longitudinal, and latitudinal
c                  components of tides at each site
c          64 bit: use IERS tide model
c
c JCT(11):  binary coded control for solid body tides on the Moon
c           1 bit: calculate the effect of Earth-induced tides on radio
c                  and laser ranging observables in radmn and fermn
c           2 bit: calculate the effect of sun-induced tides on the Moon
c           4 bit: print out love numbers and lag angle used
c           8 bit: print out radial, longitudinal and latitudinal
c                  components of tides at each site
c
c JCT(12):  controls meaning of dt parameters
c           -1:  dt(i) is constant value of a1-ut1 or wobble
c                over interval jddt(i) to jddt(i+1)
c            0:  dt(i) is value of a1-ut1 or wobble at jddt(i),
c                tabular point in piecewise linear function
c
c JCT(13):  controls epoch of inertial equatorial coordinate system
c           0: B1950.0
c           1: J2000
c
c JCT(14):  controls use of cross partials on embary and observed body
c           tapes (i.e. 406 on embary tape or 306 on mars tape)
c           does not affect target body partials of non-observed bodies
c           <=0: don't use (default)
c            >0: use
c
c JCT(15):  controls underflow interrupt during inversion of d-matrix
c           in partial prereduction
c           0: underflows generate interrupts which are intercepted by
c              the fortran error handler.  result is set to zero and
c              traceback is printed for first 15 occurrences.
c           1: underflow interrupt is disabled.  result is set to zero
c              by hardware.  there is no indication of any error, nor
c              is there any overhead to process the error.
c
c JCT(16):  controls time tag used for planetary rotation (bodies other
c           than Earth)
c           0: (default) use proper time (omly implemented for mars new
c              rotation model)
c           1: use ephemeris time
c
c JCT(17):  controls Lorentz contraction calculation for observing and
c           observed bodies
c           0: (default) apply contraction
c           1: omit contraction
c
c JCT(20):  binary coded control determines whether PEP stops or
c           continues after certain errors.  0 --> continue,
c           1 --> abort.
c            1 bit: reserved for arithmetic exception in input link
c            2 bit: errors in multiple parameter solution control
c
c JCT(21): controls precession model (also mean sidereal time)
c           -1 use pre-1976 IAU model (default if nothing specified)
c            0 use IAU 1976 (default if non-old precession constant)
c            1 use IAU 2000 (not implemented)
c            2 use IAU 2006
c
c JCT(24)=-1  print out all planetocentric elevation angles for radio
c             observations of a probe within the atmosphere of another
c             planet
c JCT(24)= 0  do not print out planetocentric elevation angles
c JCT(24)= n  print out elevation angles .le. n degrees
c
c JCT(25) =-2 do not calculate transmitter freq. drift rate corrections
c             in alsep differenced n-count observable
c JCT(25) =-1 do not include transmitter freq. drift corrections in
c             alsep differenced n-count observable
c JCT(25) = 0 include such correction
c
c JCT(26) = 0 if ilib.eq.0 get Moon physical libration from Moon or
c             nbody tape
c JCT(26) = 1 if ilib.eq.0 get Moon physical libration by evaluating
c             analytical theory
c
c JCT(27) = 0 usual ordering required in input stream
c         = 1 free ordering of input stream under control of * commands
c             with the following constraints:
c              a. title card
c              b. &nmlst1
c              c. any *object packets
c              d. any other * commands (note: *apriori must follow any
c                 input giving parameters to be adjusted.)
c             note: with * commands the usual &end or blank card
c             is not needed to end a given portion of the input stream.
c
c JCT(28): binary coded control for corrections to Earth rotation.
c          1 bit=1, include corrections to nutation to account for
c                   the Earth's elasticity and fluid core.
c                   (ercond(13-16) in prcnut)
c          2 bit=1, include corrections to a1-ut1 to restore short-
c                   period terms removed by b.i.h. smoothing.
c          4 bit=1, include corrections to diurnal polar motion.
c          8 bit=1, include free nutation corrections.
c                   (ercond(11-12) in prcnut)
c
c JCT(29) = 0 erotat con(1-6) are ad hoc small rotation
c             rates and angles for Earth
c         > 0 con(1-6) are annual wobble terms and free nutation terms
c
c JCT(30) = 0 no e,f,g (+edot,fdot,gdot) tape output
c             Earth centered rotating coordinates. e in true equator of
c             date towards greenwich meridian, g to north, f completes
c             right hand system
c JCT(30).ne.0 such tape output on data set iabs(JCT(30))
c              unless JCT(40) says there is angle output instead
c        .lt.0 no frintout of e,f,g
c        .gt.0 printout of e,f,g at end of dummy observation series
c JCT(31) = space defense center satellite number for e,f,g output
c JCT(32) = 0 no binary e,f,g output in subroutine efgout for fitting
c             polynomials in subroutine efgfit
c JCT(32).ne.0 such binary output on data set iabs(JCT(32)) if
c              JCT(30).ne.0
c        .gt.0 JCT(30) output with JCT(32) output
c        .lt.0 JCT(30) output is supressed
c
c JCT(33) = 0 use BIH UT1 and wobble if available (before 1969 Jan 3
c             must use USNO values; before 1956 Jan 17, must use IPMS
c             wobble and cannot calculate A.1-UT1)
c         =-1 use USNO values of A.1-UT1 and wobble (available from
c             1956 Jan 3 to 1975 Feb 26 for A.1-UT1, to 1975 May 7 for
c             wobble)
c       .gt.0 read A.1-UT1 and wobble tables from data sets
c             JCT(33) and JCT(33)+1, respectively.
c
c JCT(34) = 0 no polynomial correction to CT-UT.
c         = 1 such a polynomial correction, before 1956. (coded for
c             mercury transit observations only)
c             coefficients are ercond(17-19),[con(11-13)].
c             Epoch is hardwired as 1956 Jan 17.0 (see CTUTF).
c             Note: can logically be used to control the polynomial
c             correction to A.1-UT1 after 1956(also using ERCOND(17-19)).
c             However, the current code does not test JCT(34), but
c             rather CON1(1).gt.0 (see A1UT1F).  See ICT(34).
c
c JCT(35) = hour of start of dummy observations each day
c JCT(36) = minute of start
c JCT(37) = hour of  end  of dummy observations each day
c JCT(38) = minute of end
c If JCT(35 to 38) are zero, every time within day is allowed.
c Else, dummy observations done just between JCT(35&36) and JCT(37&38).
c If JCT(37&38).lt.JCT(35&36), then period of obs includes 0 hr UTC.
c These time limits could be the start and end of an optical observatory
c observing evening. This feature is not needed for radar predictions.
c If JCT(35 to 38) nonzero, best to have ICT(37) & ICT(38) .le.0.
c
c JCT(39)=-1 right ascension and declination rates instead of angles
c            themselves are calculated in subroutines of OPTIC
c            (for use in calculating covariance of rate predictions
c            because the partial derivatives of the rates are coded).
c JCT(39)= 0 photographic topcentric right ascension-declination
c            observations are referred to the mean equinox and equator
c            of the reference epoch
c JCT(39)= 1 photographic topocentric right ascension-declination
c            observations are referred to the true equinox and equator
c            of date
c JCT(39)= 2 same as JCT(39)=1 with hour angle replacing right asc.
c            (also applies to transit-circle observations)
c
c JCT(40) = 0 no angle output for optical observatories
c JCT(40).ne.0 such tape output on data set iabs(JCT(30))
c              in format compatable with data general nova computer
c JCT(40).gt.0 printout of angles, data set JCT(40) used as intermediate
c              ebcdic buffer
c JCT(40).le.0 no printout of angles
c
c JCT(41) = 0 angle output referred to true equinox and equator of date
c             for JCT(40) output
c JCT(41).gt.0 angle output referred to mean equinox and equator of
c              1900+JCT(41). for eample, JCT(41)=76 denotes output
c              referred to mean equinox and equator of 1976.0
c
c JCT(42) = 0 no refraction correction made in computing JCT(40) angles
c JCT(42) = 1 refraction correction made in computing JCT(40) angles
c             (optical frequency correction)
c JCT(42) = 2 refraction correction made in computing JCT(40) angles
c             (radio frequency correction, millstone sub.dell)
c
c JCT(43) =   file number for JCT(40) angle output (0 is file 1, 1 is
c             file 2, etc, for data general nova control word at end of
c             514 character record)
c
c JCT(44) = 0 no punch card output for haystack
c JCT(44).gt.0 80 character card image output for haystack on data set
c              JCT(44) (punch if 7) if JCT(40).gt.0
c
c JCT(45) = 0 no prediction print and tape output for Arecibo and other
c             radar sites observing asteroids and Jupiter satellites
c             and also planets.  See also JCT(70)
c JCT(45) < 0 such print output (az,el,ra,decl,delay,Doppler in ANGOT)
c JCT(45) > 0 such print and arecibo bcd site tape from ANGOT.  Output
c             file is JCT(45) and, if >50, also work file on JCT(45)+10.
c
c JCT(46) = 0 JCT(40) print angle output has range in meters and
c             range rate in cm/min
c JCT(46) = 1 JCT(40) print angle output has one way time delay in
c             seconds and one way Doppler shift in Hz
c
c JCT(47) = 0 usual JCT(40) print and tape output
c JCT(47) > 0 JCT(40) tape output has satellite radius instead
c             of Doppler rate, print output has special les-8/9
c             format with 6 frequency Doppler
c JCT(47) = 8 les-8 output
c JCT(47) = 9 les-9 output
c
c JCT(48) = 0 right ascension and hour angle rates are seconds
c             of time per second in JCT(40) print output
c JCT(48) = 1 right ascension and hour angle rates are seconds
c             of arc per second in JCT(40) print output
c
c JCT(49) = control for fluid displacement affecting sites + unit number
c           of external dataset of displacements
c           control has packed bits stored as control*100, independently
c           selecting effects to include in the calculated displacement
c           1: ocean loading
c           2: atmospheric loading
c           4: groundwater loading
c           also a print flag stored as print*10000 to control whether
c           the results are to be printed out as calculated
c
c JCT(50)= 0 usual stoptc,stangl calculations
c JCT(50)= 1 for satellite based Moon observations in stoptc & stangl,
c            do sun also
c
c JCT(51)< 0 restore normal equations after solution and use disparate
c            parameter facility to compute postfit rms for each set of
c            normal equations
c JCT(51)= 0 no multiple parameter set solutions
c JCT(51)> 0 solve for multiple parameter sets
c
c JCT(52)= 0 no summary printout for multiple parameter set (or
c            multiple epoch ) runs
c JCT(52)> 0 summary printout for multiple parameter set (or
c            multiple epoch ) runs
c
c JCT(53)= 0 no partial prereduction of normal equations
c JCT(53)= 1 partially prereduce normal equations using ne generated
c            during this run
c JCT(53)= 2 partially prereduce normal equations using ne generated
c            during an earlier run only (total sne)
c
c JCT(54)= 0 take anything from s-body data set
c JCT(54)=-2 take nothing from s-body data set (compar link)
c JCT(54)= 1 use any available body tapes in preference to s-body tape
c
c JCT(55) is a packed-bits control for the planet/system center of mass
c             correction performed in plcms.
c JCT(55) = 0 coordinates computed in pltrp are purely from tabulated
c             values for the planet itself.
c    1 bit = 1: get satellite coordinates from interpolation if possible
c    2 bit = 1: get satellite coordinates from elliptic if necessary
c    4 bit = 1: use only satellites that are on s-body tape
c    8 bit = 1: for elliptic method, obtain elements from s-body tape
c                (if possible).  warning: this option should be set as a
c                rule because the initial conditions in /empcnd/ are
c                likely to have been overwritten from tape anyway and
c                may not match the epoch in con1(1).
c   16 bit = 1: allow calculation of barycenter of planet+observed
c                satellite instead of planet (not done unless partials
c                are being calculated)
c   32 bit = 1: use precessing elliptic orbit instead of stationary
c
c JCT(56)= 0 NUMINT calls SBOUT, not SBFOUT.
c JCT(56)= 1 NUMINT calls SBFOUT, not SBOUT.
c JCT(56)> 1 NUMINT calls SBFOUT, not SBOUT; SBFOUT restarts filter
c            by reading ICONOF with ITERAT = 1.
c
c JCT(57)<=0 no expanded print of restored normal equation
c            headers by frmhed
c JCT(57)> 0 expanded rne header print by frmhed
c
c JCT(58)<=0 station-to-spacecraft l.o.s.vector not put on obslib
c JCT(58)> 0 l.o.s. vector (in planet rotating frame) output in
c            save vector
c
c JCT(59)= 0 'equator-equinox' terms applied to interferometric
c            observations are used as a pseudo-clock polynomial
c            affecting the observable (just as in radio observations)
c JCT(59)> 0 'equator-equinox' terms are used as true clock corrections
c            (affecting ct epoch + phase and frequency) for
c            interferometric observations only
c
c JCT(60)=-1 like 0, but print diagnostic information
c        = 0 copy a priori matrix + mean residual sensitivity vector
c            to sne datasets and "d" variances to pprsne
c        = 1 do not write a priori or mean residual info (for downward
c            compatibility).  this can allow restoring sne datasets
c            with imbedded a priori matrices without having assigned
c            imat2 (q.v.)
c
c JCT(61)=-1 do not skip dummy observations when sun is above horizon
c          0 skip such a dummy observation
c
c JCT(62)=-1 do not skip dummy observations when sun is above horizon
c            at the observed spot on the observed body
c          0 skip such a dummy observation
c
c JCT(67)= 0 use true velociies in compar link for observations
c        = 1 use apparent velocities in compar link (logic contained
c            in subroutine vlrtrd).  also apply light time iteration
c            corrections to delay observable partials in compar link
c            (forces computation of velocities).
c
c JCT(68)    controls the logical time direction for observation series
c            except dummy observations (which go from start time to
c            finish time)
c JCT(68)=-1 logical time direction is backward
c JCT(68)= 1 logical time direction is forward
c JCT(68)= 0 logical time direction obeys usual criteria, i.e.,
c            forward if nplnt0>30 or ncodf<4 or ncodf>9 or nspot>0
c                       or nspot2>0 or non-Earth observing site
c            otherwise, backward.
c
c JCT(69)= 0 input observation header cards follow old format(see
c            subroutine dltred for details)
c JCT(69)= 1 input observation header cards follow new format(see
c            subroutine dltred for details)
c
c JCT(70)    packed bits controls for ANGOT
c    1 abbreviated print to JOUT in ANGOT if JCT(45).ne.0 (else verbose)
c    2 use mean equator/equinox of reference epoch (else true of date)
c    4 increase reference epoch by 50 yrs (else use default)
c
c JCT(71)= 0 normal choice of source ct-at values
c JCT(71)> 0 radar link uses JCT(71) as selector for ct-at values
c            in ctatf calling sequence (variable nterm)
c JCT(71)= 2 radar link uses analytic annual ct-at term
c JCT(71)= 4 radar link uses integrated annual ct-at term
c
c JCT(78)= 0 Unnormalized gravitational potential tesseral harmonic
c            partial derivatives in SBFN (also zonals)
c JCT(78)= 1 Normalized gravitational potential tesseral harmonic
c            partial derivatives in SBFN (zonals still unnormalized)
c
c JCT(79)= 0 no reintegrations before calculating dummy observations
c            after orbit fit convergence
c JCT(79)= 1 reintegration of motion without partials after orbit fit
c            convergence to be used in calculation of dummy observations
c            if iabs(ICT(79))=1
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           initialize control constants for least squares adjustment of
c           variables in /param/ and /empcnd/ labeled commons
      call ZFILL(Lprm, zlcntr)
c
c lprm(i)=k parameter PRMTER(k) to be adjusted (i=1,100)(k=1,100). if
c           k=0 no parameter adjusted. we must have
c           (1) if j.gt.i and lprm(i)=0, then lprm(j)=0
c           (2) if j.gt.i and lprm(i),lprm(j) not 0, then
c               lprm(j).gt.lprm(i)
c
c lem(i)  =least-squares adjustment control constants for econd,i=1,30
c lmn(i)  =least-squares adjustment control constants for mcond,i=1,30
c ler(i)  =least-squares adjustment control constants for ercond,i=1,30
c lmr(i)  =least-squares adjustment control constants for mrcond,i=1,30
c lpl(i,j)=least-squares adj.control const for pcond(.,j),i=1,30,j=1,
c           numpln
c see subroutine bodred for explanation of lem,lmn,ler,lmr,lpl where
c these quantities are read in &nmlst2 namelist as l.
c
c
c           initialize solar system parameters
      Mass(1)  = 1._10/6.031916E6_10
      Mass(2)  = 1._10/4.08522E5_10
      Mass(3)  = 1._10/3.289001E5_10
      Mass(4)  = 1._10/3.0987E6_10
      Mass(5)  = 1._10/1.0473908E3_10
      Mass(6)  = 1._10/3.4992E3_10
      Mass(7)  = 1._10/2.2930E4_10
      Mass(8)  = 1._10/1.9260E4_10
      Mass(9)  = 1._10/4.0E6_10
      Mass(10) = 1._10/82.301_10
      do i = 11, 90
         prmter(i) = 0._10
      end do
      Relfct     = 1._10
      prmter(41) = 1._10
      prmter(42) = 1._10
      prmter(43) = 1._10
      prmter(44) = 1._10
      Sunpoe     = -7.8E-4_10
      Aultsc     = 499.00478_10
      Reldel     = 1._10
      prmter(62) = 1._10
      prmter(63) = 1._10
      prmter(81) = 1._10
      Ecinc      = 23.4457870675_10
      Seqinc     = 7.25_10
      Seqasc     = 75.0625_10
      Sunrad     = 6.96000E5_10
      prmter(95) = 6.96E5_10
      prmter(96) = 0._10
      prmter(97) = 2440000.5_10
      Mdstau     = 4.263529034E-5_10
      Mdstsc     = -1.E-10_10
      Ltvel      = 299792.458_10
      Gauss      = 0.017202098950000_10
c
c MASS(j) =(mass of body j)/(mass of sun)
c           1 Mercury
c           2 Venus
c           3 Earth+Moon
c           4 Mars
c           5 Jupiter
c           6 Saturn
c           7 Uranus
c           8 Neptune
c           9 Pluto
c MASS(10)=(mass of Moon)/(mass of Earth+Moon)
c MASS(j),j=11,30 masses of additional bodies (divided by mass of sun
c           if no nplnt(i)=j (i=1,u_mxpl). if nplnt(i)=j for some i, then
c           mass is divided by the total mass of the planet-satellite
c           system npcent(i) or mass of sun if npcent(i)=0.)
c     By definition, all MASS values are less than 1, and an input value
c     greater than 1 will be replaced by its inverse, except that values
c     between 1 and 500 are allowed for body 9, which may thus be
c     "hijacked" from the solar system and treated as a distant star.
c
c PRMTER(31)=RELFCT=general relativity motion factor
c PRMTER(32)=GMVARY=time variation factor for gravitational constant
c                   (unit=1/day)
c PRMTER(33)=SUNHAR=second harmonic of gravitational potential of sun
c          (unit=(sun radius)**2
c PRMTER(34-38)=densities (g/cc) of five classes of asteroids
c PRMTER(39)=SUNPOE=ratio of gravitational to total energy of the sun
c            (negative, used for beta, gamma contribution to
c            principle of equivalence violation)
c PRMTER(40)= principle of equivalence parameter deltas (see SBSET).
c PRMTER(41)=metric parameter beta, excluding principle of equivalence
c            violation contribution
c PRMTER(42)=metric parameter gamma, excluding p.o.e. violation
c PRMTER(43)=metric parameter beta, including principle of equivalence
c            violation contribution
c PRMTER(44)=metric parameter gamma, including p.o.e. violation
c PRMTER(45)=mass of 2nd extra asteroid/comet belt (solar masses)
c PRMTER(46)=mass of 1st extra asteroid/comet belt (solar masses)
c PRMTER(47)=right ascension of the ascending node of the asteroid belt
c            on the mean equator of the reference epoch measured from
c            the mean equinox (degrees)
c PRMTER(48)=inclination of the asteroid belt to the mean equator of
c            the reference epoch (degrees)
c PRMTER(49)=distance of circular asteroid belt from sun (au)
c PRMTER(50)=mass of asteroid belt (solar masses)
c PRMTER(51)=AULTSC=value of the astronomical unit in light seconds
c                   (inverse of the speed of light at a large distance
c                   from the sun)
c PRMTER(52)=LTVARY=time variation factor for speed of light
c            (unit=1/day)
c PRMTER(53)=RELDEL=general relativity time delay factor
c PRMTER(54)=RELDOP=general relativity Doppler shift factor (not
c                   programmed) just used as a switch. if greater than
c                   zero, subroutine radar sets ndop for general rel.
c                   effect on Doppler with factor reldel like time delay
c PRMTER(60)= constant plasma factor for interplanetary media
c PRMTER(61)= time variable plasma factor for interplanetary media
c PRMTER(62)= Earth atmosphere effect factor for delay and Doppler
c PRMTER(63)= Earth ionosphere effect factor for delay and Doppler
c PRMTER(72)=CTVARY=coefficient of scale variation for time varying
c            gravitatonal constant
c PRMTER(73)=radius of 1st extra asteroid belt (au). see PRMTER(47-49)
c PRMTER(74)=inclination of 1st extra asteroid belt (deg).    ditto
c PRMTER(75)=ascending node of 1st extra asteroid belt (deg). ditto
c PRMTER(76)=radius of 2nd extra asteroid belt (au).          ditto
c PRMTER(77)=inclination of 2nd extra asteroid belt (deg).    ditto
c PRMTER(78)=ascending node of 2nd extra asteroid belt (deg). ditto
c PRMTER(81)= scale factor for atomic time to coordinate time conversion
c PRMTER(82)= phase for atomic time to coordinate time conversion
c ECINC   = mean inclination of ecliptic on equator at epoch (degrees)
c SEQINC  = inclination of suns  equator on ecliptic (degrees)
c SEQASC  = longitude of ascending node of suns  equator on ecliptic
c           measured from the mean equinox at reference epoch (degrees)
c SUNRAD  = radius of sun in kilometers
c PRMTER(95)= radius of sun in kilometers for transit observations
c PRMTER(97)=2440000.5=epoch for time varying gravitational constant
c MDSTAU  = distance unit for lunar ephemeris on perturbing planet tape
c           in astronomical units (if 0, assumed to be 1)
c MDSTSC  = distance unit for lunar ephemeris on Moon tape in light-
c           seconds (if 0, assumed to be aultsc)
c LTVEL   = velocity of light in kilometers per second
c
c
c           initialize constants for n-body integration
      call ZFILL(Ibody, zbdctr)
      Ibody  = 10
      Epsbdy = 1.0E-16
      Intbdy = 2
      do i = 1, 40
         Kbdy(i) = -1
      end do
      jdpad = 0
      Kbdy(3)  = 0
      Kbdy(21) = 0
      Kbdy(22) = -1
      Kbdy(23) = -1
      Kbdy(27) = -1
      Kbdy(28) = 3
      Kbdy(29) = 11
      Kbdy(30) = 54
      Kbdy(31) = -12
      Kbdy(32) = -14
      Kbdy(36)=0
      Kbdy(39) = 0
      do i = 1, 9
         Nplbdy(i) = i
      end do
c
c nbody  = 0 no n-body integration or n-body data set
c nbody  = n.gt.0, n-body integration and/or n-body data set with n
c                  bodies plus Moon if not one of n=9 bodies
c nplbdy(j),j=1,nbody, are the planet numbers of the bodies in the
c              n-body integration in the order which they are integrated
c              the initial conditions for the integration come from
c              /empcnd/ labeled common
c nplbdy(j)= 3 body is Earth-Moon barycenter relative to the sun with
c              initial osculating elliptic orbital elements econd(1-6)
c nplbdy(j)=10 body is Moon relative to Earth with initial osculating
c              elliptic orbital elements mcond(1-6)
c nplbdy(j)=nplnt(k) for some k=1,...u_mxpl, body is planet nplnt(k)
c              relative to central body ncentr(k) (the sun if ncentr(k)
c              is 0) with initial osculating elliptic orbital
c              elements pcond(1-6,k) (nplnt and ncentr are in /namtim/
c              labeled common)
c
c ibody  =   fortran logical data set number for n-body output
c
c jdbdy0=initial epoch of n-body integration is midnight beginning of
c        day of day with julian day number jdbdy0 (if negative there is
c        checkpoint restart at the epoch -jdbdy0 unless
c        jdbdy0.le.-3000000. if jdbdy0=-1 there is checkpoint restart
c        at record just before end of file of n-body ephemeris tape.)
c        (if jdbdy0=0 or jdbdy0=-3000000, there is no
c        numerical integration, comparison of theory and observation
c        assumes motion is already on tape, but if initial conditions
c        are adjusted, jdbdy0 is set to positive value on second record
c        of tape so there will be reintegration on next least squares
c        iteration, whereas if no initial conditions are adjusted,
c        jdbdy0 is set to zero in the compar link so there will be no
c        reintegration on subsequent least squares iterations)
c        (jdbdy0.le.-3000000 is signal to compar link on first
c        iteration to override end time on tape by input jdbdy2)
c jdbdy1=julian day number of first data record on output tape
c jdbdy2 julian day number of last data record on output tape
c jdpad if positive is used to extend the range of jdpby1-jdbdy2 at both
c        ends (provided those are also both positive). in other words,
c        if jdbdy2>jdbdy1, then jdbdy2 is replaced by jdbdy2+jdpad, and
c        jdbdy1 is replaced by jdbdy1-jdpad.  if jdbdy2<jdbyd1, then
c        jdbdy2 is replaced by jdbdy2-jdpad, and jdbdy1 by jdbdy1+jdpad.
c        this same logic applies to all jd1/jd2 pairs in the same run.
c must have jdbdy0 lying between jdbdy1 and jdbdy2, unless jdbdy0 is
c        negative indicating checkpoint restart, in which case -jdbdy0
c        must be between jdbdy1 and jdbdy2, or unless jdbdy0 is zero
c        indicating no integration or less than or equal to -3000000
c        indicating no integration and input jdbdy2 to everride that on
c        second record of input  n-body data set in compar link (needed
c        because such data set might have been produced by check point
c        restart going beyound point indicated on second record)
c        (if jdbdy0=-1, consistency check for -jdbdy0 between
c        jdbdy1,jdbdy2 is not made)
c time direction of output tape is from jdbdy1 to jdbdy2, forward in
c        time if jdbdy1.lt.jdbdy2, backward in time if jdbdy1.gt.jdbdy2
c if jdbdy0.lt.0 & .gt.-3000000, output tape from a previous integration
c        is read to the time -jdbdy0 (or to the record just before end
c        of file if jdbdy0=-1) and the coordinates on the tape at
c        that epoch are used to restart that integration, proceeding
c        from -jdbdy0 to jdbdy2. integration cannot go in both
c        directions from epoch -jdbdy0 in restart mode and a
c        previous integration cannot be restarted at a time in the
c        first direction of the previous integration if it was in the
c        two direction mode, since output tape is not written until it
c        gets to second direction. in the checkpoint restart mode with
c        jdbdy0 negative, the only data which can be input for the
c        integration are jdbdy1,jdbdy0,jdbdy2,KBDY(39) because the rest
c        of the integration data is taken from the second record of the
c        tape of the previous integration.  if KBDY(39).lt.0, the tape
c        of the previous integration is rewound , whereas if
c        KBDY(39).ge.0, the restarted integration continues writing
c        the tape of the previous integration.
c if jdbdy0 is positive between jdbdy1 and jdbdy2, we have
c       (a) if jdbdy0=jdbdy1, numerical integration starts at jdbdy0
c           and proceeds to jdbdy2, writing output tape as integration
c           proceeds if KBDY(39).ge.0
c       (b) if jdbdy0 is strictly between jdbdy1 and jdbdy2, integration
c           starts at jdbdy0 and goes in direction of jdbdy1 (forward
c           in time if jdbdy0.lt.jdbdy1, backward in time if jdbdy0.gt.
c           jdbdy1). the integration polynomial coefficients determined
c
c           at jdbdy0 by the starting proceedure are saved and the
c           numerical integration proceeds towards jdbdy1, writing onto
c           disk buffer ibuf if KBDY(39).ge.0. when jdbdy1 is reached,
c           ibuf is backspaced and read to write output tape ibody in
c           time direction from jdbdy1 to jdbdy0 if KBDY(39).ge.0. then
c           using the saved starting polynomial coefficients, the
c           integration goes from jdbdy0 to jdbdy2, writing output
c           tape ibody as the integration proceeds if KBDY(39).ge.0.
c           when the integration is completed, we have an output tape
c           from jdbdy1 to jdbdy2, even though integration went from
c           jdbdy0 to jdbdy1 and then from jdbdy0 to jdbdy2. these
c           contortions are gone through in order to simplify the logic
c           needed to read the tape in the compar link.
c
c JVLBDY = 0 velocity as well as position are on output tape
c JVLBDY = m (m.gt.0) velocity as well as position are on output tape
c            forward from the initial epoch and backward from the epoch
c            for m days, but only position data is on tape backward
c            from initial epoch more than m days
c
c EPSBDY =   accuracy constant for controlling step size in n-body
c            integration
c
c INTBDY =-m, m.ge.0, tabular interval is 2**-m days for n-body integr.
c INTBDY = m, m.gt.0, tabular interval is m days for n-body integration
c           This variable is arbitrarily defined to refer to Mercury.
c           There are 5 tabular points per record on output tape for
c           Venus through Pluto, but 10 points for Mercury and 40 for
c           the Moon.  Any value other than 2 will cause problems.
c
c If the Earth-Moon barycenter is integrated in the n-body routine, then
c kbdy(j) is examined for all j=1,nbody. If, for any j, we have
c KBDY(j)>= 0, then the exact equations of motion are integrated,
c              with the coordinates of the Moon relative to the Earth
c              being determined from the integration if the Moon is
c              one of the integrated bodies and from the perturbing
c              planet data set if it is not.  Otherwise, if, for all j,
c KBDY(j)< 0, then the Earth-Moon barycenter is treated as a point particle
c              unless the Moon is also integrated.
c
c If the Moon is integrated here, then all kbdy(j) are treated as 1.
c Further, all planets in the integration are forced to be included
c as perturbers of the lunar orbit, superseding the values specified
c in K(31-39) for the Moon. (see BODRED)
c
c KBDY(19) =-1 Cowell's method used for the Moon if it is one of the
c              integrated bodies (default)
c KBDY(19) = 0 Encke's method with elliptic reference orbit is used for
c              the Moon if it is one of the integrated bodies
c KBDY(19) = 1 mean orbit method is used for the Moon if it is one of
c              the integrated bodies
c
c KBDY(21) =-1 no relativity effect in n-body integration
c KBDY(21) = 0 effect of general relativity included in n-body integr.
c              (with same relativity motion factor relfct used for
c              all bodies)
c KBDY(21) = 1 same as KBDY(21)=0 except that each body has its
c              individual relativity motion factor
c
c KBDY(22) =-1 no effect of time variation of gravitational constant
c              in n-body integration
c KBDY(22) =>0 effect of time variation of gravitational constant
c              included in n-body integration
c
c KBDY(23) =-1 no effect of second harmonic of gravitational potential
c              of sun included in n-body integration
c KBDY(23) =>0 effect of second harmonic of gravitational potential
c              of sun included in n-body integration
c
c KBDY(25)=-1  effect of any input limited asteroids not included in
c              n-body integration
c KBDY(25)=>0  effect of all input limited asteroids (if any) included in
c              n-body integration
c
c KBDY(27) = positive, step size for Adams-Moulton or second sum
c            numerical integration is KBDY(27) days
c KBDY(27).le.zero, interval size for Adams-Moulton or second sum
c            numerical integration is 2**KBDY(27) days
c
c KBDY(28) = 1 Nordsieck method used for n-body numerical integration
c KBDY(28) = 2 Adams-Moulton method used for n-body numerical integ.
c KBDY(28) = 3 royal road (second sum) method used for n-body num.int.
c
c KBDY(29) = number of predictor and corrector terms in adams-moulton
c            or royal road (second sum) n-body numerical integration
c
c KBDY(30) = number of equations controlling step size for  n-body
c            integration (must equal the number of equations)
c            interval size controlled throughout nordsieck method, but
c            only during nordsieck starting proceedure for adams-
c            moulton method.  This input variable is currently ignored,
c            and the actual number of equations is used instead.
c
c KBDY(31) = interval in starting proceedure for n-body integration
c            is 2**KBDY(31)
c
c KBDY(32) = minimum interval in n-body integration is 2**KBDY(32)
c
c KBDY(36) =-1 use simplified G.R. formulation.  See K(61) in BODRED.
c KBDY(36) = 0 use full General Relativity formulation.
c
c KBDY(37) = 0 initial conditions are osculating elliptic orbital elem.
c KBDY(37) = 1 initial conditions are position and velocity
c
c KBDY(38) = 0 printout every tabular point if KBDY(39).le.0 for n-body
c              integration
c KBDY(38) = n.gt.0, printout every n tabular points if KBDY(39).le.0
c              for n-body integration
c
c KBDY(39) = -1 printout,no tape for n-body integration
c KBDY(39) = 0 printout and tape for n-body integration
c KBDY(39) = 1 tape, no printout for n-body integration
c
c KKBDY(j) =0 asteroid with planet number j+10 does not perturb n-body
c             integraton (j=1,20)
c KKBDY(j) =1 such individual asteroid does perturb n-body integration
c
c KKBDY(70)   copy of JCT(13), not independently settable
c
c KKBDY(71) = 0 do not print accelerations
c KKBDY(71) = 1 print accelerations every tabular point on Kout, if any
c
c KKBDY(80)=0 do not include asteroid ring in n-body integration
c KKBDY(80)=1 effect of asteroid ring included in n-body integration
c
c
c           initialize dt parameters
c     1956 Jan 17.0 - beginning of A.1 time
      call ZFILL(Dt0, zdtprm)
      Jddt0 = 2435490
      Dt0   = 31.3669E0
c
c jddt0 =-1  dt table refers to wobble only
c jddt0  = 0  dt table refers to A.1-UT1 and wobble (used after 1956)
c jddt0  = 1  dt table refers to A.1-UT1 (used after 1956)
c jddt0  > 1  dt table refers to ET-UT2 (used before 1956)
c numdt  =  number of tabular points in each segment (A.1-UT1, xwob, and
c           ywob, or just A.1-UT1 or just ET-UT2) of adjustable dt table
c           if A.1-UT1, xwob, and ywob are respresented by the table,
c           the maximum number of each is 200.  if just A.1-UT1 or
c           ET-UT2 are represented, the maximum number is 600.
c dt0    = value of A.1-UT2 at Julian date jddt0-0.5
c dt(i)  =  value of A.1-UT1 or ET-UT2 at Julian date jddt(i)
c dt(i+200)=value of xwob(in arcsecs) at Julian date jddt(i), i=1,numdt
c dt(i+400)=value of ywob(in arcsecs) at Julian date jddt(i), i=1,numdt
c ldt(i) = 0 dt(i) not adjusted in least squares analysis, i=1,numdt
c ldt(i) = 1 dt(i) adjusted in least squares analysis, i=1,numdt
c see also JCT(12)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           initialize normal equation weighting factors
      do i = 1, 10
         Wgtobs(i) = 1._10
      end do
      do i = 1, maxmt0
         Wgtmt0(i) = 1._10
      end do
c initialize global deletions -- default is no deletions
      Tdlt0=0._10
      Tdlton=0._10
      Tdltof=0._10
c initialize mascons -- default is no mascons
      call ZFILL(Masmsc,zmascn)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      repeat = 0
c repeat=0  job ends at end of run with given input data (default)
c repeat=1  job repeats at end of run with new input data
c
      noprnt = 0
c noprnt =0 print out all data and control constants in sub.prntot
c noprnt =1 suppress such printout
c
      iseq = 0
c     iseq   =0 check tape and series sequencing numbers for increasing
c               order for overriding error weighting and dummy
c               observtion cards in subroutine prnobs and stop program
c               if there is a sequencing error
c     iseq   =1 do not stop program if there is such a sequencing error
c
c        extprc =-1 extended precision arithmetic is done double prec.
c        extprc = 0 extended precision arithmetic uses hardware
c        extprc = 1 extended precision arithmetic uses software
      Extprc = 0
c
c typout = 0 no messages typed to operator in subroutine typout
c typout = 1 messages are typed to operator in subroutine typout
      Typout = 0

      maxpln = u_mxpl
c     maxpln= maximum number of &nmlst2 with nplnt not 3,10,-3,-10
c
c
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c     setup for jcal default - - - - default = 00
      do i = 1, 100
         Jcal(i) = 0
      end do
c
c     jcal controls the calculation and use of propagation media
c     corrections in the radar and interferometry links
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
c     tens digit =
c                 1 : do not use the ith correction
c                 2 : use standard logic for correction use
c                 3 : use the correction regardless of rank
c     - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           spool &nmlst1 from in to in0 with a-format printout
      call PEPTIC(In, Iout, in0, 8, 'SOLAR SYSTEM PARAMETERS &NMLST1 ',
     .            nstop, 0)
c
c read &nmlst1 namelist from in0
      read(in0, NMLST1)
      rewind in0

      call EXTFLG(Extprc)
      sepeat = repeat
      moprnt = noprnt
      iseq1  = iseq
      Ipert0 = Ipert
      Jpert0 = Jpert
      Kpert0 = Kpert
      Kkbdy(70) = Jct(13)
      if(Nbody.gt.0) Kbdy(30) = 6*Nbody
      if(Nbody.gt.0 .and. Mdstsc.lt.0._10) Mdstsc = 0._10
c
c modify integration range if requested
      if(Jdbdy1.gt.0 .and. Jdbdy2.gt.0 .and. jdpad.gt.0) then
         if(Jdbdy2.ge.Jdbdy1) then
            Jdbdy2=Jdbdy2+jdpad
            Jdbdy1=Jdbdy1-jdpad
         else
            Jdbdy1=Jdbdy1+jdpad
            Jdbdy2=Jdbdy2-jdpad
         endif
      endif
c
c invert inverse masses
      do i = 1, 30
         if(Mass(i).gt.500._10 .or. (Mass(i).gt.1._10.and.i.gt.9))
     .    Mass(i) = 1._10/Mass(i)
      end do
c
c set up global data deletions
      Tdltp=Tdlton+Tdltof
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c set up propco controls
      call USCOSU
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      return
      end
