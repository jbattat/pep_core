      PROGRAM ABC
      IMPLICIT NONE
C
C         Prints, plots and performs Fourier analysis on quantities in
C         PEP obslib tape; copies obslib tape deleting series and points
C         as specified.
C
C         Input parameters fall into several groups:
C
C         Parameters describing program activity
C             REPEAT - allows multiple input namelists
C                 0 - (default) no more input
C                     (must be explicitly set on last of multiple inputs
C                 1 - read another input, no reinitialization
C                 2 - read another input, first reinitializing
C             TASK - where task=n implies all tasks .le. n - int*4
C                 0 - do nothing (default)
C                 2 - select and display data for analysis
C                 3 - do Fourier transform on data
C                 4 - do autocorrelation on power spectrum
C             FFT - if true, use fast Fourier transform of binned data
C                     (default if TASK is 3 or more, required otherwise)
C                   if false, use unbinned data and forbid data plots
C             WTFT - if true and if FFT is false, use data weights in
C                      estimating Fourier transfor
C                    if false of if FFT is true, use unweighted Fourier
C                      transform (default)
C             LOOK - binary coded option for just looking at data - I*4
C                 separate pass through tape is made and dependent
C                 variable can be printed and/or plotted
C                    (ignore if .le.0)
C                 1 - print
C                 2 - plot
C                16 - print all partial derivatives
C             KBYTES - storage array size in kilobytes (def=64)
C             NOPRNT - flag for printout of input-stream values
C                 0 - print all input constants (default)
C                 1 - do not print the input constants
C             ISEED - seed for random-number generator
C             RNDMDT - if true, then replace input data with zero-mean
C                      gaussian-distributed random values
C
C         Data sets and related parameters
C             CONTRL - input parameters - default value 5
C             PRINT - printer output - default value 6
C             OBSLIB(10) - observation library tapes - default 1,0,...,0
C             OBSOUT - copy of OBSLIB(1) with selected series
C             SPOOL - spooling dataset - default 99
C             IPRINT - binary coded printing option - default 15
C                 1 - data read (for selection; see also LOOK)
C                 2 - data selected
C                 4 - data transformed
C                 8 - data autocorrelated
C                16 - all partial derivatives
C
C          Parameters for skipping data
C             NOPART - (default=.false.)
C                    =.true., no partials copied to obsout
C             NTH - accept only every nth observation in each series
C                 = 0 (default)  no skipping of points
C                 > 0 - every nth point accepted according to
C                        the following qualifiers
C             NTHBEG - start the sequence of accepted observations at
C                      number NTHBEG (default 1)
C             DTMIN - time interval in minutes within which
C                     subsequent observations are skipped
C                   = 0.0 (default)
C             DTMAX - do not skip observation if time interval
C                     from last record exceeds dtmax minutes
C                   = 1.0e10 (default)
C
C         Parameters to be checked on obslib tape (blank value or -999
C         implies no check - see initialization).  50 series allowed.
C             NSEQ(50) - series sequences
C             NCODF(50) - type of observation (type 3 record)
C             NPLNT(50) - planet number
C             RSITE(50) - receiving site
C             SERIES(50) - series name
C             SSITE(50) - sending site
C             NCENTB(50) - central body number
C             SPOT(50) - observed spot
C             SPOT2(50) - 2nd observed spot
C             FREQ(50) - observing frequency (tolerance 1 ppm)
C
C         Parameters describing range of time accepted
C             START - lower bound on range of values considered - REAL*8
C                 (if time used then may have also)
C                 STHR, STMIN, STSEC - REAL*4
C                 (which determine a fraction of start, a Julian day)
C             STOP - upper bound - REAL*8
C                 SPHR, SPMIN, SPSEC - REAL*4
C             TSPAN - maximum time span treated as a single series in
C                     "looker" plots (days) - REAL*8
C
C         Parameters describing range of independent variable
C             NIV - selects independent variable, default 0 (time)
C                 - .lt. 0 - index into save array
C                 - .eq. 1 - uses result (observable)
C                 - .gt. 1 - index into deriv
C                   1000+n - index (n) into cal array
C             TWOOBS - (default false) if true then independent variable
C                      is selected from the 1st observable if 0<niv<1000
C                      and dependent variable from the 2nd if 0<ndv<1000
C                      and only observations with both are selected
C             NORMIV - if true, then divide the independent variable by
C                      the corresponding error.  Default is false.
C             MINIV,MAXIV - if NIV ne 0,
C                           Min and max value for independent variable
C             PERIV(2) - declares periodic boundary conditions with
C                        period PERIV(1) and forced into the span
C                        beginning at PERIV(2).  Takes effect after
C                        MINIV/MAXIV selection.  Default: 0,0 (none).
C
C         Parameters describing sampling configuration in terms of
C         independent variable.  (Two and only two of the following
C         three parameters must be non-zero.  The third is determined
C         from the relation:  SIZSPN = NPOINT*DELTA.)
C             SIZSPN - range of one span - REAL*8
C                 SZHR, SZMIN, SZSEC - REAL*4
C             NPOINT - number of sample points in one span - INT*4
C                 (If NPOINT<0, NPOINT = 2**(-NPOINT) and
C                 FOUR1 is used instead of FOURG.)
C             DELTA - interval between points - REAL*8
C                 DEHR, DEMIN, DESEC - REAL*4
C             NSPAN - number of spans - INT*4
C             NSERIE - number of series to be read - INT*4
C             STRTVL(20) - span starting values (optional) - REAL*8
C             AVGSPN - .true. if power spectra to be averaged over spans
C             PADZ - if true then npoint zeroes are tacked on to end
C                 of data as selected and a 2*NPOINT span is analyzed.
C
C         Parameters describing dependent variable
C             NCODE - type of observation (selects records in RDLIB)
C                 1 - accept only 1st observable
C                 2 - accept either observable (choose 1st if both)
C                 3 - accept only 2nd observable
C             NDV - selects dependent variable - default 2 (O-C)
C                <1 - index into SAVE array
C                 1 - uses RESULT (observable)
C                >1 - index into DERIV array
C                 1000+n - index (n) into CAL array
C             NORMDV - if true, then divide the dependent variable by
C                      the corresponding error.  Default is false.
C             PERDV(2) - declares periodic boundary conditions with
C                        period PERDV(1) and forced into the span
C                        beginning at PERDV(2).  Default: 0,0 (none).
C             FNAME & LNAME - first & last names of a PEP parameter
C                 If NDV .ge. 3 and FNAME .ne. blanks, the dependent
C                 variable becomes the partial wrt parameter with given
C                 names, if found - default blanks
C                 Note: UT & wobble parameters may not be specified
C             NAMES - if .true. then names of parameters wrt which
C                 there are partials are printed for each series
C             CUT(2) - cutoff for residuals (NDV = 2) - REAL*4
C                 (if abs(O-C) .ge. abs(CUT*ERROR), residual ignored)
C                 CUT(1) for time delay, (2) for Doppler.
C                 If the corresponding CUT(i).eq.0 then do not cut.
C
C         Parameters describing smoothing filter (local) and weighting.
C             FILTER - choice of filter - INT*4
C                 1 - nearest point
C                 2 - rectangular
C
C         Sort options
C             SORT=.true.    (default=false)
C                 This will cause a sort to be done on data for "looker"
C                 plots.  Sort-Oneaxe plots will self-scale like
C                 non-Oneaxe plots unless EXTRMX, EXTRMY are specified.
C                 (if EXTRMX is non-zero, x data will be scaled like
C                 normal Oneaxe plots.  If EXTRMY is nonzero, y data will
C                 be scaled like normal Oneaxe plots.)  Default is
C                 self-scaling.
C             BGAP - when set equal to nonzero(days) in sort mode,
C                    successive points with separation bigger than
C                    gap will signal a new graph.
C             BSPAN - when set equal to nonzero(days) in sort mode,
C                     the independent axes will be limited by BSPAN
C                     plus (BSHR(hours) + BSMIN(min) + BSSEC(sec))
C                     if they are specified.
C             BURSTL - when set nonzero, it is used together with
C                      BGAP/BSPAN to find the scale of the x axis.
C
C         Parameters describing plotting
C             XFRAME - maximum extent allowed in x direction
C             YFRAME - maximum extent allowed in y direction
C                       (defaults depend on hardware)
C             HIPREC - if .true., use high precision labelling
C                      option for the x-axis (default=f)
C             IPLOT - binary coded plotting option - default 0
C                 0 - no plot
C                 2 - data selected
C                 4 - power spectra of transforms
C                 8 - autocorrelation
C             XDEN - density of points in x-direction - REAL*4
C                    (default = 10.0)
C             XLEN - maximum length of x axis - REAL*4
C                    default = XFRAME
C                    For ONEAXE=t, XLEN is used.
C                    For ONEAXE=f, smaller of XLEN,N/XDEN
C                         is used where N= no. of points.
C             YLEN - maximum length of y axis - REAL*4
C                 Default 10. (in.)  If less than 10., successive frames
C                 in looker plots (REPEAT>0) are stacked.
C             TMODE - selector for plot title configuration.  Default
C                     depends on environment (generally 1 or 4).  1 and
C                     3 are interchangeable as input for looker plots.
C                  0: no titles
C                  1: all titles run vertically up the paper (only if
C                     non-Oneaxe).
C                  2: tape titles run up, plot titles across top
C                  3: all titles run across at left (only if Oneaxe).
C                  4: no tape titles, plot titles across top
C             FACT - scale factor for Calcomp plot - REAL*4
C             Oneaxe - if true then all looker plots done on one set of
C                 axes.  Need EXTRMX,EXTRMY unless sort specified.
C             JDREF - time origin for looker plots with NIV=0
C                      default=2440000
C             SYMOPT - if .true., points having the same
C                      series name and the same sending and receiving
C                      sites but belonging to different sequences will
C                      plotted with the same symbol. (default=f)
C             PLTSYM - If 999, choose plotting symbol by algorithm (def).
C                      Else use code 'pltsym' for all points, and, if
C                      ZLINE is false, assume this run is part of
C                      a pair with identical indep.var range -- plot
C                      vertical axis either left or right accordingly.
C             OUTIN - if true, dependent variable falling outside
C                     boundary of plot is plotted at boundary
C             EXTRMX(3) - min and max for x axis on Oneaxe plot - R*4
C                           (in JD, relative to 2440000)
C                         non-zero values overides automatic x scaling
C                         for sort-Oneaxe plots.
C                         3rd element modifies automatic x scaling by
C                         increasing the range on each end.
C             EXTRMY(3) - min and max for y axis on Oneaxe plot - R*4
C                         non-zero values overides automatic y scaling
C                         for sort-Oneaxe plots.
C                         3rd element modifies automatic y scaling by
C                         increasing the range on each end.
C             ASPECT - if true, scale looker plots of RA vs. DEC to
C                      give the true aspect ratio, default false.
C             INDEX(20) - pairs of indices into power spect. - INT*4
C                 (alternate terms comprise a pair which tell a part of
C                 the power spectrum to be plotted.)
C             ZLINE - if .true., draw zero line for dependent var. and
C                     also draw ticks on right and top of plot
C             SYMBEL - if .true., plot all points with symbols, changing
C                 symbols for each series.  If .false, plot first point
C                 of series with symbol, connect rest with line.
C                 (default=t)
C             DOTPLT - if true, then use dot instead of + as preferred
C                     sole plotting symbol (must have symbel=t), default
C             NDEC - number of decimals after the decimal point
C                    in axis numbering (default=2)
C                  - same convention as Calcomp axis routine
C                    (i.e.,ndec=-1 - no decimal point written)
C             NBOX - number of boxes for histograms
C             SZBOX - size of each box in histograms (ignored unless
C                     nbox is .le.1).  If either of these is non-zero,
C                     then LOOKER calls HSTGRM instead of PLOTER.
C             TLAT(2) - min and max latitude values of subradar point
C                       allowed if NIV= -1 (plot longitude)
C             LPLEMS - length in bytes of special message to display
C                      after each plot frame (default 0)
C             PLEMS  - special message to display after each plot frame
C
C         COMMON
      include 'burst.inc'
      include 'bursu.inc'
      include 'checks.inc'
      include 'graphs.inc'
      include 'hrmnsc.inc'
      REAL*4 HMS(12)
      EQUIVALENCE (HMS(1),STHR)
      include 'inodtabc.inc'
      include 'misc.inc'
      include 'namen.inc'
      include 'pltfrm.inc'
      include 'pltlab.inc'
      include 'plttyp.inc'
      include 'skip.inc'
      include 'sortv.inc'
      include 'span.inc'
C
      INTEGER*4 I,NOPRNT
      REAL*4 XFRAMD,YFRAMD
      REAL*8 RANDOM,TRAN
C
C         NAMELIST
      NAMELIST /INPUT/ START,STOP,SIZSPN,DELTA,NSERIE,NSPAN,
     1 NPOINT,NDV,AVGSPN,PADZ,IPRINT,CONTRL,PRINT,OBSLIB,OBSOUT,
     2 RSITE,SSITE,SERIES,NSEQ,NCODF,NPLNT,NCENTB,NCODE,
     3 STHR,STMIN,STSEC,SPHR,SPMIN,SPSEC,SZHR,SZMIN,SZSEC,
     4 DEHR,DEMIN,DESEC,LOOK,REPEAT,ZLINE,SYMBEL,KTH,FFT,WTFT,
     5 FILTER,IPLOT,STRTVL,INDEX,XDEN,XLEN,YLEN,FACT,TASK,CUT,
     6 EXTRMX,EXTRMY,ONEAXE,FNAME,LNAME,NAMES,PROBNO,PROGNO,RNDMDT,
     7 DTMIN,DTMAX,NTH,NTHBEG,OUTIN,NOPART,NDEC,NIV,TLAT,MINIV,MAXIV,
     8 SPOT,SPOT2,NBOX,SZBOX,SORT,HIPREC,SYMOPT,XFRAME,YFRAME,TMODE,
     9 BGAP,BSPAN,BSHR,BSMIN,BSSEC,BURSTL,JDREF,KBYTES,NOPRNT,DOTPLT,
     A LPLEMS,PLEMS,NORMIV,NORMDV,TWOOBS,ASPECT,TSPAN,PERIV,PERDV,
     B PLTSYM,ISEED,FREQ
C
C* START=100
C
C         LOCAL
      CHARACTER*8 BLANK /'        '/
      CHARACTER*4 BLNK4
      EQUIVALENCE (BLANK,BLNK4)
      CHARACTER*8 PROBNX /'ABC'/, PROGNX /'OBSLIB'/
      CHARACTER*80 BUFF
      INTEGER*4 TMODD
C
C           INITIALIZE LOCALE VARIABLES
      IPC=.FALSE.
      XMRG=0.67
      YMRG=0.67
      TMODD=1
      GOTO (210,220,230,240,250,260,270,280,290), LOCALE
C           DRAPER CALCOMP
  210 XFRAMD=240.
      YFRAMD=11.
      GOTO 300
C           DRAPER PRINTER
  220 XFRAMD=9.6
      YFRAMD=7.5
      TMODD=4
      GOTO 300
C           DRAPER VERSATEC
  230 XFRAMD=60.
      YFRAMD=10.55
      GOTO 300
C           IPS CALCOMP
  240 IPC=.TRUE.
      XFRAMD=100.
      YFRAMD=11.
      GOTO 300
C           HP TERMINAL
  250 XFRAMD=14.4
      YFRAMD=7.2
      TMODD=4
      GOTO 300
C           OIT CALCOMP
  260 XFRAMD=150.
      YFRAMD=11.
      GOTO 300
C           LCG PLOTTER
  270 XFRAMD=150.
      YFRAMD=30.
      YMRG=5.0
      GOTO 300
C           APPLE LASERWRITER
  280 XFRAMD=10.6
      YFRAMD=8.0
      TMODD=4
      GOTO 300
C           TEKTRONIX TERMINAL (MS-KERMIT)
  290 XFRAMD=13.65
      YFRAMD=10.24
      TMODD=4
  300 CONTINUE
C
C           ITEMS INITIALIZED JUST ONCE
      THGT=0.105
C
      DUAL=0
      HLEN=PROPOR*THGT*88+0.5
      PFIRST=.TRUE.
      REPEAT=0
      TFIRST=.TRUE.
      XLNOM=0.
C
C         INITIALIZE
C*  START=1000
 1000 CONTINUE
      ASPECT=.FALSE.
      AVGSPN = .TRUE.
      BGAP=0.0
      BSPAN=0.0
      BSHR=0.0
      BSMIN=0.0
      BSSEC=0.0
      BURSTL=0.0
      BURSTS=0.0
      CONTRL = 5
      CUT(1) = 0.0
      CUT(2) = 0.0
      DELTA=0D0
      DOTPLT=.FALSE.
      DTMIN = 0.0D0
      DTMAX = 1.0D10
      EXTRMX(1) = 0.0
      EXTRMX(2) = 0.0
      EXTRMX(3) = 0.0
      EXTRMY(1) = 0.0
      EXTRMY(2) = 0.0
      EXTRMY(3) = 0.0
      FACT = 1.0
      FFT=.TRUE.
      FILTER = 1
      FNAME = BLANK
      HIPREC=.FALSE.
      IPLOT = 0
      IPRINT = 15
      ISEED=0
      ISYMNO=0
      ISYMNX=0
      JDREF=2440000
      KBYTES=64
      KTH = 1
      LOOK = 0
      LNAME = BLANK
      LPLEMS=0
      MINIV= -1.0E+20
      MAXIV=  1.0E+20
      NAMES = .FALSE.
      NBOX=0
      NCODE = 2
      NDEC= 2
      NDV = 2
      NIV= 0
      NOPART= .FALSE.
      NOPRNT=0
      NORMDV=.FALSE.
      NORMIV=.FALSE.
      NPOINT = 0
      NSERIE = 1
      NSPAN = 1
      NSYM=-1
      NTH = 0
      NTHBEG = 1
      OBSLIB(1) = 1
      DO 1020 I = 2, 10
      OBSLIB(I) = 0
 1020 CONTINUE
      OBSOUT = 0
      ONEAXE = .FALSE.
      OUTIN= .FALSE.
      PADZ = .FALSE.
      PERIV(1)=0D0
      PERIV(2)=0D0
      PERDV(1)=0D0
      PERDV(2)=0D0
      PLTSYM=999
      PRINT = 6
      PROBNO=PROBNX
      PROGNO=PROGNX
      RNDMDT=.FALSE.
      SIZSPN=0D0
      SORT=.FALSE.
      SPOOL = 99
      START=0D0
      STOP=1D10
      SYMBEL = .TRUE.
      SYMOPT=.FALSE.
      SZBOX=0.
      TAPENO = 0
      TASK = 0
      TLAT(1)= -30.0
      TLAT(2)=  30.0
      TMODE=TMODD
      TSPAN=0D0
      TWOOBS=.FALSE.
      USEF1 = .FALSE.
      WTFT=.FALSE.
      XDEN = 10.0
      XFRAME=XFRAMD
      XLEN=XFRAMD
      YFRAME=YFRAMD
      YH=YFRAMD-1.0
      YLEN = 10.0
      ZLINE = .FALSE.
C
      DO 1030 I = 1, 50
      RSITE(I) = BLANK
      SSITE(I) = BLANK
      SERIES(I) = BLNK4
      NSEQ(I) = -999
      NCODF(I) = -999
      NPLNT(I) = -999
      NCENTB(I) = -999
      SPOT (I)=BLNK4
      SPOT2(I)=BLNK4
      FREQ(I)=-999._10
 1030 CONTINUE
      DO 1040 I=1,20
      INDEX(I) = 0
      STRTVL(I) = 0.0D0
 1040 CONTINUE
      DO 1060 I=1,12
 1060 HMS(I)=0.
      IF (REPEAT .GT. 0) GO TO 3000
C
C*  START=2000
C         WRITE DAY AND TIME
      CALL DAYTIM
C
C         SPOOL INPUT STREAM
      WRITE(PRINT,100)
  100 FORMAT('0*** INPUT STREAM ***',T74,'XXXXXXXX'/)
      DO WHILE (.TRUE.)
         READ(CONTRL,110,END=2010) BUFF
  110    FORMAT(A80)
         WRITE(PRINT,120) BUFF
  120    FORMAT(' ', A80)
         BUFF(73:80)=BLANK
         I=1
         DO WHILE (I.LT.73)
            IF(BUFF(I:I).EQ.'$') THEN
               BUFF(I:72)=BLANK
               I=73
            ENDIF
            I=I+1
         END DO
         WRITE (SPOOL,110) BUFF
      END DO
 2010 END FILE SPOOL
      REWIND SPOOL
C
C*  START=3000
C         READ INPUT PARAMETERS AND ECHO OUT
 3000 CONTINUE
      READ(SPOOL, INPUT)
      IBURST=0
      IF(.NOT.SORT) GOTO 3002
      IF((BSPAN+BSHR+BSMIN+BSSEC).EQ.0.) GOTO 3001
      IBURST=3
      BURSTS=BSPAN+((BSSEC/60.+BSMIN)/60.+BSHR)/24.
 3001 CONTINUE
      IF(BGAP.EQ.0.0)GOTO 3002
      IBURST=2
      BURSTS=BGAP
 3002 CONTINUE
      HISTOG= NBOX.GT.0.OR.SZBOX.GT.0.
      IF(MOD(LOOK/2,2).EQ.1 .AND. IPLOT.EQ.0) IPLOT=1
      NCHKS=0
      IF(RSITE(1).EQ.BLANK .AND. SSITE(1).EQ.BLANK .AND.
     1 SERIES(1).EQ.BLNK4 .AND. NSEQ(1).EQ.-999 .AND.
     2 NCODF(1).EQ.-999 .AND. NPLNT(1).EQ.-999 .AND.
     3 NCENTB(1).EQ.-999 .AND. SPOT(1).EQ.BLNK4 .AND.
     4 SPOT2(1).EQ.BLNK4 .AND. FREQ(1).EQ.-999._10)  GOTO 3020
      DO 3010 I=1,50
      NCHKS=51-I
      IF(RSITE(NCHKS).EQ.BLANK .AND. SSITE(NCHKS).EQ.BLANK .AND.
     1 SERIES(NCHKS).EQ.BLNK4 .AND. NSEQ(NCHKS).EQ.-999 .AND.
     2 NCODF(NCHKS).EQ.-999 .AND. NPLNT(NCHKS).EQ.-999 .AND.
     3 NCENTB(NCHKS).EQ.-999 .AND. SPOT(NCHKS).EQ.BLNK4 .AND.
     4 SPOT2(NCHKS).EQ.BLNK4 .AND. FREQ(NCHKS).EQ.-999._10)  GOTO 3010
C        FOUND HIGHEST NUMBERED CHECK ITEM
      GOTO 3020
 3010 CONTINUE
 3020 CONTINUE
C
      IF(YLEN.GT.YFRAME-YMRG) YLEN=YFRAME-YMRG
      IF(YLEN.LT.1.0) YLEN=1.0
      YLEN=AINT(YLEN)
      YL=YLEN
C
      IF(PLTSYM.EQ.999 .OR. ZLINE) THEN
         DUAL=0
      ELSE
         DUAL=DUAL+1
         IF(DUAL.GT.2) DUAL=1
      ENDIF
C
C        PRINT OUT INPUT CONSTANTS
      IF(NOPRNT.EQ.0) CALL PRNTCN
C            IF TIME RANGE IS GIVEN, CONVERT HR,MIN,SEC TO FRAC OF DAY
      IF(START.GT.0.0D0)  START= START +
     1    ((STSEC/60.0D0+STMIN)/60.0D0 + STHR)/24.0D0
      IF(STOP.LT.1.0D10)  STOP = STOP +
     1    ((SPSEC/60.0D0+SPMIN)/60.0D0 + SPHR)/24.0D0
C        EXERCISE RANDOM-NUMBER GENERATOR
      TRAN= RANDOM(ISEED)
      TRAN= RANDOM(ISEED)
      TRAN= RANDOM(ISEED)
      TRAN= RANDOM(ISEED)
      TRAN= RANDOM(ISEED)
      TRAN= RANDOM(ISEED)
C
      CALL WRKARR(KBYTES)
C
C         GET ANOTHER INPUT STREAM
      IF (REPEAT .GE. 1) CALL AGAIN
      IF (REPEAT .EQ. 1) GO TO 3000
      IF (REPEAT .EQ. 2) GO TO 1000
C
C         CLOSE PLOT TAPE IF USED
      IF(.NOT.PFIRST) CALL PLFINL
      STOP
      END
