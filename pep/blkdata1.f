      block data BLKDATA1
 
      implicit none

c
c|blkdata1  m.e.ash   march 1976   setup constants for main link
c

c array dimensions
      include 'globdefs.inc'
 
      include 'funcon.inc'
      data Gauss/0.017202098950000_10/
      data Pi/3.14159265358979323846_10/,
     . Twopi/6.28318530717958647693_10/,
     . Pitwo/1.57079632679489661923_10/,
     . Convd/1.74532925199432957692E-2_10/,
     . Convds/4.84813681109535993590E-6_10/,
     . Convhs/7.27220521664303990385E-5_10/
c         gauss = gaussian gravitational constant
c         pi    = ratio between the circumference and diameter of a
c                   circle
c         twopi = 2*pi
c         pitwo = pi/2
c         convd = pi/180
c         convds= convd/3600
c         convhs= twopi/86400
c
      include 'namtim.inc'
      data Aplnt(-3)/' EMBARY '/
      data Aplnt(-2)/'  MOON  '/
      data Aplnt(-1)/' EROTAT '/
      data Aplnt(0) /' MROTAT '/
      data Aplnt(u_mxpl+1)/'  MOON  '/
      data Aplnt(u_mxpl+2)/'  SUN   '/
 
      include 'loadid.inc'
 
      data Lnkdat/'20170713'/
      data Lnklvl/' 790'/
      data Lnkdsn/'PEP.PEPLOAD.PEP790'/
c --- this to be updated as needed ---
c
      include 'zeroes.inc'
      data Izero/50*0/
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c   the length of some common blocks is 'hard-wired' into code -
c
c common       routine(s)             stand-alone routine(s)
c ------       ---------------------  ----------------------
c anctrl       anset
c bdctrl       prmred
c bdydta       prtcrd
c bernum       prdobs
c dtparm       prmred
c empcnd       bdyred                 addmoon
c eqnphs       cmpar2
c ethhar       bdyred
c france       bdyred
c ktrap        prdobs
c lcntrl       prmred
c leon         comset                 libang
c lfix         partl
c ltrap        cmpar1, comset         plastest
c mascon       prmred
c mnsprt       cmpar2                 plastest
c monhar       bdyred
c mtrap        cmpar1
c namtim                              addmoon, intrpp, orbdif, plastest
c                                     pmoon
c number       cmpar1, comset         plastest
c plndta                              pmoon
c plnhar       bdyred
c psrstf       bdyred
c scoef4       bdyred
c skymap       cmpar2
c sitcrd       cmpar1                 plastest
c smlbdy       bdyred
c stats        comset                 plastest
c tidal        fermtr
c tidal        radar
c wrkcom       bodred
c xprcom       morset, sbset
c
c   in addition, some routines have dimension limits from include
c   files with names beginning with 'max' - whenever an array size
c   is changed, the corresponding include, if any, should be also.
c
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c documentation of pep 770 object module produced 1990 november
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c        principal overlay links among the total of 86 overlay links
c
c        (a purist would say 'segment' instead of 'link'. however,
c        this program on the ibm 370 computer, now 120000 cards long,
c        grew out of a 6000 card long program on the ibm 7094
c        computer. the overlay terminology in the 7094 operating
c        system used 'link', whereas the terminology in os/360 is
c        'segment'. the word 'link' is imbedded at a number of
c        places in the program in printout formats.)
c
c  link   main    origin    purpose
c        program
c
c    1    main              always in core storage. controls iterative
c                           least squares process, provides service
c                           routines.
c
c    2    punch   subtnlnk  punch out input data at end of run, print
c                           total time elapsed.
c
c    3    julday  subtnlnk  subroutines which are needed in integration
c                           and compar segments of pep
c
c    6    input   comonlnk  initilize, read and printout input data
c                           occupies less than 150 k bytes
c
c   16    numint  integlnk  numerical intergration, perturbing planet,
c                           elliptic orbit routines brought into storage
c                           when planet,moon or rotat is called
c
c   17    planet  secndlnk  main link for n-body, earth-moon barycenter,
c                           planet-asteroid, or satellite-probe
c                           numerical integrations.
c
c   18    bodfn   plantlnk  perform n-body numerical integration
c                           (motion), write results on magnetic tape or
c                           disk.
c
c   21    fn      plantlnk  perform earth-moon barycenter and planet-
c                           asteroid numerical integrations (motion and
c                           partial derivatives), write results on
c                           magnetic tape or disk (one data set per
c                           body).
c
c   26    sbfn    satpblnk  perform satellite-probe numerical
c                           integrations (motion and partial
c                           derivatives), write results on magnetic tape
c                           or disk (one data set per body).
c
c   32    moon    secndlnk  perform moon numerical integration (motion
c                           and partial derivatives), write results
c                           on magnetic tape or disk.
c
c   36    rotat   secndlnk  dummy link for planet rotation numerical
c                           integrations
c
c   40    compar  extralnk  process observations of sun, moon, planets,
c                           asteroids, satellites, artificial space
c                           probes, stars
c                           (compute observed minus theoretical values
c                           and partial derivatives of observations,
c                           write results on magnetic tape or disk)
c
c   41    comset  comprlnk  setup for comparison of theory and obser-
c                           vation, read first two records of input
c                           observation library tapes
c
c   42    cmpar1  stobslnk  logic for setting up processing
c                           of observation series
c
c   43    cmpar2  cmparlnk  logic for setting up processing
c                           of observation series
c
c   44    bdrd1   cmparlnk  read routines for first five records of
c                           n-body,moon,earth-moon barycenter,planet-
c                           asteroid,satellite-probe data sets
c
c   45    stradr  stobslnk  process radar observations of one
c                           space craft by another
c
c   46    stoptc  stobslnk  process look angle observations of
c                           the sun or a space craft by a space craft
c
c   47    stmrdr  stobslnk  process time delay measurement: ground to
c                           spacecraft to another spacecraft to ground
c
c   48    ctutf   stobslnk  ephemeris time minus universal time routine
c
c   49    optic   opticlnk  process optical observations (meridian
c                           circle,photographic) of sun,moon,planets,
c                           satellites
c
c   54    trnsit  opticlnk  process transit and occultation observations
c
c   55    media   comprlnk  routines for interplanetary media and
c                           earth atmosphere and ionosphere
c
c   56    radar   rdfrmlnk  process radar observations of planets
c
c   68    fermtr  rdfrmlnk  process radio interferometer and
c                           pulsar observations
c
c   72    analiz  subtnlnk  labeled commons with right side,coefficient
c                           matrix and control integers for normal
c                           equations, and main routine for solution
c
c   73    alaliz  nrmeqlnk  main routine to form and solve normal
c                           equations to get least squares adjustments
c                           to parameters
c
c   77    nrmfrm  analzlnk  restore saved solution of normal equations
c                           or form normal equations from subset of
c                           saved normal equations. solve by inverting
c                           matrix with double precision arithmetic
c                           using fact that matrix is symmetric.
c                           (extended precision intermediate
c                           operations if desired)
c
c   81    nrmict  analzlnk  form normal equations from input
c                           observation library tapes
c
c   83    adjust  adjstlnk  adjust parameters with solutions to normal
c                           equations.
c
c   85    prterb  nrmeqlnk  replace ephemerides on perturbing planet
c                           tape with results of numerical integrations.
c                           (dummy, function served by n-body integ.)
c
c   86    prdict  nrmeqlnk  predict observed minus theory
c
c
c
      end
