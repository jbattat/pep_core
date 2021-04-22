      subroutine COMSET
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 oblq, oblq0, oblqty, omeg, prcn, prcn0, qq, vmag
      integer   i, j
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   jan 1967     subroutine comset
c setup for comparison of theory and observation
c
 
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'bdctrl.inc'
      include 'comdateq.inc'
      include 'empcnd.inc'
      real*10 eflat
      equivalence (Econd(8),eflat)
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'leon.inc'
      integer*4 zleon/2208/   !r8=1136,r10=2208
      include 'ltrapx.inc'
      integer*4 zltrap/16238/   !r8=11502,r10=16238
      include 'namtim.inc'
      include 'obscrd.inc'
      include 'number.inc'
      include 'param.inc'
      include 'plndta.inc'
      include 'stats.inc'
      integer*4 zstats/1058/   !r8=666,r10=1058
      include 'timstf.inc'
c
c        secday   = 8.64E4_10 number of seconds in a day
c        secdy2   = 4.32E4_10 number of seconds in half a day
c        aultvl   = (astronomical unit in light seconds)/8.64E4_10
c        auacc    = aultvl/8.64E4_10
c        mnfact   = (mass of moon)/(mass of earth+moon) times the value
c                   of distance unit of moon ephemeris on earth-moon
c                   tape in light-seconds
c        mnvel    = mnfact/8.64E4_10
c        mnacc    = mnvel/4.32E4_10
c        mnfct    = (mass of moon)/(mass of earth+moon) times the value
c                   of distance unit of moon ephemeris on earth-moon
c                   tape in astronomical units
c        mnau     = moon distance unit in au
c        velipt(1-3) = constant component of earth velocity for
c                      elliptic aberration
c        aukm     = 1 au in km
c
c         lnotm = 1 if any l vector is on for which corresponding
c             m vector is off, else = 0
c
      logical*4 j2000
 
      j2000 = (Jct(13).gt.0)
c
c setup errtot constants
      do i = 17, 25
         Erquan(i) = 0.0_10
      end do
      Nit(1) = Npage
      Nit(2) = Npage
      do i = 3, 12
         Nit(i) = 0
      end do
      Nit(13) = Ireal0
      Nit(14) = Itotsk
      Nit(15) = Ireal0
      Nit(16) = Itotsk
c
c set up constants
      Aukm   = Aultsc*Ltvel
      Aultvl = Aultsc/Secday
      Auacc  = Aultvl/Secday
      if(Mdstsc.gt.0.0_10) then
         Mnfact = Mass(10)*Mdstsc
         Mnau   = Mdstsc/Aultsc
         Mnltsc = Mdstsc
      else
         Mnfact = Mass(10)*Aultsc
         Mnau   = 1.0_10
         Mnltsc = Aultsc
      endif
      Mnvel     = Mnfact/Secday
      Mnacc     = Mnvel/Secday
      Mnfct     = Mnfact/Aultsc
      Mfctwo    = Mnfct*2.0_10
      Mnfcte    = Mnau*(1.0_10 - Mass(10))
      omeg      = Convd*102.08055555555555556_10
      oblqty    = Convd*23.445833333333333333_10
      vmag      = Gauss*0.01673_10*Aultvl
      Velipt(1) = -vmag*SIN(omeg)
      Velipt(2) = vmag*COS(omeg)*COS(oblqty)
      Velipt(3) = vmag*COS(omeg)*SIN(oblqty)
      Lnotm     = 0
c
c scale factors for partial derivatives w.r.t. prmter(31-50)
      do i = 1, 20
         Hippo(i) = 1._10
      end do
      Hippo(1)  = (Gauss*Aultvl)**2
      Gmc2      = Hippo(1)*Aultsc
      Hippo(2)  = 365.25_10
      Hippo(3)  = (Sunrad/Aukm)**2
      Hippo(11) = Hippo(1)*Aultsc
      Hippo(12) = Hippo(11)
      Gmerth    = Gmc2*Mass(3)*(1._10 - Mass(10))
      Gmmoon    = Gmc2*Mass(3)*Mass(10)
c
c ss(j),j=1,4= geodetic flattening sine coefficients
c cc(j),j=1,4= geodetic flattening cosine coefficients
      qq    = -0.5_10*eflat
      Ss(1) = 1._10 + qq*(3._10 - 0.125_10*eflat*(5._10+1.5_10*eflat))
      Cc(1) = 1._10 - qq*(1._10 + 0.125_10*eflat*(5._10+3.5_10*eflat))
      Ss(2) = qq*(1._10 - eflat*(1._10+0.15625_10*eflat))
      Cc(2) = qq*(1._10 + eflat*(1._10+0.84375_10*eflat))
      qq    = qq*eflat
      Ss(4) = qq*eflat*0.15625_10
      Cc(4) = Ss(4)
      qq    = -qq*0.375_10
      Ss(3) = qq*(1._10 - 0.5_10*eflat)
      Cc(3) = qq*(1._10 + 1.5_10*eflat)
c
c   read con1(1-12) for earth, moon, earth rotation, and moon rotation.
c   (these may be read again in emrd1, mnrd1, or rtrd1 if tape data set
c   exists.)  pcom, con1 for a planet, is read in plrd1.  pcomc(sbcom),
c   con1 for a satellite or probe, is read in sbrd1.  sccom in common
c   scdta, con1 for a second probe, is read in scrd1.
      read(Iplcon) Ecom
      read(Iplcon) Mcom
      read(Iplcon) Ercom
      read(Iplcon) Mrcom
      rewind Iplcon
c
c zero all precession, nutation quantities
      call ZFILL(Dhprec,zleon)
c constants for precession matrix
c exact quantities for new iau precession constant at epoch 1950.0
c are coded in preces (r.w.king 10 jun 80)
      if(j2000) then
         prcn0 = 5029.0966_10*Convds
         oblq0 = 84381.448_10*Convds
         if(Jct(21).eq.1) then
            prcn0 = 5028.79695_10*Convds
         else if(Jct(21).eq.2) then
            prcn0 = 5028.796195_10*Convds
            oblq0 = 84381.406_10*Convds
         endif
      else
         if(Jct(21).lt.0) then
            prcn0 = 5026.75_10*Convds
            oblq0 = 8.440484E4_10*Convds
         else if(Jct(21).eq.0) then
            prcn0 = 5027.878_10*Convds
            oblq0 = 84404.85522_10*Convds
         else if(Jct(21).eq.1) then
            prcn0 = 5027.57835_10*Convds
            oblq0 = 84404.85522_10*Convds
         else
            prcn0 = 5027.690726_10*Convds
            oblq0 = 84404.824187_10*Convds
         endif
      endif

      prcn = prcn0
      if(Ercond(28).gt.0) prcn = Ercond(28)*Convds
      Coblq0 = COS(oblq0)
      Soblq0 = SIN(oblq0)
      oblq   = oblq0
      if(Ercond(29).gt.0) oblq = Ercond(29)*Convds
c     partials of moon observables w.r.t. prcn and oblq do not
c     effect of moon rotation, so must use nominal values if
c     a cassini-type m.o.d. ecliptic libration model used
c     (e.g., jpl llb-5).  if euler angles w.r.t. 1950.0 equatorial
c     system used, there is no problem.
c
      Fact1 = COS(oblq)
      Fact2 = SIN(oblq)
      Fact3 = prcn*Fact2
      Fact4 = prcn*Fact1
      Hxz1  = (Fact4 - prcn0*Coblq0)/Convds*0.5_10
      Hxz2  = (Fact3 - prcn0*Soblq0)/Convds
      Fact1 = Fact1*0.5_10*Convds
      Fact3 = Fact3*0.5_10*Convds
      Fact2 = Fact2*Convds
      Fact4 = Fact4*Convds
c
c setup nutation constant
      Nutn0 = 9.21_10*Convds
      Nutn  = Nutn0
      if(Ercond(30).gt.0) Nutn = Ercond(30)*Convds
c
c setup for small rotation matrix which multiplies precession
      if(Jct(29).le.0) then
         Kddprc = 1
         do i = 1, 6
            Dprang(i) = Ercond(i + 6)
         end do
      endif
 
c further initialize quantities
      call PRECES(0._10)
c
c           zero out initial times so that equations of motion will not
c           be reintegrated unless corresponding initial conditions
c           are being adjusted in least squares analysis (in this latter
c           case, the initial times and initial conditions are set equal
c           to the initial times and initial conditions on the corres-
c           ponding peripheral data sets in subroutines emrd1,mnrd1,
c           plrd1,sbrd1,bdrd1)
      do j = 1, 20
         Jdpl9(j - 4) = Jdpl0(j - 4)
         Jdpl0(j - 4) = 0
      end do
      Jdbdy9 = Jdbdy0
      Jdbdy0 = 0
 
      call ZFILL(Nsite1,2*50)
      call ZFILL(Stat,zstats)
      do i = 1, 4
         Ntaps(i) = -999999
      end do
      Ntaps(2) = 0
      Ntapsb   = -999999
      Klam = -999
c
c
c zero all "merged" l-vectors to start
      call ZFILL(Deriv,zltrap)
c
c write title page for comparison of theory and observation
      call OPRMSG('COMPARISON OF THEORY AND OBSERVATION', 9)
      call PAGSET('COMPARISON OF THEORY AND OBS', 7)
      call NEWPG
      write(Iout,100)
      write(Iout,200)
  100 format('-'/'-'/'-'/'-'//23x,
     .     ' **********    *********    **            **   ********** '
     .     , '    *********    ********** '/23x,
     .     '***********   ***********   ***          ***   ***********'
     .     , '   ***********   ***********'/23x,
     .     '**            **       **   ****        ****   **       **'
     .     , '   **       **   **       **'/23x,
     .     '**            **       **   ** **      ** **   **       **'
     .     , '   **       **   **       **'/23x,
     .     '**            **       **   **  **    **  **   **       **'
     .     , '   **       **   **       **'/23x,
     .     '**            **       **   **   **  **   **   **       **'
     .     , '   **       **   **       **'/23x,
     .     '**            **       **   **    ****    **   ***********'
     .     , '   ***********   ***********')
  200 format(23x,
     .     '**            **       **   **     **     **   ********** '
     .     , '   ***********   ********** '/23x,
     .     '**            **       **   **            **   **         '
     .     , '   **       **   **  **     '/23x,
     .     '**            **       **   **            **   **         '
     .     , '   **       **   **   **    '/23x,
     .     '**            **       **   **            **   **         '
     .     , '   **       **   **    **   '/23x,
     .     '**            **       **   **            **   **         '
     .     , '   **       **   **     **  '/23x,
     .     '***********   ***********   **            **   **         '
     .     , '   **       **   **      ** '/23x,
     .     ' **********    *********    **            **   **         '
     .     , '   **       **   **       **')
      call PLINK
 
      return
      end
