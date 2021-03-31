      subroutine SITRED(in0,nstop,init)
 
      implicit none

c m.e.ash   march 1969    subroutine sitred
c observing site names and coordinates are initialized and read
c (geocentric, geodetic or cylindrical)

c parameters
      integer*4 in0,nstop
      logical   init
c
c modified for *command july 1978  r.b. goldstein

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'empcnd.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'stcord.inc'
      integer*2 k
      equivalence (Numsit,k)
c
c numsit      = number of observing sites (between 1 and 100)
c for k=1,...,numsit we have
c
c site(1-2,k) = 8 character observing site name
c
c kscrd(k)    =-1 cylindrical observing site coordinates relative to
c                 center of mass of earth
c kscrd(k)    = 0 spherical observing site coordinates relative to
c                 center of mass of earth
c kscrd(k)    = n.gt.0, spherical observing site coordinates relative
c               to center of mass of body with planet number n
c               (not programmed yet, program can only handle earth based
c               observing sites for now)
c all coordinates are in body fixed reference systems
c
c for spherical site coordinates
c scord(1,k)  = radius in kilometers
c scord(2,k)  = west longitude in degrees
c               (west of greenwich for geocentric)
c scord(3,k)  = north latitude in degrees
c
c for cylindrical site coordinates
c scord(1,k)  = radius in x-y plane in kilometers
c scord(2,k)  = west longitude in degrees
c               (west of greenwich for geocentric)
c scord(3,k)  = z distance in kilometers
c
c for cartesian site coordinates (for conversion to cylindrical)
c x direction is the 0,0 of latitude and longitude
c scord(1,k)  = x distance in kilometers
c scord(2,k)  = y distance in kilometers
c scord(3,k)  = z distance in kilometers
c
c for all sites
c scord(4,k)  = upward velocity (mm/yr)
c scord(5,k)  = westward velocity (mm/yr)
c scord(6,k)  = northward velocity (mm/yr)
c t0site(k)   = reference epoch for site coordinates: julian day number
c
c lscrd(i,k)  = 0 scord(i,j) not adjusted (i=1,6)
c lscrd(i,k)  = 1 scord(i,j) adjusted in least squares analysis (i=1,6)
c
c                 iv observing sites             (sitred)
c  for each site
c card 1a
c  columns
c    1- 8  site name (only first 4 char used operationally)  (1a8)
c    9-24  first  site coordinate                            (f16.9)
c            or  9-26   (f18.11) ##
c   25-40  second site coordinate                            (f16.9)
c            or 27-44   (f18.11) ##
c   41-56  third  site coordinate                            (f16.9)
c            or 45-62   (f18.11) ##
c      58  velocity flag: if '6', then another card follows  (1x,a1)
c            or    64   (1x,a1) ##
c   59-64  not used
c            or 65-70 ##
c   65-70  ls= 1  adjust each of three site coordinates      (3i2)
c          ls= 0  do not adjust each of three site coordinates
c            or 71-76   (3i2) ##
c   71-72  ks  type of input site coordinates                (i2)
c            or 77-78   (i2) ##
c          ks=  n  (n.gt.0)  spherical coords. on body  n
c          ks=  0            spherical coords. on earth
c          ks= -1            cylindrical coords. on earth
c          ks= -2            geodetic (spheroidal) coords. on earth
c                            with radius and flattening read in
c                            or setup before sitred called
c          ks= -3            geodetic (spheroidal) coords. on earth
c                            same as ks=-2 except
c                            radius and flattening individually input
c          ks= -4            XYZ coords. on earth (convert to sph)
c                            also applies to velocities, if any, which
c                            will then be converted to topocentric
c if the second site card is omitted, then the velocity is assumed to
c be zero, and the components are not to be adjusted
c ## if columns 79-80 are '##", use the alternate (new) format
c
c card 1b
c    1- 8  not used, but must be non-blank
c    9-24  upward site velocity                              (f16.9)
c            or  9-26   (f18.11) ##
c   25-40  westward site velocity                            (f16.9)
c            or 27-44   (f18.11) ##
c   41-56  northward site velocity                           (f16.9)
c            or 45-62   (f18.11) ##
c   57-64  reference epoch                                   (f8.0)
c            or 63-70   (f8.0) ##
c   65-70  ls= 1  adjust each of three site coordinates      (3i2)
c          ls= 0  do not adjust each of three site coordinates
c            or 71-76   (3i2) ##
c
c card 2     (for ks=-3 only)
c    1-19  earth radius   (kilometers)
c   20-39  flattening
c
c external functions
      real*10 DOT

c local variables
      real*10 eflat,erad,q2,q3,q4,qq,rc,rs
      integer i,j,ks,ns
      integer*2 ls(6)
      character*8 blnk8/'        '/,st8
      character*4 st(2)
      equivalence (st8,st)
      real*10 sc(6),ss(4),cc(4),latr,t0
      character*1 vflg
      real*10 jd2000/2451545._10/
      character*80 card
c
c standard sites, default values
      integer*4 numstd/36/
      integer*2 kscst(36)/0,0,0,0,0,0,0,0,-1,-1,-1,-1,-1,
     .          -1,-1,-1,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,
     .          0,0,0,0,0,0/
      real*10 stdcrd(3,36)
      real*10 stdcr1(3,18)/
     . 6368.551653028_10,  71.48866666666667_10,  42.431518383_10,
     . 6368.485020842_10,  71.266693300_10,       42.267879061_10,
     . 6368.475524712_10,  71.266093280_10,       42.268159595_10,
     . 6368.478870405_10,  71.265773370_10,       42.2683426048_10,
     . 6368.475794561_10,  71.265765270_10,       42.2682041262_10,
     . 6368.563831130_10,  71.49138888888889_10,  42.425660969_10,
     . 6376.560245971_10,  66.75302777777778_10,  18.2287613852_10,
     . 6372.1770000_10,   116.7940075_10,         35.066598100_10,
     . 5206.350322378_10, 116.84977459871_10,   3673.785176_10,
     . 5212.0508000_10,  -243.19463_10,         3665.646800_10,
     . 5203.9974000_10,  -243.11059_10,         3677.063000_10,
     . 5450.1978000_10,  -136.88759_10,        -3302.326200_10,
     . 5205.3610282_10,  -148.980957988_10,    -3674.612989_10,
     . 5742.9380000_10,   -27.68546_10,        -2768.719300_10,
     . 4862.6044000_10,  -355.75109_10,         4114.851800_10,
     . 4860.8114000_10,  -355.63229_10,         4116.966000_10,
     . 5391.8270000_10,   110.7244167_10,       3400.679000_10,
     . 6374.717549422_10, 104.02225_10,           30.503081043_10/
      real*10 stdcr2(3,18)/
     . 6378._10,           77.0660375_10,          0._10,
     . 6378._10,           77.06554166666667_10,   0._10,
     . 6378._10,           77.0654625_10,          0._10,
     . 6378._10,           77.06554166666667_10,   0._10,
     . 6378._10,          -18.47658333333333_10,   0._10,
     . 6378._10,            0.0_10,                0._10,
     . 6378._10,             -9.47916652E-2_10,      0._10,
     . 6378._10,            1.251666666666667_10,  0._10,
     . 6378._10,           75.716458320_10,        0._10,
     . 6378._10,           -2.3371249980_10,       0._10,
     . 6378._10,           -1.46249996_10,         0._10,
     . 6378._10,           -7.3004166650_10,       0._10,
     . 6378._10,           -5.98924996_10,         0._10,
     . 6378._10,           -4.358208333333333_10,  0._10,
     . 6378._10,         -139.54075_10,            0._10,
     . 6378._10,           -7.768333333333330_10,  0._10,
     . 6378._10,          -13.1066666666666_10,    0._10,
     . 6378._10,            3.18333333333333_10,   0._10/
 
      equivalence (stdcrd(1,1),stdcr1(1,1)), (stdcrd(1,19),stdcr2(1,1))
 
      character*8 stdnam(36) /
     .   'HAYSTACK', 'B30LNCLN', 'B10LNCLN', 'B04LNCLN', 'B01LNCLN',
     .   'MILLSTON', 'ARECIBO ', '85JPLVNS', '11DSPION', '12DSECHO',
     .   '14DSMARS', '41DSWOOM', '42DSCANB', '51DSJOHA', '61DSMADR',
     .   '62DSCEBR', 'AFLASER ', 'MCDONALD', '6USNAVAL', '8USNAVAL',
     .   '9USNAVAL', 'MUSNAVAL', 'CAPETOWN', 'GRENWICH', 'CAMBRIDG',
     .   'RADCLIFF', 'OTTAWA  ', 'PARIS   ', 'TOULOUSE', 'NICE    ',
     .   'BESANCON', 'UCCLE   ', 'GTOKYO  ', 'STRASBRG', 'BERLIN  ',
     .   'EDINBRG ' /
c
c initialize site data
      do j = 1,u_mxsit
         Sitd(j)  = blnk8
         Kscrd(j) = 0
         do i = 1,6
            Scord(i,j) = 0._10
            Lscrd(i,j) = 0
         end do
         T0site(j) = 0._10
      end do
c
c setup standard site names and coordinates
      do j = 1,numstd
         Kscrd(j) = kscst(j)
         do i = 1,3
            Scord(i,j) = stdcrd(i,j)
         end do
         Sitd(j) = stdnam(j)
      end do
      Numsit = numstd
 
      if(init) return
c
c*  start=100
c spool site cards from in to in0 with a-format printout
      call PEPTIC(In,Iout,in0,3,'SITE CARDS  ',nstop,1)
      do while(.true.)
c
c read site data
         read(in0,30) card
   30    format(a80)
         if(card(79:80).eq."##") then
            read(card,40) st8,(sc(i),i=1,3),vflg,(ls(i),i=1,3),ks
   40       format(a8,3f18.11,1x,a1,6x,4i2)
         else
            read(card,50) st8,(sc(i),i=1,3),vflg,(ls(i),i=1,3),ks
   50       format(a8,3f16.9,1x,a1,6x,4i2)
         endif
         if(st8.eq.blnk8) goto 200
         if(vflg.eq.'6') then
            if(card(79:80).eq."##") then
               read(in0,555) (sc(i),i=4,6),t0,(ls(i),i=4,6)
  555          format(8x,3f18.11,f8.0,3i2)
            else
               read(in0,565) (sc(i),i=4,6),t0,(ls(i),i=4,6)
  565          format(8x,3f16.9,f8.0,3i2)
            endif
         else
            t0=jd2000
            do i=4,6
               sc(i)=0._10
               ls(i)=0
            end do
         endif
c
c convert cartesian to spherical coordinates
c and XYZ velocity to up/west/north
         if(ks.eq.-4) then
            qq = SQRT(DOT(sc,sc))
            if(vflg.eq.'6') then
c use ss array as temporary for converted velocities and cc as direction
               ss(1)=DOT(sc,sc(4))/qq
c longitudinal direction
               rc=SQRT(sc(1)**2+sc(2)**2)
               cc(1)=sc(2)
               cc(2)=-sc(1)
               cc(3)=0._10
               ss(2)=DOT(cc,sc(4))/rc
c meridional direction
               cc(1)=-sc(1)*sc(3)/rc
               cc(2)=-sc(2)*sc(3)/rc
               cc(3)=rc
               ss(3)=DOT(cc,sc(4))/qq
               do i=4,6
                  sc(i)=ss(i-3)
               end do
            endif
            latr=ASIN(sc(3)/qq)/Convd
            sc(2)=-ATAN2(sc(2),sc(1))/Convd
            sc(1)=qq
            sc(3)=latr
            ks=0
         endif
c
c     convert input height above mean sea level sc(1) and spheroidal
c     latitude sc(3) to distance from center of earth and geocentric
c     latitude in the same storage locations
c
c           ks =positive spherical coordinates on body ks
c           ks =  0   geocentric spherical coordinates on earth
c           ks = -1   cylindrical coordinates on earth
c           ks = -2   geodetic (spheroidal) coordinates on earth with
c                     radius and flattening as input to program as whole
c           ks = -3 same as ks=-2 except that radius and flattening
c                     to be specially read in
         if(ks.lt.-1) then
            if(ks.ge.-2) then
               erad  = Econd(7)
               eflat = Econd(8)
            else
               read(in0,60) erad,eflat
   60          format(2F20.10)
            endif
c ss(j),j=1,4= geodetic flattening sine coefficients
c cc(j),j=1,4= geodetic flattening cosine coefficients
            qq    = -0.5_10*eflat
            ss(1) = 1._10+qq*(3._10-0.125_10*eflat*(5._10+1.5_10*eflat))
            cc(1) = 1._10-qq*(1._10+0.125_10*eflat*(5._10+3.5_10*eflat))
            ss(2) = qq*(1._10-eflat*(1._10+0.15625_10*eflat))
            cc(2) = qq*(1._10+eflat*(1._10+0.84375_10*eflat))
            qq    = qq*eflat
            ss(4) = qq*eflat*0.15625_10
            cc(4) = ss(4)
            qq    = -qq*0.375_10
            ss(3) = qq*(1._10-0.5_10*eflat)
            cc(3) = qq*(1._10+1.5_10*eflat)
            latr  = Convd*sc(3)
            q2    = COS(2._10*latr)
            q3    = COS(4._10*latr)
            q4    = COS(6._10*latr)
            rc    = (sc(1)*1E-3_10 + erad*(cc(1)+cc(2)*q2+cc(3)*q3+
     .              cc(4)*q4))*COS(latr)
            rs    = (sc(1)*1E-3_10 + erad*(ss(1)+ss(2)*q2+ss(3)*q3+
     .              ss(4)*q4))*SIN(latr)
            sc(1) = SQRT(rc**2 + rs**2)
            sc(3) = ATAN2(rs,rc)/Convd
            ks    = 0
         endif
c
c*  start=200
c check whether this site already set up
         if(Numsit.gt.0) then
            do i=1,Numsit
               if(st(1).eq.Site(1,i)) then
                  ns = i
                  goto 100
               endif
            end do
         endif
         if(Numsit.lt.u_mxsit) then
            Numsit = Numsit + 1
            Sitd(Numsit) = st8
 
c fill in new coordinates
            ns = Numsit
         else
 
c too many input sites, set error flag
            write(Iout,80) u_mxsit
   80       format('0 *** MORE THAN',i4,
     .             ' INPUT SITES, ERROR IN SITRED ***')
            goto 200
         endif
  100    do i = 1,6
            Scord(i,ns) = sc(i)
            Lscrd(i,ns) = ls(i)
         end do
         Kscrd(ns) = ks
         T0site(ns)= t0
      end do
c
c*  start=900
  200 rewind in0
      return
      end
