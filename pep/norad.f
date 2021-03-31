      subroutine NORAD(jdq,ihr,imin,sec,cond,in0)
 
      implicit none
c
c     m.e.ash   may 1974   subroutine norad
c     read norad 2 card mean elements and convert to osculating
c     cartesian coordinates by applying kozai corrections
c     norad = north american air defense, satellite part of which is now
c     called the space defense center (sdc)
c
c arguments
      integer*4 jdq,in0
      real*10 sec,cond(6)
      integer*2 ihr,imin

c array dimensions
      include 'globdefs.inc'

c common 
      include 'crdbuf.inc'
      include 'ebloc.inc'
      real*10 thing(38)
      equivalence (Time,thing)
      include 'funcon.inc'
      include 'inodta.inc'
      include 'param.inc'
c
c format of space defense center 2 card classical mean elements
c card 1
c col   name  description                            units  field format
c 1     linno line number (always 1 for card 1)      none              x
c 2     blank
c 3-7   satno satellite number                       none          xxxxx
c 8     clasi element classification(u=unclassified, none              x
c             c=confidential, s=secret)
c 9     blank
c 10-11 idyr  international designator (last two     launch           xx
c             digits of launch year)                 yr
c 12-14 idln0 international designator (launch       none            xxx
c             number of the year
c 15-17 idpn0 international designator (piece of     none            xxx
c             launch)
c 18    blank
c 19-20 epyr  epoch year (last two digits of year)   epoch yr         xx
cc 21-32
c 21-32 epoch epoch (day and fractional days of the  utc    xxx.xxxxxxxx
c             year (jan 1 is day 1)                  days
c 33    blank
c 34-43 ndot2 first time derivative of the mean motion
c 34-43 ndot2 first time derivative of the mean      revolu- +-.xxxxxxxx
c        or   motion or ballistic coefficient        tions   (if ndot2
c       bterm (depending on ephemeris type)          /day**2 .gt.1, +
c                                                    or m**2 w/o sign)
c                                                    /kg
c 44    blank
c 45-52 ndot6 second time derivative of mean motion  revolu-   +-xxxxx-x
c             (field blank if ndot6 not applicable)  tions    (decimal
c                                                    /day**3  pt assumed
c                                                             between
c                                                             cols 45&46
c 53    blank
c 54-61 bstar bstar drag term if gp4 general pertur-           +-xxxxx-x
c        or   bations theory was used. otherwise will
c       agom  be the radiation pressure coefficient.
c 62    blank
c 63    ephtyp ephemeris type (specifies ephemeris   none              x
c             theory used to produce elements)
c 64    blank
c 65-68 elno  element number                         number         xxxx
c 69    cksum check sum (modulo 10)                  none              x
c card 2
c col   name  description                            units  field format
c 1     linno line number (always 2 for card 2)      none              x
c 2     blank
c 3-7   satno satellite number                       none          xxxxx
c 8     blank
c 9-16  ii    inclination                            degrees    xxx.xxxx
c 17    blank
c 18-25 node  right ascension of ascending node      degrees    xxx.xxxx
c 26    blank
c 27-33 ee    eccentricity (decimal pt assumed)      none       .xxxxxxx
c 34    blank
c 35-42 omega argument of perigee                    degrees    xxx.xxxx
c 43    blank
c 44-51 mm    mean anomaly                           degrees    xxx.xxxx
c 52    blank
c 53-63 nn    mean motion                            revolu- xx.xxxxxxxx
c                                                    tions
c                                                    /day
c 64-68 revno revolution number at epoch             revolutions   xxxxx
c 69    cksum check sum (modulo 10)                  none              x
c
c only epyr,epoch,ii,node,ee,omega,mm,nn are read
      integer*2 epyr,ione/1/,izer2/0/
      real*10 epoch,ii,node,ee,omega,mm,nn
      real*10 aukm,ctutc,fract
      integer   i,idum,jd0,JULDAY,LEG,ncard,nnb
      character*1 dollar/'$'/
c
c spool 2 card mean element set to disk
      nnb   = 0
      ncard = 0
      do while( .true. )
         call PEPTIN(In,Iout,idum)
   50    format(a80)
         if(nnb.gt.0) then
            write(Iout,60) Card80
   60       format(1x,a80)
         else
            nnb = 1
            write(Iout,80) Card80
   80       format(1x,a80,' MEAN ELEMENTS FOR ABOVE &NMLST2')
         endif
         if(LEG(1,1,Card8,1,dollar).ne.0) then
            write(in0,50) Card80
            ncard = ncard + 1
            if(ncard.ge.2) then
 
c end file in0
               rewind in0
c
c read 2 card mean element set
               read(in0,90) epyr,epoch
   90          format(18x,i2,f12.8)
               read(in0,100) ii,node,ee,omega,mm,nn
  100          format(8x,f8.4,1x,f8.4,1x,f7.7,1x,f8.4,1x,f8.4,
     .                1x,f11.8)
               rewind in0
c
c zero out common /ebloc/ for sgp
               do i = 1,38
                  thing(i) = 0.0_10
               end do
               Year = 0
c     time=tepoch=0. so sgp will apply kozai corrections at epoch only
c     and will not propigate orbit ahead in time
c     unless input jdq.gt.0 (see below)
c
c           determine julian day number and coordinate time of epoch
               jd0   = epoch
               fract = epoch - jd0
               jd0   = JULDAY(ione,izer2,epyr) + jd0
               ctutc = 45.1843817_10 + (epyr - 74)
 
c above formula for ct-utc subject to revision
               fract = fract + ctutc/8.64E4_10
               if(fract.ge.1._10) then
                  fract = fract - 1._10
                  jd0   = jd0 + 1
               endif
               if(jdq.gt.0) then
c
c propagate mean elements from epoch on card to epoch of
c numerical integration
                  Time = ihr*3600+imin*60
                  Time = (Time + sec)/8.64E4_10 - fract
                  Time = Time + (jdq - jd0)
               else
c
c epoch of integration is on mean element card
c no propagation of mean elements
                  fract = fract*8.64E4_10
                  ihr   = fract/3600._10
                  fract = fract - 3600._10*ihr
                  imin  = fract/60._10
                  sec   = fract - 60._10*imin
                  jdq   = jd0
               endif
c
c setup right units for sgp (radians and days)
               Eo     = ee
               Aio    = ii*Convd
               Omegao = omega*Convd
               Anodeo = node*Convd
               Amo    = mm*Convd
               Tno    = nn*Twopi
               Al     = Amo + Omegao + Anodeo
c
c apply kozai corrections to get osculating cartesian coord
               call SGP
c
c convert from km,km/sec to au,au/day
               aukm = prmter(51)*prmter(100)
               do i = 1,3
                  cond(i)     = R(i)/aukm
                  cond(i + 3) = Rdot(i)/aukm*8.64E4_10
               end do
 
               return
            endif
         endif
      end do
      end
