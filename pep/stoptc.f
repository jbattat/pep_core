      subroutine STOPTC(mocpar)
 
      implicit none

c     m.e.ash   june 1970   subroutine stoptc
c     process angle observations from body other than earth

c arguments
      integer*4 mocpar
c     mocpar = 0 compar called in midst of least squares iteration
c     read all three data sets iobcon, iobs, iabs1 which become in
c     subroutines of compar iiobcn, iiobs, iiabs1
c     mocpar = 1 compar called at end of least squares iteration to
c                calculate dummy observations
c     read only iobcon which become in subroutines of compar iiobcn
c
c array dimensions
      include 'globdefs.inc'

c        commons
      include 'comdat.inc'
      include 'coord.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'ltrapx.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'obscrd.inc'
      real*4 dltqq
      equivalence (Estf(10),dltqq)
      real*10 ctrecf, fdev, tmdly0, dop0
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),
     . (Result,tmdly0),(Result(2),dop0)
      include 'obstap.inc'
      include 'param.inc'
      include 'redobs.inc'
      include 'sitcrd.inc'
      character*4    sit1
      equivalence (Sitf,sit1)
      include 'statsrad.inc'
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
      real*10 DOT
      real*4 STORNE
      integer*4 JULDAY

c local variables
      character*2 astrik(2)/'  ', '* '/
      character*4 nmles8/'LES8'/
      character*12 colhd1(2,2)/' PITCH ANGLE','  ROLL ANGLE',
     .                         'RT.ASCENSION',' DECLINATION'/
      real*10 xsun(13),psun,de(3),dde(3),angsun,saz,sel,secbx
      character*1 umbra(3)/' ', 'P', 'U'/
      character*1 confus(3)/' ', '*', '#'/
      integer*4 i,icnfus,iumbra,ntype1
      integer*2 iyr19
c
c initialization
      ntype1 = 1
      if(Ncodf.gt.4) ntype1 = 2
      Ncod(1) = 1
      Ncod(2) = 2
c
c determine receiving site quantities
      call ESHAPE(1)
c
c write first page heading
      call OBSBEG(' E-T LOOK ANGLE ', Aplnt(Klans1),colhd1(1,ntype1),
     .            3, colhd1(2,ntype1),3)
      if(Nplnt0.ne.10 .or. (Jct(50).le.0)) then
         write(Intern,50) (colhd1(i,ntype1),i = 1,2)
   50    format(
     .         '  GRNWCH  JULIAN  REC UTC (WWV)  UT1-UTC AT-UTC  CT-AT'
     .         , 6x, A12, ' (DEGREES)', 10x, A12, ' (DEGREES)    '/
     .         '   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC',
     .         5x, 'OBSERVED     ERROR     OBS-TH', 4x,
     .         'OBSERVED     ERROR     OBS-TH')
      else
         if(sit1.ne.nmles8) then
            write(Intern,660)
  660       format(
     .'0IDEAL REFERENCE FRAME, PITCH TO NORTH OF ORBIT PLANE, ROLL ALONG
     .POSITIVE VELOCITY')
         else
            write(Intern,640)
  640       format(
     .'0UPSIDE DOWN REFERENCE FRAME, PITCH TO SOUTH OF ORBIT PLANE, ROLL
     .ALONG NEGATIVE VELOCITY')
         endif
         write(Intern,680)
  680    format('0 GRNWCH  JULIAN REC UTC(WWV)', 3x, 'MOON', 9x,
     .    'MOON', 8x, 'MOON', 4x, 'MOON', 3x, 'MOON', 3x,
     .    'MOON', 7x, 'SUN', 9x, 'SUN', 4x, 'ECLI-', 1x,
     .    'CON-'/'   DATE   DAY NUM ', 'HR MN SEC', 4x,
     .    'DIST(KM)', 3x, 'AZIMUTH(D)    ELEV(D)', 1x,
     .    'ILLUM', 3x, 'PHASE', 2x, 'DIAM', 4x,
     .    'AZIMUTH(D)    ELEV(D)', 2x, 'PSE', 2x, 'FUSE')
      endif
      call OBSBGT
      goto 300
c
c read observation data card
c
c if err on iiobs, dummy read necessary
  100 write(Iout,200) Iiobs
  200 format(' **** ERROR ON DATA SET IIOBS=', i2,
     . ', OBSERVATION CARD RECORD SKIPPED IN STOPTC *********')
      Line = Line + 1
      read(Iiobs,500)
      Nerrra = Nerrra + 1
      goto 400
  300 if(mocpar.gt.0) goto 600
      if(Niobs.ne.0 .or. Iiobs.le.0) goto 600
  400 read(Iiobs,500,err=100) Ncodeb(1,1),Ihrb(1,1),Iminb(1,1),secbx,
     . Resltb(1,1),Errorb(1,1),Resltb(2,1),Errorb(2,1),
     . Atutsb(1,1),Ututsb(1,1),Imnthb(1,1),Idayb(1,1),Iyearb(1,1)
  500 format(i1,2I3,f8.4,f13.7,e7.0,f15.8,e7.0,f8.4,f7.4,3I2)
      Niobs = -1
      Secb(1,1) = STORNE(secbx)
      if(Ncodeb(1,1).gt.0) then
         Fdsb(1,1) = Ihrb(1,1)*3600._10 + Iminb(1,1)*60._10 + Secb(1,1)
         iyr19=Iyearb(1,1)+ctime*100
         Jdsb(1,1) = JULDAY(Imnthb(1,1),Idayb(1,1),iyr19)
         Jdb(1,1)  = Jdsb(1,1)
         Niobs = 1
      endif
  600 call OBSRED(mocpar)
      if(Nice.lt.-1) then
c
c printout at end of observation series
         call OBSEND
 
         return
      endif
c
c calculation of theoretical value of look angles
      call STANGL(xsun)
      if(Jd.le.-10) goto 300
      if(Jd.gt.0) then
 
c idumob=1  dummy mode            idumob=-1  not dummy mode
         if(Idumob.eq.1) then
 
c for dummy mode only
            Result(1) = Deriv(2,1)
            Result(2) = Deriv(2,2)
         endif
c
c change pitch error if obs.not from obs.library tape
         if(Idumob.ne.-2) then
            dltqq = Deriv(2,2)*Convd
            dltqq = cos(dltqq)
            Deriv(1,1) = Deriv(1,1)/dltqq
            Error(1)    = Error(1)/dltqq
         endif
      endif
 
      call OBSCNT
      if(Jd.le.0) goto 300
c
c special sun-moon mode
c     xsun(1-3)=position of satellite relative to sun with aberration
c                   corrections
c     xsun(7)=r**2
c     xsun(8)=r
c     xsun(9)=pitch angle or azimuth positive from yaw to roll
c     xsun(10)=roll angle or elevation.
c     these angles in ideal satellite reference frame pitch to north,
c     roll along satellite velocity vector
c     xsun(11)=fraction of moon's disk illuminated by sun
c              as seen from satellite
c     xsun(12)=moon phase angle as seen from satellite
c     xsun(13)=moon angular diameter as seen from satellite
c
      if(Nplnt0.eq.10 .and. (Jct(50).gt.0)) then
c           les-8 is upside down
         if(sit1.eq.nmles8) then
            Result(1) = -Result(1)
            Result(2) = -Result(2)
            xsun(9)   = -xsun(9)
            xsun(10)  = -xsun(10)
         endif
c
c is satellite in eclipse
         iumbra = 1
         psun   = -DOT(xsun,Xscsun)
         if(psun.lt.0._10) then
            do i = 1, 3
               de(i)  = -xsun(i)*psun/xsun(7)
               dde(i) = Xscsun(i,1) - de(i)
            end do
            psun   = SQRT(DOT(dde,dde))*Ltvel
            angsun = 6.95E5_10/xsun(8)*SQRT(DOT(de,de))
            if(psun.le.6378.16_10 + angsun) then
               iumbra = 2
               if(psun.le.6378.16_10 - angsun) iumbra = 3
            endif
         endif
c
c will moon confuse attitude control sensors
         icnfus = 1
         sel    = ABS(Result(2))
         saz    = ABS(Result(1))
         if((sel.ge.7.4_10) .and. (sel.le.14.3_10)) then
            if(saz.lt.5._10) icnfus = 2
         else if(sel.le.1.3_10) then
            if(saz.gt.5._10 .and. saz.lt.16._10) icnfus = 3
         endif
c
c write out moon and sun angles
         write(Iout,610) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,
     .    Rsitp(2),Result(1),Result(2),xsun(11),xsun(12),xsun(13),
     .    xsun(9),xsun(10),umbra(iumbra),confus(icnfus)
  610    format(i3,'/', i2, '/', i2, i8, i3, ':', i2, ':',
     .    f3.0, f12.1, f13.5, f11.5, f6.2, f8.2, f7.4,
     .    f13.5, f11.5, 3x, a1, 4x, a1)
         Line = Line + 1
c
c printout observed minus theory
      else if(Ict(2).eq.0 .or. mocpar.gt.0) then
         if(Nice.lt.0) then
 
c pitch only
            write(Iout,620) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,Ututs,
     .       Atuts,Ctat,Result(1),(Deriv(i,1),i = 1,2),astrik(Nast11)
         else if(Nice.eq.0) then
 
c pitch and roll
            write(Iout,620) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,Ututs,
     .       Atuts,Ctat,Result(1),(Deriv(i,1),i = 1,2),astrik(Nast11),
     .       Result(2),(Deriv(i,2),i = 1,2),astrik(Nast12)
  620       format(i3,'/', i2, '/', i2, i8, 2I3, 4F8.4,
     .       2(f12.5,2F10.5,a1))
         else
 
c roll only
            write(Iout,630) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,Ututs,
     .       Atuts,Ctat,Result(2),(Deriv(i,2),i = 1,2),astrik(Nast12)
  630       format(i3,'/', i2, '/', i2, i8, 2I3, 4F8.4, 33x,
     .       f12.5, 2F10.5, a1)
         endif
 
c increment line counter after printing observation results
         Line = Line + 1
      endif
c
c calculate partial derivatives w.r.t. quantities to be adj.
      if(Ict(1).lt.0) then
      else if(Ict(1).eq.0) then
c
c write(obs-th),partials buffer
         call COMRIT(0)
      else if(Idumob.ne.1 .or. (Ict(3).ge.0)) then
         call PARTL(2)
         call COMRIT(0)
      endif
 
      goto 300
      end
