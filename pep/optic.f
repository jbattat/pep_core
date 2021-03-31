      subroutine OPTIC(mocpar)
 
      implicit none

c
c           m.e. ash  feb 1968  subroutine optic
c      ash/becker     oct 1967    subroutine optic
c     computation of theoretical value of right ascension and
c     declination for a meridian (transit) circle or photographic
c     observation of the sun, moon or a planet
c
c parameter
      integer*4 mocpar
c     mocpar = 0 compar called in midst of least squares iteration
c     read all three data sets iobcon, iobs, iabs1 which become in
c     subroutines of compar iiobcn, iiobs, iiabs1
c     mocpar = 1 compar called at end of least squares iteration to
c                calculate dummy observations
c     read only iobcon which become in subroutines of compar iiobcn

c array dimensions
      include 'globdefs.inc'
c
c        commons
      include 'comdateq.inc'
      include 'coord.inc'
      include 'empcnd.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'ltrapobs.inc'
      include 'number.inc'
      include 'obscrd.inc'
      real*4 acctim, accdst, accprc, dltqq
      real*10 fdev
      equivalence (acctim, Estf), (accdst, Estf(2)),
     .            (accprc, Estf(3)), (dltqq, Estf(10)),
     .            (fdev, Dstf(10))
c     ctime = centuries from 1900 (needed because only last two
c               digits of year are read in)
      include 'obstap.inc'
      include 'param.inc'
      include 'redobs.inc'
      include 'sitcrd.inc'
      include 'statsrad.inc'

c external functions
      real*10 DMS2S
      integer*4 JULDAY

c local variables
      real*10 alse, convf, cycf, cycle, dese, dra
      integer   i, ialm, idem, khrdg, klap2, llab1,
     .          llab2, nkodf, nlab1, nlab2, nrahr, ntype
      character*3 ialh, ided
      character*2 astrik(2)/'  ', '* '/
      character*16 typobs(4)/' AZIMUTH-ELEVAT.', ' MERIDIAN CIRCLE',
     .                       '  ASTROMETRIC   ', ' ASTROGRAPHIC   '/
      character*6 cepch(2)/'1950.0','2000.0'/
c
c savb(i,1), i=1,numsav,  saved precession-nutation,etc. from input
c                         observation library tape for radar, optical
c savb(i,2), i=1,numsav,  the same quantities read from iiobs
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
c           set up observation data set quantities
c           niobc,niobs,niabs1  indicate status of data sets  iiobc,
c           iiobs, iiabs1
c              -1  not to be read
c              0  ready to be read
c              1  already read
c           if  iiobc, iiobs, or iiabs1  are initially zero,  data set
c           is not read in any case.
c
      integer*2 iyr19
      character*16 rahr(2)
      character*4 obslab(8)
      character*8 obsl8(8)/' AZIMUTH',' ELEVATI','ON  DECL','INATION ',
     .          ' RIGHT A', 'SCENSION', '     HOU', 'R ANGLE '/
      equivalence (obsl8(1),obslab(1)), (obsl8(5),rahr(1))
      character*4 qsec(2)/'SEC',''''' '/
      character*8 qhrmn(2)/'HR MIN  ','DEG ''   '/
c*  start=100
c
c set up indices
      nrahr = 1
      if(Jct(39).gt.1) nrahr = 2
      nkodf = Ncodf - 4
      if(Ncodf.gt.20) nkodf = nkodf - 20
      if(Ncodf.eq.20) nkodf = -1
      if(Klanb.gt.0) then
         if(Ncp0.eq.3 .or. Ncp0.eq.10) then
            if(nkodf.le.0) nkodf = -1
         endif
      endif
 
c set up for measurement type labels
      if(nkodf.ge.0) then
         nlab1 = 5 + 4*nrahr
         llab1 = 4
         nlab2 = 6
         llab2 = 3
      else
         nlab1 = 1
         llab1 = 2
         nlab2 = 3
         llab2 = 3
      endif
c
c read data about observation series
      accdst = 1.0E-2
      if(nkodf.ne.0) accdst = acctim
      if(Klap.lt.u_mxpl+1) then
         D = Pcond(7, Klap)
         if(Ncodf.gt.20) then
            klap2 = Klan
            if(Klans1.gt.0) klap2 = Klans1
            if(klap2.gt.0 .and. klap2.le.u_mxpl) D = D - Pcond(7, klap2)
         endif
      else if(Klap.eq.u_mxpl+1) then
         D = Mcond(7)
      else
         D = Sunrad
      endif
      D     = 2._10*D/Aukm
      dltqq = 1.
      Dstf(5) = 0._10
      Sclplt  = 1._10 + Pnox
      cycle   = Secday
      if(nkodf.lt.0) cycle = 360._10
      cycf = cycle*1E-2_10
c
c determine receiving site quantities
      call ESHAPE(1)
      call REFRC1(1)
c
c*  start=300
c write first page heading
      ntype = nkodf + 2
      call OBSBEG(typobs(ntype), Sitf(1), obslab(nlab1), llab1,
     .            obslab(nlab2), llab2)
      khrdg = 1
      if(ntype.gt.2) then
         if(Ncodf.gt.20) then
            khrdg = 2
            if(Ncodf.ne.26) go to 100
         else if(Jct(39).lt.1) then
            if(Jct(39).lt.0) then
               write(Iout, 10) cepch(Jct(13)+1)
   10          format(
     .'0TOPOCENTRIC RIGHT ASCENSION-DECLINATION RATES (SEC/SEC AND "/SEC
     .) REFERRED TO THE MEAN EQUINOX AND EQUATOR OF ', a6)
               go to 300
            endif
         else if(Jct(39).eq.1) then
            go to 100
         else
            write(Iout, 20)
   20       format(
     .'0TOPOCENTRIC HOUR ANGLE-DECLINATION FOR DUMMY OBSERVATION PREDICT
     .IONS')
            go to 300
         endif
         write(Iout, 50) cepch(Jct(13)+1)
   50    format(
     .'0TOPOCENTRIC RIGHT ASCENSION-DECLINATION REFERRED TO THE MEAN EQU
     .INOX AND EQUATOR OF ', a6)
      endif
      go to 300
  100 write(Iout, 200)
  200 format(
     .'0TOPOCENTRIC RIGHT ASCENSION-DECLINATION REFERRED TO THE TRUE EQU
     .ATOR AND EQUINOX OF DATE')
  300 Line = Line + 2
      if(Ncodf.gt.20) then
         write(Iout, 350)
  350    format(' --- DIFFERENCES ARE (DELTA RA)(COS DEC) AND (DELTA ',
     .          'DEC)  -- BOTH IN " OF ARC')
         Line = Line + 1
      endif
      if(nkodf.lt.0) then
         Line = Line + 2
         if(fdev.le.(1._10-900E-10_10)) then
            write(Iout, 360)
  360       format('0RADIO FREQUENCY REFRACTION CORRECTION MADE TO ',
     .  'THEORETICAL VALUE OF ELEVATION BEFORE COMPARING WITH ',
     .  'OBSERVATION')
         else
            write(Iout, 380)
  380       format(
     .'0NO REFRACTION CORRECTION IN THEORETICAL VALUE OF ELEVATION, ASSU
     .MED REMOVED FROM DATA')
         endif
      endif
c*  start=600
c decide what kind of observed body
      if(Nplnt0.gt.0) then
         Kst  = 4
         Kst1 = 3
         if(Klanb.gt.0) then
c satellite/probe
c check for cis-lunar object
            if(Ncp0.eq.3 .or. Ncp0.eq.10) Kst1 = 2
 
c moon or any satellite
            Kst = 1
         else if(Nplnt0.eq.10) then
 
c moon
            Kst1 = 2
            Kst  = 1
         endif
      else
 
c sun (or star)
         Kst  = 2
         Kst1 = 1
         if(Nplnt0.eq.-4) Kst=3
      endif
      if(nkodf.lt.0) then
         write(Intern, 400)
  400    format(/
     .         '  GRNWCH  JULIAN  REC UT(SIGNAL) UT1-UTS AT-UTS  CT-AT'
     .         , 8x, 'AZIMUTH (DEGREES)', 14x, 'ELEVATION (DEGREES)'/
     .         '   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC',
     .         2x, 2(3x,'OBSERVED',4x,'ERROR',4x,'OBS-TH',3x))
      else
         write(Intern, 450) rahr(nrahr), qhrmn(khrdg),
     .                      (qsec(khrdg), i = 1, 3)
  450    format('0 GRNWCH   JULIAN   UT  TIME OF   AT-UT    CT-AT', 10x,
     .          a16, 20x, 'DECLINATION'/ '   DATE     DAY', 5x,
     .          'OBSERVATION', 22x, 'OBSERVED     ERROR    OBS-TH',
     .          6x, 'OBSERVED     ERROR   OBS-TH'/
     .          11x, 'NUMBER  HR MIN  SEC', 5x, 'SEC', 6x,
     .          'SEC', 4x, a8, a3, 5x, a3, 6x, a3, 6x,
     .          'DEG  ''   ''''', 5x, '''''', 6x, '''''', 4x,
     .          'CLMP LIMB OBS')
      endif
      call OBSBGT
      go to 700
c
c*  start=1000
c           observation cards are read
c
c           if err on iiobs, dummy read necessary
  500 write(Iout, 600) Iiobs
  600 format(' **** ERROR ON DATA SET IIOBS=', i2,
     .', OBSERVATION CARD RECORD SKIPPED IN OPTIC **********')
      Line = Line + 1
      read(Iiobs, 850)
      Nerrra = Nerrra + 1
      go to 800
  700 if(mocpar.gt.0 .or. Niobs.ne.0 .or. Iiobs.le.0) goto 1000
  800 if(nkodf.ge.0) then
c
c read right ascension (hr,min,sec) and declination (deg,','')
         read(Iiobs, 850, err=500) Ncodeb(1, 1), Ihrb(1, 1),
     .                               Iminb(1, 1), Secb(1, 1), ialh,
     .                               ialm, alse, Errorb(1, 1), ided,
     .                               idem, dese, Errorb(2, 1),
     .                               Atutsb(1, 1), Clampb(1, 1),
     .                               (Limbb(i,1), i = 1, 2),
     .                               Obsrvb(1, 1), Imnthb(1, 1),
     .                               Idayb(1, 1), Iyearb(1, 1)
  850    format(i1, 2I3, f8.4, 1x, a3, i3, f8.4, f7.4, a3, i3, 2F7.3,
     .          f8.4, 1x, 3A1, a2, 1x, 3I2)
      else
c
c read azimuth and elevation (degrees)
         read(Iiobs, 900, err=500) Ncodeb(1, 1), Ihrb(1, 1),
     .                               Iminb(1, 1), Secb(1, 1),
     .                               (Resltb(i,1), Errorb(i,1), i = 1,
     .                               2), Atutsb(1, 1), Ututsb(1, 1),
     .                               Imnthb(1, 1), Idayb(1, 1),
     .                               Iyearb(1, 1)
  900    format(i1, 2I3, f8.4, f13.7, f7.4, f15.7, f7.4, f8.4, f7.4,
     .          3I2)
      endif
c
c calculate quantities from card
      Niobs = -1
      if(Ncodeb(1,1).gt.0) then
         Niobs = 1
         Fdsb(1, 1) = Ihrb(1, 1)*3600._10 + Iminb(1, 1)
     .                *60._10 + Secb(1, 1)
         iyr19 = Iyearb(1,1) + ctime*100
         Jdsb(1, 1) = JULDAY(Imnthb(1,1), Idayb(1,1), iyr19)
         Jdb(1, 1)  = Jdsb(1, 1)
         if(nkodf.ge.0) then
            Resltb(1, 1) = DMS2S(ialh, ialm, alse)
            Resltb(2, 1) = DMS2S(ided, idem, dese)
         endif
      endif
c
c*  start=1500
c read observation library tape
 1000 call OBSRED(mocpar)
      if(Nice.lt.-1) then
c
c*  start=9000
c printout at end of observation series
         call OBSEND
 
         return
      else
c
c calculation of theoretical value of right ascension and/or
c declination  or azimuth-elevation as the case may be
         call ANGCTL(nkodf)
         if(Jd.gt.-10) then
            if(Jd.gt.0) then
c
c modify error of first observable if not from obs.lib.tape
               if(Idumob.ne.-2) then
                  if(nkodf.le.0 .or. Jct(39).ge.0) then
                     if(Ncodf.le.20) then
                        convf = Convd
                        if(nkodf.ge.0) convf = Convds
                        dltqq = COS(Deriv(2,2)*convf)
                        Deriv(1, 1) = Deriv(1, 1)/dltqq
                        Error(1)    = Error(1)/dltqq
                     endif
                  endif
               endif
c
c idumob=1  dummy mode            idumob=-1  not dummy mode
               if(Idumob.eq.1) then
 
c for dummy mode only
                  Result(1) = Deriv(2, 1)
                  Result(2) = Deriv(2, 2)
               endif
c
c see if conversion is required for observable
               if(Idumob.ne.-1) then
                  if(nkodf.ge.0) then
 
c convert result(1) to ialh, ialm, alse
                     if(Nice.le.0)
     .                   call S2DMS(ialh, ialm, alse, Result(1), 4)
 
c convert result(2) to ided, idem, dese
                     if(Nice.ge.0)
     .                   call S2DMS(ided, idem, dese, Result(2), 3)
                  endif
               endif
            endif
c
c test for ra residual of 24 hrs (real obs only)
            if(Nice.le.0 .and. Idumob.ne.1) then
               if(nkodf.le.0 .or. Jct(39).ge.0) then
                  dra = Result(1) - Deriv(2, 1)
                  if(ABS(ABS(dra)-cycle) .le. cycf) Deriv(2, 1)
     .                = Deriv(2, 1) + SIGN(cycle, dra)
               endif
            endif
 
            call OBSCNT
            if(Jd.gt.0) then
               if(Ict(2).eq.0 .or. mocpar.gt.0) then
c
c*  start=3000
c printout observed minus theory
                  if(Nice.lt.0) then
 
                     if(nkodf.lt.0) then
                        write(Iout, 1010) Imonth, Iday, Iyear, Jds,
     .                        Ihr, Imin, Sec, Ututs, Atuts, Ctat,
     .                        Result(1), (Deriv(i,1), i = 1, 2),
     .                        astrik(Nast11)
 1010                   format(i3, '/', i2, '/', i2, i8, 2I3, 4F8.4,
     .                         2(f12.5,2F10.5,a1))
                     else
                        write(Iout, 1020) Imonth, Iday, Iyear, Jds,
     .                        Ihr, Imin, Sec, Atuts, Ctat, ialh, ialm,
     .                        alse, (Deriv(i,1), i = 1, 2),
     .                        astrik(Nast11), Clamp, Limb, Observ
 1020                   format(i3, 2('/',i2), i8, i4, i3, 2F8.4, f9.4,
     .                         1x, a3, i3, 2F8.4, f10.5, a1, 33x, a1,
     .                         3x, 1A1, 1x, 1A1, 3x, 1A2)
                     endif
                  else if(Nice.eq.0) then
 
                     if(nkodf.lt.0) then
                        write(Iout, 1010) Imonth, Iday, Iyear, Jds,
     .                        Ihr, Imin, Sec, Ututs, Atuts, Ctat,
     .                        Result(1), (Deriv(i,1), i = 1, 2),
     .                        astrik(Nast11), Result(2),
     .                        (Deriv(i,2), i = 1, 2), astrik(Nast12)
                     else
                        write(Iout, 1030) Imonth, Iday, Iyear, Jds,
     .                        Ihr, Imin, Sec, Atuts, Ctat, ialh, ialm,
     .                        alse, (Deriv(i,1), i = 1, 2),
     .                        astrik(Nast11), ided, idem, dese,
     .                        (Deriv(i,2), i = 1, 2), astrik(Nast12),
     .                        Clamp, Limb, Observ
 1030                   format(i3, 2('/',i2), i8, i4, i3, 2F8.4, f9.4,
     .                         1x, a3, i3, 2F8.4, f10.5, a1, 2x, a3,
     .                         i3, 2F7.3, f9.4, 1A1, 1x, a1, 3x, 1A1,
     .                         1x, 1A1, 3x, 1A2)
                     endif
 
                  else if(nkodf.lt.0) then
                     write(Iout, 1040) Imonth, Iday, Iyear, Jds, Ihr,
     .                                 Imin, Sec, Ututs, Atuts, Ctat,
     .                                 Result(2),
     .                                 (Deriv(i,2), i = 1, 2),
     .                                 astrik(Nast12)
 1040                format(i3, '/', i2, '/', i2, i8, 2I3, 4F8.4, 32x,
     .                      f12.5, 2F10.5, a1)
                  else
                     write(Iout, 1050) Imonth, Iday, Iyear, Jds, Ihr,
     .                                 Imin, Sec, Atuts, Ctat, ided,
     .                                 idem, dese,
     .                                 (Deriv(i,2), i = 1, 2),
     .                                 astrik(Nast12), Clamp, Limb,
     .                                 Observ
 1050                format(i3, 2('/',i2), i8, i4, i3, 2F8.4, f9.4,
     .                      36x, a3, i3, 2F7.3, f9.4, 1A1, 1x, a1, 3x,
     .                      1A1, 1x, 1A1, 3x, 1A2)
                  endif
c
c*  start=4000
                  Line = Line + 1
               endif
c
c calculate partial derivatives w.r.t. quantities to be adj.
               if(Ict(1).lt.0) then
               else if(Ict(1).eq.0) then
c
c write(obs-th),partials buffer
                  call COMRIT(0)
               else if(Idumob.ne.1 .or. Ict(3).ge.0) then
                  call PARTL(2)
                  call COMRIT(0)
               endif
            endif
         endif
         go to 700
      endif
      end
