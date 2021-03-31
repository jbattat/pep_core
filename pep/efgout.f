      subroutine EFGOUT(sidtm1)
 
      implicit none
c
c     m.ash   march 1974   subroutine efgout
c     write output efg satellite coordinates for cloudcroft 48 inch
c     optical observatory (+ others)
c     binary efg output for least squares fit of polynomials in efgfit
c     option to output angles instead of e,f,g
c     called from mdeldp. entry points efg1,efg2 called from radar.
c
c calling parameter
      real*10 sidtm1
c     sidtm1 = greenwich true sideral time in radians at receive time
c              epoch. satellite at receive time if delit=1.e6, otherwise
c              at retarded time.

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdat.inc'
      character*8 satnam
      equivalence (Comcon(127),satnam)
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      include 'param.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'
c
c jct(30) = 0 no e,f,g (+edot,fdot,gdot) tape output
c             earth centered rotating coordinates. e in true equator of
c             date towards greenwich meridian, g to north, f completes
c             right hand system
c jct(30).ne.0 such tape output on data set iabs(jct(30))
c              unless jct(40) says there is angle output instead
c        .lt.0 no frintout of e,f,g
c        .gt.0 printout of e,f,g at end of dummy observation series
c jct(31) = space defense center satellite number for e,f,g output
c jct(32) = 0 no binary e,f,g output in subroutine efgout for fitting
c             polynomials in subroutine efgfit
c jct(32).ne.0 such binary output on data set iabs(jct(32)) if
c              jct(30).ne.0
c        .gt.0 jct(30) output with jct(32) output
c        .lt.0 jct(30) output is supressed
c
c jct(40) = 0 no angle output for optical observatories
c jct(40).ne.0 such tape output on data set iabs(jct(30))
c              in format compatable with data general nova computer
c jct(40).gt.0 printout of angles, data set jct(40) used as intermediate
c              ebcdic buffer
c jct(40).le.0 no printout of angles
c jct(41) = 0 angle output referred to true equinox and equator of date
c             for jct(40) output
c jct(41).gt.0 angle output referred to mean equinox and equator of
c              jct(41). for eample, jct(41)=76 denotes output referred
c              to mean equinox and equator of 1976.0
c jct(42) = 0 no refraction correction made in computing jct(40) angles
c jct(42) = 1 refraction correction made in computing jct(42) angles
c jct(43) =   file number for jct(40) angle output (0 is file 1, 1 is
c             file 2, etc, for data general nova control word at end of
c             514 character record)
c jct(44) = 0 no punch card output for haystack
c jct(44).gt.0 80 character card image output for haystack on data set
c              jct(44) (punch if 7) if jct(40).gt.0
c
c jct(46) = 0 jct(40) print angle output has range in meters and
c             range rate in cm/min
c jct(46) = 1 jct(40) print angle output has one way time delay in
c             seconds and one way doppler shift in hz
c
c jct(47) = 0 usual jct(40) print and tape output
c jct(47) > 0 jct(40) tape output has satellite radius instead
c             of doppler rate, print output has special les-8/9
c             format with 6 frequency doppler
c jct(47) = 8 les-8 output
c jct(47) = 9 les-9 output
c
c jct(48) = 0 right ascension and hour angle rates are seconds
c             of time per second in jct(40) print output
c jct(48) = 1 right ascension and hour angle rates are seconds
c             of arc per second in jct(40) print output
c
c     xplsc     = position and velocity of earth satellite, lunar
c                 orbiter or spot on moon relative to center of earth
c                 referred to mean equinox and equator of the
c                 reference epoch at retarded
c                 time (units are light seconds and light seconds per
c                 second) (calculated in subroutine radmn)
c
c local variables
      real*10 az,ctim,dd8,ddecl,decls,dha,dra,el,has,hax,
     .          oneddp,onedel,onedop,ras,sec8,stim,
     .          subln8,sublt8
      integer   i,i89,icross,iefg,ifrst8,ihr8,isat,isec8,
     .          iumbra,j,jang,jct44,jdy8,jdyr,jdyr8,jefg,jjj,
     .          jmn8,jmon,jyr8,k,min8,nerr,nrec,nref
      character*4 aa(33),bb(32)
      character*1 echar/'E'/,iaa(132)
      equivalence (iaa, aa)
      real*10 efg(6), xsat(3), vsat(3), xrot(3), vrot(3)
      integer*2 angbuf(257), iangbf
      real*10 bessel(3, 3)
      character*2 head(51),minus(5)
      integer*4 jjj1(23, 5)
      character*8
     1    refct(2,2)/' NO REFR','ACTION  ',' REFRACT','ION COR '/
      character*4 iblank/'    '/,iminus/'----'/, ihsgn
      integer*2 imon(12)/0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
     .          304, 334/
      real*10 freq89(6,2)/37.98E9_10, 38.04E9_10, 37.4400E9_10,
     . 2.24E9_10, 248.86E6_10, 302.7E6_10, 36.90E9_10, 36.84E9_10,
     . 37.4400E9_10, 2.25E9_10, 249.36E6_10, 303.4E6_10/
      real*10 radius, wlong, nlat, dopler(6), twodel
      character*1 umbra(3)/' ','U','P'/
      character*2 minus1/'- '/
      character*8
     1    ascnod(2,2)/'ASCENDIN','G NODE  ','DESCENDI','NG NODE '/
c
c compute quantities for angle tape
      if(Jct(40).eq.0) then
c
c
c transform from mean equinox and equator of reference epoch to
c true equinox and equator of date
         call CORCHN(xsat, Xplsc(1))
         call CORCHN(vsat, Xplsc(4))
c
c transform from true equinox and equator of date to
c rotating efg system
         stim    = SIN(sidtm1)
         ctim    = COS(sidtm1)
         xrot(1) = xsat(1)*ctim + xsat(2)*stim
         xrot(2) = -xsat(1)*stim + xsat(2)*ctim
         xrot(3) = xsat(3)
         vrot(1) = vsat(1)*ctim + vsat(2)*stim + xrot(2)*Sidvel
         vrot(2) = -vsat(1)*stim + vsat(2)*ctim - xrot(1)*Sidvel
         vrot(3) = vsat(3)
c     add derivative of precession-nutation to velocity formula (later)
c     (if statement on kind)
c
c           include effect of wobble
c     if(iwob.le.0) go to 71
c
c           change units to km, km/sec
         do i = 1, 3
            efg(i)     = xrot(i)*Ltvel
            efg(i + 3) = vrot(i)*Ltvel
         end do
c
c write out e,f,g,edot,fdot,gdot
         if(jefg.ne.0) then
            write(jefg) isat, satnam, Jds, Iyear, Imonth, Iday, Ihr,
     .                  Imin, Sec, efg
            if(Jct(32).lt.0) go to 100
         endif
         write(Intern, 50) isat, satnam, Iyear, Imonth, Iday, Ihr,
     .                     Imin, Sec, efg
   50    format(i5, 1x, a8, 5I3, f7.3, 1p, 6D16.9)
         rewind Intern
         read(Intern, 1300) aa
         rewind Intern
         do i = 49, 129, 16
            iaa(i) = echar
         end do
         write(iefg, 1300) aa
      else
         call ANGOUT(iangbf, angbuf, bessel, sidtm1)
         if(iangbf.ge.5) then
            if(Jct(40).gt.0) write(jang, 60) angbuf
            write(iefg, 60) angbuf
   60       format(255A2, 2A2)
            iangbf = 0
         endif
      endif
 
  100 return
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c initialization, write header record on tape
      entry EFG1
      iefg = Jct(30)
      iefg = iabs(iefg)
      isat = Jct(31)
      if(Jct(40).eq.0) then
c
c for e,f,g output
         jefg = Jct(32)
         jefg = iabs(jefg)
         if(Jct(32).ge.0) then
            write(iefg, 120) isat, satnam, Heding, Date, Sitf(1)
  120       format(i5, 1x, a8, 1x, 'E,F,G,EDOT,FDOT,GDOT', 1x, 18A4,
     .             1x, 2A4, ' (', a8, ')')
         endif
      else
c
c for angle output
         nref = 1
         if(Jct(42).gt.0) nref = 2
         if(Jct(41).gt.0) then
            write(Intern, 140) isat, satnam, Sitf(1), 1900+Jct(41),
     .                         (refct(i,nref), i = 1, 2), Date
  140       format(i5, 1x, a8, 1x, a8,
     .             ' UTC,R.A.,DECL,H.A.,AZ,EL,RANGE,RDOT,RDOT2', 1x,
     .             'MEAN ', i4, '.0', 1x, 2A8, 2A4)
         else
            write(Intern, 160) isat, satnam, Sitf(1),
     .                         (refct(i,nref), i = 1, 2), Date
  160       format(i5, 1x, a8, 1x, a8,
     .             ' UTC,R.A.,DECL,H.A.,AZ,EL,RANGE,RDOT,RDOT2', 1x,
     .             'TRUE OF DATE', 2A8, 2A4)
         endif
         rewind Intern
         read(Intern, 200) (angbuf(i), i = 1, 52)
  200    format(51A2)
         rewind Intern
         iangbf = 1
         angbuf(256) = Jct(43)
         angbuf(257) = Jct(43)
         jang = Jct(40)
         if(Jct(44).gt.0) then
            jct44 = Jct(44)
            write(jct44, 220) isat, satnam, Sitf(1), Date
  220       format(i5, 1x, a8, 1x, a8,
     .             ' ADD 46.184 SEC TO UTC TO GET ET', 3x, 2A4)
         endif
         if(Jct(41).gt.0) then
            call PRECES(2442778.721_10 + (Jct(41)-76)*365.2421988_10)
            do i = 1, 3
               do j = 1, 3
                  bessel(i, j) = Prec(i, j)
               end do
            end do
         endif
      endif
      return
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c rewind tape, printout results
      entry EFG2
c
c for angle output
      if(Jct(40).eq.0) then
c
c for e,f,g output
         if(jefg.eq.0) go to 800
 
c end file jefg
         rewind jefg
         write(Iout, 250) jefg
  250    format('0E,F,G,EDOT,FDOT,GDOT OUTPUT ON BINARY DATA SET', i3)
         Line = Line + 2
         if(Jct(32).lt.0) return
         go to 800
      else
 
c end file iefg
         rewind iefg
         write(Iout, 300) iefg
  300    format(/
     .      ' ANGLE AND RANGE INFORMATION OUTPUT ON CHARACTER DATA SET'
     .      , i3)
         if(Jct(40).lt.0) return
         sublt8 = -1E10_10
 
c end file jang
         rewind jang
         read(jang, 350) head,
     .                   ((jjj1(i,k),i=1,10), minus(k), (jjj1(i,k),
     .                   i=11,23), k = 2, 5)
  350    format(51A2, 4(6I2,i2,i2,i5,i6,a1,i2,i2,i4,i6,i2,i2,i5,i6,2I8,
     .          i10,i10,i9))
         Line = 66
         jjj  = 2
         i89  = 1
         if(Jct(47).eq.9) i89 = 2
      endif
      do while(.true.)
         do j = jjj, 5
            ifrst8 = 0
            ras    = jjj1(9, j)*1E-3_10
            dra    = jjj1(10, j)*1E-4_10
            decls  = jjj1(13, j)*1E-2_10
            ddecl  = jjj1(14, j)*1E-3_10
            has    = jjj1(17, j)*1E-3_10
            dha    = jjj1(18, j)*1E-4_10
            az     = jjj1(19, j)*1E-5_10
            el     = jjj1(20, j)*1E-5_10
            jmon   = jjj1(1, j)
            jdyr   = jjj1(2, j) + imon(jmon)
            if((jjj1(3,j)/4)*4 .eq. (jjj1(3,j)) ) then
               if(jmon.gt.2) jdyr = jdyr + 1
            endif
            do while(.true.)
c
c page heading
               if(Line.ge.59) then
                  if(Jct(47).le.0) then
                     write(Iout, 360) head, Npage
  360                format('1', 51A2, 20x, 'PAGE', i5)
                     if(Jct(46).gt.0) then
                        write(Iout, 370)
  370                   format(/3x, 'DATE', 5x, 'UTC', 4x,
     .                         'RT.ASCENSION', 2x, 'DRA/DT', 1x,
     .                         'DECLINATION', 2x, 'DDEC/DT', 1x,
     .                         'LOC.HR. ANGLE', 2x, 'DHA/DT', 2x,
     .                         'AZIMUTH', 2x, 'ELEVATION', 1x,
     .                         '1-WAY DELAY', 2x, 'DOPPLER', 3x,
     .                         'DOPDOT'/1x, 'MO DY YR', 1x, 'HR MN SC',
     .                         1x, 'HR MN  SEC', 3x, 'SEC/SEC', 1x,
     .                         'DEG', 2x, '''', 2x, '''''', 3x,
     .                         '''''/SEC', 3x, 'HR MN  SEC', 3x,
     .                         'SEC/SEC', 3x, 'DEG', 7x, 'DEG', 6x,
     .                         'SECONDS', 5x, 'HERTZ ', 4x, 'HZ/MIN')
                     else if(Jct(48).le.0) then
                        write(Iout, 380)
  380                   format(/3x, 'DATE', 5x, 'UTC', 4x,
     .                         'RT.ASCENSION', 2x, 'DRA/DT', 1x,
     .                         'DECLINATION', 2x, 'DDEC/DT', 1x,
     .                         'LOC.HR. ANGLE', 2x, 'DHA/DT', 2x,
     .                         'AZIMUTH', 2x, 'ELEVATION', 3x, 'RANGE',
     .                         7x, 'RDOT', 6x, 'RDOT2'/1x, 'MO DY YR',
     .                         1x, 'HR MN SC', 1x, 'HR MN  SEC', 3x,
     .                         'SEC/SEC', 1x, 'DEG', 2x, '''', 2x,
     .                         '''''', 3x, '''''/SEC', 3x,
     .                         'HR MN  SEC', 3x, 'SEC/SEC', 3x, 'DEG',
     .                         7x, 'DEG', 7x, 'METERS', 5x, 'CM/MIN',
     .                         2x, 'CM/MIN**2')
                     else
                        write(Iout, 390)
  390                   format(/3x, 'DATE', 5x, 'UTC', 4x,
     .                         'RT.ASCENSION', 2x, 'DRA/DT', 1x,
     .                         'DECLINATION', 2x, 'DDEC/DT', 1x,
     .                         'LOC.HR. ANGLE', 2x, 'DHA/DT', 2x,
     .                         'AZIMUTH', 2x, 'ELEVATION', 3x, 'RANGE',
     .                         3x, 'DDEC/DT', 1x, 'DHA/DT', 2x,
     .                         'UTC'/1x, 'MO DY YR', 1x, 'HR MN SC',
     .                         1x, 'HR MN  SEC', 3x, ' ''''/SEC', 1x,
     .                         'DEG', 2x, '''', 2x, '''''', 3x,
     .                         '''''/SEC', 3x, 'HR MN  SEC', 3x,
     .                         ' ''''/SEC', 3x, 'DEG', 7x, 'DEG', 7x,
     .                         'METERS', 2x, ' ''''/SEC', ' ''''/SEC',
     .                         1x, 'HR MN')
                     endif
                  else
                     write(Iout, 400) Heding(2), (head(i), i = 4, 12),
     .                                (freq89(i,i89), i = 1, 6), Date,
     .                                Npage
  400                format('1', a4, 1x, 9A2, ' FREQ(GHZ) K-HORN=',
     .                      -9pf5.2, ' K-DISH=', -9pf5.2, ' K-UP=',
     .                      -9pf7.4, ' S-BAND=', -9pf4.2, ' UHF-DN=',
     .                      -9pf6.5, ' UHF-UP=', -9pf6.5, 1x, 2A4,
     .                      ' PAGE', i5)
                     write(Iout, 410)
  410                format(/3x, 'DATE', 2x, 'DY OF', 2x, 'UTC', 4x,
     .                      'AZIMUTH', 2x, 'ELEV. 1-WAY DLY', 15x,
     .                      '1-WAY DOPPLER (HERTZ)', 18x, '2-WAY DLY',
     .                      1x, 'RADIUS  W.LONG  N.LAT')
                     write(Iout, 420)
  420                format(1x, 'MO DY YR', 2x, 'YR HR MN SC', 3x,
     .                      'DEG', 5x, 'DEG', 2x, 'MICROSEC', 3x,
     .                      'K-HORN', 4x, 'K-DISH', 5x, 'K-UP', 5x,
     .                      'S-BAND', 2x, 'UHF-DN', 2x,
     .                      'UHF-UP MICROSEC', 3x, 'KM', 5x, 'DEG', 5x,
     .                      'DEG')
                  endif
                  Line  = 4
                  Npage = Npage + 1
               endif
c
c angle output
               ihsgn = iblank
               hax   = jjj1(15, j)*3600._10 + jjj1(16, j)*60._10 + has
               if(Jct(47).gt.0) then
c
c calculate subsatellite longitude and latitude
                  if(ifrst8.le.0) then
                     wlong = hax/240._10+Coords(2,1)
                     wlong = MOD(wlong,360._10)
                     if(wlong.gt.180._10) wlong = wlong-360._10
                     if(wlong.lt.-180._10) wlong = wlong+360._10
                     nlat = jjj1(11,j)+jjj1(12,j)/60._10+decls/3600._10
                     if(minus(j).eq.minus1) nlat = -nlat
                  endif
c
c nodal crossing printout
                  if(sublt8.ge.-1E9_10) then
                     if(sublt8*nlat.le.0._10) then
                        dd8  = ABS(sublt8)
                        dd8  = dd8/(dd8 + ABS(nlat))
                        sec8 = sec8 +
     .                         dd8*((jjj1(6,j)-sec8) + 60._10*(jjj1(5,j)
     .                         -min8) + 3600._10*(jjj1(4,j)-ihr8)
     .                         + 8.64E4_10*(jdyr-jdyr8)) + min8*60._10 +
     .                         ihr8*3600._10
                        i    = sec8/8.64E4_10
                        if(sec8.lt.0._10) i = i - 1
                        jdyr8  = jdyr8 + i
                        jdy8   = jdy8 + i
                        sec8   = sec8 - i*8.64E4_10
                        ihr8   = sec8/3600._10
                        sec8   = sec8 - ihr8*3600._10
                        min8   = sec8/60._10
                        sec8   = sec8 - min8*60._10
                        isec8  = sec8 + 0.500000000000001_10
                        icross = 1
                        if(sublt8.gt.0._10) icross = 2
                        subln8 = (wlong - subln8)*dd8 + subln8
                        write(Iout, 430) jmn8, jdy8, jyr8, jdyr8, ihr8,
     .                        min8, isec8, (ascnod(i,icross), i = 1, 2)
     .                        , subln8
  430                   format(i3, '/', i2, '/', i2, i4, i3, ':', i2,
     .                         ':', i2, 1x, 2A8, 79x, f8.3, '  0.00')
                        sublt8 = -1E10_10
                        ifrst8 = 1
                        Line   = Line + 1
                        go to 500
                     endif
                  endif
c
c les-8/9 output
                  onedel = jjj1(21,j)*1E3_10/Ltvel
                  iumbra = 1
                  if(onedel.lt.0._10) iumbra = 2
                  onedel = ABS(onedel)
                  twodel = 2._10*onedel
                  onedop = -jjj1(22,j)/60E5_10/Ltvel
                  do i = 1, 6
                     dopler(i) = freq89(i, i89)*onedop
                  end do
                  radius = jjj1(23,j)*1E-3_10
                  if(radius.lt.0) iumbra = 3
                  radius = ABS(radius)
                  write(Iout, 440) (jjj1(i,j), i = 1, 3), jdyr,
     .                             (jjj1(i,j), i = 4, 6), az, el,
     .                             onedel, (dopler(i), i = 1, 6),
     .                             twodel, radius, wlong, nlat,
     .                             umbra(iumbra)
                  Line = Line + 1
  440             format(i3, '/', i2, '/', i2, i4, i3, ':', i2, ':',
     .                   i2, f8.2, f7.2, f9.1, f10.1, f10.1, f10.1,
     .                   f9.1, f8.1, f8.1, f9.1, f8.1, f7.2, f7.2, a1)
                  sublt8 = nlat
                  subln8 = wlong
                  jmn8   = jjj1(1, j)
                  jdy8   = jjj1(2, j)
                  jyr8   = jjj1(3, j)
                  jdyr8  = jdyr
                  ihr8   = jjj1(4, j)
                  min8   = jjj1(5, j)
                  sec8   = jjj1(6, j)
               else
                  if(hax.gt.4.32E4_10) then
                     ihsgn = iminus
                     hax   = ABS(hax-8.64E4_10)
                     jjj1(15, j) = hax/3600._10
                     hax = hax - jjj1(15, j)*3600._10
                     jjj1(16, j) = hax/60._10
                     has = hax - jjj1(16, j)*60._10
                  endif
                  if(Jct(46).gt.0) then
                     onedel = jjj1(21,j)*1E-3_10/Ltvel
                     onedop = -jjj1(22,j)/60E5_10/Ltvel*Freq
                     oneddp = -jjj1(23,j)/60E5_10/Ltvel*Freq
                     write(Iout, 450) (jjj1(i,j), i = 1, 8), ras, dra,
     .                                minus(j), jjj1(11, j),
     .                                jjj1(12, j), decls, ddecl, ihsgn,
     .                                jjj1(15, j), jjj1(16, j), has,
     .                                dha, az, el, onedel, onedop,
     .                                oneddp
  450                format(i3, '/', i2, '/', i2, i3, ':', i2, ':', i2,
     .                      i3, i3, f7.3, f8.4, 1x, a1, i2, i3, f6.2,
     .                      f8.3, 1x, a1, i2, i3, f7.3, f8.4, f10.5,
     .                      f9.5, f11.8, f11.2, f10.2)
                  else if(Jct(48).le.0) then
                     write(Iout, 460) (jjj1(i,j), i = 1, 8), ras, dra,
     .                                minus(j), jjj1(11, j),
     .                                jjj1(12, j), decls, ddecl, ihsgn,
     .                                jjj1(15, j), jjj1(16, j), has,
     .                                dha, az, el,
     .                                (jjj1(i,j), i = 21, 23)
  460                format(i3, '/', i2, '/', i2, i3, ':', i2, ':', i2,
     .                      i3, i3, f7.3, f8.4, 1x, a1, i2, i3, f6.2,
     .                      f8.3, 1x, a1, i2, i3, f7.3, f8.4, f10.5,
     .                      f9.5, i11, i11, i10)
                  else
                     dra = dra*15._10
                     dha = dha*15._10
                     write(Iout, 470) (jjj1(i,j), i = 1, 8), ras, dra,
     .                                minus(j), jjj1(11, j),
     .                                jjj1(12, j), decls, ddecl, ihsgn,
     .                                jjj1(15, j), jjj1(16, j), has,
     .                                dha, az, el, jjj1(21, j), ddecl,
     .                                dha, jjj1(4, j), jjj1(5, j)
  470                format(i3, '/', i2, '/', i2, i3, ':', i2, ':', i2,
     .                      i3, i3, f7.3, f8.3, 1x, a1, i2, i3, f6.2,
     .                      f8.3, 1x, a1, i2, i3, f7.3, f8.3, f10.5,
     .                      f9.5, i11, f8.3, f7.3, i3, ':', i2)
                  endif
                  Line = Line + 1
                  if(Jct(44).gt.0) then
                     write(jct44, 480) (jjj1(i,j), i = 1, 8), ras, dra,
     .                                 minus(j), jjj1(11, j),
     .                                 jjj1(12, j), decls, ddecl,
     .                                 jjj1(21, j), jjj1(22, j)
  480                format(i2, '/', i2, '/', i2, i3, ':', i2, ':', i2,
     .                      i3, i3, f7.3, f8.4, 1x, a1, i2, i3, f6.2,
     .                      f8.3, i11, i10)
                  endif
               endif
               go to 550
  500       end do
c
c end of do loop
  550    end do
         read(jang, 600, end=700) ((jjj1(i,k),i=1,10), minus(k), (
     .                              jjj1(i,k),i=11,23), k = 1, 5)
  600    format(5(6I2,i2,i2,i5,i6,a1,i2,i2,i4,i6,i2,i2,i5,i6,2I8,i10,
     .          i10,i9))
         jjj = 1
      end do
  700 rewind jang
      if((Jct(44).gt.0) .and. (jct44.ne.7) ) rewind jct44
      return
 
c end file iefg
  800 rewind iefg
      write(Iout, 900) iefg
  900 format(/' E,F,G,EDOT,FDOT,GDOT OUTPUT ON EBCDIC DATA SET', i3)
      Line = Line + 2
      if(Jct(30).lt.0) return
      Line = 66
      nerr = 0
      nrec = 0
      read(iefg, 1300, err=1000) bb
      go to 1200
 1000 read(iefg, 1300) bb
      write(Iout, 1100) iefg
 1100 format(' **** ERROR ON FIRST HEADER RECORD OF DATA SET', i3)
      nerr = 1
 1200 read(iefg, 1300, err=1400, end=1900) aa
 1300 format(33A4)
      go to 1600
 1400 read(iefg, 1300) aa
      write(Iout, 1500) iefg
 1500 format(' **** BELOW IS AN ERROR RECORD ON DATA SET', i3)
      Line = Line + 1
      nerr = nerr + 1
 1600 nrec = nrec + 1
      if(Line.ge.59) then
         write(Iout, 1650) bb, Npage
 1650    format('1', 32A4, i4)
         write(Iout, 1700)
 1700    format(/2x, 'SDC', 1x, 'SATELLITE', 1x, 'YR', 1x, 'MN', 1x,
     .          'DY', 1x, 'HR', 1x, 'MIN', 1x, 'SEC', 7x, 'E (KM)',
     .          10x, 'F (KM)', 10x, 'G (KM)', 7x, 'EDOT (KM/SEC)', 3x,
     .          'FDOT (KM/SEC)', 3x, 'GDOT (KM/SEC)')
         Npage = Npage + 1
         Line  = 3
      endif
      write(Iout, 1800) aa
 1800 format(1x, 33A4)
      Line = Line + 1
      go to 1200
 1900 rewind iefg
      write(Iout, 2000) nrec, iefg, nerr
 2000 format(/' HEADER +', i5,
     . ' INFORMATION RECORDS (132 CHARACTERS EACH) WRITTEN ON DATA SET'
     . , i3, ' (', i3, ' ERROR RECORDS)')
 
      return
      end
