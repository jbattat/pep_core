      subroutine TRNSIT(mocpar)
 
      implicit none
 
c           ash/forni/slade  1971 subroutine trnsit
c           revised 1977 to allow mutual occultations - j.f.chandler
c           main program for processing transit and occultation
c           observations

c arguments
      integer*4 mocpar
c     mocpar = 0 compar called in midst of least squares iteration
c     read all three data sets iobcon, iobs, iabs1 which become in
c     subroutines of compar iiobcn, iiobs, iiabs1
c     mocpar = 1 compar called at end of least squares iteration to
c                calculate dummy observations
c     read only iobcon which become in subroutines of compar iiobcn
c        procedural note - all coordinates are to be expressed in au,
c        au/day.  this means that values read from tape need not be
c        scaled (except moon), but that site coordinates must be.
c
c array dimensions
      include 'globdefs.inc'
 
c commons
      include 'comdateq.inc'
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'ltrapobs.inc'
      include 'namtim.inc'
      character*8 sname
      equivalence (sname,Aplnt(18))
      include 'number.inc'
      include 'obscrd.inc'
      integer*2 ngo
      equivalence (ngo,Clamp)
c usage in 'transit' -
c observ has ksite(1), clamp has sun or moon limb no. (1 or 2)
      real*4    acctim, accdst, accprc
      equivalence (acctim,Estf),(accdst,Estf(2)),
     .            (accprc,Estf(3))
c     ctime = centuries from 1900 (needed because only last two
c               digits of year are read in)
      include 'obstap.inc'
      include 'param.inc'
      include 'redobs.inc'
      integer*2 ksitt
      equivalence (ksitt,Obsrvb(1,1))
      character*8 sitt,numst
      real*10 scordt(3),starp(7),svelt(3),t0st,ra,xmotra,dec,xmotdc,
     . prlx,rvel
      equivalence (starp,ra),(starp(2),xmotra),(starp(3),dec),
     .            (starp(4),xmotdc),(starp(5),prlx),
     .            (starp(6),rvel),(starp(7),numst)
      equivalence (sitt,Savb(22,2)),(scordt,Savb(23,2)),
     . (t0st,Savb(18,2)),(svelt,Savb(19,2)),(starp,Savb(28,2))
      integer*2 iclmpb(2,2)
      equivalence (iclmpb(1,1),Clampb(1,1))
      include 'sitcrd.inc'
      include 'statsrad.inc'
      include 'stcord.inc'
      include 'trnocc.inc'
      real*10 mrad0, mrad
      equivalence (Rhos,mrad0),(Rhom,mrad)
      include 'watstf.inc'
c
c     Note: the arrays in REDOBS common for quantities to be read from
c     card or tape were all doubled to accommodate the scheme of up to 4
c     observables per card in TRNSIT.  The complete list is:
c      Resltb,Errorb,Ncodeb,Ihrb,Iminb,Secb,Atutsb,Ututsb,Clampb,Limbb,
c      Obsrvb,Imnthb,Idayb,Iyearb,Jdsb,Jdb,Ctatb,Ctrecb,Fdsb
c     However, the second halves of only two of these are ever used:
c      Errorb and Ncodeb.
c savb(i,1), i=1,numsav,  saved precession-nutation,etc. from input
c                         observation library tape
c savb(i,2) i=1,numsav   these saved quantities from input observation
c                         cards  (i.e., to be saved on output tape)
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
c clampb(1,j) = 1 or 2 depending on whether 1st or 2nd half of
c               observation card.
c
c           niobc,niobs,niabs1  indicate status of data sets  iiobc,
c           iiobs, iiabs1
c              -1  not to be read
c              0  ready to be read
c              1  already read
c           if  iiobc, iiobs, or iiabs1  are initially zero,  data set
c           is not read in any case.
c
c        type of observation specified by ibtrn
c           ibtrn=1 moon occults star
c           ibtrn=2 planet occults star
c           ibtrn=3 moon occults planet
c           ibtrn=4 planet transits sun
c           ibtrn=5 (not used)
c           ibtrn=6 planet occults planet (mutual event)
c           ibtrn=7 (reserved for black dot photographic obs)
c           ibtrn=8 planet occults planet (mid-time)
c           ibtrn=9 planet eclipses planet (mid-time)
c
c
c local
      real*10 ascdc, ccse, ctut,rhom0, secra
      integer*4 i, i4, ij, imndc, imra, iobn2, itst1, itst2,
     . jdso, k, klap2, nkodf, npair, ntype
      character*1 astrik(2)/' ','*'/
      character*8 typobs(2,5)/'  OCCULT','ATION   ', '    TRAN',
     1 'SIT     ', '  SPOT T','RANSIT  ', 'TRANSIT ','MID-TIME',
     2 'ECLIPSE ','MID-TIME'/
      character*4
     .    obsnam(3,3)/'1ST ','CROS','SING', '2ND ','CROS','SING',
     .                            'SEPA', 'RATI','ON  '/
      character*8 eclocc(2)/' CROSSES','ECLPS BY'/,ecoc
      character*4 oclmd/'OMID'/,eclmd/'EMID'/
      character*8 sitfo,sitco8,astar
      character*4 sitco
      character*8 ecentr/'ECENTER '/,ord1(2)/'  FIRST ','  SECOND'/
      character*8 occmes(7)/'********',' IS NOT ','CONCENTR','IC WITH ',
     .          '********', ', STOP I', 'N TRNSIT'/
      character*8 zname,aname
      equivalence (occmes,zname),(occmes(5),aname)
      character*1 lbf(17)
      character*8 blnk8/'        '/,bf(3),dtblnk/'   /  / '/
      character*4 blank,amper/'&&&&'/,sitec
      equivalence (blank,blnk8),(lbf,bf)
      character*1 vflg
      integer*2 ih(4),im(4)
      real*4 sc(4)
      integer*2 iyr19
      real*10 jd2000/2451545._10/
 
c constant for pre-rounding hms (.5*least significant digit)
      real*10 round/5E-4_10/
 
c value of century in days
      real*10 centry/36524.21988_10/
 
c times for printout
      integer*2 ihot(2),imot(2),ialh,ialm,ided,idem
      real*4    scot(2),alse,dese
      equivalence (ihot(1),ialh),(ihot(2),ided),
     .            (imot(1),ialm),(imot(2),idem),
     .            (scot(1),alse),(scot(2),dese)

c used for reading astronomical notation ra, dec
      character*3 hra,ddc

c external functions
      real*10 DMS2D
      integer*4 JULDAY
c
c set up indices
      nkodf = Ncodf - 8
c
c setup data about observation series
      sitco = blank
      sitfo = blnk8
      jdso  = 0
 
c fdsbo=-1._10
      Dstf(5) = 0._10
c
c calculate rhos and rhom in au
      Rhos  = Sunrtr/Aukm
      zname = sname
      rhom0 = Pcond(29,Klap)/Aukm
      aname = Pname
 
c set rhos here for satellite mutual occultation
      if(Nps1.gt.0) then
         klap2 = Klan
         if(Klans1.gt.0) klap2 = Klans1
         Rhos   = Pcond(29,klap2)/Aukm
         zname  = Aplnt(klap2)
         Nplnz  = Nps1
         Npzcnt = Ncs1
 
c make sure that only one planet tape needed for mutual event
         if(Mnplnt.ne.0 .and. Ncs1.gt.0) then
c they don't match
            if(Ncs1.ne.Nplnt(Klan)) call SUICID(occmes,14)
         endif
      endif
      Ibtrn = 0
      iobn2 = 2
 
c distinguish mid-time events
      if(Sita2(Ntape).eq.oclmd) Ibtrn = 8
      if(Sita2(Ntape).eq.eclmd) Ibtrn = 9
      if(Ibtrn.eq.0) then
         Ibtrn = nkodf*3 + 4
         if(Nps1.gt.0) Ibtrn = Ibtrn + 2
         if(Mnplnt.ne.0 .and. Ibtrn.eq.1) Ibtrn = 2
      else
         iobn2 = 3
         ecoc  = eclocc(Ibtrn - 7)
      endif
      if(Ibtrn.le.3) then
 
c set occultation radius
         mrad0 = Mcond(29)/Aukm
         if(Mnplnt.ne.0) mrad0 = rhom0
      endif
c
c write first page heading
      ntype = nkodf + 2
      if(Ibtrn.gt.7) ntype = Ibtrn - 4
      call OBSBEG(typobs(1,ntype),blnk8,obsnam,3,obsnam(1,iobn2),3)
      if(Ibtrn.le.3) then
         write(Intern,50)
   50    format('0 GRNWCH   JULIAN CT-UT        IMMERSION  TIME',13x,
     .    'EMMERSION  TIME      OCCULTED    SITE        SITE   ',
     .    'COORDINATES'/
     .    '   DATE     DAY            OBSERVED   ERROR  OBS-TH    ',
     .    'OBSERVED   ERROR  OBS-TH    STAR'/
     .    11x,'NUMBER  SEC   HR MIN  SEC   SEC    SEC    HR MIN  SEC',
     .    2('   SEC '),5x,'Z.C.#   NAME      RADIUS  LONGITUDE  ',
     .    'LATITUDE')
      else if(Ibtrn.lt.7) then
         write(Intern,100) zname
  100    format('0 GRNWCH   JULIAN CT-UT  FIRST  LIMB  OF  PLANET', 5x,
     .    'SECOND  LIMB  OF  PLANET    CROSSES    SITE        SITE   ',
     .    'COORDINATES'/
     .    '   DATE     DAY            OBSERVED   ERROR  OBS-TH    ',
     .    'OBSERVED   ERROR  OBS-TH  ', a8/
     .    11x,'NUMBER  SEC   HR MIN  SEC   SEC    SEC    HR MIN  SEC',
     .    2('   SEC '),4x,'LIMB    NAME      RADIUS  LONGITUDE  ',
     .    'LATITUDE')
      else
         write(Intern,150) ecoc,zname
  150    format('0 GRNWCH   JULIAN CT-UT     CENTER  OF  PLANET', 34x,
     .    a8, '   SITE        SITE   COORDINATES'/
     .    '   DATE     DAY',12X,'OBSERVED   ERROR  OBS-TH  SEPARA',
     .    'TION  ERROR  OBS-TH   ', a8/
     .    11x,'NUMBER  SEC   HR MIN  SEC   SEC    SEC',8x,
     .    '/SUM OF RADII',21x,'NAME      RADIUS  LONGITUDE  LATITUDE')
      endif
      call OBSBGT
      goto 300
c*  start=300
c
c err on iiobs   dummy read to accomodate fortran release 11
  200 read(Iiobs,500) (Ncodeb(i,1),i = 1,2)
      Nerrra = Nerrra + 1
      Niobs  = -1
      goto 400
c
c get observation data into storage
  300 if(mocpar.gt.0) goto 1000
      if(Niobs.ne.0 .or. Iiobs.eq.0) goto 1000
 
c read observation card
      Nsavb2 = 0
  400 if(Ibtrn.le.3) then
c read card for stellar or planetary occultation
c ncodeb(1,1) can only be 0 or 1 (=>  1st limb of moon only)
c ncodeb(2,1) can only be 0 or 3 (=>  2nd limb of moon only)
         do i = 1, 4
            ih(i) = 0
            im(i) = 0
            sc(i) = 0.
            Errorb(i,1) = 0.
         end do
         Ncodeb(1,1) = 0
         Ncodeb(2,1) = 0
         read(Iiobs,450,err=200) itst1,itst2,sitec,ih(1),im(1),sc(1),
     .    Errorb(1,1),ih(4),im(4),sc(4),Errorb(4,1),Ang1,Ang2,
     .    Posang,Corect,Derpos,Imnthb(1,1),Idayb(1,1),Iyearb(1,1)
  450    format(2I1,1x,a4,2(1x,2I2,f6.3,f4.1),2F7.3,f8.3,f5.2,
     .          f5.2, 3x, 3I2)
         if(itst1.eq.1) Ncodeb(1,1) = 1
         if(itst1.eq.2) Ncodeb(2,1) = 3
         if(itst2.eq.1) Ncodeb(2,1) = 3
      else
c
c           ncodeb(i,1)=1  transit of first limb of planet only
c           ncodeb(i,1)=2  transit of first and second limb
c           ncodeb(i,1)=3  transit of second limb only
c           i= 1   transit is across first limb of sun
c           i= 2   transit is across second limb of sun
         read(Iiobs,500,err=200) (Ncodeb(i,1),i = 1,2),sitec,
     .    (ih(i),im(i),sc(i),Errorb(i,1),i=1,4),
     .    Imnthb(1,1),Idayb(1,1),Iyearb(1,1)
  500    format(2I1,1x,1A4,4(1x,2I2,f6.3,f4.1),5x,3I2)
      endif
 
      Niobs = -1
      if(Ncodeb(1,1).le.0 .and. Ncodeb(2,1).le.0) goto 1000
      Niobs = 1
c
c read site card if necessary
      if(sitec.ne.blank) then
         if(sitec.eq.amper) then
 
c site is at center of earth
            sitt = ecentr
            do i = 1, 3
               scordt(i) = 0._10
               svelt(i) = 0._10
            end do
            t0st=0._10
            ksitt=0
            goto 700
         else
c
c see if same site as last card
            if(sitec.eq.sitco) goto 700
c different site from last card
c search list of 'standard' (adjustable) sites
            do k = 1, Numsit
               if(sitec.eq.Site(1,k)) then
c found proper input site
c move name and coordinates to waiting area (must interleave
c data with input obslib tape in obsred)
                  ksitt = Kscrd(k)
                  t0st=T0Site(k)
                  sitt=sitd(k)
                  do i=1,3
                     scordt(i) = Scord(i,k)
                     svelt(i) = Scord(i+3,k)
                  end do
                  goto 700
               endif
            end do
 
c site not found
            write(Iout,520) sitec,Imnthb(1,1),Idayb(1,1),Iyearb(1,1)
  520       format('0  NO MATCH FOUND FOR SITEC=',a4,'  DATE=',i2,2I3)
            call SUICID('NO INPUT SITE FOUND, STOP IN TRNSIT ', 9)
         endif
      endif
c*  start=500
c read card with site coordinates
      read(Iiobs,600) sitt,scordt,vflg,ksitt
  600 format(a8,3F16.11,1x,a1,12x,i2)
      if(vflg.eq.'6') then
         read(Iiobs,610) svelt,t0st
  610    format(8x,3f16.9,f8.0)
      else
         t0st=jd2000
         do i=1,3
            svelt(i)=0._10
         end do
      endif
 
c remember this site name for shortcut next time
  700 sitco8 = sitt
      Nsavb2 = 25
      if(Ibtrn.le.2) then
c
c read star card for stellar occultation
         read(Iiobs,750) numst,hra,imra,secra,xmotra,ddc,imndc,
     .                    ascdc, xmotdc, prlx, rvel
  750    format(a8,2(a3,i3,f11.6,f10.5),2F9.4)
c
c           input units: ra,xmotra - hms, s/century
c                       dec,xmotdc - dms, s/century  (arc)
c                        prlx - arc-sec
c                        rvel - km/s
c           for storage, ra and dec are converted to radians,
c           xmotra and xmotdc to radians/day
c           prlx to radians (1/au)
c           rvel to au/day
c
         ra     = DMS2D(hra,imra,secra)*15._10*Convd
         xmotra = xmotra*Convhs/centry
         dec    = DMS2D(ddc,imndc,ascdc)*Convd
         xmotdc = xmotdc*Convds/centry
         prlx   = prlx*Convds
         rvel   = rvel*Secday/Aukm
         Nsavb2 = 34
      endif
c*  start=600
c
c conversions of data from observation card
      iyr19 = Iyearb(1,1) + ctime*100
      Jdsb(1,1)   = JULDAY(Imnthb(1,1),Idayb(1,1),iyr19)
      Jdb(1,1)    = Jdsb(1,1)
      iclmpb(1,1) = 1
      if(Ncodeb(1,1).gt.0) goto 900
c card has only 2nd half of observation left to do
c move 2nd half quantities to 1st half storage locations
  800 iclmpb(1,1) = 2
      Ncodeb(1,1) = Ncodeb(2,1)
      if(Ncodeb(1,1).le.0) goto 300
 
c signal that card contains more info
      Niobs = 1
 
      errorb(1,1)=errorb(3,1)
      errorb(2,1)=errorb(4,1)
      ih(1) = ih(3)
      ih(2) = ih(4)
      im(1) = im(3)
      im(2) = im(4)
      sc(1) = sc(3)
      sc(2) = sc(4)
c fill 'time of observation' with observed time of event
c note that the time is the observable
  900 ij = 1
      if(Ncodeb(1,1).ge.3) ij = 2
      Ihrb(1,1)  = ih(ij)
      Iminb(1,1) = im(ij)
      Secb(1,1)  = sc(ij)
      do i = 1,2
         Resltb(i,1) = (60*ih(i) + im(i))*60
         Resltb(i,1) = Resltb(i,1) + sc(i)
      end do
      Fdsb(1,1) = Resltb(ij,1)
c
c read observation library tape and select next data point
 1000 call OBSRED(mocpar)
      if(Nice.lt.-1) then
c*  start=9000
c
c printout at end of observation series
         call OBSEND
 
         return
      else
c
c setup observing site quantities if different site
         if(Sitf(1).ne.sitfo) then
            call ESHAPE(1)
            sitfo = Sitf(1)
c decide whether to abbreviate 2nd line
c require: 1) same page, 2) 2nd half, 3) same date as last one
         else if(Line.lt.58) then
            if(ngo.eq.2) then
               if(Jds.eq.jdso) then
 
c 2nd half of same observation card
                  npair = 1
                  bf(1) = blnk8
                  bf(2) = blnk8
                  bf(3) = blnk8
                  goto 1100
               endif
            endif
         endif
c not the same date, site, etc.
c encode date and jdate into array
         npair = 0
         bf(1) = dtblnk
         i4    = Imonth
         call EBCDI(i4,lbf(1),3)
         i4 = Iday
         call EBCDI(i4,lbf(5),2)
         i4 = Iyear
         call EBCDI(i4,lbf(8),2)
         call EBCDI(Jds,lbf(10),8)
         jdso = -10
      endif
c
c*  start=800
c start of first and second limb of sun logic
 1100 if(ngo.eq.2) then
         Rhom     = -rhom0
         Trnrd(1) = Aukm
         Trnrd(2) = -Aukm
      else
         Rhom     = rhom0
         Trnrd(1) = -Aukm
         Trnrd(2) = Aukm
      endif
c
c set up deriv instead of result as first guess
      if(Idumob.ne.-2) then
 
c data from card input
         Deriv(2,1) = Result(1)
         Deriv(2,2) = Result(2)
      else
 
c data from obslib tape
         if(Nice.le.0) Deriv(2,1) = Result(1) - Deriv(2,1)
         if(Nice.ge.0) Deriv(2,2) = Result(2) - Deriv(2,2)
      endif
 
      call OCCULT(nkodf)
      if(Jd.le.-10) goto 300
      if(Jd.gt.0) then
 
c idumob=1  dummy mode            idumob=-1  not dummy mode
         if(Idumob.eq.1) then
 
c for dummy mode only
            Result(1) = Deriv(2,1)
            Result(2) = Deriv(2,2)
c
c for processing real observations
         else if(Idumob.eq.-1) then
 
c observation from card, use time from card
            ihot(1) = ih(1)
            ihot(2) = ih(2)
            imot(1) = im(1)
            imot(2) = im(2)
            scot(1) = sc(1)
            scot(2) = sc(2)
            goto 1200
         endif
 
c convert results to h,m,s (dummy or data from obslib)
         if(Nice.le.0) then
            ccse = Result(1) + round
            ialm = ccse/60._10
            alse = ccse - ialm*60 - round
            if(alse.lt.0.) alse = 0.
            ialh = ialm/60
            ialm = ialm - ialh*60
            if(Nice.lt.0) goto 1200
         endif
         ccse = Result(2) + round
         idem = ccse/60._10
         dese = ccse - idem*60 - round
         if(dese.lt.0.) dese = 0.
         ided = idem/60
         idem = idem - ided*60
      endif
c*  start=1000
c count cards and test page logic
 1200 call OBSCNT
      if(Jd.le.0) goto 300
      if(Ict(2).eq.0 .or. mocpar.gt.0) then
 
c printout observed minus theory
         ctut = Ctat + Atuts
 
         if(Ibtrn.gt.3) astar = ord1(ngo)
         if(Ibtrn.gt.7) astar = blnk8
         if(Nice.lt.0) then
 
c ngo=1 write first line   ngo=2 write second line
            if(npair.eq.1) then
               write(Iout,1210) bf,ctut,ialh,ialm,alse,
     .          Deriv(1,1),Deriv(2,1),astrik(Nast11),astar
            else
               write(Iout,1210) bf,ctut,ialh,ialm,alse,
     .          Deriv(1,1),Deriv(2,1),astrik(Nast11),astar,Sitf(1),
     .          (Coords(i,1),i = 1,3)
 1210          format(2A8,a1,f7.3,2I3,f7.3,f5.1,f8.2,a1,29x,
     .                a8, 1x, a8, f11.4, 2F10.5)
            endif
         else if(Nice.eq.0) then
            if(Ibtrn.gt.7) then
               write(Iout,1220) bf,ctut,ialh,ialm,alse,
     .          Deriv(1,1),Deriv(2,1),astrik(Nast11),
     .          Result(2),Deriv(1,2),Deriv(2,2),astrik(Nast12),astar,
     .          Sitf(1),(Coords(i,1),i = 1,3)
 1220          format(2A8,a1,f7.3,2I3,f7.3,f5.1,f8.2,a1,f12.7,
     .          f6.3, f9.5, a1, 1x, a8, 1x, a8, f11.4, 2F10.5)
            else if(npair.eq.1) then
               write(Iout,1230) bf,ctut,ialh,ialm,alse,
     .          Deriv(1,1),Deriv(2,1),astrik(Nast11),ided,idem,dese,
     .          Deriv(1,2),Deriv(2,2),astrik(Nast12),astar
            else
               write(Iout,1230) bf,ctut,ialh,ialm,alse,
     .          Deriv(1,1),Deriv(2,1),astrik(Nast11),ided,idem,dese,
     .          Deriv(1,2),Deriv(2,2),astrik(Nast12),astar,Sitf(1),
     .          (Coords(i,1),i = 1,3)
 1230          format(2A8,a1,f7.3,2(2I3,f7.3,f5.1,f8.2,a1,1x),a8,
     .                1x, a8, f11.4, 2F10.5)
            endif
         else if(npair.eq.1) then
            write(Iout,1240) bf,ctut,ided,idem,dese,Deriv(1,2),
     .       Deriv(2,2),astrik(Nast12),astar
         else
 
c no mid-time observations without ncode.le.2
            write(Iout,1240) bf,ctut,ided,idem,dese,Deriv(1,2),
     .       Deriv(2,2),astrik(Nast12),astar,Sitf(1),(Coords(i,1),i=1,3)
 1240       format(2A8,a1,f7.3,28x,2I3,f7.3,f5.1,f8.2,a1,1x,
     .             a8, 1x, a8, f11.4, 2F10.5)
         endif
 
c*   start=1500
         Line = Line + 1
      endif
c
c*  start=2000
      jdso = Jds
 
c fdsbo=fdsb(1,1)
      if(Ict(1).ge.0) then
c
c calculate partial derivatives w.r.t. quantities to be adj.
         if(Ict(1).ne.0) call PARTL(3)
c
c write(obs-th),partials buffer
         call COMRIT(0)
      endif
c
c termination of first and second limb of sun logic
      if(ngo.ge.2) goto 300
 
c if data were from card, look at 2nd half
      if(Idumob.eq.-1) goto 800
      goto 300
      end
