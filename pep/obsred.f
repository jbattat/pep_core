      subroutine OBSRED(mocpar)
 
      implicit none
c
c     t.forni  december 1967  subroutine obsred
c           radar, optic, trnsit are the calling routines
c           read one record from iiobcn and / or iiabs1 if applicable
c           set up obscrd labeled common
c           set up switches for rest of run
c
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
      include 'comdateq.inc'
      include 'empcnd.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
c ict(79)=-2 dummy observations not calculated each iteration nor
c            after convergence
c ict(79)=-1 dummy observations not calculated each iteration, but are
c            after convergence
c ict(79)= 0 dummy observations calculated each iteration, but not
c            after convergence
c ict(79)= 1 dummy observations calculated each iteration and after
c            convergence
c
c ict(80)=-1 input observation library tapes not used
c ict(80)= 0 input observation library tapes are used, program eats its
c             tale from iteration to iteration
c ict(80)= 1 input observation library tapes only are used to form the
c            normal equations for first iteration with no computation
c            of observed minus theory.
c            subsequent iterations same as ict(80)=0
c
      include 'funcon.inc'
      include 'inodta.inc'
      include 'ltrapobs.inc'
      include 'mtrapx.inc'
      integer*2 imdt(6)
      equivalence (imdt(1),Imdt1)
      include 'namtim.inc'
      include 'number.inc'
      include 'obscrd.inc'
      include 'kobequiv.inc'
      integer*2 iobsrv
      equivalence (iobsrv,Observ)
      character*8 sitt
      equivalence (Save(22),sitt)
      include 'obsdta.inc'
      include 'obstap.inc'
      include 'param.inc'
      include 'prpgat.inc'
      include 'redobs.inc'
      integer*2 nsavb(2)
      equivalence (nsavb,Nsavb1)
      integer*4 jdsb1,jdsb2
      equivalence (Jdsb(1,1),jdsb1),(Jdsb(1,2),jdsb2)
      real*4    erobs1,erobs2
      equivalence (erobs1,Erobsc),(erobs2,Erobsc(2))
      real*10 fdsb1,fdsb2
      equivalence (Fdsb(1,1),fdsb1),(Fdsb(1,2),fdsb2)
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'statsrad.inc'
      include 'trpcom.inc'
c
c savb(i,1), i=1,numsav,  saved precession-nutation,etc. from input
c                         observation library tape for radar, optical
c                         and first part of transit observation
c savb(i,2), i=1,numsav,  the same quantities read from iiobs
c     savb(1-8) are set from dstf(2-9) in comrit
c     savb(1) =dstf(2)= longitude of subradar point
c     savb(2) =dstf(3)= latitude of subradar point
c     savb(3) =dstf(4)= time delay from receive to reflection (sec)
c     savb(4) =dstf(5)= longitude of observed spot away from
c                       subradar spot
c                       -or- time delay to 2nd object
c                       -or- 2nd frequency (hz)
c     savb(5) =dstf(6)= latitude ofobserved spot away from
c                       subradar point
c                       -or- cumulative obsrvd psr phase (cycles)
c     savb(6) =dstf(7)= xmitter frequency for vlbi observations (hz)
c                       -or- interstellar freq for psr observations
c     savb(7) =dstf(8)= topography (2-way elevation in sec) of
c                       observed point
c                       -or- 1-way delay attributable to pulsar orbit
c                       (in sec) to be subtracted from the calculated
c                       delay from the pulsar system center of mass
c                       -or- 2nd xmitter frequency for vlbi (hz)
c     savb(8) =dstf(9)= instantaneous psr pulse period at reception
c     savb(9-26)are set from nutprc(1-18)
c     savb(27)= pc(1) = dpsi*cos(obliq)
c     for meridian-circle observations (ncodf.ge.4) savb(18-21) have a
c     different meaning:
c     savb(18)= salph
c     savb(19)= calph
c     savb(20)= tdelt
c     savb(21)= estf(10)
c     for transit/occultation data, savb(18-34) are treated differently
c     savb(18)= plate motion reference epoch
c     savb(19-21)= coords(4-6,1)
c     savb(22)= sitf(1)
c     savb(23-25)= coords(1-3,1)
c     savb(28-34)= occultation star data
c     savb(28-39) are used for radio observations
c     savb(28)= doppler counting interval
c     savb(28)= biased doppler or ddr counting interval
c     savb(29)= reference frequency for site1
c     savb(30)= reference frequency for site2
c     savb(31)= received frequency at site1
c     savb(32)= received frequency at site2
c     savb(33)= counted-cycle vlbi multiplication factor
c     savb(34)= ddr output bias freq (=xfactr*syn.freq./16)
c     savb(35)= synthesizer offset  between site1 and site2
c     savb(36)= half the modulo number for ranging
c     savb(37)= raw modded range datum
c     savb(38)= received freq. for second object at site1
c     savb(39)= received freq. for second object at site2
c     savb(40)= double precision 'seconds' part of observation epoch
c     savb(41-46)= meteorological data - should be set to -1 if unknown
c                  defaults are 273.15, 1000, and 0.0
c     savb(41)= surface temperature in deg k. at site1
c     savb(42)= surface pressure in millibars at site1
c     savb(43)= surface relative humidity in percent/100 at site1
c     savb(44)= surface temperature in deg k. at site2
c     savb(45)= surface pressure in millibars at site2
c     savb(46)= surface relative humidity in percent/100 at site2
c     savb(47)= predicted range residual on prdict tape
c     savb(48)= receive time offset between sites
c     savb(48)= height of spacecraft from planet center in km
c     savb(49-51)= receiving station to spacecraft line of sight
c                  vector in planet rotating frame
c     savb(52-55)= zenith angles at 1st & 2nd sites of 1st & 2nd
c                  objects
c     savb(56)= reference epoch for s/c state vector (sec past J2000)
c     savb(57-62)= s/c state vector in km, km/s relative to system c-o-m
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
c resltb(3,j) = time that first limb of planet cuts second limb of sun
c                 for transit observation (second part) (seconds)
c resltb(4,j) = time that second limb of planet cuts second limb of sun
c                 for transit observation (second part) (seconds)
c
c errorb(i,j) = error in resltb(i,j) in same units (i=1,4)
c        (i=3,4  not used)
c
c
c           niobc,niobs,niabs1  indicate status of data sets  iiobc,
c           iiobs, iiabs1
c              -1  not to be read
c              0  ready to be read
c              1  already read
c           if  iiobc, iiobs, or iiabs1  are initially zero,  data set
c           is not read in any case.
c
c
c external functions
      integer*4 ITYPOB
c
c           internal to this subroutine only
      real*10 ape,app,dth,dthd,fdedit,fds,fdysc1,
     .          fdysc2,fdysct,ome,omp,sdyscp,tdlt,tobs,tpe,tpp
      integer   i,ib,ibt,ict41,iict41,isw700,j,jbkfor,
     .          jdc1,jdc2,jdcp,jdct,jdedit,klang1,
     .          mocpar,nbkfor,nfirst
      integer   njdc,npe,npp
      real*10 dvx(2,2)/4*0._10/
c dvx is assumed not needed here, so  not fully dimensioned
c if correlated partials should ever be needed, then
c dvx would be expanded and must be correctly read
c at label=1030
      integer*2 nmp2/1/
      real*4    erwgtb(2)
      integer*2 itim
      real*10 reslt0(2),savb28
      character*2 blnk/'  '/
      character*1 warn1/'0'/, warn2/'N'/
c WARN1 tracks state of each series:
c  '0'=>initial or in midst of real data,
c  '1'=>first observation,
c  '+'=>subsequent observation, but previous data were all suppressed
c
c*  start=100
c
      if(Nk1.le.-1) then
 
c setup for first observation of series (for counting interval diff)
         reslt0(1) = 0._10
         reslt0(2) = 0._10
 
c set up time direction logic:  first, assume backwards in time
         jbkfor = 1
c forward in time if planet number greater than 30,
c for modern observation types, for spot on moon, or
c if observing site on non-earth body
         if(Jct(68).gt.0 .or. (Jct(68).eq.0 .and.
     .      (Nplnt0.gt.30 .or. (Ncodf.ge.1.and.Ncodf.lt.4) .or.
     .       Ncodf.gt.9 .or. Nspot.gt.0 .or. Nspot2.gt.0 .or.
     .       (Ksite(1).gt.0 .and. Ksite(1).ne.3))))
     .     jbkfor=2
 
         if(Idumob.ne.1) nbkfor = jbkfor
         if(warn1.eq.'0') then
            warn1='1'
            warn2='N'
         else
            warn1='+'
         endif
      else
         warn1='0'
      endif
c
c read iiobcn for one record of obs series (from disk)
c containing changes in observation if iiobcn.gt.0 and niobc=0
      isw700 = 0
 
  100 if(Iiobcn.le.0) goto 900
      if(Niobc.ne.0) goto 900
      goto 400
 
c err on iiobcn  skip record  no message
  200 if(Idumob.ne.1) call SUICID(
     .           'ERROR ON IOBCON IN NON-DUMMY MODE, STOP IN OBSRED   '
     .           , 13)
c err on iiobcn  dummy read to accomodate fortran release 11
c what if zero record --- program will hang up later
      write(Iout,300) Iiobcn
  300 format(' **** ERROR ON DATA SET IIOBCN=', i2,
     .' IN DUMMY MODE, DUMMY OBSERVATION INTERVAL RECORD SKIPPED IN OBSR
     .ED LABEL=130 ******')
      Line = Line + 1
      read(Iiobcn)
      Neribc = Neribc + 1
  400 read(Iiobcn,err=200) (Jdc(i),Ihrc(i),Iminc(i),Secc(i),
     .                        i = 1, 2), (Erobsc(i),i = 1,2),Intdyc,
     .                        Intscc
      call HRZSWT
 
c test for zero record
      if(Jdc(1).le.0) then
         Niobc  = -1
         nfirst = 3
         goto 900
      else
         Niobc = 1
c big question   to check or not to check input on iiobcn
c less than 60 secs, 60 mins, 24 hrs
c maybe check intday greater than o  intsec between o and 8.64E4_10
         Fdyscp = Ihrc(1)*3600._10 + Iminc(1)*60._10 + Secc(1)
         fdysct = Ihrc(2)*3600._10 + Iminc(2)*60._10 + Secc(2)
         jdcp   = Jdc(1)
         jdct   = Jdc(2)
 
c if dummy mode, direction given by input times
         if(Idumob.ne.1) then
            if(jbkfor.eq.1) goto 500
            if(jbkfor.eq.2) goto 600
         endif
         if(jdct.gt.jdcp .or.
     .      (jdct.eq.jdcp .and. fdysct.ge.Fdyscp)) goto 600
      endif
c*  start=300
c backwards in time: see if consistent
  500 nbkfor = 1
      if(jdct.gt.jdcp .or.
     .   (jdct.eq.jdcp .and. fdysct.gt.Fdyscp)) goto 700
 
c dates on iiobcn ok order
      njdc = 1
      goto 800
 
c forwards in time: see if consistent
  600 nbkfor = 2
      if(jdct.lt.jdcp .or.
     .   (jdct.eq.jdcp .and. fdysct.lt.Fdyscp)) goto 700
 
c dates on iiobcn ok order
      njdc = 1
      goto 800
 
c to fix not backwards or not forwards
  700 njdc   = 2
      sdyscp = Fdyscp
      Fdyscp = fdysct
      fdysct = sdyscp
      jdcp   = jdct
      jdct   = Jdc(1)
 
  800 nfirst = 1
      jdc1   = jdcp
      jdc2   = jdct
      fdysc1 = Fdyscp
      fdysc2 = fdysct
c
c*  start=1000
c read iiabs1  for one record of obs series from observation
c library tape
  900 if(isw700.gt.0) goto 2200
      if(Iiabs1.le.0) goto 1700
      if(Niabs1.ne.0) goto 1700
      savb28 = 0._10
 
      goto 1200
c err on iiabs1  skip record with message and statistics
c err on iiabs1  dummy read to accomodate fortran release 11
 1000 write(Iout,1100) Iiabs1
 1100 format(' **** ERROR ON DATA SET IIABS1=', i2,
     .', INPUT OBSERVATION LIBRARY RECORD SKIPPED IN OBSRED *******')
      Line = Line + 1
      read(Iiabs1)
      Nerrob = Nerrob + 1
 1200 do while( .true. )
         read(Iiabs1,err=1000,end=1300) Ncodeb(1,2),Ihrb(1,2),
     .        Iminb(1,2),Secb(1,2),
     .        (Resltb(i,2),Errorb(i,2),i = 1,2),Atutsb(1,2),
     .        Ututsb(1,2),Clampb(1,2),(Limbb(i,2),i = 1,2),
     .        Obsrvb(1,2),Imnthb(1,2),Idayb(1,2),Iyearb(1,2),
     .        Jdsb(1,2),Jdb(1,2),Ctatb(1,2),Ctrecb(1,2),Nsavb1,
     .        (Savb(i,1),i = 1,Nsavb1),Mun1,Munpar,
     .        ((Vired(i,j),i=1,Munpar),j = 1,Mun1),imdt,nmp2,
     .        ((dvx(1,j),i=1,nmp2),j = 1,Mun1),Ncal,
     .        (Cal(i),Scal(i),Ical(i),i = 1,Ncal),
     .        (Sumcor(i),i = 1,2),Mnshob,Mixshp,
     .        (Mshobs(i),i = 1,Mnshob)
         Mun2 = 2
         if(Mun1.eq.1) Mun2 = 1
 
c test for zero record
         if(Ncodeb(1,2).le.0) goto 1500
         call DTCKI(imdt)
         if(Ict(41).le.0) goto 1600
         if(Nk1.eq.-1) ict41 = -1
         ict41  = ict41 + 1
         iict41 = Ict(41)
         if(mod(ict41,iict41).eq.0) then
            Savb(28,1) = Savb(28,1) + savb28
            goto 1600
         else if(ITYPOB(Ncodf).eq.4 .and. nintrf.ge.0) then
 
c for accumulated n-count need count time from last observation
            savb28 = savb28 + Savb(28,1)
         endif
      end do
 1300 if(Line.ge.57) call OBSPAG
      write(Iout,1400) Iiabs1
 1400 format('0* * * END OF FILE ACCEPTED ON IABS1=', i2, ' * * *')
      Line   = Line + 2
      Iiabs1 = -1
 1500 Niabs1 = -1
      goto 1700
 1600 Niabs1 = 1
      fdsb2  = Ihrb(1,2)*3600._10 + Iminb(1,2)*60._10 + Secb(1,2)
c
c*  start=2000
c test for dummy mode
 1700 if(Idumob.eq.1) then
 
c dummy mode
         Numsav = 0
 
         if(nfirst.eq.1) then
 
c place initial values into /obscrd/
            nfirst = 2
            Jds    = Jdc(njdc)
            Jd     = Jds
            Sec    = Secc(njdc)
            Ihr    = Ihrc(njdc)
            Imin   = Iminc(njdc)
c erobs1= erobsc(1)
c erobs2= erobsc(2)
            Error(1) = erobs1
            Error(2) = erobs2
            do i = 1, 2
               Deriv(1,i) = Error(i)*Erwgta(i,Ntape)
            end do
c set ncode   for nonzero values of:
c       1         erobs1
c       2         erobs1 erobs2
c       3                erobs2
            if(Ncodf.le.0 .or. (Ncodf.ge.7 .and. Ncodf.le.9))
     .           call SUICID(' NCODF=0,7,8 OR 9,STOP IN OBSRED', 8)
            if(erobs1.eq.0.0 .and. erobs2.eq.0.0)
     .           call SUICID('BOTH EROBS1&2=0, STOP OBSRED',7)
            Ncode = 2
            if(erobs2.eq.0.) Ncode = 1
            if(erobs1.eq.0.) Ncode = 3
 
c set iabs1=0  in order that partials are calculated
            Iabs1 = 0
 
c initialize quantities not read from iiobcn
            Ctat    = 32.15_10
            Atuts   = 0.
            Ututs   = 0.
            Clamp   = blnk
            Limb(1) = blnk
            Limb(2) = blnk
            Observ  = blnk
         else if(nfirst.eq.2) then
c*  start=2300
c dummy mode with nfirst=2 (within given read of iiobcn)
c reset first error because of division by cos(decl) in optic
            Error(1)    = erobs1
            Deriv(1,1) = Error(1)*Erwgta(1,Ntape)
 
c place incremented values into /obscrd/
            if(nbkfor.eq.2) then
               Fdyscp = Fdyscp + Intscc
               jdcp   = jdcp + Intdyc
            else
               Fdyscp = Fdyscp - Intscc
               jdcp   = jdcp - Intdyc
            endif
            do while( Fdyscp.lt.0._10 )
               Fdyscp = Fdyscp + Secday
               jdcp   = jdcp - 1
            end do
            do while( Fdyscp.ge.Secday )
               Fdyscp = Fdyscp - Secday
               jdcp   = jdcp + 1
            end do
c finish test for dummy mode of interval being processed
c from iiobcn
            if((nbkfor.eq.2 .and.
     .          (jdct.lt.jdcp .or. (jdct.eq.jdcp.and.fdysct.lt.Fdyscp)))
     .         .or.
     .         (nbkfor.eq.1 .and.
     .          (jdct.gt.jdcp .or. (jdct.eq.jdcp.and.fdysct.gt.Fdyscp)))
     .     ) then
 
c read iiobcn for next time interval
               Niobc = 0
               goto 100
            endif
            Jds  = jdcp
            Jd   = Jds
            Imin = (Fdyscp + 1E-3_10)/60._10
            Sec  = Fdyscp - Imin*60
            if(Sec.lt.0.) Sec = 0.
            Ihr  = Imin/60
            Imin = Imin - Ihr*60
 
c nfirst=3  dummy mode completed  set switches and return
         else
            Ncode = 0
            goto 2400
         endif
c
c*  start=2500
c note if jds altered in chkhms,  jd(1) has not been altered
         call CHKHMS(Jds,Ihr,Imin,Sec)
         call MDYJUL(Imonth,Iday,Iyear,itim,Jds)
         fds = Fdyscp
         if(itim.ne.ctime .and. warn2.eq.'N') then
            warn2='Y'
            call SUICID('DUMMY DATA CENTURY DOES NOT MATCH ITIME ',-10)
         endif
         goto 2400
c
c
c           not dummy mode  iiobs and/or iiabs1 being processed
c           mtape=0, niabs1=1, niobs=1  compare iiobs,iiabs1
c                                                niobs=-1 store iiabs1
c           mtape=0, niabs1=-1, niobs=1  store iiobs
c                                                niobs=-1  return
c           mtape=1, niobs=1  store iiobs        niobs=-1  return
c           mtape=2, niabs1=1 store iiabs1      niabs1=-1  return
c           maybe check mtape = 0,1,2
c*  start=3000
      else if(Mtape.lt.1) then
         if(Niabs1.lt.0) then
            if(Niobs.lt.0) then
c*  start=8000
c completion of modes being processed by obsred
c set ncode=0  so that radar or optic or transit will return
c to compar
               Ncode = 0
               goto 2400
            else if(Niobs.eq.0) then
               goto 2600
            else
c
c*  start=5000
c partials to be calculated for process of iiobs alone
               Iabs1 = 0
               goto 1900
            endif
         else if(Niabs1.eq.0) then
            goto 2600
         else
            if(Niobs.lt.0) goto 2000
            if(Niobs.eq.0) goto 2600
         endif
      else if(Mtape.eq.1) then
         if(Niobs.lt.0) then
            Ncode = 0
            goto 2400
         else if(Niobs.eq.0) then
            goto 2600
         else
            Iabs1 = 0
            goto 1900
         endif
      else if(Niabs1.lt.0) then
         Ncode = 0
         goto 2400
      else if(Niabs1.eq.0) then
         goto 2600
      else
         goto 2000
      endif
 
c compare observation time from iiobs and iiabs1
      if(nbkfor.eq.2) then
         if(jdsb1.lt.jdsb2) then
            Iabs1 = 0
            goto 1900
         else if(jdsb1.eq.jdsb2) then
            if(ABS(fdsb2-fdsb1).ge.2._10) then
               if(fdsb1.gt.fdsb2) goto 2000
               Iabs1 = 0
               goto 1900
            endif
         else
            goto 2000
         endif
      else if(jdsb2.lt.jdsb1) then
         Iabs1 = 0
         goto 1900
      else if(jdsb2.eq.jdsb1) then
         if(ABS(fdsb2-fdsb1).ge.2._10) then
            if(fdsb2.gt.fdsb1) goto 2000
            Iabs1 = 0
            goto 1900
         endif
      else
         goto 2000
      endif
c set constants for process of iiabs1 as well as iiobs
c partials from iiabs1
      Iabs1 = Iiabs1
      if(Ict(4).gt.0) Iabs1 = 0
 
c set niabs1 to read iiabs1 again
      Niabs1 = 0
 
c process iiobs alone or with iiabs1
 1900 ib  = 1
      fds = fdsb1
 
c set niobs to read iiobs again
      Niobs = 0
      goto 2100
 
c process iiabs1
 2000 ib  = 2
      fds = fdsb2
 
c partials from iiabs1
      Iabs1 = Iiabs1
      if(Ict(4).gt.0) Iabs1 = 0
 
c set niabs1 to read iiabs1 again
      Niabs1 = 0
c ib=1 transfer information from iiobs  to /obscrd/
c ib=2 transfer information from iiabs1 to /obscrd/
 2100 do i = 1, 2
         Result(i)   = Resltb(i,ib)
         Error(i)    = Errorb(i,ib)
         Limb(i)     = Limbb(i,ib)
         Deriv(1,i) = Error(i)*Erwgta(i,Ntape)
      end do
      Ctat   = Ctatb(1,ib)
      Atuts  = Atutsb(1,ib)
      Ututs  = Ututsb(1,ib)
      Jds    = Jdsb(1,ib)
      Jd     = Jdb(1,ib)
      Sec    = Secb(1,ib)
      Ihr    = Ihrb(1,ib)
      Imin   = Iminb(1,ib)
      Imonth = Imnthb(1,ib)
      Iday   = Idayb(1,ib)
      Iyear  = Iyearb(1,ib)
      if(Iyear.gt.99) Iyear=Iyear-100
      Clamp  = Clampb(1,ib)
      Observ = Obsrvb(1,ib)
      Ncode  = Ncodeb(1,ib)
      Idumob = -ib
      ibt    = 3 - ib
      Numsav = nsavb(ibt)
      if(Numsav.gt.0) then
 
c copy savb into save vector
         do i = 1, Numsav
            Save(i) = Savb(i,ibt)
         end do
      endif
c
c*  start=5500
c copy particular items from observation record
      if(ib.eq.2) then
         Deriv(2,1) = Vired(2,1)
         Deriv(2,2) = Vired(2,Mun2)
         Dstf(1)     = Ctrecb(1,2)
         do i = 1, 8
            Dstf(i + 1) = Save(i)
         end do
c to be inserted
c contents of savb vector to be transfered to other cells
         if(Ncodf.ne.20) then
            if(Ncodf.lt.4 .or. Ncodf.gt.9) goto 2200
         endif
         Salph    = Save(18)
         Calph    = Save(19)
         Tdelt    = Save(20)
         Estf(10) = Save(21)
      endif
      if(Ncodf.ge.7 .and. Ncodf.le.9) then
 
c pull transit/occultation stuff out
         if(iobsrv.gt.-9 .and. iobsrv.lt.9) then
            Ksite(1)=iobsrv
         else
            j=1
            call DECODI(observ,j,2,i)
            Ksite(1)=i
         endif
c Save(22) has site name
         Sitf(1)     = sitt
         Coords(1,1) = Save(23)
         Coords(2,1) = Save(24)
         Coords(3,1) = Save(25)
         if(Save(18).lt.2.E7_10 .or. Save(18).gt.3.E7_10) then
            do i=18,21
               Save(i)=0._10
            end do
         endif
         T0sit(1)   = Save(18)
         Coords(4,1) = Save(19)
         Coords(5,1) = Save(20)
         Coords(6,1) = Save(21)
      endif
c
c*  start=7000
c           ntape not 1    return
c           ntape=1  niobc=1  iiobcn is read     niobc=-1  return
c           compare iiobcn with date in /obscrd/
 2200 isw700 = 0
      if(Ntape.eq.1) then
         if(Niobc.ge.0) then
            if(Niobc.eq.0) then
               goto 2600
            else
c is date in /obscrd/ within interval of iiobcn
c force 0 <= seconds < 86400
               jdedit = Jds
               fdedit = fds
               if(fdedit.lt.0._10) then
                  fdedit = fdedit + Secday
                  jdedit = jdedit - 1
               else if(fdedit.ge.Secday) then
                  fdedit = fdedit - Secday
                  jdedit = jdedit + 1
               endif
               if(nbkfor.eq.2) then
                  if(jdc1.lt.jdedit) then
                  else if(jdc1.eq.jdedit) then
                     if((fdedit-1E-3_10).lt.fdysc1) goto 2220
                  else
                     goto 2220
                  endif
                  if(jdc2.lt.jdedit) then
                  else if(jdc2.eq.jdedit) then
                     if((fdedit+1E-3_10).le.fdysc2) goto 2300
                  else
                     goto 2300
                  endif
               else
                  if(jdedit.lt.jdc1) then
                  else if(jdedit.eq.jdc1) then
                     if((fdedit-1E-3_10).gt.fdysc1) goto 2220
                  else
                     goto 2220
                  endif
                  if(jdedit.lt.jdc2) then
                  else if(jdedit.eq.jdc2) then
                     if((fdedit+1E-3_10).ge.fdysc2) goto 2300
                  else
                     goto 2300
                  endif
               endif
c*  start=7200
c outside interval read another record from iiobcn
c return to label=7000
               Niobc = 0
 
c set another switch so that return to label=7000 will occur
               isw700 = 1
               goto 100
            endif
 
c outside interval  don't read another record of iiobcn
 2220       erwgtb(1) = 1.0
            erwgtb(2) = 1.0
         endif
      endif
      goto 2400
 
c between time interval of iiobcn
 2300 erwgtb(1) = erobs1
      erwgtb(2) = erobs2
c reading from iiabs1 only within iiobcn interval.
c shall this record of iiabs1 be skipped?
      if(Ntape.ne.1 .or. Mtape.ne.2 .or.
     . erwgtb(1).gt.-1.E3 .or. erwgtb(2).gt.-1.E3) goto 2350
c
c skip record within series on iiabs1
c
c if observable is accumulated n-count must restart count with
c first included observation (convenient also for phase obs.)
 2320 if(ITYPOB(Ncodf).eq.4 .and. nintrf.ge.0 .and.
     .          Ncodeb(1,2).le.1) then
         reslt0(1) = Result(1)
         reslt0(2) = Result(2)
      endif
      isw700 = 0
      goto 100
c
c*  start=7500
c calculate derivatives
 2350 do i = 1, 2
         Deriv(1,i) = Error(i)*Erwgta(i,Ntape)*erwgtb(i)
 
c set switches
      end do
 2400 if(Ncode.gt.0) then
 
c check for global deletions
         if(Tdlt0.gt.0._10 .and. Tdltp.gt.0._10) then
            tdlt = MOD(Jds + fds/86400._10 - 0.5_10 - Tdlt0,Tdltp)
            if(tdlt.LT.0._10) tdlt = tdlt+Tdltp
            if(tdlt.gt.Tdlton) then
 
c skip if not in proper range (unweight if IOBS)
               if(Idumob.ne.-1) goto 2320
               Deriv(1,1)=0._10
               Deriv(1,2)=0._10
            endif
         endif
 
         if(Idumob.ne.1) then
            Result(1) = Result(1) - reslt0(1)
            Result(2) = Result(2) - reslt0(2)
         endif
 
c reslt0 non-zero only for n-count when changing initial point
c
c
c set tguess for 1st point of series (if possible)
         npe    = 3
         klang1 = Klan
         if(Klan.le.0) klang1 = Klanb
         npp = Nplnt(klang1)
         if(Nk1.lt.0 .and. (Ksite(1).le.0 .or. Ksite(1).eq.3) .and.
     .    (klang1.gt.0 .and. klang1.le.u_mxpl) .and. Icnd(klang1).eq.0
     .    .and. Pcond(2,klang1).lt.0.9_10) then
            app  = Pcond(1,klang1)
            tobs = Jds + fds/Secday
            tpp  = T0svpl
            if(klang1.eq.Klanb) tpp = T0svsb
            tpp = tobs - tpp
            ape = Econd(1)
            tpe = tobs - T0svem
            omp = 1._10
            if(npp.le.30) omp = 1._10 + Mass(npp)
            omp  = SQRT(omp/app**3)
            ome  = SQRT((1._10+Mass(npe))/ape**3)
            dthd = 0._10
            do i = 4, 6
               dthd = dthd + Pcond(i,klang1) - Econd(i)
            end do
            dth    = Gauss*(omp*tpp - ome*tpe) + Convd*dthd
            Tguess = SQRT(ape*(ape-2._10*app*COS(dth))+app**2)*Aultsc
            if(mod(Jct(6)/64,2).ne.0) then
               write(Iout,2410) Tguess,dth
 2410          format(' OBSRED: INITIAL DISTANCE APPROX.',1pd20.10,
     .          ' SEC  AT LONGITUDE DIFFERENCE',0pf10.4, ' RAD')
               Line = Line + 1
            endif
         endif
c
c*  start=9500
c flush iiobcn to end of series if ncode=0
      else if(Iiobcn.gt.0) then
         if(Niobc.ge.0) then
            do while( .true. )
               read(Iiobcn,err=200) Jdc(1)
               if(Jdc(1).le.0) goto 2500
            end do
         endif
      endif
c
c*  start=9900
c set 'nice' flag
 2500 Nice = Ncode - 2
      if(Nice.gt.1) Nice = Nice - 3
      if(Ncode.le.0 .or. Nice.gt.1) then
 
c end of observation series, finish writing obslib series
         if(Ict(1).ge.0) call COMRIT(1)
 
c test for bad ncode
         if(Ncode.lt.0 .or. Nice.gt.1) then
            write(Iout,2520) Ncode
 2520       format('0*** BAD OBSERVATION CODE:', i7)
            call SUICID('BAD OBSERVATION CODE, STOP IN OBSRED', 9)
         endif
      endif
      return
 
2600  call SUICID(' NIOBS, NIABS1, NIOBC ARE ZERO INCORRECTLY IN '//
     .            'OBSRED', 13)
 
      end
