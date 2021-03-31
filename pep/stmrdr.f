      subroutine STMRDR(mocpar)
 
      implicit none

c     m.e.ash   aug 1970    subroutine stmrdr
c     main program for processing delay and doppler for signal
c     sent from a ground station to a satellite over to another
c     satellite and then down to another ground station

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
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'ltrapx.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,tmdly0,dop0
      equivalence (ctrecf,Dstf),(fdev,Dstf(10))
      equivalence (Result,tmdly0),(Result(2),dop0)
      include 'kobequiv.inc'
      include 'obstap.inc'
      include 'param.inc'
      include 'redobs.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'statsrad.inc'
c
c        nprec  = 0 nutation-precession determined midway between
c                   receiving and sending times
c        nprec  = 1 nutation-precession determined at both receiving and
c                   sending times
c
c        ndprec = 0 derivative of nutation-precession matrix not used in
c                   determining velocities
c        ndprec = 1 derivative of nutation-precession matrix used in
c                   determining velocities
c
c        neatm  =-1 effect of earth atmosphere on radar signal is
c                   ignored
c        neatm  = 0 effect of earth atmosphere on time delay only is
c                   included
c        neatm  = 1 effect of earth atmosphere on time delay and doppler
c                   shift is included
c
c        neion  =-1 effect of earth ionosphere on probe signal ignored
c        neion  = 0 effect of earth ionosphere on probe delay is calc.
c        neion  = 1 effect of earth ionosphere on probe dopler is calc.
c
c        iwob   = 0 wobble effect not included
c        iwob   = 1 wobble effect is included
c
c        kob(j) = possible additional control integers, j=11,16)
c
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
      real*4 STORNE
      integer*4 JULDAY

c local variables
      real*10 secbx
      character*2 astrik(2)/'  ', '* '/
      logical*4 cislun
      integer*4 i,ictst
      integer*2 iyr19
c
c initialization
      cislun = Ncp0.eq.3 .or. Ncp0.eq.10
      Tguess = 100._10
      if(Ncp0.eq.3) Tguess  = 0.2_10
      if(Ncp0.eq.10) Tguess = 1.5_10
c
c determine data about observation series
      nprec  = 0
      ndprec = 0
      ictst  = 1
      if(Ncodf.eq.2) ictst = 0
      if(Ict(29).ge.ictst) ndprec = 1
      if(Ict(28).ge.ictst .and.  .not. cislun) nprec = 1
      Iwob = 1
c
c determine receiving and sending site quantities
      call ESHAPE(2)
c
c write first page heading
      call OBSBEG(' 2 SPCRFT DELAY ', Sitf(1),'  TIME DELAY', 3,
     .            'DOPPLER SHIFT   ', 4)
      if(Ncodf.ne.14) then
         write(Intern,740) Sitf(2),Aplnt(Klanb),Aplnt(Klans1),Sitf(1)
  740    format('0RADIO TRANSPONDER OBSERVATIONS FROM SEND SITE ',a8,
     .    ' TO SPCRFT ',a8,' TO SPCRFT ',a8,' TO RCV SITE ',a8)
      else
         write(Intern,730) Sitf(2),Aplnt(Klans1),Aplnt(Klanb),Sitf(1)
  730    format('0RADIO TRANSPONDER OBSERVATIONS FROM SEND SITE ',a8,
     .    ' TO SPCRFT ',a8,' TO & BACK FROM SPCRFT ',a8,
     .    ' TO RCV SITE ',a8)
      endif
      write(Intern,100)
  100 format('0 GRNWCH  JULIAN  REC UTC (WWV)  UT1-UTC AT-UTC  CT-AT',
     . 11x,'TIME DELAY (SECONDS)',13x,'DOPPLER SHIFT (CYCLES/SECOND)'/
     . '   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC',
     . 7x, 'OBSERVED     ERROR     OBS-TH', 9x,
     .     'OBSERVED     ERROR     OBS-TH')
      call OBSBGT
      goto 400
c
c read observation data card
c
c if err on iiobs, dummy read necessary
  200 write(Iout,300) Iiobs
  300 format(' **** ERROR ON DATA SET IIOBS=', i2,
     .', OBSERVATION CARD RECORD SKIPPED IN STMRDR **********')
      Line = Line + 1
      read(Iiobs,600)
      Nerrra = Nerrra + 1
      goto 500
  400 if(mocpar.gt.0) goto 700
      if(Niobs.ne.0 .or. Iiobs.le.0) goto 700
  500 read(Iiobs,600,err=200) Ncodeb(1,1),Ihrb(1,1),Iminb(1,1),secbx,
     . Resltb(1,1),Errorb(1,1),Resltb(2,1),Errorb(2,1),
     . Atutsb(1,1),Ututsb(1,1),Imnthb(1,1),Idayb(1,1),Iyearb(1,1)
  600 format(i1,2I3,f8.4,f13.7,e7.0,d15.8,e7.0,f8.4,f7.4,3I2)
      Niobs = -1
      Secb(1,1) = STORNE(secbx)
      if(Ncodeb(1,1).gt.0) then
         Fdsb(1,1) = Ihrb(1,1)*3600._10 + Iminb(1,1)*60._10 + Secb(1,1)
         iyr19=Iyearb(1,1)+ctime*100
         Jdsb(1,1) = JULDAY(Imnthb(1,1),Idayb(1,1),iyr19)
         Jdb(1,1)  = Jdsb(1,1)
         Niobs = 1
      endif
  700 call OBSRED(mocpar)
      if(Nice.lt.-1) then
c
c printout at end of observation series
         call OBSEND
 
         return
      endif
c
c calculation of theoretical value of time delay and/or
c doppler shift
c set tmdly0  used by deldop
      if(Idumob.eq.1 .or. Nice.gt.0) tmdly0 = Tguess*2.01_10
      call STMDLD
      if(Jd.le.-10) goto 400
 
c idumob=1  dummy mode            idumob=-1  not dummy mode
      if(Idumob.eq.1) then
 
c for dummy mode only
         Result(1) = Deriv(2,1)
         Result(2) = Deriv(2,2)
      endif
      call OBSCNT
      if(Jd.le.0) goto 400
      if(Ict(2).le.0 .or. mocpar.gt.0) then
c
c printout observed minus theory
         if(Nice.lt.0) then
            write(Iout,710) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,
     .       Ututs,Atuts,Ctat,tmdly0,(Deriv(i,1),i=1,2),astrik(Nast11)
         else if(Nice.eq.0) then
            write(Iout,710) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,
     .       Ututs,Atuts,Ctat,tmdly0,(Deriv(i,1),i=1,2),astrik(Nast11),
     .       dop0,(Deriv(i,2),i=1,2),astrik(Nast12)
  710       format(i3,'/',i2,'/',i2,i8,2I3,4F8.4,
     .       f16.10,1pe8.1,1pe13.5,1A1,1pd17.9,1pe8.1,1pe13.5,1A1)
         else
            write(Iout,720) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,
     .       Ututs,Atuts,Ctat,dop0,(Deriv(i,2),i=1,2),astrik(Nast12)
  720       format(i3,'/',i2,'/',i2,i8,2I3,4F8.4,38x,
     .       1pd17.9,1pe8.1,1pe13.5,1A1)
         endif
         Line = Line + 1
      endif
c
c calculate partial derivatives w.r.t. quantities to be adj.
      if(Ict(1).lt.0) goto 400
      if(Ict(1).eq.0) then
c
c write(obs-th),partials buffer
         call COMRIT(0)
      else if(Idumob.ne.1 .or. (Ict(3).ge.0)) then
         call PARTL(1)
         call COMRIT(0)
      endif
      goto 400
      end
