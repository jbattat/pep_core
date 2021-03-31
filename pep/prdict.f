      subroutine PRDICT(jda,jdb,mode)
 
      implicit none

c m.ash/l.friedman  july 1969   subroutine prdict
c main program for prediction of new observed minus theory
c residuals with error statistics and/or for performing harmonic
c analysis of residuals

c arguments
      integer*4 mode
      real*10 jda,jdb

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'bernum.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'ktrap.inc'
      integer*2 numobs
      equivalence (Num2,numobs)
      include 'obsdta.inc'
      include 'rtsidesl.inc'
      include 'prdmov.inc'
      integer*4 ist1,ist2
      equivalence (Istrik,ist1),(Istrik(2),ist2)
      real*10 ersav(5,2)
      include 'wrkcompr.inc'
c
c local
      real*10 fnobs,jdobs,tfnmo1,tfnmo2
      integer   i,j,jdu,jdu1,jline,k,keasm1,keasm2,keasmt,
     .          measm9,neasm9,nfirst,ns1,ntyp9,nvmobt
      character*16 typobs(14)/'     RADAR      ', '      RADIO     ',
     .    '      RADAR     ', '*** optical ****', '  PHOTOGRAPHIC  ',
     .    '  PHOTOGRAPHIC  ', '  OCCULTATION   ', '    TRANSIT     ',
     .    '  SPOT TRANSIT  ', '  PULSAR        ', ' INTERFEROMETER ',
     .    '   ??           ', ' 2 SPCRFT DELAY ', ' DOPPLER COUNT  '/
      character*16 azel(3)/   'MERIDIAN CIRCLE ', 'AZIMUTH-ELEVATN ',
     .                        '  LOOK ANGLE    '/
      real*10 a(2, 3), enobs(2)
      character*2 astrik(2)/'  ', '* '/
      character*8 prdnam/' PRDICT '/
      character*20 colhed(2,6)/
     .         'RIGHT ASCENSION (S) ', '  DECLINATION ('''')  ',
     .         '   AZIMUTH (DEG)    ', '  ELEVATION (DEG)   ',
     .         ' PITCH ANGLE (DEG)  ', '  ROLL ANGLE (DEG)  ',
     .         '  TIME DELAY (SEC)  ', 'DOPPLER SHIFT (C/S) ',
     .         ' FIRST MEASUREMENT  ', ' SECOND MEASUREMENT ',
     .         'DIFFERENTL DELAY (S)', 'DFRNTL DLY RATE  S/S'/
      real*4    hr, xnout(3, 2)
c
c     iabs1=input obs-theory and partial derivitive tapes
c     iabs2=output if desired
c
c           write prdict link title page
      call PAGSET('PREDICTED OBSERVED-THEORY   ',7)
      call NEWPG
      write(Iout,100)
      write(Iout,200)
  100 format('-'/'-'/'-'/'-'//30x,
     .    '**********    **********    **********    **    ********** '
     .    ,'  ************'/30x,
     .    '***********   ***********   ***********   **   *********** '
     .    ,'  ************'/30x,
     .    '**       **   **       **   **       **   **   **          '
     .    ,'       **     '/30x,
     .    '**       **   **       **   **       **   **   **          '
     .    ,'       **     '/30x,
     .    '**       **   **       **   **       **   **   **          '
     .    ,'       **     '/30x,
     .    '**       **   **       **   **       **   **   **          '
     .    ,'       **     '/30x,
     .    '***********   ***********   **       **   **   **          '
     .    ,'       **     ')
  200 format(30x,
     .    '**********    **********    **       **   **   **          '
     .    ,'       **     '/30x,
     .    '**            **  **        **       **   **   **          '
     .    ,'       **     '/30x,
     .    '**            **   **       **       **   **   **          '
     .    ,'       **     '/30x,
     .    '**            **    **      **       **   **   **          '
     .    ,'       **     '/30x,
     .    '**            **     **     **       **   **   **          '
     .    ,'       **     '/30x,
     .    '**            **      **    ***********   **   *********** '
     .    ,'       **     '/30x,
     .    '**            **       **   **********    **    ********** '
     .    ,'       **     ')
      call PLINK
      call OPRMSG('PREDICTION OF OBSERVED MINUS THEORY OR DE-CORRELAT'//
     .  'ION OF PARTIAL DERIVATIVES',19)
      if(Ict(11).lt.-2) call DECORR
      if(Ict(11).le.-2) return
      if(Mout.gt.0 .and. Ict(11).eq.-1) write(Mout,300)
  300 format(' PRDICT WITH PRINT ONLY')
c
c initialization
      jline  = 60
      Neasmt = 0
      Measmt = 0
      do j = 1,2
         Ucert(j) = 0._10
         do i = 1,3
            Ermesp(i,j) = 0._10
         end do
      end do
      Negu  = 0
      Niobc = 0
 
c niobc =-1 do not read iobcon
c niobc = 0 read first record of series but not interior records
c niobc = 1 read interior records of series but not first record
      if(Iobcon.le.0) Niobc = -1
c
c increment tape counter
      nvmobt = Numobt
      if(Ict(16).gt.0) Numobt = Ict(16)
      Jtape = 0
  400 call PRDOBS(2,1)
      if(Iabs1.le.0) then
c
c error analysis for old & new observed minus theory residuals
         call PAGSET(-1,Jout)
         call NEWPG
         call PAGSET(-1,0)
         if(Ict(11).gt.-2) then
            fnobs = Measmt
            if(Measmt.eq.0) fnobs = 1._10
            do j = 1,2
               a(j,1)  = Ermesp(1,j)/fnobs
               a(j,2)  = Ermesp(2,j)/fnobs
               enobs(j)= Ermesp(3,j)/fnobs
               a(j,3)  = SQRT(enobs(j))
            end do
            write(Iout,420) Measmt,Nparam,a,enobs,
     .                       (Ermesp(3,j),j = 1,2)
  420       format('0ERROR ANALYSIS FOR THE',i8,
     .   ' MEASUREMENTS USED TO FORM THE NORMAL EQUATIONS OF ORDER',i4/
     .       35x,' FORMER       PREDICTED'/
     .       35x,'RESIDUALS     RESIDUALS'/
     .       '          AVERAGE (OBS-TH)/ERROR',1p,2e14.5/
     .       '       AVERAGE ABS(OBS-TH)/ERROR',2e14.5/
     .       ' ROOT MEAN SQUARE (OBS-TH)/ERROR',2e14.5/
     .       '     AVERAGE ((OBS-TH)/ERROR)**2',2e14.5/
     .       '         SUM ((OBS-TH)/ERROR)**2',2E14.5)
            keasmt = Neasmt - Measmt
            write(Iout,440) keasmt,Neasmt
  440       format(/i8,' OF THE TOTAL',i8,
     .             ' MEASUREMENTS WERE DELETED *'/)
            Line = Line + 12
            if(Jout.gt.0) then
 
c prdict statistics output on jout
               write(Jout,420) Measmt,Nparam,a,enobs,
     .                          (Ermesp(3,j),j = 1,2)
               write(Jout,440) keasmt,Neasmt
            endif
 
            if(Niobc.ge.0) then
               rewind Iobcon
               Itrwnd(Iobcon) = 0
            endif
         endif
 
         call TIMRIT('PREDICTION OF OBSERVED MINUS THEORY ',9)
 
         Numobt = nvmobt
c at this point file jout used to be rewound,
c just as at the end of each iteration in the analiz link.
c now jout is left open to the end of the job.
         return
      endif
 
 
      if(Mout.gt.0 .and. Ict(11).ge.0) write(Mout,450) Iabs2
  450 format(' PRDICT WITH IABS2=',i3)
 
c
c write first two records of output observation library tape
      Tapnam = prdnam
      if(Ict(11).ge.0) call PRDOBS(-2,1)
 
c
c read first record of observation series
  500 call PRDOBS(2,3)
c
c write first record of observation series
      if(Ict(11).ge.0) call PRDOBS(-2,3)
      if(Ncodf.le.0) goto 400
      nfirst = 0
      if(Jtypob.eq.2 .or. Jtypob.eq.3) then
         if(Ict(10).lt.0) goto 1000
      else if(Ict(10).gt.0) then
         goto 1000
      endif
c
c setup error statistics for this series
      call ZFILL(Erm22,16*28)
      do i=1,2
         Meas22(i)=0
         Neas22(i)=0
      end do

      Line = 60
c
c get correct observation type name
      if(Ncodg.ne.20) then
         if(Ncodg.eq.14) Ncodg = 13
         if(Ncodg.eq.19) Ncodg = 14
         if(Ncodg.eq.10 .or. Ncodg.eq.12) Ncodg = 11
         if(Ncodg.eq.18) Ncodg = 10
         if(Ncodg.gt.14) Ncodg = 12
         ntyp9 = 1
         if(Ncodg.eq.4) then
            if(Klans1.gt.0) then
               ntyp9 = 3
            else if(Klanb.gt.0) then
               if(Ncp0.eq.3 .or. Ncp0.eq.10) ntyp9 = 2
            endif
         endif
      else
         Ncodg = 4
         ntyp9 = 2
      endif
      typobs(4) = azel(ntyp9)
      if(Ncodg.lt.4) then
         ntyp9 = 4
      else if(Ncodg.gt.6) then
         if(Ncodg.eq.13) then
            ntyp9 = 4
         else
            ntyp9 = 5
            if(Ncodg.eq.11) ntyp9 = 6
         endif
      endif
c
c write header on auxillary data set
      if(Nout.gt.0) write(Nout,600)
  600 format(//)
      if(Mout.gt.0) write(Mout,700) Ntape,Nseq,Series,Ncodf,
     .                   Sitf(1,1),Sitf(1,2),Nplnt0,Spota
  700 format(' NTAPE=',i3,'  NSEQ=',i4,' SER=',a4,'  NCODF=',i2,
     .       '  SITES=',a4,1x,a4,'  NPLNT0=',i2,' SPOT =',a4)
      if(Mout.gt.0 .and. Nplnt2.ne.0) write(Mout,800) Nplnt2,Spota2
  800 format(59x,'NPLNT2=',i2,' SPOT2=',a4)
      if(Nout.gt.0) write(Nout,900)
  900 format(' JD HR MIN SEC   HR    OLD O-C    NEW O-C   ERROR   OLD ',
     . 'O-C  NEW O-C ERROR')
 1000 do while( .true. )
c
c read observation record
         call PRDOBS(2,4)
         if(Ncode.le.0) goto 1200
c
c        decide if a predict is desired for this observation and get
c        applicable solutions into core.
c             mode=1  implies solution is in core already
c             mode=2  implies that we must check
c
         jdobs = Jd + Fdysc0/8.64E4_10
         if(jdobs.ge.jda .and. jdobs.le.jdb) then
 
            if(mode.eq.2) call GETSOL(jdobs)
 
            if(Jtypob.eq.2 .or. Jtypob.eq.3) then
               if(Ict(10).lt.0) goto 1100
            else if(Ict(10).gt.0) then
               goto 1100
            endif
c
c calculate iptr vector
            call PRDEQS(nfirst)
c
c predict
            call PRDSID
c
c printout new observed minus theory
            if(Ict(11).le.0) then
               if(Line.gt.57) then
                  call NEWPG
                  write(Iout,1010) typobs(Ncodg),Plnnam,Nplnt0,Spota,
     .                              (Sitf(i,1),i = 1,2),
     .                              (Sitf(i,2),i = 1,2),Series,
     .                              Ntape,Nseq,
     .                              (colhed(i,ntyp9),i = 1,2)
 1010             format('0',a16,'OBS-TH RESIDUALS OF ',2A4,
     .                   ' (NPLNT0=',i2,',SPOT=',a4,') MADE AT ',
     .                   2(2A4,1x),' SERIES=',a4,' NTAPE=',i3,
     .                   ' NSEQ= ',i5/
     .                   '0 GRNWCH  JULIAN  REC UT(SIGNAL)',12x,a20,
     .                   23x,a20/ '   DATE   DAY NUM HR MIN  SEC  ',
     .                   2(3x,'OLD',7x,'NEW',5x,' NEW-OLD ',3x,'ERROR',
     .                   1x,'UNCERTAINTY'))
                  Line = Line + 4
               endif
               tfnmo1 = Obsth(1) - Deriv(2,1)
               if(Nice.lt.0) then
                  write(Iout,1020) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,
     .             Deriv(2,1),Obsth(1),tfnmo1,Deriv(1,1),astrik(ist1),
     .             Ucert(1)
 1020             format(i3,'/',i2,'/',i2,i8,2I3,f8.4,
     .                   2(1p,2D10.3,1pd11.3,1pd9.2,1A1,1pd9.2))
               else if(Nice.eq.0) then
                  tfnmo2 = Obsth(2) - Deriv(2,2)
                  write(Iout,1020) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,
     .             Deriv(2,1),Obsth(1),tfnmo1,Deriv(1,1),astrik(ist1),
     .             Ucert(1),Deriv(2,2),Obsth(2),tfnmo2,Deriv(1,2),
     .             astrik(ist2),Ucert(2)
               else
                  write(Iout,1030) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,
     .             Deriv(2,1),Obsth(1),tfnmo1,Deriv(1,1),astrik(ist1),
     .             Ucert(1)
 1030             format(i3,'/',i2,'/',i2,i8,2I3,f8.4,50x,1p,
     .                   2D10.3,1pd11.3,1pd9.2,1A1,1pd9.2)
               endif
               Line = Line + 1
            endif
            if(Nice.lt.0) then
 
c save predicted range residual
               ns1 = Numsav + 1
               if(ns1.le.46) then
                  do i=ns1,46
                     Save(i)=0._10
                     if(i.ge.41) Save(i)=-1._10
                  end do
               endif
               Save(47) = Obsth(1)
               if(Numsav.lt.47) Numsav = 47
            endif
            if(Ict(11).le.0) then
c
c write residuals on auxilliary data set
               if(Nout.gt.0) then
                  jdu = Jds - 2440000
                  if(Meas22(1).le.1 .and. Meas22(2).le.1) then
                     do j = 1,2
                        do i = 1,3
                           xnout(i,j) = 0.
                        end do
                     end do
                     jdu1 = jdu
                  endif
                  if((jdu-jdu1).gt.1 ) jdu1 = jdu
                  hr = (Ihr*3600+Imin*60+Sec)/36E2_10+(jdu-jdu1)*24
                  do j = 1,2
                     xnout(1,j) = Deriv(2,j)
                     xnout(2,j) = Obsth(j)
                     xnout(3,j) = Deriv(1,j)
                  end do
                  if(Nice.lt.0) then
                     write(Nout,1040) jdu,Ihr,Imin,Sec,hr,
     .                                 (xnout(i,1),i = 1,3)
 1040                format(i4,2I3,f5.1,f9.5,1p,2D10.2,d8.1)
                  else if(Nice.eq.0) then
                     write(Nout,1050) jdu,Ihr,Imin,Sec,hr,
     .                                 ((xnout(i,j),i=1,3),j = 1,2)
 1050                format(i4,2I3,f5.1,f9.5,1p,2D10.2,d8.1,
     .                      2D10.2,d8.1)
                  else
                     write(Nout,1060) jdu,Ihr,Imin,Sec,hr,
     .                                 (xnout(i,2),i = 1,3)
 1060                format(i4,2I3,f5.1,f9.5,20x,1p,2D10.2,d8.1)
                  endif
               endif
            endif
c
c replace old observed minus theory with new
            do j = 1,numobs
               Deriv(2,j) = Obsth(j)
            end do
            goto 1200
         endif
 1100 end do
c
c write new observed minus theory tape
 1200 if(Ict(11).ge.0) call PRDOBS(-2,4)
      if(Ncode.gt.0) goto 1000
c
c accumulate total error statistics
      measm9 = Meas22(1) + Meas22(2)
      Measmt = Measmt + measm9
      neasm9 = Neas22(1) + Neas22(2)
      Neasmt = Neasmt + neasm9
      do j = 1,2
         do i = 1,3
            ersav(i,j)  = Erm22(i,j,1) + Erm22(i,j,2)
            Ermesp(i,j) = Ermesp(i,j) + ersav(i,j)
         end do
      end do
c
c printout error statistics for this observation series
      if(Ict(11).gt.-2) then
         fnobs = measm9
         if(measm9.le.0) fnobs = 1._10
         do j = 1,2
            ersav(5,j) = ersav(3,j)
            do i = 1,3
               ersav(i,j) = ersav(i,j)/fnobs
            end do
            ersav(4,j) = ersav(3,j)
            ersav(3,j) = SQRT(ersav(3,j))
         end do
         do k = 1,2
            fnobs = Meas22(k)
            if(Meas22(k).le.0) fnobs = 1._10
            do j = 1,2
               if(Meas22(k).le.0) Erm22(7,j,k) = 1._10
               Erm22(5,j,k) = Erm22(3,j,k)
               Erm22(7,j,k) = SQRT(Erm22(3,j,k)/Erm22(7,j,k))
               Erm22(4,j,k) = Erm22(3,j,k)/fnobs
               Erm22(3,j,k) = SQRT(Erm22(4,j,k))
               Erm22(2,j,k) = Erm22(2,j,k)/fnobs
               Erm22(1,j,k) = Erm22(1,j,k)/fnobs
               Erm22(6,j,k) = SQRT(Erm22(6,j,k)/fnobs)
            end do
         end do
         if(Line.gt.45) then
            call NEWPG
            write(Iout,1010) typobs(Ncodg),Plnnam,Nplnt0,Spota,
     .                        (Sitf(i,1),i = 1,2),
     .                        (Sitf(i,2),i = 1,2),Series,Ntape,
     .                        Nseq,(colhed(i,ntyp9),i = 1,2)
            Line = Line + 4
         endif
         write(Iout,1250) measm9,Meas22,(colhed(i,ntyp9),i=1,2),
     .                    (ersav(i,1),ersav(i,2),
     .                     (Erm22(i,1,k),Erm22(i,2,k),k=1,2),i = 1,5),
     .                    ((Erm22(i,1,k),Erm22(i,2,k),k=1,2),i = 6,7)
 1250    format('0ERROR STATISTICS FOR THE',i8,' =',i7,' +',i7,
     .' RESIDUALS INCLUDED IN THE NORMAL EQUATIONS FROM THE ABOVE ',
     . 'SERIES'/41x,'TOTAL',11x,2(10x,a20)/
     . 35x,3(' FORMER       PREDICTED',7x)/
     . 35x,3('RESIDUALS     RESIDUALS',7x)/
     . '          AVERAGE (OBS-TH)/ERROR',3(1p,2E14.5,2x)/
     . '       AVERAGE ABS(OBS-TH)/ERROR',3(2E14.5,2x)/
     . ' ROOT MEAN SQUARE (OBS-TH)/ERROR',3(2E14.5,2x)/
     . '     AVERAGE ((OBS-TH)/ERROR)**2',3(2E14.5,2x)/
     . '         SUM ((OBS-TH)/ERROR)**2',3(2E14.5,2x)/
     . '       ROOT MEAN SQUARE (OBS-TH)',30x,2(2E14.5,2x)/
     . ' WEIGHTED ROOT MEAN SQUARE (OBS-TH)',27x,2(2E14.5,2x))
         keasmt = neasm9 - measm9
         keasm1 = Neas22(1) - Meas22(1)
         keasm2 = Neas22(2) - Meas22(2)
         write(Iout,1300) keasmt,keasm1,keasm2,neasm9,Neas22
 1300    format(/i8,' =',i7,' +',i7,' OF THE TOTAL',i8,' =',i7,
     .          ' +',i7,' MEASUREMENTS WERE DELETED *')
         Line = Line + 13
         if(Jout.gt.0) then
 
c prdict statistics output on jout
            if(jline.gt.43) then
               write(Jout,1310) Iterat,Heding,Date
 1310          format('1PREDICTED OBSERVED-THEORY  ITERAT=',i3,3x,
     .                18A4,1x,2A4)
               jline = 1
            endif
            write(Jout,1320) typobs(Ncodg),Plnnam,
     .                        Nplnt0,Nspot,(Sitf(i,1),i = 1,2),
     .                        (Sitf(i,2),i = 1,2),Series,Ntape,
     .                        Nseq,Npage
 1320       format(// 1x,a16,'OBS-TH RESID.OF ',2A4,' (NPLNT0=',
     .             i2,', NSPOT=',i2,') MADE AT ',2(2A4,1x),
     .             ' SERIES=',a4,' NTAPE=',i3,' NSEQ=',i5,
     .             ' PAGE',i4)
            jline = jline + 3
            write(Jout,1250) measm9,Meas22,(colhed(i,ntyp9),i = 1,2),
     .                    (ersav(i,1),ersav(i,2),
     .                     (Erm22(i,1,k),Erm22(i,2,k),k=1,2),i = 1,5),
     .                    ((Erm22(i,1,k),Erm22(i,2,k),k=1,2),i = 6,7)
            write(Jout,1300) keasmt,keasm1,keasm2,neasm9,Neas22
            jline = jline + 13
         endif
      endif
      goto 500
      end
