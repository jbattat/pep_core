      subroutine OBSBEG(typobs,site,lab1,l1,lab2,l2)
 
      implicit none

c arguments
      character*4 typobs(4),lab1(l1),lab2(l2)
      character*8 site
      integer*4 l1,l2
c
c        subroutine obsbeg - j.f.chandler, 1977 oct
c        compar link page format control routine
c        each observables main program should call obsbeg at the very
c        beginning to set up data set flags, note the current page and
c        print a header line.  header looks like:
c typetypetypetype observations of   body   made from sitename
c series=....  (spot=....)  frequency= ................
c        a second line is also printed for ncodf.gt.20
c        this header also appears on every subsequent page of the series
c
c        input to obsbeg:
c   typobs - 16-character name of observation type
c   site - 8-character site name
c
c        after setup is done, the program should call obsbgt to get the
c        elapsed time printed out and the header read back from "intern"
c        and printed.
c        no arguments for obsbgt.
c
c        for every observation record (except a below-horizon dummy) the
c        main routine should call obscnt to form residuals, accumulate
c        error quantities, and check for full page.  ordinarily, no
c        printout appears except the usual one line per record, but
c        if any other routines also print on iout, then
c        whenever the line count exceeds 57, the responsible
c        program should call entry obspag.  this will print out any
c        deletions on the current page and start a new one.  the main
c        program must have called obsbgt already to set up the
c        observation data header.  it is also printed.
c
c   lab1 - name of first observable
c   l1 - length (in words) of lab1.  should be 2-4
c   lab2 - name of second observable
c   l2 - length of lab2.
c
c        at the end of a series the main program should call obsend.
c                  this will print out any
c        deletions on the last page and then print out the usual trailer
c        detailing the type of observation, the pages used, the average
c        number of iterations needed, and a header for the error
c        analysis summary which appears at the end of the series.
c        no arguments for obscnt, obsend or obspag.
c
c array dimensions
      include 'globdefs.inc'

c commons
      include 'comdat.inc'
      character*8 pname
      equivalence (Comcon(127),pname)
      include 'eqnphs.inc'
      real*4    eqnx(3)
      equivalence (eqnx,Pnox)
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'ltrapx.inc'
      include 'namtim.inc'
      include 'number.inc'
      integer*2 nsite(2)
      equivalence (nsite,Nsite1)
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*4    acctim, accdst, accprc
      equivalence (acctim,Estf),(accdst,Estf(2)),
     .            (accprc,Estf(3))
      real*10 ctrecf,fdev,freq2
      equivalence (Dstf,ctrecf),(Dstf(10),fdev)
      equivalence (Dstf(5),freq2)
      real*10 th0(2),th1(2),th2(2),erob1(2),erob2(2),erth0(2),
     .          erth1(2),erth2(2)
      equivalence (Erquan,th0),(Erquan(3),th1),
     .            (Erquan(5),th2),(Erquan(7),erob1),
     .            (Erquan(9),erob2),(Erquan(11),erth0),
     .            (Erquan(13),erth1),(Erquan(15),erth2)
      real*10 erstf(3,2)
      equivalence (Erquan(17),erstf(1,1))
      include 'obstap.inc'
      include 'obstuf.inc'
      include 'redobs.inc'
      include 'sitcrd.inc'
      include 'skymap.inc'
      include 'statsrad.inc'

c local variables
      integer*4 mast(2)
      character*4 and/' AND'/
      character*2 bl2
      character*4 n1/'    '/,n2/'    '/,blnk
      character*8 blank/'        '/, rcsnd(4)/' RECEIVE','    SEND',
     1 '   FIRST','  SECOND'/, point/'. . .   '/, sit1,sit2,pnam2
      character*1 space,sitsep
      equivalence (blank,blnk,space,bl2)
 
      character*4 setlin(12)/'  SE','TUP ','FOR ',4*' ',' OBS','ERVA',
     .          'TION', ' SER', 'IES '/, forobs(10)
      character*4 stpobs(7)
      equivalence (setlin(4),forobs(2),stpobs)
      character*1 chdr(128,4),skip/'0'/
      equivalence (Hdr(1,1),chdr(1,1))
      real*10 fnobs,plsps,ww
      integer*4 i,icall,il1,il2,ix,j,jtypob,kit,klap2,
     . linf,mpage,n,nf,nn,npart,nprtl,ns1,ns2,nsv
 
c allow up to 4 iteration averages
      real*4    f1, f2, avgit(4)
      character*1 rpar/')'/,comma/','/,sep(4)
      real*4    a(2,9)
c external functions
      integer*4 ITYPOB
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c initialize all three data sets
      Niobc  = 0
      Niobs  = 0
      Niabs1 = 0
 
c initialize quantities not always read from iiobs
      Ctatb(1,1)  = 32.15_10
      Ututsb(1,1) = 0.
      Clampb(1,1) = bl2
      Limbb(1,1)  = bl2
      Limbb(2,1)  = bl2
      Obsrvb(1,1) = bl2
      Nsavb2  = 0
      mast(1) = 0
      mast(2) = 0
      mpage   = Npage - 1
      do i = 1, 4
         Svtpbs(i) = typobs(i)
         stpobs(i) = typobs(i)
      end do
      do i = 1, 10
         Labsv(i) = blnk
      end do
      il1 = min0(l1,4)
      ix = max0(0,4 - il1)
      do i = 1, il1
         Labsv(i + ix) = lab1(i)
      end do
      il2 = min0(l2,4)
      do i = 1, il2
         Labsv(i + 5) = lab2(i)
      end do
      sit1 = site
      if(sit1.eq.blank) sit1 = point
      sit2   = Sitf(2)
      sitsep = comma
      if(sit2.eq.blank) sitsep = space
      pnam2 = pname
      klap2 = Klap
      if(Nplnt2.eq.Ncp0 .and. Klan.gt.0) klap2 = Klan
      if(Klans1.gt.0) klap2 = Klans1
      if(klap2.gt.0) pnam2  = Aplnt(klap2)
      icall = 1
c*   start=100
c
c now write standard series header line
  100 write(Iout,200) stpobs,pname,sit1,sitsep,sit2,Series,
     .                 Spotf, Freq
  200 format('0', 7A4, 'S OF ', a8, ' MADE FROM ', a8, a1, a8,
     .       '  SERIES=', a4, '  (SPOT=', a4, ')  FREQUENCY=',
     .       1pd20.13)
 
c increment line counter  for this header + observation header
      Line = Line + 2
      if(Nplnt2.ne.0) then
         write(Iout,250) pnam2,Nspot2,Spotf2,freq2
  250    format(37x,'SECOND OBSERVED OBJECT IS ', a8, '  NSPOT2=', i3,
     .          '  (SPOT=', a4, ')  FREQUENCY=', 1pd20.13)
         Line = Line + 1
      endif
  300 if(icall.eq.2 .or. icall.eq.4) goto 1100
      if(icall.eq.3) goto 1300
      if(icall.eq.5) return
c
c write site information, if any
      jtypob = ITYPOB(Ncodf)
      if(jtypob.ne.3) then
         ns1 = 2
         if(Sitf(1).ne.blank) ns1 = 1
         ns2 = 1
         if(Sitf(2).ne.blank) ns2 = 2
         if(ns1.le.ns2) then
            nf = 0
            if(jtypob.eq.4) nf= 2
            if(ns2.eq.1) nf= 0
            nsv=0
            do j=ns1,ns2
               do i=4,6
                  if(Coords(i,j).ne.0._10 .or. Lsite(i,j).gt.0) nsv=1
               end do
            end do
            if(nsv.eq.0) then
               write(Iout,320) (rcsnd(j+nf),Sitf(j),nsite(j),
     .          (Lsite(i,j),i=1,3),(Coords(i,j),i=1,3),Ksite(j),
     .          j=ns1,ns2)
  320          format('0', 17x, 'NAME  NSITE LSITE', 5x, 'RADIUS', 7x,
     .          'LONGITUDE',6x,'LATITUDE  KSITE'/
     .          (2x,a8,' SITE=',a8,i4,1x,3I2,f15.9,f14.8,f15.8,i3))
            else
               write(Iout,325) (rcsnd(j+nf),Sitf(j),nsite(j),
     .          (Lsite(i,j),i=1,6),(Coords(i,j),i=1,6),Ksite(j),
     .          j=ns1,ns2)
  325          format('0',17x,'NAME  NSITE   LSITE',9x,'RADIUS',7x,
     .          'LONGITUDE',6x,'LATITUDE',8x,'VELOCITY',16x,'KSITE'/
     .          (2x,a8,' SITE=',a8,i4,1x,6I2,f15.9,f14.8,f15.8,
     .          3f10.4,i3))
            endif
            Line = Line + 3 + ns2 - ns1
         endif
      endif
      write(Iout,400) accprc,accdst,acctim,
     .                 (i,Erwgta(i,Ntape),i = 1,2),fdev
  400 format('0ACCPRC=', 1pe10.3, '  ACCDST=', e10.3, '  ACCTIM=',
     .       e10.3, 2('  ERWGT',i1,'=',e10.3), '  FDEV=', 3pd12.3)
      fdev = fdev + 1._10
      Line = Line + 2
c write series biases, if any
c*   start=300
      if(Neqnox.gt.0) then
         write(Iout,450) Neqnox,(i,Leqnx(i),i = 1,3),eqnx
  450    format('0NEQNOX=', i4, 3(' LEQNX(',i1,')=',i2), ' DEQUINOX=',
     .          1pe12.5, ' DEQUATOR=', e12.5, ' DLATITUDE=', e12.5)
         Line = Line + 2
      endif
      if(Nstar.gt.0) then
         nn = Nskyc/4
         write(Iout,500) Nstar,Ctlg,nn,(Sky(i),i = 1,Nskyc)
  500    format('0NSTAR=', i3, ' (', a8, ') N=', i3, ' SKY=', 1p,
     .          10D10.2/(1x,1p,13D10.2))
         write(Iout,550) (Lsky(i),i = 1,Nskyc)
  550    format(' LSKY= ', 80I1)
         Line = Line + 3 + (Nskyc + 3)/13
      endif
      if(Nplsr.gt.0) then
         plsps = Plsper + Psrprm(6)
         write(Iout,600) Nplsr,Jdps0,plsps,Lpsrx
  600    format('0NPLSR=', i3, ' JD0=', i8, ' PER=', 1pd22.15, '  L=',
     .          24I3)
         Line = Line + 2
      endif
      if(Nphase.gt.0) then
         write(Iout,650) Nphase,(i,Lphsx(i),i = 1,9),
     .                    (blnk,i,Aphs(i),i = 1,Ncph)
  650    format('0NPHASE=', i4, 9(' LPHSX(',i1,')=',i2)
     .          /(12x,5(a1,'APHASE(',i1,')=',1pe12.5)))
         Line = Line + 3 + Ncph/6
      endif
      if(Nrbias.gt.0) then
         write(Iout,700) Nrbias,Lrbsx,Rbsx
  700    format('0NRBIAS=', i4, '  LRBIAS=', 2I2, '   RBIAS=', 1p,
     .          2E13.5)
         Line = Line + 2
      endif
 
c write control integers for series
      write(Iout,800) Nspot,Kob,Iwob
  800 format('0  NSPOT =', i3, '   NPREC =', i3, '   NDPREC=', i3,
     .       '   NTMDLY=', i3, '   NDOP  =', i3, '   NMEDIA=', i3,
     .       '   NEATM =', i3, '   NPATM =', i3, '   NPSHP =', i3,
     .       '   NPLNG =', i3/'   NEION =', i3, '  KOB(11)=', i3,
     .       '   LOPTRF=', i3, '   NDDIFF=', i3, '   NLIBPR=', i3,
     .       '  KOB(15)=', i3, '   CALCVL=', i3, '   ITIME =', i3,
     .       '   NTIME =', i3, '    IWOB =', i3)
      Line = Line + 3
 
c write abbreviated output on auxilliary data set
      if(Mout.gt.0) then
         write(Mout,850) Ntapa(Ntape),Ncodf,Nplnt0,pname,Freq,
     .                    Nseqa(Ntape),Series,sit1,sit2,Spotf
  850    format(' NTAPE =', i3, '   NCODF =', i3, '   NPLNT0 =', i2,
     .          2x, a8, '   FREQ =', 1pd22.15/' NSEQ =', i5,
     .          '   SER = ', a4, '   SITE1 = ', a8, '   SITE2 = ', a8,
     .          '   SPOT = ', a4)
         if(Nplnt2.ne.0) then
            write(Mout,860) Nplnt2,pnam2,Spotf2,freq2
  860       format(' NPLNT2=', i2, 2x, a8, '   SPOT2=', a4,
     .             '   FREQ2=', 1pd22.15)
         endif
      endif
      if(Nout.gt.0) write(Nout,900)
  900 format(/
     .' JD  HR MIN SEC    HR    RESULT(1)  ERROR   OBS-TH  RESULT(2)  ER
     .ROR   OBS-TH')
 
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=500
c write elapsed time and data header
      entry OBSBGT
      icall = 5
      rewind Intern
      Nlhd = 0
      do while( .true. )
         read(Intern,950,end=1000) (Hdr(i,Nlhd+1),i = 1,16)
  950    format(16A8)
         Nlhd = Nlhd + 1
         if(Nlhd.ge.6) goto 1000
      end do
 
c found end of header array
 1000 chdr(1,1) = skip
      nprtl = Nlhd
 
c get actual number of print lines
      do j = 1, Nlhd
 
         if(chdr(1,j).eq.skip) nprtl = nprtl + 1
      end do
      rewind Intern
 
c now print out elapsed time
      call TIMRIT(setlin,12)
c
c now print out observation data header
 1100 write(Iout,950) ((Hdr(i,j),i=1,16),j = 1,Nlhd)
      Line = Line + nprtl
      if(icall.eq.1 .or. icall.eq.2 .or. icall.eq.5) return
      if(icall.eq.3) goto 1300
      if(icall.eq.4) goto 2100
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*   start=600
c entry obspag for just page footing
      entry OBSPAG
      icall = 2
 
c do page logic
 1200 if(Nast(1) + Nast(2).gt.0) then
         call EBCDI(Nast(1),n1,3)
         call EBCDI(Nast(2),n2,3)
         write(Iout,1250) n1,(Labsv(ix+i),i=1,il1),and,
     .                     n2, (Labsv(i+5),i=1,il2)
 1250    format(
     .' * MEASUREMENTS ON THIS PAGE DELETED IN LEAST SQUARES ANALYSIS:'
     ., 16A4)
         Line = Line + 1
         do i = 1, 2
            mast(i) = mast(i) + Nast(i)
            Nast(i) = 0
         end do
      endif
      if(Line.lt.37 .and. icall.eq.3) goto 300
      call NEWPG
      goto 100
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*   start=700
      entry OBSEND
      icall = 3
      goto 1200
 
c print info for end of series
 1300 linf    = Npage - mpage
      Nast(1) = mast(1)
      Nast(2) = mast(2)
      do i = 1, 4
         sep(i) = comma
      end do
      f2 = Ncard(1)
      if(Ncard(1).eq.0) f2=1.
      kit = 0
      avgit(1) = 0.
c 2-spacecraft iterations are accumulated in nit(17-18), while ordinary
c observations are in nit(19-20).
      do i = 1, 4
         if(Nit(21-i).le.0) goto 1400
         kit = kit + 1
         f1  = Nit(21 - i)
         avgit(kit) = f1/f2
      end do
 1400 if(kit.eq.0) kit = 1
      sep(kit) = rpar
      write(Iout,1500) forobs,linf,mpage,Ncard(1),
     .                  (avgit(i),sep(i),i = 1,kit)
 1500 format('0ERROR ANALYSIS ', 10A4, 'EXTENDING FOR', i4,
     .       ' PAGES STARTING ON PAGE',i5/
     .       ' AVERAGE NUMBER OF ITERATIONS FOR', i6,
     .       ' OBSERVATIONS WAS (', 4(f9.5,a1))
      write(Iout,1600) Labsv
 1600 format(30x,20A4)
      Line = Line + 4
c
c printout error analysis for observation series
      write(Iout,1700) Nast,Nobs
 1700 format(' NUMBER OF MEASUREMENTS DELETED ', i10,
     .       i14/' NUMBER OF MEASUREMENTS INCLUDED', i10, i14)
      do i = 1, 2
         if(Nobs(i).gt.0) then
            fnobs = Nobs(i)
            if(Nobs(i).le.0) fnobs = 1._10
            a(i,1) = erob1(i)/fnobs
            a(i,2) = SQRT(erob2(i)/fnobs)
            a(i,3) = erth0(i)/fnobs
            a(i,4) = erth1(i)/fnobs
            ww = erth2(i)/fnobs
            a(i,6) = ww
            a(i,5) = SQRT(ww)
            a(i,7) = th0(i)/fnobs
            a(i,8) = th1(i)/fnobs
            a(i,9) = SQRT(th2(i)/fnobs)
         else
            do j = 1, 9
               a(i,j) = 0.
            end do
         endif
      end do
      write(Iout,1800) ((a(i,j),i=1,2),j = 7,9),
     .                  ((a(i,j),i=1,2),j = 1,6),erth2
 1800 format(16x,'AVERAGE (OBS-TH)', 1p, 2E14.5/
     . 13x,'AVERAGE ABS(OBS-TH)', 2E14.5/
     . 6x,' ROOT MEAN SQUARE (OBS-TH)', 2E14.5/
     . 19x,'AVERAGE ERROR', 2E14.5/
     . 9x,' ROOT MEAN SQUARE ERROR', 2E14.5/
     . 10x,'AVERAGE (OBS-TH)/ERROR', 2E14.5/
     . 7x,'AVERAGE ABS(OBS-TH)/ERROR',2E14.5/
     . ' ROOT MEAN SQUARE (OBS-TH)/ERROR', 2E14.5/
     . 5x,'AVERAGE ((OBS-TH)/ERROR)**2', 2E14.5/
     . 5x,'    SUM ((OBS-TH)/ERROR)**2', 2E14.5)
      Line = Line + 12
c
c increment error analysis quantities for observations of a
c given type of a given body
      Nit(5) = Nit(5) + Ncard(1)
      Nit(6) = Nit(6) + Ncard(2)
      do i = 1, 2
         Nit(i + 8)  = Nit(i + 8) + Nast(i)
         Nit(i + 10) = Nit(i + 10) + Nobs(i)
         erstf(1,i) = erstf(1,i) + erth0(i)
         erstf(2,i) = erstf(2,i) + erth1(i)
         erstf(3,i) = erstf(3,i) + erth2(i)
      end do
c
c statistics for error reads
      write(Iout,1900) Iiobcn,Iiobs,Iiabs1,(stat(i),i = 1,3)
 1900 format(' NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON=', i2,
     .       ', IOBS=', i2, ', IABS1=', i2,
     .       ') FOR THIS OBSERVATION SERIES WERE (', i5, ',', i5, ',',
     .       i5, ')')
      do i = 1, 3
         do j = 3, 9, 3
            stat(j + i) = stat(j + i) + stat(i)
         end do
         stat(i) = 0
      end do
c
c printout timer information for observation series
      npart = Numpar - 2
      write(Iout,2000) Ncard,npart
 2000 format('0', i6, ' +', i4,
     .   ' OBSERVATION CARDS PROCESSED (NUMBER OF PARTIAL DERIVATIVES='
     .   , i3, ')')
      Line = Line + 3
      call TIMRIT('  PROCESSING OBSERVATION SERIES ', 8)
      Nast(1) = 0
      Nast(2) = 0
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=1000
c keep count of observation records and set up numobs
      entry OBSCNT
      icall = 4
      if(Jd.le.-10) return
      if(Jd.le.0) then
 
c not processed, but message printed
         Ncard(2) = Ncard(2) + 1
         Line     = Line + 1
      endif
 
c page logic
      if(Line.ge.56) goto 1200
c
c*  start=1100
c error calculations
 2100 if(Jd.le.0) return
      Ncard(1) = Ncard(1) + 1
c
c calculate error quantities for observations
c this code is derived from old subroutine obserr
c (m.e.ash    oct 1966)
      n = 1
      if(Nice.gt.0) goto 2300
 2200 Deriv(2,n) = Result(n) - Deriv(2,n)
      nast1(n)    = 1
      j = Ncod(n)
      if(ABS(Deriv(2,n)).lt.Eps(j)*Deriv(1,n)) then
         Nobs(n)  = Nobs(n) + 1
         th0(n)   = th0(n) + Deriv(2,n)
         th1(n)   = th1(n) + ABS(Deriv(2,n))
         th2(n)   = th2(n) + Deriv(2,n)**2
         erob1(n) = erob1(n) + Deriv(1,n)
         erob2(n) = erob2(n) + Deriv(1,n)**2
         ww = Deriv(2,n)/Deriv(1,n)
         erth0(n) = erth0(n) + ww
         erth1(n) = erth1(n) + ABS(ww)
         erth2(n) = erth2(n) + ww*ww
      else
         nast1(n) = 2
         Nast(n)  = Nast(n) + 1
      endif
      if(n.eq.2) goto 2400
 2300 n = 2
      if(Nice.ge.0) goto 2200
c
c*  start=1200
c set up numobs for writing output tape
 2400 if(Nice.lt.0) then
         Num1   = 1
         Num2   = 1
         Numobs = 1
      else if(Nice.eq.0) then
         Num1   = 1
         Num2   = 2
         Numobs = 2
      else
         Num1   = 2
         Num2   = 2
         Numobs = 1
      endif
c
c*  start=1300
c done 1st point in series,
c set nk1=1 if partl not to be called
      if(Ict(1).le.0) then
         if(Nk1.le.0) Nk1 = 1
      else if(Idumob.eq.1 .and. Ict(3).lt.0) then
         if(Nk1.le.0) Nk1 = 1
      endif
c
c
c*   start=9000
      return
      end
