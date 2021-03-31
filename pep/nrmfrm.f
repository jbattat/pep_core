      subroutine NRMFRM(icv,iptr,b)
 
      implicit none

c     m.ash  feb 1970   r.reasenberg may 1973   subroutine nrmfrm
c     main program for restoring normal equations from saved normal eqs.
c       changed to use qread and iptr oct 1970  r.reasenberg
c     modified for filter oct 1970  r.reasonberg
c     further changes jan 1974 by d.white
c
c         parameters
c     iptr = integer vector to hold order of rows of normal eqs
c     b    = lower diagonal half of symmetric coefficient matrix of
c            restored normal equations (passed through calling sequence
c            rather than through common /nrmmat/ for filter)
c     icv  = input/output control integers (see below for explanation)
c        icv(1)    >0        jtape<=1
c                            imats = imat1 (or imat2, see icv(5))
c                  =0        0<jtape<=nummt0
c                            imats = imat0(jtape)
c                  <0        apriori format(title+ptrs+equations)
c                            imats = icv(4)
c        icv(2)=niobc
c        niobc =-1 do not read iobcon
c        niobc = 1 iobcon has been read
c        niobc = 0 ready to read iobcon
c        icv(3)   not used
c        icv(4)    imats + 100 * jmata
c                  on return, set to last-used imats
c        icv(5) = 0 normal call
c               = 1 set imats to imat2 if sne already formed on imat2
c              .ge.2 call only to get iptr from eqns on imat2
c                      (imats = imat2)
c                   also return mesmt1,ermes1 from imat2
c        icv(6) = 0 old logic, no DPSNEC
c               > 0 binary coded DPSNEC control
c                  1 bit: copy all nominals from imat0(1)
c                  2 bit: print free-parameter differences
c                  4 bit: print fixed-parameter differences
c                 64 bit: suppress rhs correction
c        icv(7) = binary flags for data/a priori restoration
c                  1: restoring a priori normal eqns
c                  2: print out summaries
c                  4: saved a prioris already encountered
c        jmata
c        0   restore n.e.
c        1   skip r.h.s. of n.e.  (increment left side)
c                   (also skip error statistics)
c        2   skip l.h.s. of n.e.  (increment right side)
      real*10 b(10)
      integer*2 iptr(10),icv(12)
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'aprtbf.inc'
      include 'ciptrx.inc'
      include 'dltflg.inc'
      include 'dpsnec.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'iptrst.inc'
      include 'mcnfrm.inc'
      include 'nrmwgt.inc'
      include 'obsdta.inc'
      include 'restor.inc'
      include 'rtside.inc'
c
c local
      integer*4 jmata, jmatb/-1/
      real*10 ai,ctlg1,dum,freq1,freq2,rms,wgtfct,wgtk
      integer*4 i,ia,iaprio,ic,icprm,icv74,ihedp,imats,iobk,
     . ippr,iprtap,iuseap,ivectk,jdccc,jtape,kall,key,m2,n1,niobc,
     . niptr,nparm,nseq1,nseq2,nsize,ntape1,ntape2
      equivalence (ai,ia)
      real*4    erwgt1(2),erwgt2(2)
      character*4 site1,ser1,site2,spot1,spot2
      real*4    fdev1,acctm1
      integer*2 ncodf1,nplnt1,itime1,nrewn1,nplnt2
      character*40 mesg/'  RESTORING NORMAL EQNS OF ORDER*****   '/
      character*4 parens/'(  )'/, subscr, blank/'    '/
      character*8 sidtyp(3) /'BOTH','LEFT','RIGHT'/
      character*4 iob(2) /'NO','YES'/
      character*5 tmat(4) /'IMATS','IMAT0','IMAT1','IMAT2'/, qmat
      character*120 phdr/'        NTAPE1    NSEQ1     MESMT1     SUM O-C
     ./ERR  SUMABS O-C/ERR  SUMSQS O-C/ERR         WEIGHT  IOBCON   SIDE
     .    RMS '/
      equivalence (phdr,qmat)
c
c initialize
      kall   = icv(1)
      jtape  = 0
      wgtfct = 1._10
      subscr = blank
 
c set up print logic
      iuseap = icv(7)
      iuseap = mod(iuseap,2)
      iprtap = mod(icv(7)/2,2)
      if(Jct(60).lt.0) iprtap = 1
      icv74 = icv(7) - 4*mod(icv(7)/4,2) + 4
      if(iuseap.eq.0 .and. kall.eq.0 .and. Nummt0.ne.0) then
         do i = 1,Nummt0
            Iaprsv(i) = 0
         end do
      endif
 
c set up DPSNEC logic
      icprm = 0
      Icdif = icv(6)
      if(Icdif.gt.0) Icdif = Icdif + 4096
      if(icv(5).gt.1 .or. icv(1).ne.0) Icdif = 0
      Itdif = 0
  100 do while( .true. )
c
c bump jtape
         jtape = jtape + 1
         niptr = 0
         ippr  = 0
c
c switch on kall to get imats
         if(kall.lt.0) then
            imats = icv(4)
            qmat  = tmat(1)
            goto 200
         else if(kall.eq.0) then
            if(jtape.gt.Nummt0) goto 1400
            if(iuseap.le.0 .or. Iaprsv(jtape).ne.0) then
               qmat   = tmat(2)
               subscr = parens
               call EBCDIX(jtape,subscr,2,2)
               imats  = Imat0(jtape)
               wgtfct = Wgtmt0(jtape)
               if(wgtfct.ne.0._10) goto 200
               if(iprtap.ne.0) then
                  call PAGCHK(60,2,0)
                  write(Iout,110) qmat,subscr,imats,wgtfct
  110             format('0INPUT FROM ',a5,a4,'=',i3,
     .                   ' SKIPPED: WEIGHT=',f9.4)
               endif
               goto 600
            endif
         else
            qmat  = tmat(4)
            imats = Imat2
            if(icv(5).le.0) then
               qmat  = tmat(3)
               imats = Imat1
            endif
            goto 200
         endif
      end do
  200 call PAGSET(phdr,-30)
      if(imats.gt.0) then
c
c strip jmata off
         jmata  = imats/100
         imats  = imats-jmata*100
         icv(4) = imats
c
c if restoring same side on subsequent tape, do not rewind iob
         if(jmata.ne.jmatb) then
            jmatb = jmata
c
c check status of iobcon and imats
            niobc = icv(2)
            if(niobc.ge.0) then
               rewind Iobcon
               Itrwnd(Iobcon) = 0
            endif
         endif
         Itrwnd(imats) = 1
c
c read header records  (assign weight for bz)
         ihedp = 0
 
c if read apriori, use newhed
         if(kall.lt.0) then
            call NEWHED(imats,jtape,0,iptr,iprtap,Apssq1)
            Weight = 1._10
            Itdif  = 0
            Znsqsn = 0._10
            ntape1 = 0
            nseq1  = 0
            Mesmt1 = 0
            do i = 1,3
               Ermes1(i) = 0._10
            end do
            Zsmsn  = 0._10
            iaprio = 0
            ivectk = 0
            goto 1300
         else
c
c if making second or later pass in multiple parameter set mode,
c skip reading of m-vector
            if(.not. (Hedskp)) then
               if(mod(Icdif,2).eq.1 .and. jmata.eq.0) icprm = 1
 
c skip i.c. copy if DPSNEC operating
               key = 1
               if(icv(6).gt.0) key = 2
               call FRMHED(imats,qmat,subscr,key,ippr,iprtap)
               call PAGSET(phdr,-30)
 
c compare parameter values
               Itdif = 0
               if(Icdif.gt.0 .and. jmata.eq.0)
     .             call DIFNOM(Icdif,Itdif,Prmdif)
               if(icprm.gt.0) Icdif = (Icdif/2)*2
               if(mod(Icdif/64,2).eq.1 .or. Ibuf1.le.0) Itdif = 0
            else
               call FRMHED(imats,qmat,subscr,0,ippr,iprtap)
               call PAGSET(phdr,-30)
            endif
            goto 400
         endif
      else
         call PAGCHK(60,2,0)
         write(Iout,250)
         if(Mout.gt.0) write(Mout,250)
  250    format('0ERROR IN NRMFRM, NEGATIVE DATA SET')
         goto 600
      endif
c
c* start=300
c finished reading one set/series of equations
  300 if(kall.lt.0) goto 500
c
c read first record of series for saved normal equations
  400 read(imats,end=500) Mesmt1,ntape1,nseq1,
     . (erwgt1(i),i = 1,2),(Ermes1(i),i = 1,3),Znsqsn,Zsmsn,
     . iaprio,ivectk,Apssq1
      if(ntape1.ne.0) goto 700
      if(Mesmt1.gt.0) goto 700
  500 rewind imats
      Itrwnd(imats) = 0
  600 if(iprtap.ne.0) then
         call PAGCHK(60,2,0)
         write(Iout,650) qmat,subscr,imats
  650    format(' FINISHED RESTORING NORMAL EQUATIONS FROM ',a5,a4,
     .          '=',i3/)
      endif
      if(kall.eq.0) goto 100
      goto 1400
  700 iobk = 1
c
c*  start=500
c read iobcon for overridding error weighting
      if(niobc.lt.0) then
c
c*  start=700
c no overriding of error weighting for this series
         Weight = wgtfct
         goto 1100
      else if(niobc.ne.0) then
         goto 900
      endif
  800 read(Iobcon) ncodf1,nplnt1,site1,ser1,site2,spot1,
     .             (erwgt2(i),i = 1,2),acctm1,itime1,fdev1,freq1,
     .             nrewn1,ntape2,nseq2,nplnt2,spot2,freq2,ctlg1,
     .             Gncode,Gnplt1,Gsite1,Gseres,Gsite2,Gspot,
     .             Gerwgt,Gacctm,Gitime,Gfdev,Gfreq,Gnrwnd,
     .             Gnplt2,Gspot2,Gfreq2,Gctlg
      Itrwnd(Iobcon) = 1
      if(ncodf1.le.0 .and. nseq2.eq.0) then
c
c end of iobcon data set
         niobc = -1
         rewind Iobcon
         Itrwnd(Iobcon) = 0
         Weight = wgtfct
         goto 1100
      else
         do while( .true. )
            read(Iobcon) jdccc
            if(jdccc.le.0) then
               niobc = 1
               goto 900
            endif
         end do
      endif
  900 if(ntape2.lt.ntape1) goto 800
      if(ntape2.eq.ntape1) then
         if(nseq2.lt.nseq1) goto 800
         if(nseq2.eq.nseq1) then
            if((ntape2.le.0) .and. (Ict(7).le.0)) goto 800
            niobc = 0
            if(.not. Gerwgt(1)) erwgt2(1) = erwgt1(1)
            if(.not. Gerwgt(2)) erwgt2(2) = erwgt1(2)
c
c override error weighting for this series
            if((erwgt2(1).gt.0.E0) .or. (erwgt2(2).gt.0.E0)) then
               n1 = 1
               if(erwgt1(1).le.0.) n1 = 2
               do i = n1,2
                  if(erwgt2(i).gt.0.) then
                     if(abs(erwgt1(i)-erwgt2(i))/erwgt2(i).le.1E-3)
     .                   then
                        Weight = wgtfct
                     else
                        Weight = erwgt1(i)
                        Weight = wgtfct*(Weight/erwgt2(i))
                        iobk   = 2
                     endif
                     goto 1100
                  endif
               end do
            endif
         else
            Weight = wgtfct
            goto 1100
         endif
      else
         Weight = wgtfct
         goto 1100
      endif
c
c skip this series entirely
 1000 call BSKIP(imats,Mparam)
      if(iaprio.gt.0) call BSKIP(imats,Mparam)
      if(iprtap.ne.0) then
         call PAGCHK(60,1,1)
         write(Iout,1050) qmat,subscr,imats,nseq1,Mesmt1
 1050    format(1x,a5,a4,'=',i3,'  NSEQ1=',i5,'  MESMT1=',i8,
     .          '  SERIES SKIPPED')
      endif
      goto 300
c
c apply error weight, if any, to statistics
 1100 Ermes1(1) = Ermes1(1)*Weight
      Ermes1(2) = Ermes1(2)*Weight
      Ermes1(3) = Ermes1(3)*Weight*ABS(Weight)
      Apssq1    = Apssq1*Weight*ABS(Weight)
      if(Weight.lt.0._10) Mesmt1 = -iabs(Mesmt1)
      rms = 0._10
      if(Mesmt1.ne.0) rms = SQRT(ABS(Ermes1(3)/Mesmt1))
c
c printout restoring normal equations message
      if(iprtap.ne.0) then
         if(ihedp.ne.1) then
            if(Line.ge.50) call NEWPG
            call PAGHED(0)
            ihedp = 1
         endif
         call PAGCHK(60,1,1)
         if((erwgt1(1).le.0.E0) .and. (erwgt1(2).le.0.E0))
     .      goto 1000
         if(jmata.ne.1) then
c accumulating statistics and at least right-hand side
            write(Iout,1150) imats,ntape1,nseq1,
     .       Mesmt1,Ermes1,Weight,iob(iobk),sidtyp(jmata+1),rms
 1150       format(i6,i9,i9,i11,1p,3D16.5,1pd15.5,5x,a3,3x,a5,g11.4)
            if(Apssq1.ne.0._10) then
               call PAGCHK(60,1,1)
               write(Iout,1160) Apssq1
 1160          format(' A PRIORI CONTRIBUTION:',t68,1pd16.5)
            endif
         else
c skipping right-hand side and statistics
            write(Iout,1200) imats,ntape1,nseq1,
     .       Mesmt1,Weight,iob(iobk),sidtyp(jmata+1)
 1200       format(i6,i9,i9,i11,48x,1pd15.5,5x,a3,3x,a5)
         endif
      endif
 
      wgtk   = Weight
      Weight = Weight*ABS(Weight)
c
c read pointer group for pre-reduced sne
      m2    = Mparam*2
      Zsmpp = 0._10
      if(ippr.ne.0) then
         if(ivectk.eq.0) read(imats) nparm,Ncparm,Ndparm,
     .                            Znsqpp,Nserpp,Nauxpp,
     .                            (Icrest(i),i = 1,Ncparm),
     .                            (Idrest(i),i = 1,Ndparm)
         if(ivectk.gt.0) read(imats) nparm,Ncparm,Ndparm,
     .                            Znsqpp,Nserpp,Nauxpp,
     .                            (Icrest(i),i = 1,Ncparm),
     .                            (Idrest(i),i = 1,Ndparm),
     .                            (dum,i = 1,m2),Zsmpp
      endif
c
c if making second or later pass in multiple parameter set mode,
c copy pseudo-iptr from iptrx
      if(Hedskp) then
         do i = 1,Mparam
            iptr(i) = Iptrx(i)
         end do
c
c           setup iptr vector
c                 iptr(ia) is the location in side (or in b) of the
c                 element ia of the vector sav
c     niptr = 0 iptr not yet setup for this tape
c     niptr = 1 iptr already setup for this tape
c                 mparam found in frmhead read number 2
c                 nparam in /fcntrl/
      else if(niptr.le.0) then
         if(kall.ge.0) niptr = 1
         call ZFILL(Sigma,16*Nparam)
         call ZFILL(iptr,2*Mparam)
         do i = 1,Mparam
            ia     = i
            Sav(i) = ai
         end do
         call FRMMVE(Sigma,Nparam)
         do i = 1,Nparam
            ai = Sigma(i)
            if(ia.ne.0) iptr(ia) = i
         end do
         if(Itdif.gt.0) then
c
c scale parameter differences and setup for rhs
            read(Ibuf1)
            read(Ibuf1) (ai,Sigma(i),i = 1,Nparam)
            rewind Ibuf1
            do i = 1,Nparam
               Sigma(i) = Prmdif(i)/Sigma(i)
            end do
            nsize = Mparam
            if(ippr.ne.0) nsize = Ncparm
            call ZFILL(Prmdif,16*nsize)
            do ic = 1,nsize
               ia = ic
               if(ippr.ne.0) ia = Icrest(ic)
               i = iptr(ia)
               if(i.gt.0) Prmdif(ic) = Sigma(i)
            end do
         endif
c
c if multiple parameter set mode is to be invoked in this run,
c pass iptr to subroutine analix via common ciptrx
         if(Jct(51).gt.0) then
            do i = 1,Mparam
               Iptrx(i) = iptr(i)
            end do
         endif
      endif
c
c*  start=1000
c increment equations
 1300 nsize = Mparam
      if(ippr.ne.0) then
         nsize = Ncparm
         do i = 1,nsize
            Icrest(i) = iptr(Icrest(i))
         end do
      endif
      if(iaprio.gt.0) icv(7) = icv74
      if(iuseap.eq.0 .and. kall.eq.0) Iaprsv(jtape) = iaprio
      if(icv(5).le.1) then
         if(iaprio.gt.0 .and. Imat2.eq.0 .and. Jct(60).le.0)
     .       call SUICID(
     .'RESTORING IMBEDDED A PRIORI, BUT IMAT2=0. STOP IN NRMFRM',14)
         if(ippr.ne.0) then
            if(iuseap.le.0 .or. iaprio.ne.0)
     .       call NRMADD(b,Sav,Icrest,Ncparm,imats,Mesmt1,Weight,
     .       ntape1,nseq1,jmata,Itdif,Prmdif,Znsqpp+Znsqsn,
     .       Side,Ermes1,Sumzns,Ermeas,Measmt,iuseap,iaprio,
     .       ivectk,Vectk,wgtk,Sumzsm,Zsmsn+Zsmpp,Apssq1,Sumaps)
            goto 500
         else
            call NRMADD(b,Sav,iptr,Mparam,imats,Mesmt1,Weight,
     .       ntape1,nseq1,jmata,Itdif,Prmdif,Znsqsn,
     .       Side,Ermes1,Sumzns,Ermeas,Measmt,iuseap,iaprio,
     .       ivectk,Vectk,wgtk,Sumzsm,Zsmsn,Apssq1,Sumaps)
            goto 300
         endif
      else
         rewind imats
         Itrwnd(imats) = 0
         call PAGCHK(60,1,0)
         write(Iout,1350)
 1350    format(' M-VECTORS RESTORED TO GET IPTR; B, SIDE NOT RESTORED')
         return
      endif
c
c* start=2000
c print out timer information
 1400 if(niobc.ge.0) then
         rewind Iobcon
         Itrwnd(Iobcon) = 0
      endif
      call EBCDIX(Nparam,mesg,33,5)
      if(iprtap.eq.0) return
      call TIMRIT(mesg,10)
      return
      end
