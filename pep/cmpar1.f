      subroutine CMPAR1(jtape,mocpar)
 
      implicit none
 
c amuchastegui/ash - april 1970 - subroutine cmpar1

c array dimensions
      include 'globdefs.inc'
c common
      include 'comdat.inc'
      include 'dltflg.inc'
      include 'dtparm.inc'
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'loadid.inc'
      include 'lcntrl.inc'
      include 'leon.inc'
      include 'mnsprt.inc'
      include 'ltrapx.inc'
      include 'mtrapx.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'obscrd.inc'
      real*10 freq2
      equivalence (Dstf(5),freq2)
      real*4    acctim
      equivalence (Estf,acctim)
      include 'kobequiv.inc'
      include 'obsdta.inc'
      include 'obstap.inc'
      include 'param.inc'
      include 'sitcrd.inc'
      integer*4 zsitcr/3308/   !r8=1684,r10=3308
      include 'skymap.inc'
      include 'spqind.inc'
      include 'stats.inc'
 
c temproary storage
      common/YVECT/ Uucond(u_nmbod),Vvcond(u_nmbod),Vvcone(u_nmbod,2),
     . Spcdxx(6,2),Nmobst(2,2),Uname(2),Rbsu(2),Eqnu(3),Phsu(9),Jevv0(2)
      character*4 Nmobst, Uname
      real*10 Uucond,Vvcond,Vvcone,Spcdxx
      real*4    Rbsu,Eqnu,Phsu
      integer*4 Jevv0
 
      include 'zeroes.inc'
      character*8 czer8
      character*4 czer4
      equivalence (Zero(1),czer8,czer4)
c
c local
      integer*2 mczon1, mctes1, mumtr1, mszon1, mstes1, numdt1
      integer*4 lhprct(8)/1,2,3,4,5,6,22,23/
      character*8 namobs/' PLNOBS '/
      integer*2 ncentu, nrbius, neqnux
      integer*2 mcobs1, mcobs2, ntest, i2d,i2e,i2f, nsitcr/6/,nsptcr/6/
      real*10 jd2000/2451545._10/
      character*8 blank8/'        '/
      character*4 blank
      equivalence (blank8,blank)
      character*8 exnams
      integer*4 ione/1/
      character*1 htype, colon/':'/, vflg
      character*4 al15/'AL15'/, ap15/'AP15'/
      integer*4 i, j, jdcc, jduu0, jdvv0, jobf, jobs, jobt,
     . jtape, jtypob, k, mocpar, ndt9, next, niabs9, niobc9,
     . niobs9, npagea, npshpt, ntape0, ntest2
c external functions
      integer*4 ITYPOB, LEG
c
c*  start=1000
c read first record of obs.series on observation library tape
  100 if(Niabs1.ne.0) goto 400
      if(Iabs1s.le.0) goto 400
      Niabs1 = -1
 
c test for end of file in obsred
      if(Iiabs1.eq.-1) goto 300
 
c clear all undefined m-vectors (Mprmx, Memx, Mmnx, Merx, Mmrx, Mplx,
c Mdtx, Mumdtx, Msitcr, and Msptcr already set)
      call ZFILL(Mplx,2*2558)
 
      read(Iabs1s,end=300) Nseqa(3),Ncoda(3),Nplnta(3),Sita1(3),
     . Sita1(4),Sera(3),Sita2(3),Sita2(4),Spota(3),(Erwgta(i,3),i=1,2),
     . Acctma(3),Fdeva(3),Freqa(3),Itima(3),Nrewna(3),npagea,Uname,
     . ncentu,jduu0,(Uucond(i),i=1,Ncnmo),(Mplx(i),i=1,Ncnmo),
     . jdvv0,(Vvcond(i),i=1,Ncnmo),(Msbx(i),i=1,Ncnmo),
     . ((Mscrdx(i,j),i=1,Msitcr),j=1,2),nrbius,Rbsu,
     . Mrbsx,neqnux,Eqnu,Meqnx,Ncphu,Phsu,Mphsx,Nobcon,(Obscon(i),i=1,
     . MAX(1,0+Nobcon)),((Scord1(i,j,3),i=1,Msitcr),j=1,2),
     . Kscrd1(1,3),Kscrd1(2,3),(Spcdxx(i,1),i=1,Msptcr),
     . (Mspcdx(i,1),i=1,Msptcr),
     . Mczone,mczon1,(Mczhar(i),i=1,mczon1),Mctess,mctes1,(Mcchar(i),
     . i=1,mctes1),(Mcshar(i),i=1,mctes1),Mumtar,mumtr1,(Mtrg(i),
     . (Mtbod(j,i),j=1,Ncnmo),Mtzone(i),(Mtzhar(j,i),j=1,4),Mttess(i),
     . (Mtchar(j,i),j=1,5),(Mtshar(j,i),j=1,5),i=1,mumtr1),mcobs1,
     . mcobs2,(Nmobst(1,j),Nmobst(2,j),(Mscx(i,j),i=1,mcobs2),Jevv0(j),
     . (Vvcone(i,j),i=1,mcobs2),j=1,mcobs1),Mszone,mszon1,(Mszhar(i),
     . i=1,mszon1),Mstess,mstes1,(Mschar(i),i=1,mstes1),(Msshar(i),
     . i=1,mstes1),Nplta2(3),Spota2(3),(Spcdxx(i,2),i=1,Msptcr),
     . (Mspcdx(i,2),i=1,Msptcr),Freqa2(3),Memmn,
     . next,(exnams,exnams,i=1,MAX(1,0+next)),Mngd,Mnfour,Ctlga(3),
     . Mskyc,(Msky(i),i=1,MAX(1,0+Mskyc)),i2d,(Mpsrx(i),i=1,i2d),
     . Mmpex,i2e,i2f,(Mplex(j),(Mpex(i,j),i=1,i2f),j=1,i2e),
     . T0st1(1,3),T0st1(2,3),T0sp1
 
      if(Mumtar.gt.i_mxtrg)
     .     call SUICID('TOO MANY TARGET BODIES ON INPUT OBSLIB, '//
     .    'STOP IN CMPAR1  ', 14)
 
      if(spota2(3).eq.czer4) spota2(3)=blank
      if(Mskyc.le.1) Mskyc = 0
c fill in site velocity if missing
      if(Msitcr.lt.6) then
         do j=1,2
            T0st1(j,3)=jd2000
            do i=Msitcr+1,6
               Mscrdx(i,j)=0
               Scord1(i,j,3)=0._10
            end do
         end do
      endif
c fill in spot velocity if missing
      if(Msptcr.lt.6) then
         do j=1,2
            T0sp1(j)=jd2000
            do i=Msptcr+1,6
               Mspcdx(i,j)=0
               Spcdxx(i,j)=0._10
            end do
         end do
      endif
      if(Ncnmo.lt.u_nmbod) then
         do i=Ncnmo+1,u_nmbod
            Mplx(i)=0
            Msbx(i)=0
            do j=1,mumtr1
               Mtbod(i,j)=0
            end do
         end do
      endif
      if(mcobs1.gt.0 .and. mcobs2.lt.u_nmbod) then
         do j=1,mcobs1
            do i=mcobs2+1,u_nmbod
               Mscx(i,j)=0
            end do
         end do
      endif
      if(i2e.gt.0 .and. i2f.lt.u_nmbod) then
         do j=1,i2e
            do i=i2f+1,u_nmbod
               Mpex(i,j)=0
            end do
         end do
      endif
      if(Ncoda(3).le.0) goto 300
      call SEQCHK(3,.false.,i)
      if(i.le.0) then
         Niabs1 = 1
         goto 400
      endif
  200 do while( .true. )
         read(Iabs1s,err=200) ntest
         if(ntest.le.0) goto 100
      end do
  300 rewind Iabs1s
      Itrwnd(Iabs1s) = 0
      Iabs1s = 0
 
c reset obsred tape number
      Iiabs1 = 0
c
c*  start=2000
c read first card of observation series
  400 if(Niobs.ne.0 .or. Iobs.le.0) goto 700
      Niobs    = -1
      Fdeva(2) = 0.
      Ctlga(2) = blank8
      if(Jct(69).gt.0) then
 
c new format header card a
         read(Iobs,450,end=600) Ncoda(2),Nplnta(2),htype,
     .    Sita1(2),Sera(2),Sita2(2),Spota(2),(Erwgta(i,2),i=1,2),
     .    Acctma(2),Itima(2),Freqa(2),Nrewna(2),Ntapa(2),Nseqa(2)
  450    format(i2,i3,a1,3(a4,1x),a4,3E6.0,i2,d18.11,i2,i3,i4)
         if(LEG(1,1,htype,1,colon).eq.0)
     .     read(Iobs,460,end=600) htype,Fdeva(2),Ctlga(2)
  460    format(a1,f7.2,a8)
      else
 
c old format header card a
         read(Iobs,500,end=600) Ncoda(2),Nplnta(2),Sita1(2),
     .    Sera(2),Sita2(2),Spota(2),(Erwgta(i,2),i=1,2),Acctma(2),
     .    Itima(2),Fdeva(2),Freqa(2),Nrewna(2),Ntapa(2),Nseqa(2)
  500    format(i2,i3,4(1x,a4),3E6.0,i2,f7.2,d18.11,i2,i3,i5)
      endif
      Nplta2(2) = 0
      Spota2(2) = blank
      if(Ncoda(2).gt.20)
     .  read(Iobs,550) Nplta2(2),Spota2(2),Freqa2(2)
  550 format(i3,a4,d18.11)
      Itrwnd(Iobs) = 1
      if(Ncoda(2).gt.0) then
         call SEQCHK(2,Ict(80).lt.0,i)
         Niobs = 1
 
c read site data from iobs if so indicated
         if(Nrewna(2).lt.0 .or. Nrewna(2).gt.1) then
            do j = 1, 2
               read(Iobs,560) Site1(j,2),(Scord1(i,j,2),i = 1,3),vflg,
     .          Kscrd1(j,2)
  560          format(a8,3f16.9,1x,a1,12x,i2)
               if(vflg.eq.'6') then
                  read(Iobs,565) (Scord1(i,j,2),i=4,6),T0st1(j,2)
  565             format(8x,3f16.9,f8.0)
               else
                  T0st1(j,2)=jd2000
                  do i=4,6
                     Scord1(i,j,2)=0._10
                  end do
               endif
               if(Ncoda(2).gt.3) goto 700
            end do
         endif
         goto 700
      endif

c If IOBS is a separate file, keep reading it to the end.
c This guards against the truncation of input by an extra blank line.
      if(Iobs.ne.In .and. (Nplnta(2).eq.0 .and. Sita1(2).eq.blank 
     . .and. Sera(2).eq.blank .and. Sita2(2).eq.blank)) then
         Niobs=0
         goto 400
      endif
 
c end of observations on iobs
  600 if(Iobs.ne.In) rewind Iobs
      Itrwnd(Iobs) = 0
      if(Ict(80).ge.0) Iobs = 0
c
c*  start=3000
c read corrections to observation series
  700 if(Niobc.ne.0) goto 1200
 
c negative tape sequence numbers are skipped at label=3020
  800 Niobc = -1
      read(Iobcon,end=1100) Ncoda(1),Nplnta(1),Sita1(1),Sera(1),
     . Sita2(1),Spota(1),(Erwgta(i,1),i=1,2),Acctma(1),Itima(1),
     . Fdeva(1),Freqa(1),Nrewna(1),Ntapa(1),Nseqa(1),Nplta2(1),
     . Spota2(1),Freqa2(1),Ctlga(1),Gncode,Gnplt1,Gsite1,Gseres,
     . Gsite2,Gspot,Gerwgt,Gacctm,Gitime,Gfdev,Gfreq,Gnrwnd,
     . Gnplt2,Gspot2,Gfreq2,Gctlg
      Itrwnd(Iobcon) = 1
      if(Ncoda(1).le.0 .and. Nseqa(1).eq.0) goto 1100
      if(Ntapa(1).gt.0) goto 1000
      if(mocpar.gt.0) goto 1000
      if(Ict(79).ge.0) goto 1000
 
c skip observation series records on iobcon
      if(Nrewna(1).lt.0 .or. Nrewna(1).gt.1) read(Iobcon)
  900 do while( .true. )
         read(Iobcon,err=900) jdcc
         if(jdcc.le.0) goto 800
      end do
 
c consistency checks for iobcon
 1000 call SEQCHK(1,.false.,i)
      Niobc = 1
 
c read site data from iobcon if so indicated
      if(Nrewna(1).lt.0 .or. Nrewna(1).gt.1) read(Iobcon)
     . (Site1(j,1),(Scord1(i,j,1),i=1,6),Kscrd1(j,1),T0st1(j,1),j=1,2)
      goto 1200
 1100 rewind Iobcon
      Itrwnd(Iobcon) = 0
c
c
c*  start=4000
c           logic to determine ntape,mtape
c ntape=1  data read from iobcon overrides that read from iobs,iabs1
c          (whether iobs,iabs1 are read depends on mtape)
c ntape=2  data read from iobs overrides that read from iabs1  (iobcon
c          is not read, whether iabs1 is read depends on mtape)
c ntape=3  data read from iabs1, but not from iobcon,iobs
c mtape=-1 dummy observations generated from iobcon (ntape must be 1)
c mtape=0  data read from both iobs and iabs1 (ntape can be 1 or 2)
c mtape=1  data read from iobs, not from iabs1 (ntape can be 1 or 2)
c mtape=2  data read from iabs1, not from iobs (ntape can be 1 or 3)
c niobc,niobs,and/or niabs1 are set to zero if mtape,ntape indicate
c that iobcon,iobs and/or iabs1, respectively, are to be read beyond
c first record of observation series
 1200 Iiobcn = 0
      Iiobs  = 0
      Iiabs1 = 0
      if(mocpar.le.0) then
c
c look for next pair (tape,seq)
         jobf   = 0
         jobt   = 999999
         jobs   = 999999
         Idumob = -1
         if(Niobc.le.0) goto 1300
 
c iobcon exists, check if real or dummy
         if(Ntapa(1).gt.0) then
c
c only real observations from here to label=5000
c
c iobcon best so far
            jobt = Ntapa(1)
            jobs = Nseqa(1)
            jobf = 1
            goto 1300
         endif
      else
         Ncodf = 0
         if(Niobc.lt.0) return
         if(Ntapa(1).gt.0) return
      endif
 
c dummy observations
      Mtape  = -1
      Ntape  = 1
      Niobc  = 0
      Iiobcn = Iobcon
      Idumob = 1
      goto 1600
 1300 if(Niobs.gt.0) then
 
c iobs exists, test it
         if(Ntapa(2).gt.jobt .or.
     .       (Ntapa(2).eq.jobt .and. Nseqa(2).gt.jobs)) goto 1400
         if(Ntapa(2).lt.jobt .or.
     .       (Ntapa(2).eq.jobt .and. Nseqa(2).lt.jobs)) goto 1350
 
c iobs ok, include it
         jobf = jobf + 2
         goto 1400
 
c iobs best so far, reject any other match
 1350    jobt = Ntapa(2)
         jobs = Nseqa(2)
         jobf = 0
         jobf = jobf + 2
      endif
 1400 if(Niabs1.gt.0) then
 
c iabs1 exists, test it
         if(Ntapa(3).gt.jobt .or.
     .       (Ntapa(3).eq.jobt .and. Nseqa(3).gt.jobs)) goto 1500
         if(Ntapa(3).lt.jobt .or.
     .       (Ntapa(3).eq.jobt .and. Nseqa(3).lt.jobs)) goto 1450
 
c iabs1 ok, include it
         jobf = jobf + 4
         goto 1500
 
c iabs1 best so far, reject any other match
 1450    jobt = Ntapa(3)
         jobs = Nseqa(3)
         jobf = 0
         jobf = jobf + 4
      endif
c
c
c interpret jobf; binary coded flags
c first, check if any input at all
 1500 if(jobf.eq.0) goto 2300
 
c check for tape sequence jump
      if(Niabs1.le.0 .and. Ict(80).ge.0 .and.
     .    jobt.gt.Ntapsb .and. Ntapsb.gt.0) goto 2300
      if(jobf.gt.1) then
c
c
c           all series on a single obslib tape have, by definition,
c           the same value of ntape.  thus, unlesss ict(80) is
c           negative (no input obslib) or the series is a dummy
c           (already checked), a jump in ntape from iobs should
c           force a new output tape.
c
c           determine tape numbers:
c           iiobcn,iiobs,iiabs1 are iobcon,iobs,iabs1 or 0
c           according to whether they are read or not
c           set ntape,mtape from jobf
         Mtape = mod(jobf/2,3)
         if(jobf.ge.4) then
 
c iabs1 included
            Ntape  = 3
            Iiabs1 = Iabs1s
            Niabs1 = 0
         endif
         if(mod(jobf,4).ge.2) then
 
c iobs included, overrides iabs1
            Ntape = 2
            Iiobs = Iobs
            Niobs = 0
         endif
         if(mod(jobf,2).ne.0) then
 
c iobcon included, overrides all others
            ntape0 = Ntape
            Ntape  = 1
            Iiobcn = Iobcon
            Niobc  = 0
 
c fill blank iobcon fields with values from obslib (or cards)
            if(.not. Gncode) Ncoda(1)  = Ncoda(ntape0)
            if(.not. Gnplt1) Nplnta(1) = Nplnta(ntape0)
            if(.not. Gsite1) Sita1(1)  = Sita1(ntape0)
            if(.not. Gseres) Sera(1)   = Sera(ntape0)
            if(.not. Gsite2) Sita2(1)  = Sita2(ntape0)
            if(.not. Gspot) Spota(1)   = Spota(ntape0)
            if(.not. Gerwgt(1)) Erwgta(1,1) = Erwgta(1,ntape0)
            if(.not. Gerwgt(2)) Erwgta(2,1) = Erwgta(2,ntape0)
            if(.not. Gacctm) Acctma(1) = Acctma(ntape0)
            if(.not. Gitime) Itima(1)  = Itima(ntape0)
            if(.not. Gfdev) Fdeva(1)   = Fdeva(ntape0)
            if(.not. Gfreq) Freqa(1)   = Freqa(ntape0)
            if(.not. Gnrwnd) Nrewna(1) = Nrewna(ntape0)
            if(.not. Gnplt2) Nplta2(1) = Nplta2(ntape0)
            if(.not. Gspot2) Spota2(1) = Spota2(ntape0)
            if(.not. Gfreq2) Freqa2(1) = Freqa2(ntape0)
         endif
      else
 
c extra iobcon cards, skip these
         if(Iterat.le.1) then
            if(Line.gt.50) call NEWPG
            Line = Line + 3
            write(Iout,1520) (i,Ntapa(i),i,Nseqa(i),i = 1,3)
 1520       format('0 SERIES OUT OF ORDER OR EXTRA ON IOBCON, SKIP.'/3
     .             ('   NTAPE(',i1,')=',i4,' NSEQ(',i1,')=',i5))
         endif
         goto 900
      endif
c if(.not.gctlg) ctlga(1)=ctlga(ntape0)
c
c*  start=5000
c initialize quantities for this series
 1600 npshpt = Npshp1
      call ZFILL(Nsite1,2*50)
      Npshp1 = npshpt
 
c clear all undefined l-vectors (Lprmx, Lemx, Lmnx, Lerx, Lmrx, Ldtx,
c and Numdtx already set)
      call ZFILL(Lplx,2*2558)
 
c clear freq and site quantities
      call ZFILL(Xsite,zsitcr)
 
      Nplnt0 = Nplnta(Ntape)
      Freq   = Freqa(Ntape)
      Nrewnd = Nrewna(Ntape)
      Ncodf  = Ncoda(Ntape)
      Nplnt2 = Nplta2(Ntape)
      freq2  = Freqa2(Ntape)
      Ctlg   = Ctlga(Ntape)
      if(Ctlg.eq.czer8) Ctlg=blank8
      Series = Sera(Ntape)
      if(Spota(Ntape).eq.czer4) Spota(Ntape)=blank
      if(Spota2(Ntape).eq.czer4) Spota2(Ntape)=blank
c alsep 15 and the apollo 15 laser retroreflector are the same
c logical spot. (see also coding in cmpar2, cmpar3, nrmict, prdict.)
      if(Spota(Ntape).eq.al15) Spota(Ntape)   = ap15
      if(Spota2(Ntape).eq.al15) Spota2(Ntape) = ap15
      Spotf  = Spota(Ntape)
      Spotf2 = Spota2(Ntape)
 
      if(Ncodf.le.0) return
 
      do i = 1, 18
         Kob(i) = -1
      end do
      ctime=Itima(Ntape)
      itime=Itima(Ntape)
      if(ctime.ge.0) then
         ctime=ctime/10
         itime=itime-ctime*10
      else
         jtypob=ITYPOB(Ncodf)
         if(jtypob.ne.2 .and. jtypob.ne.3) then
            ctime=0
            itime=2
         else
            itime=0
         endif
      endif
      acctim = Acctma(Ntape)
c
c check on output observation library tape sequence numbers
      if(Ict(80).ge.0) then
         i = Ntape
         call SEQCK1(i,4)
         Ntapsb = Ntaps(4)
      endif
c*  start=5100
c
c           calculate lx vectors for parameters affecting every observ.
c       use of 'life' for status of iabs2
c        -2 - checkpoint restart just begun
c        -1 - nothing done yet
c         0 - l-vectors merged for header (if needed)
c         1 - header written
      if(Life.lt.0) then
         if(Life.lt.-1) Life = 1
         if(Life.eq.-1) Life = 0
         if(Ict(1).lt.0) goto 1700
         if(Ict(1).ne.0) then
            Iabs1 = Iabs1s
            if(Ict(4).gt.0) Iabs1 = 0
            call LVTPRM(Lprmx,Lprm,Mprmx,100)
            call LVTBDY(Lemx,Lem,Memx,u_nmbod)
            call LVTBDY(Lmnx,Lmn,Mmnx,u_nmbod)
            call LVTBDY(Lerx,Ler,Merx,u_nmbod)
            do i = 1, 8
               Lhprc(i) = 0
               do j = 7, u_nmbod
                  if(Lerx(j).le.0) goto 1620
                  if(Lerx(j).eq.lhprct(i)) then
                     if(Iabs1.gt.0) then
                        do k = 7, u_nmbod
                           if(Merx(k).eq.lhprct(i)) goto 1620
                        end do
                     endif
                     Lhprc(i) = 1
                     goto 1620
                  endif
               end do
 1620       end do
            call LVTBDY(Lmrx,Lmr,Mmrx,u_nmbod)
            Numdtx = Numdt
            if(Numdtx.gt.0) then
               if(Jddt0.le.0) Numdtx = Numdt + 400
               ndt9 = -Numdtx
               call LVTBDY(Ldtx,Ldt,Mdtx,ndt9)
            endif
         endif
      else if(Life.ne.0) then
         goto 1700
      endif
c
c*  start=5300
c write first two records of output observation library tape
      if(Life.le.0) then
         if(Iabs2.gt.0) then
            if(Mtape.ge.0 .or. Ict(3).gt.0) then
               Life   = 1
               numdt1 = Numdt
               if(numdt1.gt.0 .and. Jddt0.le.0) numdt1 = Numdt + 400
               if(numdt1.le.0) numdt1 = 1
               write(Iabs2) namobs,Heding,Date,Lnklvl
               write(Iabs2) Ntapa(Ntape),Npage,Iterat,u_nmprm,u_nmbod,
     .          Jdem0,Jdmn0,Jder0,Jdmr0,prmter,Econd,Mcond,Ercond,
     .          Mrcond,Lprmx,Lemx,Lmnx,Lerx,Lmrx,
     .          Numdt,numdt1,(Jddt(i),i=1,numdt1),(Dt(i),i=1,numdt1),
     .          (Ldtx(i),i=1,numdt1),ione,izero(1),
     .          Jddt0,Lnklvl,nsitcr,nsptcr,izero,izero,izero,
     .          izero,izero,izero,izero,izero
            endif
         endif
      endif
c
c*  start=5500
c set flag for source of partials
 1700 Iabs1 = Iiabs1
      if(Ict(4).gt.0) Iabs1 = 0
c
c determine if observation series is to be skipped
      if(Ntape.gt.1) return
      if(Mtape.lt.0) return
      if(Erwgta(1,1).gt.-1E3) return
      if(Erwgta(2,1).gt.-1E3) return
      niobc9 = 0
      niobs9 = 0
      niabs9 = 0
      if(Mtape.eq.1) niabs9 = 1
      if(Mtape.eq.2) niobs9 = 1
 1800 if(niobc9.le.0) then
         read(Iobcon) jdcc
         if(jdcc.le.0) niobc9 = 1
      endif
 1900 if(niobs9.le.0) then
         read(Iobs,1950) ntest,ntest2
 1950    format(2I1)
         if(Ncoda(2).ge.7 .and. Ncoda(2).le.9) ntest = ntest + ntest2
         if(ntest.le.0) niobs9 = 1
      endif
 2000 if(niabs9.le.0) then
         read(Iabs1s,err=2100) ntest
         if(ntest.le.0) niabs9 = 1
      endif
      goto 2200
 2100 read(Iabs1s)
 2200 if(niobc9.le.0) goto 1800
      if(niobs9.le.0) goto 1900
      if(niabs9.le.0) goto 2000
      goto 100
c
c*  start=9950
c skip to next obslib tape
 2300 jtape  = 0
      Ntapsb = -999999
c
c*  start=9990
      return
      end
