      subroutine PRDOBS(nrmprd,itype)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 fdedit,freq1,frq2,tdlt
      integer   i,ijend,itype,ITYPOB,j,jbcc,jdedit,jterat,
     .          match,natch,nbkfor,np,npaga,nrmprd,nseq1,
     .          ntape1
 
c*** end of declarations inserted by spag
 
 
c           subr. prdobs - j.f.chandler - 1983 feb
c     read obslibs and iobcon for nrmict and prdict.
c     handle logic for input (and output) tapes, overriding errors,
c     and data editing.
c
c  nrmprd= code for the caller: 1-nrmict, 2-prdict
c  itype = type of record desired:
c          1-get next tape and read the type 1 and type 2 records
c            (iabs1=0 if none)
c          2-(no effect)
c          3-read the type 3 record, set up pointers and match iobcon,
c            (ncodf=0 if none)
c          4-read type 4 record and set error weights (ncode=0 if none)

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'bernum.inc'
      include 'dltflg.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'loadid.inc'
      include 'ktrap.inc'
      integer*2 numobs,mspcx1(3),mspcx2(3),imdt(6)
      equivalence (Num2,numobs),(Mspcrd(1,1),mspcx1),
     .            (Mspcrd(1,2),mspcx2),(imdt(1),Imdt1)
      include 'namobs.inc'
      include 'nrmwgt.inc'
      include 'obsdta.inc'
      include 'prdmov.inc'
      include 'wrkcompr.inc'
      character*8 scnam(2),sitf8(2),pname
      equivalence (scnam(1),Nmobst(1,1)),(sitf8(1),Sitf(1,1)),
     .            (pname,Plnnam)
      include 'zeroes.inc'
      character*4 czer4,lnklvm
      character*8 czer8
      equivalence (zero(1),czer4,czer8)
c
c data from records of data set iobcon
      integer*2 ncodf1,nplnt1,itime1,nrewn1,nqlnt2,ihrc(2),
     .          iminc(2)
      integer*4 jdc(2)
      character*8 ctlg1
      character*4 site1,ser1,site2,spot1,spot2
      real*4 acctm1,fdev1
      real*4    erwgt1(2),secc(2),erobsc(2)
      real*10 fdysc1(2)
 
      character*8 blank8/'        '/
      character*4 blank
      equivalence (blank8,blank)
      real*4    cut(2)/0.,-1E3/
      integer*2 i2d,i2e,i2f,m2p,mmpex1,mmpex2
      integer*2 n16/16/,n6/6/
      real*10 jd2000/2451545._10/
 
      if(nrmprd.lt.0) then
         if(itype.eq.1) then
c
c write first two records of output observation library tape
            if(Iabs2.gt.0) then
               write(Iabs2) Tapnam,Heding,Date,Lnklvl
               write(Iabs2) Ntape,Npage,Iterat,u_nmprm,u_nmbod,
     .          Jdem9,Jdmn9,Jder9,Jdmr9,Prmtr9,Econd9,Mcond9,Ercnd9,
     .          Mrcnd9,Mprm,Mem,Mmn,Mer,Mmr,Mumdt,
     .          Mumdt1,(Jddtm(i),i=1,Mumdt1),
     .          (Dtm(i),i=1,Mumdt1),(Mdt(i),i=1,Mumdt1),
     .          Nlabel,(Label(j),j=1,Nlabel),Jddtm0,Lnklvl,n6,
     .          izero,izero,izero,izero,izero,izero,
     .          izero,izero
               Itrwnd(Iabs2) = 1
            endif
         else if(itype.eq.2) then
         else if(itype.eq.3) then
c
c write first record of observation series
            if(Iabs2.gt.0) then
               np = Next
               if(Next.le.0) np = 1
               m2p = Mskyc
               if(Mskyc.le.0) m2p = 1
               if(Mmpex.gt.0) then
                  mmpex1=Mmpex
                  mmpex2=u_nmbod
               else
                  mmpex1=1
                  mmpex2=1
               endif
               write(Iabs2) Nseq,Ncodf,Nplnt0,sitf8(1),Series,sitf8(2),
     .          Spota,Ugta,Acctim,Fdev,Freq,Itime,Nrewnd,Npage,pname,
     .          Ncentb,Jdpp0,Ppcond,Mpl,Jdbb0,Bbcond,Msb,Mscrd,Nrbias,
     .         Rbs,Mrbs,Neqnox,Eqn,Meqn,Ncph,Pas,Mphs,Nobcon,(Obscon(i),
     .          i=1,Nobcon),Scoord,Ksite,Spcdx,mspcx1,Mczone,Mczon1,
     .          (Mczhar(i),i=1,Mczon1),Mctess,Mctes1,(Mcchar(i),i=1,
     .          Mctes1),(Mcshar(i),i=1,Mctes1),Mumtar,Mumtr1,(Mtrg(i),
     .          (Mtbod(j,i),j=1,u_nmbod),Mtzone(i),(Mtzhar(j,i),j=1,4),
     .          Mttess(i),(Mtchar(j,i),j=1,5),(Mtshar(j,i),j=1,5),i=1,
     .          Mumtr1),Mcobs1,Mcobs2,(scnam(j),(Msc(i,j),i=1,Mcobs2),
     .          Jevv0(j),(Vvcone(i,j),i=1,Mcobs2),j=1,Mcobs1),Mszone,
     .          Mszon1,(Mszhar(i),i=1,Mszon1),Mstess,Mstes1,(Mschar(i),
     .          i=1,Mstes1),(Msshar(i),i=1,Mstes1),Nplnt2,Spota2,Spcdx2,
     .          mspcx2,Freq2,Memmn,np,(Exnams(1,j),Exnams(2,j),j=1,np),
     .          Mngd,Mnfour,Ctlgm,m2p,(Msky(i),i=1,m2p),n16,Mpsrx,
     .          Mmpex,mmpex1,mmpex2,(Mplex(j),(Mpex(i,j),i=1,mmpex2),
     .          j=1,mmpex1),T0st1,
     .         izero,izero,izero,izero,izero,izero,izero,izero
 
               if(Ncodf.le.0) then
c
c end of input tape - close output
                  rewind Iabs2
                  Itrwnd(Iabs2) = 0
               endif
            endif
         else if(itype.eq.4) then
c
c write new observed minus theory tape
            if(Iabs2.gt.0) write(Iabs2) Ncode,Ihr,Imin,Sec,
     .       (Result(j),Error(j),j=1,2),Atuts,Ututs,Clamp,Limb,
     .       Observ,Imonth,Iday,Iyear,Jds,Jd,Ctat,Ctrecf,Numsav,
     .       (Save(i),i=1,Numsav),numobs,Numpar,
     .       ((Deriv(i,j),i=1,Numpar),j=1,numobs),imdt,Nmp2,
     .       ((Dvx(i,j),i=1,Nmp2),j=1,numobs),Ncal,
     .       (Cal(i),Scal(i),Ical(i),i=1,Ncal),(Sumcor(i),i=1,2),
     .       Mnshob,Mixshp,(Mshobs(i),i=1,Mnshob),
     .       (izero(i),i = 1,50)
         else
            goto 100
         endif
         return
      endif
  100 if(itype.eq.2) return
      if(itype.eq.3) goto 500
      if(itype.eq.4) goto 1900
      do while(.true.)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c*  start=1000
c           type 1 and 2 records of next tape, if any
         Jtape = Jtape + 1
         ijend = 0
         Iabs1 = 0
         if(Jtape.gt.Numobt) return
         if(nrmprd.ne.1 .or. Wgtobs(Jtape).gt.0._10) then
            Iabs1 = Iobs1(Jtape)
            Iabs2 = Iobs2(Jtape)
            if(Ict(80).ge.0 .and. Iterat.eq.2*(Iterat/2)) then
               Iabs1 = Iobs2(Jtape)
               Iabs2 = Iobs1(Jtape)
            endif
            if(Iterat.le.1) then
               if(nrmprd.eq.2) then
                  if(Jtape.le.Ict(13)) Iabs1 = Iobs0(Jtape)
               else
                  if(Ict(80).gt.0) Iabs1 = Iobs0(Jtape)
               endif
            endif
            if(Iabs1.gt.0) then
               if(Ict(11).eq.-1 .or. nrmprd.eq.1) Iabs2 = 0
c
c read first two records of input observation library tape
               Itrwnd(Iabs1) = 1
               read(Iabs1,err=200) Title
               goto 300
            endif
         endif
      end do
  200 read(Iabs1) Title
  300 call ZFILL(Mdt,2*600)
      read(Iabs1) Ntape,npaga,jterat,Nprmo,Ncnmo,
     . Jdem9,Jdmn9,Jder9,Jdmr9,(Prmtr9(i),i=1,Nprmo),
     . (Econd9(i),i=1,Ncnmo),(Mcond9(i),i=1,Ncnmo),
     . (Ercnd9(i),i=1,Ncnmo),(Mrcnd9(i),i=1,Ncnmo),
     . (Mprm(i),i=1,Nprmo),(Mem(i),i=1,Ncnmo),(Mmn(i),i=1,Ncnmo),
     . (Mer(i),i=1,Ncnmo),(Mmr(i),i=1,Ncnmo),
     . Mumdt,Mumdt1,(Jddtm(i),i=1,Mumdt1),
     . (Dtm(i),i=1,Mumdt1),(Mdt(i),i=1,Mumdt1),Nlabel,
     . (Label(j),j=1,Nlabel),Jddtm0,lnklvm,Msitcr
      if(Nprmo.gt.u_nmprm .or. Ncnmo.gt.u_nmbod) call SUICID(
     . 'BAD NUMBER OF PARAMETERS, STOP IN PRDOBS',10)
c in principle, should set default parameter values if the input tape has
c smaller set of defined parameters
      if(Nprmo.lt.u_nmprm) then
         do i=Nprmo+1,u_nmprm
            Mprm(i)=0
         end do
      endif
      if(Ncnmo.lt.u_nmbod) then
         do i=Ncnmo+1,u_nmbod
            Mem(i)=0
            Mmn(i)=0
            Mer(i)=0
            Mmr(i)=0
         end do
      endif
      call DTCHCK(Jddtm0,Mumdt,Mumdt1,Dtm,Mdt)
      if(Nlabel.le.0) Nlabel = 1
      if(Ntape.lt.0 .and. Ict(3).gt.1) Ntape = -Ntape
      if(Msitcr.ne.6) Msitcr=3
 
      call PAGCHK(60,4,0)
      write(Iout,400) Iabs1,Title,lnklvm,Ntape,jterat,npaga,
     .                 Iabs2
  400 format('-INPUT OBS.LIB.TAPE IABS1=',i3,' TITLE=',a88/
     .       '  LEVEL=',a4,'  NTAPE=',i4,' JTERAT=',i3,
     .       ' NPAGE=',i5,40x,'(OUTPUT IS IABS2=',i3,')')
      if(Nlabel.gt.1) then
         call PAGCHK(60,2,0)
         write(Iout,450) Label
  450    format('0CORRELATED OBSERVATION RUN TITLE= ',20A4)
      endif
      return
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c*  start=3000
c           read first record of observation series
  500 if(ijend.le.0) then
         call ZFILL(Nsite1,2*50)
         call ZFILL(Mpl,2*1054)
         call ZFILL(Msc,2*1490)
         read(Iabs1,end=600) Nseq,Ncodf,Nplnt0,Sitf(1,1),Sitf(2,1),
     .    Series,Sitf(1,2),Sitf(2,2),Spota,Ugta,Acctim,Fdev,Freq,Itime,
     .    Nrewnd,npaga,Plnnam,Ncentb,Jdpp0,(Ppcond(i),i=1,Ncnmo),
     .    (Mpl(i),i=1,Ncnmo),Jdbb0,(Bbcond(i),i=1,Ncnmo),
     .    (Msb(i),i=1,Ncnmo),
     .    ((Mscrd(i,j),i=1,Msitcr),j=1,2),i2d,Rbs,Mrbs,i2d,Eqn,Meqn,
     .    Ncph,Pas,Mphs,Nobcon,(Obscon(i),i=1,Nobcon),
     .    ((Scoord(i,j),i=1,Msitcr),j=1,2),Ksite,Spcdx,mspcx1,
     .    Mczone,Mczon1,(Mczhar(i),i=1,Mczon1),Mctess,Mctes1,
     .    (Mcchar(i),i=1,Mctes1),(Mcshar(i),i=1,Mctes1),Mumtar,Mumtr1,
     .    (Mtrg(i),(Mtbod(j,i),j=1,Ncnmo),Mtzone(i),(Mtzhar(j,i),j=1,4),
     .    Mttess(i),(Mtchar(j,i),j=1,5),(Mtshar(j,i),j=1,5),i=1,Mumtr1),
     .    Mcobs1,Mcobs2,(Nmobst(1,j),Nmobst(2,j),(Msc(i,j),i=1,Mcobs2),
     .    Jevv0(j),(Vvcone(i,j),i=1,Mcobs2),j=1,Mcobs1),Mszone,Mszon1,
     .    (Mszhar(i),i=1,Mszon1),Mstess,Mstes1,(Mschar(i),i=1,Mstes1),
     .    (Msshar(i),i=1,Mstes1),Nplnt2,Spota2,Spcdx2,mspcx2,Freq2,
     .    Memmn,Next,(Exnams(1,j),Exnams(2,j),j=1,Next),Mngd,Mnfour,
     .    Ctlgm,Mskyc,(Msky(i),i=1,Mskyc),i2d,(Mpsrx(i),i=1,i2d),
     .    Mmpex,i2e,i2f,(Mplex(j),(Mpex(i,j),i=1,i2f),j=1,i2e),T0st1
 
         if(Spota2.eq.czer4) Spota2=blank
         if(Ctlgm.eq.czer8) Ctlgm=blank8
         if(Mskyc.le.1) Mskyc    = 0
         if(Nlabel.le.1) Next    = 0
c fill in site velocity if missing
         if(Msitcr.lt.6) then
            do j=1,2
               T0st1(j)=jd2000
               do i=Msitcr+1,6
                  Mscrd(i,j)=0
                  Scoord(i,j)=0._10
               end do
            end do
         endif
         if(Ncnmo.lt.u_nmbod) then
            do i=Ncnmo+1,u_nmbod
               Mpl(i)=0
               Msb(i)=0
               do j=1,Mumtr1
                  Mtbod(i,j)=0
               end do
            end do
         endif
         if(Mcobs1.gt.0 .and. Mcobs2.lt.u_nmbod) then
            do j=1,Mcobs1
               do i=Mcobs2+1,u_nmbod
                  Msc(i,j)=0
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
         Jtypob = ITYPOB(Ncodf)
         jbcc   = 2*(Jtypob - 1)
         if(jbcc.eq.8) jbcc = 0
         Ncodg = Ncodf
         if(Ncodg.gt.20) Ncodg = Ncodg - 20
         if(Ncodf.gt.0) then
c
c read first record for overriding error weighting of obs.ser.
            match     = 0
            natch     = 0
            Irid1     = 1
            Irid2     = 1
            Erwgt2(1) = Ugta(1)
            Erwgt2(2) = Ugta(2)
            if(Niobc.lt.0) goto 1800
            if(Niobc.ne.0) goto 900
            goto 800
         endif
      endif
  600 Ncodf  = 0
      Next   = 0
      Nobcon = 1
      rewind Iabs1
      Itrwnd(Iabs1) = 0
      return
  700 do while(.true.)
         read(Iobcon) jdc(1)
         if(jdc(1).le.0) goto 800
      end do
  800 read(Iobcon) ncodf1,nplnt1,site1,ser1,site2,spot1,erwgt1,
     .             acctm1,itime1,fdev1,freq1,nrewn1,ntape1,nseq1,
     .             nqlnt2,spot2,frq2,ctlg1,Gncode,Gnplt1,Gsite1,
     .             Gseres,Gsite2,Gspot,Gerwgt,Gacctm,Gitime,
     .             Gfdev,Gfreq,Gnrwnd,Gnplt2,Gspot2,Gfreq2,Gctlg
      Itrwnd(Iobcon) = 1
      Niobc = 1
      if(ncodf1.le.0 .and. nseq1.eq.0) then
         Niobc = -1
         rewind Iobcon
         Itrwnd(Iobcon) = 0
         goto 1800
      endif
  900 if(ntape1.lt.Ntape) goto 700
      if(ntape1.ne.Ntape) goto 1800
      if(nseq1.lt.Nseq) goto 700
      if(nseq1.ne.Nseq) goto 1800
      if(ntape1.le.0 .and. Ict(7).le.0) goto 700
      match = 1
      Irid1 = 2
      if(.not.Gncode) ncodf1 = Ncodf
      if(.not.Gnplt1) nplnt1 = Nplnt0
      if(.not.Gsite1) site1  = Sitf(1,1)
      if(.not.Gseres) ser1   = Series
      if(.not.Gspot) spot1   = Spota
      if(.not.Gnplt2) nqlnt2 = Nplnt2
      if(.not.Gspot2) spot2  = Spota2
      if(.not.Gsite2) site2  = Sitf(1,2)
      if(.not.Gerwgt(1)) erwgt1(1) = Ugta(1)
      if(.not.Gerwgt(2)) erwgt1(2) = Ugta(2)
      if(.not.Gctlg) ctlg1 = Ctlgm
 
      Erwgt2(1) = erwgt1(1)
      Erwgt2(2) = erwgt1(2)
c
c determine if observation series is to be skipped
      if(erwgt1(1).gt.cut(nrmprd) .or. erwgt1(2).gt.cut(nrmprd)) then
c
c check for consistency between iobcon and iabs1
         if(ncodf1.ne.Ncodf .or. nplnt1.ne.Nplnt0 .or.
     .   site1.ne.Sitf(1,1) .or. ser1.ne.Series .or. spot1.ne.Spota .or.
     .   Ctlgm.ne.ctlg1 .or. spot2.ne.Spota2 .or. nqlnt2.ne.Nplnt2)
     .      goto 1600
         if(Ncodf.gt.3 .or. Ncodf.le.9) goto 1800
         if(site2.eq.Sitf(1,2)) goto 1800
         goto 1600
      else
         call PAGCHK(60,2,2-nrmprd)
         write(Iout,950) Ntape,Nseq,Typobs(Ncodg),Ncodf,Plnnam,
     .                    Nplnt0,(Sitf(i,1),i=1,2),Series,
     .                    (Sitf(i,2),i=1,2),Spota
  950    format(/i4,i5,a6,i2,1x,2A4,i3,1x,2A4,1x,1A4,1x,
     .          2A4,1x,1A4,' **** SERIES SKIPPED ****')
      endif
 1000 do while(.true.)
         read(Iabs1,err=1500,end=1100) Ncode
         if(Ncode.le.0) goto 1300
      end do
 
c end of file on iabs1
 1100 ijend = 1
      call PAGCHK(60,2,0)
      write(Iout,1200) Ntape,Nseq,Iabs1
 1200 format('0NTAPE=',i4,' NSEQ=',i5,
     .       ' END OF FILE ACCEPTED IN MIDST OF SERIES ON IABS1=',i3)
 
c skip unused records on iobcon
 1300 if(Niobc.gt.0 .and. match.gt.0) then
         Niobc = 0
         match = 0
         if(natch.ge.0) then
            do while(.true.)
               read(Iobcon) jdc(1)
               if(jdc(1).le.0) then
                  natch = -1
                  goto 1400
               endif
            end do
         endif
      endif
 1400 if(itype.ne.4) goto 500
c
c end of series, make sure counters are minimum
      Ncode  = 0
      Numpar = 1
      numobs = 1
      Numsav = 1
      Ncal   = 1
      Mcobs2 = 1
      return
 
c error record, keep skipping
 1500 read(Iabs1)
      call PAGCHK(60,1,2 - nrmprd)
      write(Iout,2100) Iabs1
      goto 1000
 1600 call PAGCHK(60,2,2-nrmprd)
      write(Iout,1700) Iobcon,Iabs1
 1700 format('0*** DATA ON IOBCON=',i3,
     .       ' DOES NOT AGREE WITH THAT ON IABS1=',i3,
     .       ' WARNING IN PRDOBS LABEL=3180')
c
c set up pointers in bernum common
 1800 call PRDSET(nrmprd)
c
c decide if overriding series is forward or backward in time
      if(match.gt.0) then
         nbkfor = 1
         if(Jct(68).lt.0) then
         else if(Jct(68).eq.0) then
            if(nplnt1.gt.30 .or. (ncodf1.gt.1.and.ncodf1.lt.4) .or.
     .      ncodf1.gt.9 .or. Nspot.gt.0 .or. Nspot2.gt.0 .or.
     .      (Ksite(1).gt.0.and.Ksite(1).ne.3)) nbkfor = 2
         else
            nbkfor = 2
         endif
      endif
      return
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c*  start=4000
c          read observation record
 1900 read(Iabs1,err=2000,end=1100) Ncode,Ihr,Imin,Sec,
     .     (Result(j),Error(j),j = 1,2),Atuts,Ututs,Clamp,Limb,
     .     Observ,Imonth,Iday,Iyear,Jds,Jd,Ctat,Ctrecf,Numsav,
     .     (Save(i),i = 1,Numsav),numobs,Numpar,
     .     ((Deriv(i,j),i=1,Numpar),j = 1,numobs),imdt,Nmp2,
     .     ((Dvx(i,j),i=1,Nmp2),j = 1,numobs),Ncal,
     .     (Cal(i),Scal(i),Ical(i),i = 1,Ncal),Sumcor,Mnshob,
     .     Mixshp,(Mshobs(i),i = 1,Mnshob)
      Num1 = 1
      if(Ncode.le.0) goto 1300
      call DTCKI(imdt)
      Fdysc0 = Sec + 6E1_10*(Imin + 60*Ihr)
 
c force 0 <= seconds < 86400 from iabs1
      jdedit = Jds
      fdedit = Fdysc0
      if(fdedit.lt.0._10) then
         fdedit = fdedit + 86400._10
         jdedit = jdedit - 1
      else if(fdedit.ge.86400._10) then
         fdedit = fdedit - 86400._10
         jdedit = jdedit + 1
      endif
      if(Ncode.gt.3) Ncode = Ncode - 3
      Nice = Ncode - 2
      Jacc = jbcc
      if(Nice.gt.0) Jacc = Jacc + 1
 
c check for global deletions
      if(Tdlt0.gt.0._10 .and. Tdltp.gt.0._10) then
        tdlt=MOD(Jds + Fdysc0/86400._10 - 0.5_10 - Tdlt0,Tdltp)
        if(tdlt.LT.0._10) tdlt=tdlt+Tdltp
c skip if not in proper range
        if(tdlt.GT.Tdlton) goto 1900
        endif
 
c read record of iobcon series
      if(Niobc.lt.1) return
      if(match.le.0) return
      if(natch.lt.0) goto 2600
      if(natch.ne.0) goto 2300
      goto 2200
 2000 read(Iabs1)
      call PAGCHK(60,1,2 - nrmprd)
      write(Iout,2100) Iabs1
 2100 format(' **** ERROR RECORD SKIPPED ON IABS1=',i3)
      goto 1900
 2200 read(Iobcon) (jdc(i),ihrc(i),iminc(i),secc(i),i = 1,2),
     .             erobsc
      natch = 1
      do i = 1,2
         fdysc1(i) = secc(i) + 6E1_10*(iminc(i) + 60*ihrc(i))
      end do
      if(jdc(1).le.0) then
         natch = -1
         goto 2600
      endif
 2300 if(nbkfor.eq.2) then
c
c check for dates in forward in time series
         if(jdc(1).ge.jdc(2)) then
            if(jdc(1).gt.jdc(2) .or. fdysc1(1).gt.fdysc1(2)) then
               call SUICID('DATES ON IOBCON SERIES NOT FORWARD '//
     .                     'IN TIME, STOP IN PRDOBS  ',15)
            endif
         endif
         if(jdedit.lt.jdc(1)) goto 2600
         if(jdedit.eq.jdc(1)) then
            if(fdedit.lt.fdysc1(1)) goto 2600
         endif
         if(jdedit.lt.jdc(2)) then
         else if(jdedit.eq.jdc(2)) then
            if(fdedit.gt.fdysc1(2)) goto 2200
         else
            goto 2200
         endif
         goto 2500
      else
c
c check for dates in backward in time series
         if(jdc(1).le.jdc(2)) then
            if(jdc(1).lt.jdc(2) .or. fdysc1(1).lt.fdysc1(2))
     .         call SUICID('DATES ON IOBCON SERIES NOT BACKWARD IN '//
     .                  'TIME, STOP IN PRDOBS ',15)
         endif
      endif
 
      if(jdedit.lt.jdc(1)) then
      else if(jdedit.eq.jdc(1)) then
         if(fdedit.gt.fdysc1(1)) goto 2600
      else
         goto 2600
      endif
      if(jdedit.lt.jdc(2)) goto 2200
      if(jdedit.eq.jdc(2)) then
         if(fdedit.lt.fdysc1(2)) goto 2200
      endif
c
c alter error weighting for interval within series
 2500 Irid2 = 2
      if(Nice.le.0) Deriv(1,1)    = Error(1)*erwgt1(1)*erobsc(1)
      if(Nice.ge.0) Deriv(1,Num2) = Error(2)*erwgt1(2)*erobsc(2)
      return
c
c alter error weighting for whole series
 2600 if(Nice.le.0) Deriv(1,1)    = Error(1)*erwgt1(1)
      if(Nice.ge.0) Deriv(1,Num2) = Error(2)*erwgt1(2)
c
c
c*  start=9000
      return
      end
