      subroutine FRMHED(imats,qmat,subscr,key,ippr,iprnt)

      implicit none
c
c m.e.ash   feb 1970    subroutine frmhed
c read information records on saved normal equations tape
c modified for filter oct 1970  r.reasonberg
c expanded print paul macneil may, 1978
c
c arguments
      integer*4 imats,key,ippr,iprnt
      character*5 qmat
      character*4 subscr
c imats = tape that is read
c qmat,subscr = just used for a nice printout title
c key  = 0 skip records with control integers (assume they are same
c          as on a previous tape that was read)
c key.gt.0 read records with control integers
c key  = 1 printout heading information for input saved normal eqs
c key  = 2 skip restoring planet ic's from sne
c ippr is returned to calling subroutine
c ippr = 0 nrm eqns are not partially pre-reduced
c ippr = 1 nrm eqns are in partially pre-reduced format
c iprnt - code to govern printout: if 0 suppress print

c array dimensions
      include 'globdefs.inc'

c common
      include 'dtparm.inc'
      include 'empcnd.inc'
      include 'eqenxm.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'mcnfrm.inc'
      include 'mhrfrm.inc'
      include 'namtim.inc'
      include 'phasem.inc'
      include 'psrstm.inc'
      include 'rbiasm.inc'
      include 'skystm.inc'
      include 'stcrdm.inc'
      include 'sptcdm.inc'
      include 'wrkcomrs.inc'
      real*4    transf(u_stdsz,4,1000/u_stdsz)
      equivalence (transf(1,1,1),Pzham(1,1))

c local variables
      real*10 d8
      integer*4 i,i4,izero,j,jddtm0,jterat,
     . maxj,mgd2,mgdpts,mng4,mngd2,mpage,
     . msh,mtess2,mumdt2,mump4,mzone2,n,ncnmx,nmast,nprmx,nrec,nsh
      real*4 r4
      character*4 lnklvm
      character*4 grdmsg(9)/'PLAN','ET S','HAPE',' GRI','D MO','DEL,',
     . ' SUB','SET ','xx  '/
      real*10 jd2000/2451545._10/
c
c data on information records of saved normal equations tape
c which are not in common
      integer*2 mpspc(4),mcspc(4+u_mxpl),mcntr(u_mxpl)
      equivalence (mcspc(5),mcntr)
      integer*2 mumdt1
      integer*2 mzone1,mtess1,mplh3,mplh10

c local quantities for header skip
      character*88 ititl
      equivalence (ititl,mcspc)
      integer*2 i2,mudpln,mudphr,mudsit,mudspt,mudstr,mudpsr,
     .          mudrbs,mudeqn,mudphs,mdphr1,msitcr,mnmpsr,msptcr
      data izero/0/
      character*8 prmhdr(4)/'SOLAR SY','STEM PAR','AMETERS ',' '/,
     .          qprm/'PRMTR1 ='/, qmct/'MCT ='/, qdt1/'DT1 ='/,
     .          qjddt9/'JDDT9 ='/, qmdt/'MDT ='/,
     .    phdr(7)/' ',' INITIAL',' CONDITI','ONS AND ','PARAMETE',
     .          'RS PCOND', '1 ='/,
     .          hhdr(6)/' ',' GRAVITA','TIONAL P','OTENTIAL',
     .          ' HARMONI', 'CS MZH ='/,
     .    shdr(5)/' ',' STAR CA','TALOG ER','ROR COEF','FICIENTS'/
      character*1 hhdrc(48),cee/'C'/,em/'M'/,ess/'S'/,zee/'Z'/
      equivalence (hhdr,hhdrc)
      character*1 blank/' '/,astrik/'*'/
c
c read title record
      if(imats.le.0 .or. imats.gt.Ntptop) then
         write(Iout,50) imats
   50    format(' ERROR DETECTED IN FRMHED IMATS = ',i12)
         call SUICID(' IMATS OUT OF RANGE ',5)
      endif
      if(key.eq.0) then
c*  start=6000
c
c enter here to do nothing but skip header records
         Itrwnd(imats) = 1
         nrec = 1
         read(imats,end=180) ititl
         nrec = nrec + 1
         read(imats,end=180) nprmx,ncnmx,(d8,i=1,nprmx),
     .    (i2,i=1,nprmx+80),i4,i4,i4,
     .    mudpln,mudphr,mudsit,mudspt,mudstr,mudpsr,
     .    mudrbs,mudeqn,mudphs,i2,mumdt1,
     .    (r4,i4,i2,i=1,mumdt1),ippr,mdphr1,i2,(i2,i=1,mdphr1),
     .    i4,lnklvm,msitcr,mnmpsr,msptcr
         if(mnmpsr.eq.0) mnmpsr=16

c bypass head process
         maxj = mudpln + mudphr + 6
         if(mudsit.gt.0) maxj = maxj + 1
         if(mudspt.gt.0) maxj = maxj + 1
         if(mudstr.gt.0) maxj = maxj + 1
         if(mudpsr.gt.0) maxj = maxj + 1
         if(mudrbs.gt.0) maxj = maxj + 1
         if(mudeqn.gt.0) maxj = maxj + 1
         if(mudphs.gt.0) maxj = maxj + 1
         do j = 1,maxj
            nrec = nrec + 1
            read(imats,err=60,end=180)
            goto 100
   60       read(imats)
            call PAGCHK(60,1,0)
            write(Iout,80)
   80       format(' **** ERROR RECORD SKIPPED IN BYPASS HEAD PROCESS')
  100    end do
         if(iprnt.gt.0) then
            call PAGCHK(60,3,0)
            write(Iout,120) maxj,ititl
  120       format('0 HEAD PROCESS BYPASS OF ',i3,' RECORDS' /
     .             ' TITLE=',a88)
         endif
         return
      endif

      nrec = 1
      read(imats,err=200,end=180) Title
      goto 400
  180 write(iout,190) nrec,imats
  190 format('0UNEXPECTED END OF FILE AT RECORD',I5,' OF INPUT UNIT',
     . I3)
      call SUICID('EOF IN NRM. EQS. HEADERS, STOP IN FRMHED',10)

  200 read(imats) Title
      call PAGCHK(60,2,0)
      write(Iout,300) imats
  300 format(
     .  '0**** ERROR ACCEPTED ON TITLE RECORD OF SAVED NORMAL EQS.TAPE',
     .  i3)
  400 Itrwnd(imats) = 1
c
c solar system parameters
c initialize mdt vector to zero to prevent problems in frmbdy
      call ZFILL(Mdt,2*600)
      call ZFILL(Mshap,2*6)
      nrec = nrec + 1
      read(imats,end=180) nprmx,ncnmx,(Prmtr1(i),i=1,nprmx),
     . (Mprm(i),i=1,nprmx),Mct,Mparam,jterat,mpage,
     . Mumpln,Mumphr,Mumsit,Mumspt,Mumstr,Mumpsr,Mumrbs,
     . Mumeqn,Mumphs,Mumdt,mumdt1,
     . (Dt1(i),Jddt9(MIN(i,201)),Mdt(i),i=1,mumdt1),ippr,
     . Mphr1,Mshp2x,(Mshap(i),i=1,Mphr1),jddtm0,lnklvm,msitcr,mnmpsr,
     . msptcr
      if(nprmx.gt.u_nmprm .or. ncnmx.gt.u_nmbod) call SUICID(
     . 'BAD NUMBER OF PARAMETERS, STOP IN FRMHED',10)
      if(Nprmx.lt.u_nmprm) then
         do i=nprmx+1,u_nmprm
            Mprm(i)=0
         end do
      endif
      if(Mumdt.eq.0) Mdt(1) = 0
      call DTCHCK(jddtm0,Mumdt,mumdt1,Dt1,Mdt)
      if(mnmpsr.eq.0) mnmpsr=16
c
c convert i*2 to i*4
      mumdt2 = mumdt1
c
c*  start=1000
c
c read earth,moon,planet parameters
      mump4 = Mumpln + 4
      do i4 = 1,mump4
         i = i4 - 4
         if(i.gt.0) then

c ordinary planet
            nrec = nrec + 1
            read(imats,end=180) Mplnt(i),mcntr(i),Jdpl9(i),
     .         (Pcond1(j,i),j=1,ncnmx),(Mpl(j,i),j=1,ncnmx),Apln1(i)
         else

c earth-moon
            nrec = nrec + 1
            read(imats,end=180) mpspc(i4),mcntr(i),Jdpl9(i),
     .         (Pcond1(j,i),j=1,ncnmx),(Mpl(j,i),j=1,ncnmx),Apln1(i)
         endif
         if(ncnmx.lt.u_nmbod) then
            do j=ncnmx+1,u_nmbod
               Mpl(j,i)=0
            end do
         endif
      end do
c
c earth gravitational potential harmonics
      nrec = nrec + 1
      read(imats,end=180) mplh3,Mezone,mzone1,(Ezham(i),i=1,mzone1),
     .            (Mezhar(i),i=1,mzone1),Metess,mtess1,
     .            (Echam(i),i=1,mtess1),(Mechar(i),i=1,mtess1),
     .            (Esham(i),i=1,mtess1),
     .            (Meshar(i),i=1,mtess1)
c
c moon gravitational potential harmonics
      nrec = nrec + 1
      read(imats,end=180) mplh10,Mmzone,mzone1,(Mzham(i),i=1,mzone1),
     .            (Mmzhar(i),i=1,mzone1),Mmtess,mtess1,
     .            (Mcham(i),i=1,mtess1),(Mmchar(i),i=1,mtess1),
     .            (Msham(i),i=1,mtess1),
     .            (Mmshar(i),i=1,mtess1)
c
c planet gravitational potential harmonics
      if(Mumphr.gt.0) then
         do j = 1,Mumphr
            nsh = mnshap(j) + 1
            if(nsh.eq.2) then

c fourier
               nrec = nrec + 1
               read(imats,end=180) Mplhar(j),(Pzham(j,i),i=1,122),
     .                     (Mpzhar(j,i),i=1,122)
               Mpzone(j) = 122
               Mptess(j) = 0
            else if(nsh.eq.3) then

c grid
               nrec = nrec + 1
               read(imats,end=180) Mplhar(j),mgdpts,mgd2,
     .                     (Scontm(j,i),i=1,9),
     .                     (Pzham(j,i),i=1,mgd2),
     .                     (Mpzhar(j,i),i=1,mgdpts)
               Mngd(j)   = mgdpts
               Mpzone(j) = mgdpts
               Mptess(j) = 0
            else

c spherical harmonics
               nrec = nrec + 1
               read(imats,end=180) Mplhar(j),Mpzone(j),mzone1,
     .                     (Pzham(j,i),i=1,mzone1),
     .                     (Mpzhar(j,i),i=1,mzone1),Mptess(j),
     .                     mtess1,(Pcham(j,i),i=1,mtess1),
     .                     (Mpchar(j,i),i=1,mtess1),
     .                     (Psham(j,i),i=1,mtess1),
     .                     (Mpshar(j,i),i=1,mtess1)
            endif

         end do
      endif
c*  start=1500
c
c sites
      if(Mumsit.gt.0) then
         nrec = nrec + 1
         if(msitcr.lt.6) then
            read(imats,end=180) ((Site1(j,i),j=1,2),Kscrd1(i),
     .       (Scord1(j,i),j=1,3),(Mscrd(j,i),j=1,3),i=1,Mumsit)
            do i=1,Mumsit
               T0st1(i)=jd2000
               do j=4,6
                  Scord1(j,i)=0._10
                  Mscrd(j,i)=0
               end do
            end do
         else
            read(imats,end=180) ((Site1(j,i),j=1,2),Kscrd1(i),T0st1(i),
     .       (Scord1(j,i),j=1,msitcr),(Mscrd(j,i),j=1,msitcr),
     .       i=1,Mumsit)
            do i=1,Mumsit
               do j=4,6
                  if(Scord1(j,i).ne.0._10 .or. Mscrd(j,i).gt.0) goto 425
               end do
            end do
c no velocities actually specified or to be adjusted
            msitcr=3
  425       continue
         endif
      endif
c
c spots
      if(Mumspt.gt.0) then
         nrec = nrec + 1
         if(msptcr.lt.6) then
            read(imats,end=180) (Spot1(i),Msplnt(i),
     .        (Spcrd1(j,i),j=1,3),(Mspcrd(j,i),j=1,3),i=1,Mumspt)
            do i=1,Mumspt
               T0sp1(i)=jd2000
               do j=4,6
                  Spcrd1(j,i)=0._10
                  Mspcrd(j,i)=0
               end do
            end do
         else
            read(imats,end=180) (Spot1(i),Msplnt(i),T0sp1(i),
     .       (Spcrd1(j,i),j=1,msptcr),(Mspcrd(j,i),j=1,msptcr),
     .       i=1,Mumspt)
            do i=1,Mumspt
               do j=4,6
                  if(Spcrd1(j,i).ne.0._10.or.Mspcrd(j,i).gt.0) goto 427
               end do
            end do
c no velocities actually specified or to be adjusted
            msptcr=3
  427       continue
         endif
      endif
c
c star catalog corrections
      if(Mumstr.gt.0) then
         nrec = nrec + 1
         read(imats,end=180) (Ctlgn1(i),n,Mskyc(i),
     .              (Skycf1(j,i),j=1,n),(Msky(j,i),j=1,n),i=1,Mumstr)
      endif
c
c pulsar parameters
      if(Mumpsr.gt.0) then
         nrec = nrec + 1
         read(imats,end=180) (Sptps1(i),Jdpsr9(i),Plspr1(i),
     .    Ntyps1(i),(Psrcn1(j,i),j=1,mnmpsr),(Mpsrcn(j,i),j=1,mnmpsr),
     .    i=1,Mumpsr)
         if(mnmpsr.lt.u_nmpsr) then
            do i=1,Mumpsr
               do j=mnmpsr+1,u_nmpsr
                  Psrcn1(j,i)=0._10
                  mpsrcn(j,i)=0
               end do
            end do
         endif
      endif
c
c radar biases
      if(Mumrbs.gt.0) then
         nrec = nrec + 1
         read(imats,end=180) ((Rdbst1(j,i),j=1,2),Rdbsr1(i),
     .              (Rbias1(j,i),j=1,2),(Mrbs(j,i),j=1,2),Mplrbs(i),
     .              i=1,Mumrbs)
      endif
c
c equinox equator
      if(Mumeqn.gt.0) then
         nrec = nrec + 1
         read(imats,end=180) (Eqnst1(i),Eqnsr1(i),
     .             (Deqnx1(j,i),j=1,3),(Meqn(j,i),j=1,3),i=1,Mumeqn)
      endif
c
c phase corrections
      if(Mumphs.gt.0) then
         nrec = nrec + 1
         read(imats,end=180) (Phsit1(i),Phser1(i),
     .                (Aphs1(j,i),j=1,9),Mcphs(i),Mplphs(i),
     .                (Mphs(j,i),j=1,9),i=1,Mumphs)
      endif
c
c*  start=1800
c write out information records on saved normal equations tape
      if(iprnt.gt.0) then
         call PAGCHK(45,(nprmx-2)/33+4,0)
         write(Iout,450) qmat,subscr,imats,Title,lnklvm,jterat,
     .                    mpage,Mparam
  450    format('0SAVED NORMAL EQS.TAPE ',a5,a4,'=',i3,' TITLE=',
     .          a88/' LEVEL=',a4,'  JTERAT=',i3,'  MPAGE=',i5,
     .'    PARTIAL DERIVATIVE CONTROLS FOR THE SAVED NORMAL EQUATIONS OF
     . ORDER',i5,' ARE')
         write(Iout,500) (Mprm(i),i=1,nprmx)
  500    format(' SOLAR SYSTEM PARAMETERS MPRM=',34I3/(33x,33I3))

         if(Jct(57).gt.0) then
c
c*  start=2000
c
            prmhdr(4) = qprm
            call HDPNR8(Prmtr1,nprmx,1,prmhdr,8)

            prmhdr(4) = qmct
            call HDPNI2(Mct,80,1,prmhdr,8)

            if(Mumdt.gt.0) then
               prmhdr(4) = qdt1
               call HDPNR4(Dt1,mumdt2,1,prmhdr,8)
               prmhdr(4) = qjddt9
               call HDPNI4(Jddt9,Mumdt+0,1,prmhdr,8)
               prmhdr(4) = qmdt
               call HDPNI2(Mdt,mumdt2,1,prmhdr,8)
            endif
         else
            write(Iout,520) Jct(57)
  520       format(' EXPANDED RNE HEAD PRINT CANCELLED BY JCT(57) = ',
     .             i6)
            Line = Line + 1
         endif
c
c
         do i4 = 1,mump4
            i = i4 - 4
            phdr(1) = Apln1(i)
            if(Jct(57).gt.0) call HDPNR8(Pcond1(1,i),ncnmx,1,phdr,13)
            call PAGCHK(60,(ncnmx-7)/24+1,0)
            if(i.gt.0) then

c ordinary planet
               write(Iout,530) i,Apln1(i),Mplnt(i),mcntr(i),
     .                          Jdpl9(i),(Mpl(j,i),j=1,ncnmx)
            else

c earth-moon
               write(Iout,530) izero,Apln1(i),mpspc(i4),mcntr(i),
     .                          Jdpl9(i),(Mpl(j,i),j=1,ncnmx)
  530          format(i4,'. ',a8,' NPLNT=',i3,' NCNTR=',i2,
     .                ' JD0=',i7,' L=',6I2,(24I3))
            endif
         end do
         call PAGCHK(60,3,0)
         write(Iout,550) Mumphr,Mumsit,Mumspt,Mumstr,Mumpsr,
     .                    Mumrbs,Mumeqn,Mumphs,Mumdt,ippr,Mshap
  550    format('0  MUMPHR=',i4,'   MUMSIT=',i4,'   MUMSPT=',i4,
     .          '   MUMSTR=',i4,'   MUMPSR=',i4,'   MUMRBS=',i4,
     .          '   MUMEQN=',i4,'   MUMPHS=',i4,'   MUMDT =',i4/
     .          '     IPPR=',i4,'    NSHAP=',6I2)
         if(Jct(57).gt.0) then
c
c convert i*2 to i*4
            mzone2 = Mezone - 1
            mtess2 = Metess*(Metess + 1)/2 - 1
            if(mzone2.le.0) mzone2 = 1
            call PAGCHK(60,2,0)
            write(Iout,560) mplh3,Mezone,Metess
  560       format('0','MPLH3 =',i3,5x,'MEZONE =',i3,5x,
     .             'METESS =',i3)
            if(mzone2.gt.0 .or. mtess2.gt.0) then
               call HDPNR8(Ezham,mzone2,1,
     .                'EARTH GRAVITATIONAL POTENTIAL HARMONICS EZHAM = '
     .                ,12)
               call HDPNI2(Mezhar,mzone2,1,
     .                'EARTH GRAVITATIONAL POTENTIAL HARMONICS MEZHAR ='
     .                ,12)
               call HDPNR8(Echam,mtess2,1,
     .                'EARTH GRAVITATIONAL POTENTIAL HARMONICS ECHAM = '
     .                ,12)
               call HDPNI2(Mechar,mtess2,1,
     .                'EARTH GRAVITATIONAL POTENTIAL HARMONICS MECHAR ='
     .                ,12)
               call HDPNR8(Esham,mtess2,1,
     .                'EARTH GRAVITATIONAL POTENTIAL HARMONICS ESHAM = '
     .                ,12)
               call HDPNI2(Meshar,mtess2,1,
     .                'EARTH GRAVITATIONAL POTENTIAL HARMONICS MESHAR ='
     .                ,12)
            endif
c
c convert i*2 to i*4
            mzone2 = Mmzone - 1
            if(mzone2.le.0) mzone2 = 1
            mtess2 = Mmtess*(Mmtess + 1)/2 - 1
            call PAGCHK(60,2,0)
            write(Iout,580) mplh10,Mmzone,Mmtess
  580       format('0MPLH10 =',i3,5x,'MMZONE =',i3,5x,
     .             'MMTESS =',i3)
            if(mzone2.gt.0 .or. mtess2.gt.0) then
               call HDPNR8(Mzham,mzone2,1,
     .                'MOON GRAVITATIONAL POTENTIAL HARMONICS MZHAM =  '
     .                ,12)
               call HDPNI2(Mmzhar,mzone2,1,
     .                'MOON GRAVITATIONAL POTENTIAL HARMONICS MMZHAR = '
     .                ,12)
               call HDPNR8(Mcham,mtess2,1,
     .                'MOON GRAVITATIONAL POTENTIAL HARMONICS MCHAM =  '
     .                ,12)
               call HDPNI2(Mmchar,mtess2,1,
     .                'MOON GRAVITATIONAL POTENTIAL HARMONICS MMCHAR = '
     .                ,12)
               call HDPNR8(Msham,mtess2,1,
     .                'MOON GRAVITATIONAL POTENTIAL HARMONICS MSHAM =  '
     .                ,12)
               call HDPNI2(Mmshar,mtess2,1,
     .                'MOON GRAVITATIONAL POTENTIAL HARMONICS MMSHAR = '
     .                ,12)
            endif

            do j = 1,Mumphr
               call PAGCHK(60,3,0)
               write(Iout,590) Mplhar(j),Mpzone(j),Mptess(j)
  590          format('-NPLHAR=',i3,5x,'NZONE=',i3,5x,'NTESS=',i3)
               msh = mnshap(j) + 1
               if(msh.eq.2) then

                  call HDPNR8(Pzham(j,1),122,4,
     .                'PLANET SHAPE FOURIER COEFICIENTS..OVLYFS=   ',11)
                  call HDPNI2(Mpzhar(j,1),122,4,
     .                'PLANET SHAPE FOURIER COEFICIENTS..MSHP  =   ',11)
               else if(msh.eq.3) then
                  mng4  = Mngd(j)
                  mngd2 = (mng4+u_stdsz-1)/u_stdsz
                  do i=1,u_stdsz
                     call EBCDI(i,grdmsg(9),2)
                     call HDPNR4(transf(i,j,1),mngd2,4*u_stdsz,grdmsg,9)
                  end do
                  call HDPNI2(Mpzhar(j,1),mng4,4,
     .                        'PLANET SHAPE GRID PARAMETERS..MSHP= ',9)
               else

                  mzone2 = Mpzone(j) - 1
                  mtess2 = Mptess(j)*(Mptess(j) + 1)/2 - 1
                  do i = 1,Mumpln
                     if(Mplhar(j).eq.Mplnt(i)) then
                        hhdr(1) = Apln1(i)
                        goto 600
                     endif
                  end do
  600             if(mzone2.gt.0 .or. mtess2.gt.0) then
                     if(mzone2.gt.0) then
                        hhdrc(44) = blank
                        hhdrc(45) = zee
                        call HDPNR8(Pzham(j,1),mzone2,4,hhdr,12)
                        hhdrc(44) = em
                        call HDPNI2(Mpzhar(j,1),mzone2,4,hhdr,12)
                     endif
                     if(mtess2.gt.0) then
                        hhdrc(44) = blank
                        hhdrc(45) = cee
                        call HDPNR8(Pcham(j,1),mtess2,4,hhdr,12)
                        hhdrc(44) = em
                        call HDPNI2(Mpchar(j,1),mtess2,4,hhdr,12)
                        hhdrc(44) = blank
                        hhdrc(45) = ess
                        call HDPNR8(Psham(j,1),mtess2,4,hhdr,12)
                        hhdrc(44) = em
                        call HDPNI2(Mpshar(j,1),mtess2,4,hhdr,12)
                     endif
                  endif
               endif

            end do
c
c*  start=3300
            do i = 1,Mumsit
               if(Line.gt.59 .or. i.le.1) then
                  call PAGCHK(57,4,0)
                  if(msitcr.lt.6) then
                     write(Iout,610)
  610                format('0',10x,'SITE COORDINATES'/
     .                '0SITE      KSCRD     MSCRD',6x,'SITCRD')
                  else
                     write(Iout,615)
  615                format('0',10x,'SITE COORDINATES'/
     .                '0SITE      KSCRD     MSCRD',7x,'T0',10x,'SITCRD')
                  endif
               endif
               if(msitcr.lt.6) then
                  write(Iout,620) sitd1(i),Kscrd1(i),
     .             (Mscrd(j,i),j=1,3),(Scord1(j,i),j=1,3)
  620             format(1x,a8,i6,3x,3I3,1x,1p,3D25.15)
               else
                  write(Iout,625) sitd1(i),Kscrd1(i),
     .             (Mscrd(j,i),j=1,6),T0st1(i),(Scord1(j,i),j=1,6)
  625             format(1x,a8,i6,3x,6i2,f11.1,3f16.9,3f12.6)
               endif
               Line = Line + 1
            end do
c
c*  start=3400
            do i = 1,Mumspt
               if(Line.gt.59 .or. i.le.1) then
                  call PAGCHK(57,4,0)
                  if(msptcr.lt.6) then
                     write(Iout,630)
  630                format('0',10x,'SPOT COORDINATES'/
     .                   '0SPOT     NSPLNT    MSPCRD      SPCORD')
                  else
                     write(Iout,635)
  635                format('0',10x,'SPOT COORDINATES'/
     .                '0SPOT     NSPLNT    MSPCRD',7x,'T0',10x,'SPCORD')
                  endif
               endif
               if(msptcr.lt.6) then
                  write(Iout,640) Spot1(i),Msplnt(i),
     .             (Mspcrd(j,i),j=1,3),(Spcrd1(j,i),j=1,3)
  640             format(1x,a4,4x,i6,3x,3I3,1x,1p,3D25.15)
               else
                  write(Iout,645) Spot1(i),Msplnt(i),
     .             (Mspcrd(j,i),j=1,6),T0sp1(i),(Spcrd1(j,i),j=1,6)
  645             format(1x,a4,4x,i6,3x,6i2,f11.1,3f16.9,3f12.6)
               endif
               Line = Line + 1
            end do
c
c*  start=3500
            do i = 1,Mumstr
               n = Mskyc(i)
               if(n.le.0) n = 1
               shdr(1) = Ctlgn1(i)
               call HDPNR8(Skycf1(1,i),n,1,shdr,10)
               call HDPNI2(Msky(1,i),n,1,shdr,2)
            end do
c
c no interferometer clock offsets in header records
c
            do i = 1,Mumrbs
               if(Line.gt.59 .or. i.le.1) then
                  call PAGCHK(57,4,0)
                  write(Iout,650)
  650             format('0',10x,
     .                   'RADAR BIASES'/'0  SITES  SERIES    RBIAS',
     .                   35x,'MRBS  NPLRBS')
               endif
               write(Iout,660) Rdbst1(1,i),Rdbst1(2,i),
     .                          Rdbsr1(i),Rbias1(1,i),
     .                          Rbias1(2,i),Mrbs(1,i),Mrbs(2,i),
     .                          Mplrbs(i)
  660          format(3(1x,a4),1p,2E20.8,3x,2I3,i6)
               Line = Line + 1
            end do
c
c*  start=3600
            do i = 1,Mumeqn
               if(Line.gt.59 .or. i.le.1) then
                  call PAGCHK(57,4,0)
                  write(Iout,670)
  670             format('0',10x,'EQUINOX-EQUATOR CORRECTIONS'/
     .                   '0SITE SERIES     CORRECTIONS',46x,
     .                   'MEQN')
               endif
               write(Iout,680) Eqnst1(i),Eqnsr1(i),(Deqnx1(j,i),j=1,3),
     .                          (Meqn(j,i),j=1,3)
  680          format(1x,a4,2x,a4,1p,3E20.8,1x,3I3)
               Line = Line + 1
            end do
c
c*  start=3700
            do i = 1,Mumphs
               call PAGCHK(60,2,0)
               write(Iout,690) Phsit1(i),Phser1(i),Mcphs(i),Mplphs(i)
  690          format('0','PHSIT1 = ',a4,5x,'PHSER1 = ',a4,5x,
     .                'NCPHS =',i5,5x,'NPLPHS =',i5)
               call HDPNR4(Aphs1(1,i),9,1,
     .                     'PHASE CORRECTIONS APHS1 =   ',7)
               call HDPNI2(Mphs(1,i),9,1,
     .                     'PHASE CORRECTIONS MPHS =',6)
            end do

c*  start=3800
            do i = 1,Mumpsr
               call PAGCHK(60,2,0)
               write(Iout,700) Sptps1(i),Jdpsr9(i),Plspr1(i),
     .                          Ntyps1(i)
  700          format('0PULSAR ',a4,'  JD0=',i8,
     .                '  BASE PERIOD=',1pd22.15,'  NTYPE=',i3)
               call HDPNR8(Psrcn1(1,i),16,1,
     .                     'PULSAR PARAMETERS   ',5)
               call HDPNI2(Mpsrcn(1,i),16,1,'PULSAR MPSRX',3)
            end do
         endif
      endif
c
c*  start=4000
c
c test for parameter copying
      if(key.le.1) then
c
c update initial conditions
         if(Ict(80).gt.0 .and. Iterat.le.1) then
            call LIBCND(Jdem9,Econd1,Lem,Jdem0,Econd)
            call LIBCND(Jdmn9,Mcond1,Lmn,Jdmn0,Mcond)
            call LIBCND(Jder9,Ercnd1,Ler,Jder0,Ercond)
            call LIBCND(Jdmr9,Mrcnd1,Lmr,Jdmr0,Mrcond)
            do j = 1,Mumpln
               if(Numpln.le.0) goto 800
               do i = 1,Numpln
                  if(Nplnt(i).eq.Mplnt(j)) then
                     call LIBCND(Jdpl9(j),Pcond1(1,j),Lpl(1,i),
     .                           Jdpl0(i),Pcond(1,i))
                     goto 710
                  endif
               end do
  710       end do
         endif
      endif
c
c*  start=5000
c check consistency of et-ut2 tables
  800 if(Numdt.gt.0 .and. Mumdt.gt.0) then
         if(Numdt.ne.Mumdt) call SUICID(
     .' LENGTH OF ET-UT2 TABLE ON SAVED NORMAL EQS TAPE DOES NOT AGREE W
     .ITH INPUT, STOP IN FRMHED  ',23)
         nmast = 0
         do i = 1,Numdt
            Ast(i) = blank
            if(Jddt(i).ne.Jddt9(i)) then
               Ast(i) = astrik
               nmast  = nmast + 1
            endif
         end do
         if(nmast.gt.0) then
            write(Iout,820) nmast,
     .                     (i,Jddt(i),Jddt9(i),Ast(i),i=1,Numdt)
  820       format('0 N     JDDT   JDDT9 BAD=',i3/(i4,'.',2I8,1x,a1))

            call SUICID(' JD IN ET-UT2 TABLE ON SAVED NORMAL EQS '//
     .        'TAPE DOES NOT AGREE WITH INPUT, STOP IN FRMHED  ',22)
         endif
      endif
c
c*  start=9000
      return
      end
