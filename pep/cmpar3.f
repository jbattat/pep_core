      subroutine CMPAR3(mocpar, jtape)
 
      implicit none
 
 
c amuchastegui/ash -  june 1970 - subroutine cmpar3

c arguments
      integer*4 mocpar,jtape
c
c array dimensions
      include 'globdefs.inc'

c commons
      include 'bdctrl.inc'
      include 'comdateq.inc'
      include 'coord.inc'
      include 'empcnd.inc'
      include 'eqnphs.inc'
      real*4 eqnx(3)
      equivalence (eqnx, Pnox)
      include 'fcntrl.inc'
      include 'inodta.inc'
      integer*4 i2bod
      equivalence (i2bod,Jpert)
      include 'lcntrl.inc'
      include 'mnsprt.inc'
      real*10 spcdxx(3),spcdx2(3)
      equivalence (Spcdx(1,1),spcdxx(1)),(Spcdx(1,2),spcdx2(1))
      include 'ltrapx.inc'
      integer*2 lspcdy(3), lspcy2(3)
      equivalence (Lspcdx(1,1),lspcdy(1)), (Lspcdx(1,2),lspcy2(1))
      include 'mtrapx.inc'
      include 'namtim.inc'
      integer*2 np3,np10
      equivalence (np3,Nplnt(-3)),(np10,Nplnt(-2))
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*4 acctim, accdst, accprc
      real*10 fdev, freq2
      equivalence (acctim, Estf), (fdev, Dstf(10)), (Dstf(5), freq2)
      equivalence (Estf(2), accdst), (Estf(3), accprc)
      include 'kobequiv.inc'
      include 'obstap.inc'
      include 'pemctl.inc'
      include 'plndta.inc'
      include 'rotdta.inc'
      include 'rtrdvl.inc'
      include 'sbdta.inc'
      include 'sitcrd.inc'
      include 'skymap.inc'
      include 'spqind.inc'
      include 'stats.inc'
      include 'tapdta.inc'
 
c temporary storage
      common/YVECT/ Ppcond(u_nmbod),Bbcond(u_nmbod),Bbcone(u_nmbod,2),
     . Jdpp0,Jdbb0,Jebb0(2),Exnams
      character*8 Exnams
      real*10 Ppcond,Bbcond,Bbcone
      integer*4 Jdpp0,Jdbb0,Jebb0
      integer*4 zyvect/5104/   !r8=2552,r10=5104
      include 'zeroes.inc'
c
c local
      integer*2 nczon1,nctes1,nszon1,nstes1,numtr1,ncobs1,ncobs2,
     .          nsky1,itimo
      character*8 erthnm/' EARTH  '/, star/'  STAR  '/
      character*4 al15/'AL15'/, ap15/'AP15'/
      integer*4 i,itg,j,jtg,jtypob,jut1,jwob,next,npd
      integer*2 nmpex1,nmpex2
      integer*2 n16/16/
c external functions
      integer*4 ITYPOB,JBDTST
 
c first guesses at distance
      real*4 dist0(10)/300., 150., 0., 300., 2200., 4250., 9090.,
     .          14530., 19220., 1.28/
c
c*  start=100
c           reinitialize tapes if rewound
c           note: unvarying tapes are not rewound unless nrewnd.gt.0
c
      if(Nrewnd.gt.0) then
c
c read first five records of earth rotation tape if it is in
c a rewound state
         if(Inut.gt.0) then
            call RTRD1(Linut, 3)
            Linut = 1
         endif
c
c read first five records of n-body tape if it is in a rewound
c state
         if(Nbody.gt.0) then
            call BDRD1(Libdy)
            Libdy = 1
         endif
c
c read first five records of s-body tape if it is rewound
         if(i2bod.gt.0) then
            call B2RD1(Lib2y)
            Lib2y = 1
         endif
c
c read first five records of moon tape if it is rewound
         if(Imn.gt.0) then
            if(Itrwnd(Imn).le.0 .and. JBDTST(np10).le.0) then
               call MNRD1(Limn)
               Limn = 1
            endif
         endif
c
c read first five records of moon rotation tape if it is in a
c rewound state and observed body is not probe with earth
c as central body
         if(Ilib.gt.0) then
            call RTRD1(Lilib, 10)
            Lilib = 1
         endif
c
c read first five records of earth-moon barycenter tape if it
c is rewound
         if(Iem.gt.0) then
            if(Itrwnd(Iem).le.0 .and. JBDTST(np3).le.0) then
               call EMRD1(Liem)
               Liem = 1
            endif
         endif
      endif
c
c*  start=200
c           read first five records of planet tape if it is in a rewound
c           state and observed body is not moon or sun, so that it is
c           in fact planet nplnt(klan) or nplnt(klanb) with nplnt(klan)
c           as central body
      if(Nplnt0.eq.0) goto 300
      if(Klan.eq.0 .or. Klan.gt.u_mxpl) goto 200
      if(Klan.eq.Klam .and. Nrewnd.le.0) then
         if(Jplnt.gt.0 .and. Itrwnd(Jplnt).gt.0) goto 100
         if(JBDTST(Nplnt(Klan)).ne.0) goto 100
      endif
      call PLRD1(Lipl(Klan))
      Lipl(Klan) = 1
c
c read first five records of planet rotation tape if it is in
c a rewound state and observed body is not moon or sun
  100 if(Klanr.gt.0) then
         if(Klan.eq.Klam .and. Nrewnd.le.0 .and.
     .      Jpr.gt.0 .and. Itrwnd(Jpr).gt.0) goto 200
         if(Iplnt(Klanr).gt.0) then
            call RTRD1(Lipl(Klanr), 0+Nplnt0)
            Lipl(Klanr) = 1
         else
            Jpr = 0
         endif
      endif
c
c read first five records of satellite-probe tape if it is in
c a rewound state
  200 if(Klanb.gt.0) then
         if(Klanb.eq.Klamb .and. Nrewnd.le.0) then
            if(Jsb.gt.0 .and. Itrwnd(Jsb).gt.0) goto 300
            if(JBDTST(Nplnt(Klanb)).ne.0) goto 300
         endif
         call SBRD1(Lipl(Klanb))
         Lipl(Klanb) = 1
      endif
c
c read first five records of observing site not on earth tape
  300 if(Klans1.gt.0 .and. Klans1.lt.u_mxpl+1) then
         if(Klams1.eq.Klans1 .and. Nrewnd.le.0) goto 400
         call SCRD1(Lipl(Klans1))
         Lipl(Klans1) = 1
      endif
c
c read first two records of ut1 and wobble data sets if
c they are in a rewound state
  400 if(Jct(33).gt.0) then
         jut1 = Jct(33)
         if(Itrwnd(jut1).le.0) call UT1RD1
         jwob = Jct(33) + 1
         if(Itrwnd(jwob).le.0) call WOBRD1
      endif
c
c read first two records of libration correction dataset if needed
c for this series
      if(Libhoc.gt.0 .and. (Nplnt0.eq.10.or.Nplnt2.eq.10))
     . call HOCRED(0,0._10)
c
c initialize fluid displacement database if desired
      if(Jct(49).gt.0) call FLURD1
c
c*  start=400
c initialize interpolators
      call TRPNIT
c
c determine planet and satellite-probe partial derivative
c controls
      call ZFILL(Ppcond, 16*120+24)
      if(Nplnt0.gt.0 .and. Nplnt0.ne.10) then
         if(Klan.le.u_mxpl .and. Klan.gt.0) then
            do j = 1, u_nmbod
               Ppcond(j) = Pcond(j, Klan)
            end do
            call LVTBDY(Lplx, Lpl(1,Klan), Mplx, u_nmbod)
            Jdpp0 = Jdpl0(Klan)
         endif
         if(Klanb.gt.0) then
            do j = 1, u_nmbod
               Bbcond(j) = Pcond(j, Klanb)
            end do
            call LVTBDY(Lsbx, Lpl(1,Klanb), Msbx, u_nmbod)
            Jdbb0 = Jdpl0(Klanb)
         endif
      endif
      ncobs1 = 1
      ncobs2 = 1
      if(Klans1.gt.0) then
         ncobs2 = u_nmbod
         call LVTBDY(Lscx(1,1), Lpl(1,Klans1), Mscx(1,1), u_nmbod)
         do j = 1, u_nmbod
            Bbcone(j, 1) = Pcond(j, Klans1)
         end do
         Jebb0(1) = Jdpl0(Klans1)
      endif
c
c determine if embarycenter initial condition partials were
c calculated on moon tape for observations of the moon
      if(Klan.eq.u_mxpl+1) then
         do i = 1, 6
            Lemmn(i) = 0
         end do
         if(Kimn(1).gt.0) then
            do i = 8, Nkimn
               if(Kimn(i).eq.0) goto 450
               itg = (Kimn(i) - 1)/100
               if(itg.eq.3) then
                  jtg = Kimn(i) - itg*100
                  if(jtg.le.6) then
                     if(Lem(jtg).gt.0) Lemmn(jtg) = 1
                  endif
               endif
            end do
         endif
  450    if(Iabs1.gt.0) then
            do i = 1, 6
               if(Memmn(i).gt.0) Lemmn(i) = 1
            end do
         endif
      endif
c
c*  start=600
c determine central body harmonic and target body partial
c derivative controls
      call TRGCNT
c
c determine observed body gravitational potential or shape
c harmonic controls
c determine shape harmonic controls
      call SHPCNT
c set up planet center of mass (or solar-system barycenter) corrections
c must be after call to plrd1 and before resetting klam
      Cmfct = 1._10
      call PLCMS
c
c set up controls for exo-planets
      call PEXLIC
c
c write first record of observed minus theory tape for
c observation series
      Klam   = Klan
      Klamb  = Klanb
      Klaps  = Klap
      Klams1 = Klans1
      Klams2 = Klans2
      Klamr  = Klanr
      if(Spotf.eq.ap15 .and. Ncodf.gt.3) Spotf   = al15
      if(Spotf2.eq.ap15 .and. Ncodf.gt.3) Spotf2 = al15
      if(Nplnt0.ne.3) then
         Pname = Aplnt(Klap)
         if(Nplnt0.lt.0) Pname = star
      else
         Pname = erthnm
      endif
      if(Iabs2.gt.0) then
         Nobcn1 = Nobcon
         if(Iabs1s.le.0 .or. Nobcon.le.0) then
            Nobcn1    = 2
            Obscon(1) = 1._10
            Obscon(2) = 1._10
         endif
         nczon1 = Nczone
         if(Nczone.le.0) nczon1 = 1
         nctes1 = Nctess
         if(Nctess.le.0) nctes1 = 1
         nszon1 = Nszone
         if(Nszone.le.0) nszon1 = 1
         nstes1 = Nstess
         if(Nstess.le.0) nstes1 = 1
         numtr1 = Numtar
         if(Numtar.le.0) numtr1 = 1
         next  = 1
         nsky1 = Nskyc
         if(Nskyc.le.0) nsky1 = 1
         if(Nmpex.gt.0) then
            nmpex1=Nmpex
            nmpex2=u_nmbod
         else
            nmpex1=1
            nmpex2=1
         endif
         if(ctime.lt.0) then
            itimo=ctime
         else
            itimo=ctime*10+itime
         endif
         if(Idumob.ne.1 .or. Ict(3).gt.0) write(Iabs2)
     .    Nseqa(Ntape),Ncoda(Ntape),Nplnt0,Sitf(1),Series,Sitf(2),Spotf,
     .    (Erwgta(i,Ntape),i=1,2),acctim,Fdeva(Ntape),Freq,itimo,Nrewnd,
     .    Npage,Pname,Ncp0,Jdpp0,Ppcond,Lplx,Jdbb0,Bbcond,Lsbx,Lscrdx,
     .    Nrbias,Rbsx,Lrbsx,Neqnox,eqnx,Leqnx,Ncph,Aphs,Lphsx,Nobcn1,
     .    (Obscon(i),i=1,Nobcn1),Coords,Ksite,spcdxx,lspcdy,Nczone,
     .    nczon1,(Lczhar(i),i=1,nczon1),Nctess,nctes1,(Lcchar(i),
     .    i=1,nctes1),(Lcshar(i),i=1,nctes1),
     .    Numtar,numtr1,(Ntrg(i),(Ltbod(j,i),j=1,u_nmbod),
     .    Ntzone(i),(Ltzhar(j,i),j=1,4),Nttess(i),
     .    (Ltchar(j,i),j=1,5),(Ltshar(j,i),j=1,5),i=1,numtr1),ncobs1,
     .    ncobs2,Aplnt(Klans1),(Lscx(i,1),i=1,ncobs2),Jebb0(1),
     .    (Bbcone(i,1),i=1,ncobs2),Nszone,nszon1,(Lszhar(i),i=1,nszon1),
     .    Nstess,nstes1,(Lschar(i),i=1,nstes1),(Lsshar(i),i=1,nstes1),
     .    Nplnt2,Spotf2,spcdx2,lspcy2,freq2,Lemmn,next,(Exnams,Exnams,
     .    i=1,next),Lngd,Lnfour,Ctlg,nsky1,(Lsky(i),i=1,nsky1),
     .    n16,Lpsrx,Nmpex,nmpex1,nmpex2,(Nplex(j),(Lpex(i,j),
     .    i=1,nmpex2),j=1,nmpex1),T0sit,
     .    izero,izero,izero,izero,izero,izero,izero,izero
 
      endif
c
c*  start=700
c initialize error analysis quantities
      do i = 1, 16
         Erquan(i) = 0._10
      end do
      do i = 17, 20
         Nit(i) = 0
      end do
      do i = 1, 2
         Nobs(i)  = 0
         Nast(i)  = 0
         nast1(i) = 1
         Ncard(i) = 0
      end do
      Numpar = 2
c
c printout page heading
      call NEWPG
      write(Iout, 500) Nplnt0, Ntape, Ncodf, Ntape, Ntapa(Ntape),
     .                 Ntape, Nseqa(Ntape), Jiabs1, jtape, Iiabs1,
     .                 Jiabs2, jtape, Iabs2, Iiobs, Iiobcn, Nrewnd
  500 format(/' NPLNT0=', i2, '  NCODF=NCODA(', i1, ')=', i2,
     .       '  NTAPA(', i1, ')=', i3, '  NSEQA(', i1, ')=', i5,
     .       '  IABS1=IOBS', i1, '(', i2, ')=', i2, '  IABS2=IOBS', i1,
     .       '(', i2, ')=', i2, '  IOBS=', i2, '  IOBCON=', i2,
     .       '  NREWND=', i2)
      Line = Line + 2
      if(Iabs1.gt.0) then
         write(Iout,550) Klan,Mplx,Nsite1,Nsite2,Mscrdx,Mrbsx,
     .                    Klanb,Msbx,Nspot,Nspot2,
     .                    Mspcdx,Meqnx,Klans1,Mphsx,Mskyc,Sera(Ntape)
  550    format('  KLAN =', i2, ' MPL=', 30I2, 2x, ' NSITE=', 2I3,
     .    ' MSCRD=', 12I2, ' MRBS=', 2I2/
     .    '  KLANB=', i2, ' MSB=', 30I2, 2x, ' NSPOT=', i3,
     .    ' NSPOT2=', i3, ' MSPCD=', 6I2, ' MEQN=', 3I2/
     .    ' KLANS1=', i3, ' MPHS=', 9I1, ' NSKY=', i3, ' SERIES=',A4)
         Line = Line + 3
      endif
c
c precession,wobble integers
      Kindnp = 0
      Iwob   = 0
 
c zero all propagation corrections
      call PRPZRO
 
c set up initial guess of distance
      Tguess = 100._10
      npd    = Nplnt0
      if(Klanb.gt.0) npd = Ncp0
      if(npd.gt.0 .and. npd.le.10) Tguess = dist0(npd)
      Dstf(4) = 0._10
 
c velocity scale factors
      do i = 1, 9
         corfct(i + 9) = 0._10
         corfct(i)     = 1._10
      end do
c set up iteration accuracy constant
c default value for optical types
      if(acctim.le.0.) acctim = 1E-3
      accprc = 10.
      accdst = acctim
      jtypob = ITYPOB(Ncodf)
      if(jtypob.eq.5) jtypob = 1
      if(jtypob.eq.1) nddiff = 0
      Ncod(1) = 2*jtypob - 1
      Ncod(2) = Ncod(1) + 1
      fdev    = Fdeva(Ntape)*1E-10_10
      ntime   = 0
      if(fdev.ne.0._10) ntime = 1
      do i = 1, 10
         Raddum(i) = 0._10
         Angdum(i) = 0._10
      end do
 
c clear y vectors and other temporaries
      call ZFILL(Ppcond,zyvect)
c initialize flag for phase delay doppler &
c filter header/trailer crossing control
      Recalc = .false.
c
c*  start=9000
      return
      end
