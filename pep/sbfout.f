      subroutine SBFOUT(ovrlap,icall)

      implicit none
c
c        subroutine sbfout
c     output routine for satellite-probe numerical integration
c        paul macneil  nov.  1977
c        modification of:
c     becker/ash   oct. 1968   subroutine sbout
c        includes iteration with kalman filter (pm, nov., 1977)
c        nordsieck vtp output oddities corrected pm,dec.,1977
c
c arguments
      logical*4 ovrlap
      integer*4 icall
c
c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'adams.inc'
      include 'dumdum.inc'
      include 'ellips.inc'
      include 'fcntrl.inc'
      include 'filtim.inc'
      include 'filtds.inc'
      include 'funcon.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'lothrf.inc'
      include 'nmstrt.inc'
      include 'output.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'rstart.inc'
      include 'sbstuf.inc'
      include 'sscon.inc'
      include 'stbtst.inc'
      include 'stint.inc'
      include 'zeroes.inc'
c data from satellite - probe numerical integration to be
c written on tape
      common/WRKCOM/ Sb(6,i_mxeqn,5),Tabpt(3,8),V01(6*i_mxeqn),
     . V02(6*i_mxeqn),Ictran,Invtrn,Setp(18),Xnew(6),W(40),
     . Wzeror(6),Flags(20)
      real*10 Sb,Tabpt,V01,V02,Setp,Xnew,W,Wzeror
      logical*4 Flags
      real*10 dtabpt(3,8,i_mxeqn-1)
      equivalence (dtabpt,Sb(1,1,2))

c local variables
      real*10 juldat,frect,juldas
      real*10 djedz,dumjd,dy,fsb,fsb1,fsb1st,gmass,
     .          hc1,hc2,hcsave,hctest,rdv,rdvsv,ry,ry2,ry3,
     .          tnext,tprev
      real*10 tsave,ttab,ttest,tthis,wwt
      integer   i,iend,iep,iepoch,ifmt,ii,ij,indic,
     .          iorder,iordr1,iovlap,ivl,j,j1,je,jf,
     .          jorder,jordr1
      integer   jplntt,jx,k,katmos,kk,kk95,kkp2,kp98,kx,ll,lx,
     .          matmos,mgo,mtabx,natpnp,nepp1,nlflag,
     .          npcont,npnpx2
      integer   nrcfil,numprt,numpt1

      character*19 sbfmt/'(F13.6,xxP,6F20.16)'/
      character*8 sbfmth(17)/'( '' JED-', 'xxxxx00''', ',8X,',
     .          ''' X*10**', 'xx'',11X,', ''' Y*10**', 'xx'',11X,',
     .          ''' Z*10**', 'xx'',09X,', ''' DX/DT*', '10**xx'',',
     .          '7X,', ''' DY/DT*', '10**xx'',', '7X,', ''' DZ/DT*',
     .          '10**xx'')'/
      integer*4 jedz/24000/
      integer*4 nexpln/10000/
      integer*4 nic/6/

c local for filter
      logical*4 endflg

c variables for keplerian interpolation
      real*10   Ictran(6,6),Invtrn(6,6),epse/1.0e-55_10/,decro
      logical*4 offdon

c external functions
      real*10 DOT
c
c*start=1010
c
      icall = icall+1
      hc2   = Hc
      if(icall.gt.1) ttab = je+fsb
      if(icall.le.1) then
         ttab   = T0
         jplntt = Jplnt
         if(Iboth.eq.0) jplntt = Ibuf
         if(Iboth.eq.0 .and. Kp(99).ge.0 .and. Kp(88).gt.0) goto 1000
         nrcfil = 0
      endif

c skip if not in filter mode
      if(Ict(42).le.0) goto 90

      if(icall.ne.1) then
c
c decrement counter for subroutine eval
         if(Nnotst.gt.0) Nnotst = Nnotst-1

c was overlap detected during an earlier entry?
         if(ovrlap) then
c
c count records
            iovlap = iovlap+1

c is the extended starting procedure in use?
            if(Filflg(1) .and. iepoch.lt.Nepoch .and. T.ge.tnext .and.
     .         Iterat.le.1 .and. Jct(56).lt.2 .and. .not.offdon) then
c
c read wzero data set
               read(Wzero) dumjd,Pnprms,Flags
c
c setup for keplerian interpolation
               gmass = Gauss*SQRT(Mass(Ncentr))
               call JNITL(gmass,Pnprms(natpnp+1),Setp,0,dy)

c interpolate
               call JLIPT(T-tnext,Setp,1,Xnew,ry,ry2,ry3,dy)

c find offset and partials of state wrt ics
               do i = 1,6
                  Wzeror(i) = Xnew(i)-Y(i,3)
                  do j = 1,6
                     k = i+6*j
                     Ictran(i,j) = Y(k,3)
                  end do
               end do

c invert to get partials of ics wrt state
               call ECROUT(Ictran,Invtrn,decro,epse,6,-1,6,6,indic)

c find w (product of partials of ic wrt state with delta state)
               do i = 1,6
                  wwt = 0._10
                  if(Flags(natpnp+i)) then
                     do j = 1,6
                        wwt = wwt+Invtrn(i,j)*Wzeror(j)
                     end do
                  endif
                  W(Npnp+natpnp+i) = wwt
               end do
c
c other (non-ic) parameters
               if(natpnp.gt.0) then
c
c atmospheric density
                  if(Flags(natpnp)) W(Npnp+natpnp)
     .                = Pnprms(natpnp)-Rhoz
               endif
c
c write w
               write(Iconof) (W(i),i=1,npnpx2)
               offdon = .true.
            endif

c test for last integration point
         else if(ttab.lt.Tstop) then
c guarentee that there will be at least three points
c ahead of next filter span
            hctest = 12.1*Hc
            if(T.ge.(tnext-hctest)) then

c overlap starts
               ovrlap = .true.

c end of last epoch?
               if(iepoch.lt.Nepoch) then
c not last epoch
c do saves
c save time parameters
                  tsave  = T
                  hcsave = Hc

c save state variables and derivatives
                  do i = 1,N
                     A0(i)  = A(i)
                     B0(i)  = B(i)
                     C0(i)  = C(i)
                     D0(i)  = D(i)
                     Ee0(i) = Ee(i)
                     V01(i) = Y(i,1)
                     V02(i) = Y(i,2)
                     V0(i)  = Y(i,3)
                     F0(i)  = F(i,3)
                  end do
               endif
            endif
         endif
      else
         natpnp = Npnp-nic
         if(natpnp.gt.0) then

c atmospheric process noise parameters
            matmos = 100*Ncentr+28
            katmos = 0
            do i = 1,Kount
               if(Icntrl(i).eq.matmos) katmos = i
               if(katmos.gt.0) goto 20
            end do
            call SUICID(
     .    'PARTIALS WRT ATMOSPHERIC RHOZ NOT INTEGRATED SBFOUT ',13)
         endif

   20    npnpx2 = Npnp*2
         call ZFILL(W,16*npnpx2)
         rewind Iconof
c test input parameters
c the next 14 cards may be deleted experimentally pem nov.,78
c changed to test fep rather than timez & delta   (zmg 6/81)
         ttest = Fep(1)
         if(ttest.lt.T0) then
            write(Iout,30) ttest,T0
   30       format(' FILTER BEGINS AT',f19.8,5x,'T0 =',f19.8)
            call SUICID(
     .        'FILTER RANGE NOT WITHIN INT RANGE, STOP SBFOUT  ',12)
         endif
         ttest = Fep(Nepoch+1)
         if(ttest.gt.Tstop) then
            write(Iout,40) ttest,Tstop
   40       format(' FILTER ENDS AT',f19.8,5x,'TSTOP =',f19.8)
            call SUICID(
     .        'FILTER RANGE NOT WITHIN INT RANGE, STOP SBFOUT  ',12)
         endif
c
c initialize
         tprev  = 0._10
         tthis  = Fep(1)
         tnext  = Fep(2)
         offdon = .false.
         iovlap = 0
         iend   = 0
         endflg = .false.
         iepoch = 1

         nlflag = 0
         if(iepoch.eq.Nepoch) nlflag = 1

c write first header/trailer
         if(Kp(99).ge.0) then
            write(jplntt) tprev,tthis,tnext,iepoch,nlflag,nrcfil
            Nrec = Nrec+1
         endif
         write(Iout,60) tprev,tthis,tnext,iepoch,nlflag,nrcfil
   60    format(' FILTER HEADER/TRAILER: TPREV=',f15.6,' TTHIS=',
     .          f15.6,' TNEXT=',f15.6,' IEPOCH=',i4/' NLFLAG=',
     .          i2,' NRCFIL=',i3)
         Line = Line+2

c copy epochs to w data set
         Itrwnd(Iconof) = 1
         nepp1 = Nepoch+1
         if(Iterat.eq.1 .and. Jct(56).eq.1) write(Iconof)
     .      (Fep(i),i=1,nepp1)

c in pvo starting procedure, skip first epoch
         if(Filflg(1)) then
            rewind Wzero
            Itrwnd(Wzero) = 1
            read(Wzero)
         endif
c in later iterations skip epochs
c
         if(Iterat.gt.1 .or. Jct(56).ge.2) read(Iconof)
      endif
c
c        end of extended starting procedure code segment
c
c        continue sbout functions
c*start=1990
c
c
   90 if(Ntab.le.0) then

c insert initial conditions into first tabular point
         Ntab  = 1
         mtabx  = Ntab
         rdvsv = 0._10
         ifmt = Kkp(10)
         ifmt = max0(ifmt,-9)
         ifmt = min0(ifmt,99)
         call EBCDIX(ifmt,sbfmt,8,2)
         call EBCDIX(jedz,sbfmth,  9,5)
         call EBCDIX(ifmt,sbfmth, 33,2)
         call EBCDIX(ifmt,sbfmth, 49,2)
         call EBCDIX(ifmt,sbfmth, 65,2)
         call EBCDIX(ifmt,sbfmth, 85,2)
         call EBCDIX(ifmt,sbfmth,109,2)
         call EBCDIX(ifmt,sbfmth,133,2)
         call NEWPGS(Name,Jdp1,Jdp2,sbfmth)
         djedz      = 100._10*jedz
         if(Kout.gt.0) call SBEXPS(djedz,jedz)
         npcont     = 1000000
         kp98   = Kp(98)
         kp98   = iabs(kp98)
         iorder = 8
         jorder = 4
         ivl    = 6
         iordr1 = iorder-3
         numprt = Iparp-1
         numpt1 = numprt
         if(numprt.le.0) then
            jorder = 1
            numpt1 = 1
         endif
         jordr1 = jorder-3
         jx     = 1
         do i = 1,6
            Sb(i,1,1) = X0(i)
         end do
         lx = 6
         if(jx.lt.Iparp) then
            if(Ki(1).ne.0) then
               do j = 1,6
                  if(Ki(j+1).ge.0) then
                     jx = jx+1
                     do i = 1,6
                        if(Ki(1).gt.0) lx = lx+1
                        Sb(i,jx,1) = Dx0(i,j)
                     end do
                  endif
               end do
            endif
            jx = jx+1
            if(jx.le.Iparp) then
               do j = jx,Iparp
                  do i = 1,6
                     lx = lx+1
                     Sb(i,j,1) = V0(lx)
                  end do
               end do
            endif
         endif
         Line   = 60
         frect  = 0.5_10+Hmx
         juldas = T-frect-djedz
         fsb1st = (T-Hmx)-Jdp0
c
c setup initial time record for initial conditions
         if(Kp(99).ge.0) then
            juldat = T-Hmx
            je     = juldat
            fsb    = juldat-je
            if(Kp(88).le.0) then
               do i = 1,3
                  ii = i+3
                  Tabpt(i,1) = Sb(i,1,1)
                  Tabpt(i,2) = Sb(ii,1,1)
                  Tabpt(i,3) = ff0(ii,1)
                  if(iordr1.gt.0) then
                     do j = 1,iordr1
                        Tabpt(i,j+3) = aa0(ii,1,j)
                     end do
                  endif
                  if(numprt.gt.0) then
                     do k = 1,numprt
                        kk = k+1
                        dtabpt(i,1,k) = Sb(i,kk,1)
                        dtabpt(i,2,k) = Sb(ii,kk,1)
                        dtabpt(i,3,k) = ff0(ii,kk)
                        if(jordr1.gt.0) then
                           do j = 1,jordr1
                              dtabpt(i,j+3,k) = aa0(ii,kk,j)
                           end do
                        endif
                     end do
                  endif
               end do
               hc1 = Hmx
               if(Jdstrt.gt.0) hc1 = Hcc1
               hc2 = Hmx
            endif
         endif
      endif
c
c print out initial conditions
      if(Kp(99).gt.0) goto 600
      mgo    = 1
      npcont = npcont+1
      if(npcont.lt.kp98) goto 600
      npcont = 0
      if(Line.lt.58) goto 200
  100 call NEWPGT(Iout,Npage,jedz)
      Line = 2
      if(mgo.eq.2) goto 400
  200 write(Iout,sbfmt) juldas,(Sb(i,1,Ntab),i=1,6)
      Line = Line+1
      if(Kp(98).le.0) goto 600
      if(Iparp.le.1) goto 600
      mgo = 2
      jx  = 1
  300 jx  = jx+1
      if(Line.ge.58) goto 100
  400 write(Iout,500) jx,(Sb(i,jx,Ntab),i=1,6)
  500 format(i11,1x,1p,6D20.13)
      Line = Line+1
      if(jx.lt.Iparp) goto 300
c*start=2500
c
c write tabular point (saved during previous entry on int tape
  600 if(Kp(99).ge.0) then
         if(Kp(88).gt.0) then

            if(Ntab.lt.5) goto 700
            write(jplntt) je,fsb,ivl,
     .                    (((Sb(i,j,k),i=1,ivl),j=1,Iparp),k=1,5),
     .                    zero9
            Ntab = 0
            Nrec = Nrec+1
         else
            write(jplntt) je,fsb,iorder,jorder,numprt,numpt1,
     .                    ((Tabpt(i,j),i=1,3),j=1,iorder),
     .                    (((dtabpt(i,j,k),i=1,3),j=1,jorder),k=1,
     .                    numpt1),hc1,hc2,zero9
            Nrec   = Nrec+1
            nrcfil = nrcfil+1

c carry step forward for next record
            hc1 = hc2
         endif
         je  = T
         fsb = T-je
      endif
      if(Iboth.ge.0) then
         if(Nsign*((Jdp1-T)+Epsp(1)).le.0._10) Iboth = 1
      endif
  700 frect = 0.5_10
c
c*start=3000
      if(Kp(88).gt.0 .and. Kp(99).ge.0) Ntab = Ntab+1

c if writting last point, do not save a next point
      if(Nsign*(ttab-Tstop).ge.0._10) goto 1400

c save tabular point for print next entry
      jx     = 1
      kx     = 6
      juldat = T-frect

      juldas = juldat-djedz

      if(Kp(100).lt.0) then

c straight forward integration was employed
         do i = 1,6
            Sb(i,1,Ntab) = Y(i,3)
         end do
         if(Ki(1).lt.0) then
            Tlpt = T-Tlpt0
            goto 800
         else if(Ki(1).ne.0) then
            do j = 1,6
               if(Ki(j+1).ge.0) then
                  jx = jx+1
                  do i = 1,6
                     kx = kx+1
                     Sb(i,jx,Ntab) = Y(kx,3)
                  end do
               endif
            end do
         endif
         goto 900
      endif
c
c enckes method of integration was employed
      do i = 1,6
         Sb(i,1,Ntab) = Y(i,3)+Ylpt(i)
      end do
      if(Ki(1).gt.0) then
         do j = 1,6
            if(Ki(j+1).ge.0) then
               jx = jx+1
               do i = 1,6
                  kx = kx+1
                  Sb(i,jx,Ntab) = Y(kx,3)+Dylpt(i,j)
               end do
            endif
         end do
      endif
      if(Ki(1).ge.0) goto 900
  800 call ELIPT(-1,Tlpt)
      if(Icnd.eq.1) then
         do i = 1,6
            Dylpt(i,4) = Dylpt(i,4)-Dylpt(i,5)
            Dylpt(i,5) = Dylpt(i,5)-Dylpt(i,6)
         end do
      endif
      do j = 1,6
         if(Ki(j+1).ge.0) then
            jx = jx+1
            do i = 1,6
               Sb(i,jx,Ntab) = Dylpt(i,j)
            end do
         endif
c
c mean orbit integration was employed
      end do
c
c partials w.r.t. non-initial conditions
  900 jx = jx+1
      if(jx.le.Iparp) then
         do j = jx,Iparp
            do i = 1,6
               kx = kx+1
               Sb(i,j,Ntab) = Y(kx,3)
            end do
         end do
      endif
c
c
c
c*start=3500
c        save tabular point for tape write next entry
      if(Kp(99).ge.0) then
         if(Kp(88).le.0) then
            do i = 1,3
               ii = i+3
               Tabpt(i,1) = Sb(i,1,1)
               Tabpt(i,2) = Sb(ii,1,1)
               Tabpt(i,3) = ff(ii,1,3)
               if(iordr1.gt.0) then
                  ij = ii
                  do j = 1,iordr1
                     if(j.eq.5) ij = ii+2
                     Tabpt(i,j+3)  = aa(ij,1,j)
                  end do
               endif
               do k = 1,numprt
                  kk = k+1
                  dtabpt(i,1,k) = Sb(i,kk,1)
                  dtabpt(i,2,k) = Sb(ii,kk,1)
                  dtabpt(i,3,k) = ff(ii,kk,3)
                  if(jordr1.gt.0) then
                     ij = ii
                     do j = 1,jordr1
                        if(j.eq.5) ij  = ii+2
                        dtabpt(i,j+3,k)= aa(ij,kk,j)
                     end do
                  endif
               end do
            end do
         endif
      endif
c
c*start=5000
c
      mtabx = Ntab
      if(Iboth.le.0) goto 1300

c backspace buffer, write probe data set
      jplntt = Jplnt
      if(Kp(99).lt.0) goto 1100
      je   = Jd
      fsb1 = Fract
      if(Kp(88).le.0) then
         do while( .true. )
            backspace Ibuf
            read(Ibuf) je,fsb,iorder,jorder,numprt,numpt1,
     .                 ((Tabpt(i,j),i=1,3),j=1,iorder),
     .                 (((dtabpt(i,j,k),i=1,3),j=1,jorder),k=1,
     .                 numpt1),hc2,hc1
            hc2 = -hc2
            hc1 = -hc1
            backspace Ibuf
            write(Jplnt) je,fsb,iorder,jorder,numprt,numpt1,
     .                   ((Tabpt(i,j),i=1,3),j=1,iorder),
     .                   (((dtabpt(i,j,k),i=1,3),j=1,jorder),k=1,
     .                   numpt1),hc1,hc2,zero9
            if(fsb.eq.fsb1st) then
               if(je.eq.Jdp0) then
                  rewind Ibuf
                  goto 1100
               endif
            endif
         end do
      endif
 1000 do while( .true. )

c no copy for fixed interval
         call SUICID(
     .'NO BACKWARDS/FORWARDS INTEGRATION WITH OUTPUT TAPE FOR FIXED INTE
     .RVAL, STOP SBFOUT  ',21)
         jf  = je
         fsb = fsb1
         if(fsb.eq.fsb1st) then
            if(je.eq.Jdp0) then
               rewind Ibuf
               goto 1100
            endif
         endif
      end do
 1100 Iboth = -1
      Nsign = -Nsign
      Hc    = -Hsave
      Hmx   = -Hmx
      Pint5 = 4._10*Hmx
      Tstop = Jdp2
      Tstop = Tstop+Epsp(2)+Hmx
      if(Kp(88).gt.0) Tstop = Tstop+Pint5
      Both = ABS(T-T0)
      T    = T0
      L4   = -10
      L1   = -1
      Ntab = -1
      do i = 1,N
         A(i)  = A0(i)
         B(i)  = B0(i)
         C(i)  = C0(i)
         D(i)  = D0(i)
         Ee(i) = Ee0(i)
         Dydt(i,1) = Dydt0(i)
      end do
      write(Iout,1200)
 1200 format(
     .'0INTEGRATION COMPLETED IN ONE DIRECTION FROM EPOCH, STARTED IN OT
     .HER'/)
      Line = Line+3
c*start=6000
c
 1300 if(Kout.gt.0) then

         if(Kkp(2).lt.0) then
            kkp2 = -Kkp(2)
            rdv  = DOT(Sb(1,1,mtabx),Sb(4,1,mtabx))
            if(rdv*rdvsv.le.0._10) then
               nexpln = nexpln+1
               if(nexpln.ge.kkp2) then
                  nexpln = 0
                  call SBEXP(Sb,T,mtabx,sbfmt)
               endif
            endif
            rdvsv = rdv
         else if(Kkp(2).ne.0) then
            nexpln = nexpln+1
            if(nexpln.ge.Kkp(2)) then
               nexpln = 0
               call SBEXP(Sb,T,mtabx,sbfmt)
            endif
         endif

c printout for close approaches
         if(Kkp(6).gt.0) call SBEXPT(T)
      endif

c
c setup starting proceedure for thrust initiation or
c termination as signaled by subroutine lesthp
      Jthrst = 0
      Lthrst = 0
      T0     = T0sav
      if(Kthrst.gt.0) then
         Kthrst = 0
         Jthrst = 1
         Lthrst = 1
c do not set approximating polynomial coefficients to zero
c in starting proceedure
         T0   = T
         kk95 = Kkp(95)
         Hc   = 2._10**kk95

c station keeping thrusting assumes forward in time integration
         do i = 1,N
            V0(i) = Y(i,3)
         end do
         if(Kp(100).ge.0) call SUICID(
     .' CANNOT RESTART ENCKE INTEGRATION FOR THRUST, STOP IN SBFOUT',15)
      endif
c
c*start=7000
c
c skip if not in filter mode
 1400 if(Ict(42).gt.0) then

c last integration point? if so, write trailer
         if(ttab.ge.Tstop) then
c
c end of integration
c has last header/trailer been written?
            if(endflg) return
            endflg = .true.
            rewind Iconof
            tprev = Fep(Nepoch)
            tthis = Fep(Nepoch+1)
            tnext = 0._10
c
c make (permanent) wzero an exact copy of iconof,
c in case of 1st iteration restart (i.e., analiz rerun)
            if(Filflg(1)) then
               rewind Wzero
               read(Iconof) (Fep(i),i=1,nepp1)
               write(Wzero) (Fep(i),i=1,nepp1)
               do iep = 2,Nepoch
                  read(Iconof) (W(i),i=1,npnpx2)
                  write(Wzero) (W(i),i=1,npnpx2)
               end do
               rewind Iconof
               rewind Wzero
               Itrwnd(Iconof) = 0
               Itrwnd(Wzero)  = 0
            endif
         else
            if(ovrlap) then
               if(ttab.gt.tnext) iend = iend+1

c following test applies if last tnext.lt.tstop
               if(iepoch.lt.Nepoch) then
                  if(iend.eq.2) then
c
c two end points written out, one more ready to be
c written during next entry; do restores to shift integrator
                     do i = 1,N
                        A(i)    = A0(i)
                        B(i)    = B0(i)
                        C(i)    = C0(i)
                        D(i)    = D0(i)
                        Ee(i)   = Ee0(i)
                        Y(i,1) = V01(i)
                        Y(i,2) = V02(i)
                        Y(i,3) = V0(i)
                        Y(i,4) = V0(i)
                        F(i,3) = F0(i)
                     end do
                     T  = tsave
                     Hc = hcsave

c prevent interval halving for nine steps after stepover
                     Nnotst = 0
                     Nnotst = 9

                     if(.not.Filflg(1)) then
                        call ZFILL(W,16*npnpx2)

c no ic transformation for first iteration unless restarting
                        if(Iterat.eq.1 .and. Jct(56).lt.2) return
c read ic's and delta ic's from permanent data set
c values for next epoch (iepoch not yet incremented)
                        Itrwnd(Iconof) = 1
                        read(Iconof) (W(i),i=1,npnpx2)
                     endif

c print filter quantities for debug
                     ll = 3+(Npnp-1)/6*2
                     if(Line+ll.ge.58) then
                        call NEWPGT(Iout,Npage,jedz)
                        Line = 2
                     endif
                     write(Iout,1410) Npnp,(W(i),i=1,Npnp)
 1410                format(' -- FILTER STATE --  NPNP=',i3,'   W='/
     .                      (1x,1p,6D22.15))
                     write(Iout,1420) (W(Npnp+i),i=1,Npnp)
 1420                format(1x,1p,6D22.15)
                     Line = Line+ll
c
c ****************experiment nov.20, 1978 *********************
                     Nnotst = 0
c ***********************************************************
c
c
c        transform state and partials wrt ic's
c        for change in initial conditions
                     do i = 1,6
                        do j1 = 1,nic
                           k = i+6*j1

c parameter deltas follow adjustments; deltas used
                           j     = j1+Npnp+natpnp
                           A(i)  = A(i)+A(k)*W(j)
                           B(i)  = B(i)+B(k)*W(j)
                           C(i)  = C(i)+C(k)*W(j)
                           D(i)  = D(i)+D(k)*W(j)
                           Ee(i) = Ee(i)+Ee(k)*W(j)
                           Y(i,1) = Y(i,1)+Y(k,1)*W(j)
                           Y(i,2) = Y(i,2)+Y(k,2)*W(j)
                           Y(i,3) = Y(i,3)+Y(k,3)*W(j)
                           F(i,3) = F(i,3)+F(k,3)*W(j)
                        end do
                        Y(i,4) = Y(i,3)
                     end do
c
c transform state and partials wrt ics for change in
c atmospheric parameter(s)
                     if(natpnp.gt.0) then
                        do i = 1,6
                           do j1 = 1,natpnp
                              k     = i+6*katmos
                              j     = j1+Npnp
                              A(i)  = A(i)+A(k)*W(j)
                              B(i)  = B(i)+B(k)*W(j)
                              C(i)  = C(i)+C(k)*W(j)
                              D(i)  = D(i)+D(k)*W(j)
                              Ee(i) = Ee(i)+Ee(k)*W(j)
                              Y(i,1) = Y(i,1)+Y(k,1)*W(j)
                              Y(i,2) = Y(i,2)+Y(k,2)*W(j)
                              Y(i,3) = Y(i,3)+Y(k,3)*W(j)
                              F(i,3) = F(i,3)+F(k,3)*W(j)
                           end do
                           Y(i,4) = Y(i,3)
                        end do
c "transform" atmospheric parameter(s)
c
                        call SBRSET(W,npnpx2,natpnp)
                     endif
                  else if(iend.ge.3) then
c
c        third end point written out, first point of next span
c        (if any) ready to be written during next entry, write
c        header/trailer
c        find new (next) times
                     tprev = tthis
                     tthis = Fep(iepoch+1)
                     tnext = Fep(iepoch+2)
                     goto 1450
                  endif
               endif
            endif
            return
         endif

 1450    iepoch = iepoch+1

c set next-to-last header flag
         nlflag = 0
         if(iepoch.eq.Nepoch) nlflag = 1

c write header/trailer
         if(Kp(99).ge.0) then
            write(jplntt) tprev,tthis,tnext,iepoch,nlflag,nrcfil
            Nrec = Nrec+1
         endif
         write(Iout,60) tprev,tthis,tnext,iepoch,nlflag,nrcfil
         Line   = Line+2
         nrcfil = 0
c
c reset overlap logic
         ovrlap = .false.
         offdon = .false.
         iovlap = 0
         iend   = 0
      endif
c
c*start=9990
      return
      end
