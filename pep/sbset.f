      subroutine SBSET

      implicit none

c ash/becker/amuchastegui/friedman - 12/69 - subroutine sbset
c setup computations for satellite or probe numerical integration
c
c array dimensions
      include 'globdefs.inc'
c common
      include 'adams.inc'
      include 'astroi.inc'
      include 'bddtaint.inc'
      include 'b2dtaint.inc'
      include 'cnthar.inc'
      include 'drgprt.inc'
      include 'ellips.inc'
      include 'emmips.inc'
      include 'empcnd.inc'
      include 'ethhar.inc'
      include 'fcntrl.inc'
      include 'fmmips.inc'
      include 'funcon.inc'
      include 'gleaks.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'lothrf.inc'
      include 'monhar.inc'
      include 'namtimq.inc'
      include 'nmstrt.inc'
      include 'param.inc'
      include 'petina.inc'
      include 'petuna.inc'
      real*10 ax,ay,az,bz,ta,tb,tc,alpha,delta,beta,tzero,tzero2
      equivalence (ax,Con),(ay,Con(2)),(az,Con(3)),(bz,Con(4)),
     .            (ta,Con(5)),(tb,Con(6)),(tc,Con(7)),
     .            (alpha,Con1),(delta,Con1(2)),(beta,Con1(3)),
     .            (tzero,Con1(4)),(tzero2,Con1(6))
      include 'plnhar.inc'
      include 'prtpin.inc'
      include 'rstart.inc'
      include 'rtpars.inc'
      include 'sbstuf.inc'
 
c grav.grad control flags
      logical*1 kggfgs(16)/6*.true., 2*.false., 2*.true., 6*.false./
      logical*1 lggfgs(16)
      equivalence (Ggfgs,lggfgs)
 
      include 'scoef4.inc'
      include 'smlbdy.inc'
      include 'smlstf.inc'
      include 'sscon.inc'
      include 'stint.inc'
      include 'thratt.inc'
      include 'trghar.inc'
      include 'xprcom.inc'
      integer*4 zxprcm/2982/   !r8=1494,r10=2982

c
c local
      real*10 alphar,beta2,betar,deltar,den,dhsi,di4,dum,eta,
     . etutc,fcp0,fire0,gauss9(12),har,hars,hmxx,hsia,hsib,
     . massc,pfact,pfact1,quick,save1,save2,sma,xh,xn,xsign
      integer*4 i,iden,int,isv,it1,it2,j,j1,jdd,k,
     .          khar,kii,king,kk,kkcnt,kkk,kktrg,kong,kong10,
     .          koumt,kt,ktong,ktop,ktopa,l,ll,mum,n1,
     .          n2,nc,ndadx,nmsg,npt,npz,nrckm,nrckp,nt,ntaaa,
     .          ntt,ntz,nzaaa
      integer*2 kastr(12),kp8792(6),lcplp
      logical*1 tessfg((u_mxtes*(u_mxtes+1))/2-1)
      character*8 npmsg(11)/'  NO TAR', 'GET PLAN', 'ET FOR  ', ' ',
     .          ' HARMONI', 'C PARTIA', 'L, STOP ', 'IN SBSET',
     .          ' ZONAL  ', 'TESS.SIN', 'TESS.COS'/
      character*56 npmsgm
      real*4    epsp2, epsp3
      real*10 hsftr/1E15_10/
      character*4 mesage(2)/' ON ', ' OFF'/

c external functions
      real*10 LEGSCL
 
c
c*       start=1000
c
c
c           are partial derivatives being integrated this iteration
      if(Kkp(100).gt.0 .and. Iterat.gt.Kkp(100)) then
         Numki = 8
         do i = 1, Numki
            Ki(i) = 0
         end do
      endif
c integration will do what ki(1-n) says to do when iterat=1 no
c matter what the value of kkp(100)
c
c see if this is checkpoint restart
      Jdstrt = 0
      if(Jdp0.le.0) then
         Jdstrt = -Jdp0
         epsp2  = Epsp(2)
         epsp3  = Epsp(3)
         do i = 87, 92
            kp8792(i-86) = Kp(i)
         end do
         call RESTRT(Jplnt,Nplnt,Ncentr,Intp,Jdp1,Jdp0,Jdp2,Cond,
     .               Con,Con1,Epsp,Kp,Numki,Ki)
         Epsp(2) = epsp2
         if(Kkp(93).gt.0) then
            Epsp(3) = epsp3
            do i = 87, 92
               Kp(i) = kp8792(i-86)
            end do
         endif
      endif
c
c setup of constants

      Intxpr= Kkp(7)
      Gama  = Gauss**2
      Gamat = Gama
      Masse = 1._10 - Mass(10)
      do i = 1, 9
         Mass1(i) = Mass(i)
      end do
      Mass1(10) = Mass(3)*Mass(10)
      if(Ncentr.eq.3 .or. Ncentr.eq.10. or. Kp(40).ge.0)
     .   Mass1(3) = Mass(3)*Masse
      Massp = 0._10
      massc = 1._10
      do i = 1, 30
         Massfl(i) = 0
      end do
      if(Ncentr.gt.0 .and. Ncentr.le.10) then
         Massfl(Ncentr) = 2
         if(Ncentr.eq.10) Massfl(3) = 2
         if(Ncentr.eq.3) Massfl(10) = -1
         massc = Mass(Ncentr)
         if(Ncentr.eq.10) massc = Mass1(10)
         do i = 1, Numpln
            j = Nqlnt(i)
c subtract satellite masses from system mass
            if(j.gt.10 .and. j.le.30 .and. Npcent(i).eq.Ncentr .and. 
     .         (Kp(j+30).ge.0 .or. j.eq.Nplnt)) then
               Mass1(Ncentr) = Mass1(Ncentr) - Mass(j)*massc
               Massfl(j)     = -1
            endif
         end do
      endif
      if(Nplnt.le.30) then
         Massp = Mass(Nplnt)
         Massfl(Nplnt) = Massfl(Nplnt) + 1
         call SBZENG(' PLANET ')
      else
         call SBZENG(Name)
      endif
      Massp1 = 1._10 + Massp
      Goose  = Gauss
      if(Ncentr.gt.0 .and. Ncentr.le.10) Massp1 = Mass1(Ncentr)/massc +
     .         Massp
      Goose  = Gauss*SQRT(massc*Massp1)
      Gama3  = -Gauss*Gauss*massc*Massp1
      Gamt30 = Gama3
      Gamat3 = Gama3
      int    = Prmter(97) + 0.500001_10
      Tvary0 = int
      if(Kp(62).ge.0) then
         pfact1 = Gmvary*(Jdp0-Tvary0)
         pfact  = 1._10 + pfact1
         Gamat  = Gama*pfact
         Gamat3 = Gama3*pfact
         Gamt30 = Gamat3
      endif
c the gaussian grav const for ellipt. orb elem. has no time var.
c
c set all initial conditions equal to zero
      do j = 1, 6*i_mxeqn
         V0(j) = 0._10
      end do
c
c zero all grav. grad. matrices (default, in case not used)
      ndadx = 9*(5 + i_mxtrg)
      do j = 1, ndadx
         Dadxc(j,1) = 0._10
      end do
c
c           setup initial condition for equations of motion and
c           equations for partial derivatives w.r.t. initial osculating
c           eliptic orbital elements
c           or initial position and velocity (depending on icnd)
c*       start=1300
      if(Icnd.lt.0) then
         do j = 1, 6
            X0(j) = Cond(j)
            do i = 1, 6
               Dx0(i,j) = 0._10
            end do
            Dx0(j,j) = 1._10
         end do
         Kp(100) = -1
         if(Ki(1).lt.0) Ki(1) = 1
 
      else if(Icnd.ge.2) then
 
         call RADVAF(Cond, X0, Dx0, Goose, 1)
         Kp(100) = -1
         if(Ki(1).lt.0) Ki(1) = 1
      else
         save1 = Cond(5)
         save2 = Cond(6)
         if(Icnd.ne.0) then
            Cond(6) = Cond(6) - Cond(5)
            Cond(5) = Cond(5) - Cond(4)
         endif
         call INITL(Goose, 1._10, Cond, 0)
         Kp1 = Ki(1)
         Kp1 = iabs(Kp1)
         call ELIPT(Kp1, 0._10)
         if(Ki(1).lt.0) Kp1 = 0
         do j = 1, 6
            X0(j) = Ylpt(j)
            do i = 1, 6
               Dx0(i,j) = Dylpt(i,j)
               if(Icnd.eq.1) then
                  if(j.eq.4) Dx0(i,4) = Dylpt(i,4) - Dylpt(i,5)
                  if(j.eq.5) Dx0(i,5) = Dylpt(i,5) - Dylpt(i,6)
               endif
            end do
         end do
         if(Icnd.eq.1) then
            Cond(5) = save1
            Cond(6) = save2
            if(Ki(1).gt.0) Kp(100) = -1
         endif
      endif

c set up flags indicating which bodies are target bodies
      do i=1,30
         Npltg(i)=0
      end do
      do j=1,i_mxtrg
         nt=Kp(j)
         if(nt.gt.0) Npltg(nt)=j
      end do

      Nrec  = 0
      Ntab  = -1
      Iparp = 1
      N     = 6
      if(Kp(100).lt.0) then
c*  start=1400
c
         do j = 1, 6
            V0(j) = X0(j)
         end do
         if(Ki(1).gt.0) then
            do j = 1, 6
               if(Ki(j+1).ge.0) then
                  do i = 1, 6
                     N     = N + 1
                     V0(N) = Dx0(i,j)
                  end do
               endif
c*  start=1500
c
            end do
         endif
         goto 100
      else if(Kp(100).eq.0) then
c Encke method, osculating orbit. Starting values all zero
      endif
c mean orbit setup to be inserted here
c for now, treat as osculating orbit
c*  start=1600
c
      if(Ki(1).gt.0) then
         do i = 2, 7
            if(Ki(i).ge.0) N = N + 6
         end do
      endif
 
  100 if(Ki(1).ne.0) then
         do i = 2, 7
            if(Ki(i).ge.0) Iparp = Iparp + 1
         end do
      endif
c*  start=2000
c
c           setup for equation for partial derivations with respect to
c           non-initial conditions
c
c           setup of probe logic controls for partial derivatives
c           w.r.t. initial conditions
      Kount = 0
      if(Ki(1).gt.0) then
         do i = 2, 7
            if(Ki(i).ge.0) then
               Kount = Kount + 1
               Icntrl(Kount) = -30 - (i-1)
            endif
         end do
      endif
c*  start=2100
c
c           setup of probe logic controls for partial derivatives
c           w.r.t. other parameters
c     zero out harmonic vectors
      do i = 1,u_mxzon-1
         Iczone(i) = 0
      end do
      do i = 1, (u_mxtes*(u_mxtes+1))/2-1
         Iccos(i) = 0
         Icsin(i) = 0
      end do
      do i = 1, 4
         do j = 1, i_mxtrg
            Itzone(j,i) = 0
         end do
      end do
      do i = 1, 5
         do j = 1, i_mxtrg
            Itcos(j,i) = 0
            Itsin(j,i) = 0
         end do
      end do
 
c i.c. partial controls
      Lcnt = 0
      do i = 1, i_mxtrg
         Ltrg(i) = 0
      end do
c
c*  start=2300
c start of loop i=8,numki
      i = 7
  200 i = i + 1
      if(Ki(i).eq.0) goto 700
      kong = Ki(i)/100
      if(kong.le.0) goto 500
      king = Ki(i) - kong*100
 
c set up center or target index
      ktong = -2
      if(kong.le.30 .and. Npltg(kong).gt.0) ktong=Npltg(kong)
      if(kong.eq.Ncentr) ktong = 0
      if(king.le.6) then
         if(ktong.eq.0) Lcnt = 1
         if(ktong.gt.0) Ltrg(ktong) = 1
         if(ktong.lt.0) call SUICID(
     .    'NO TARGET PLANET FOR I.C. PARTIAL, STOP IN SBSET',12)
      endif
 
c zonal harmonics
      if(king.ne.31) then
c*  start=2400
c tesseral cosine harmonics
         if(king/10.ne.4) then
c*  start=2600
c tesseral sine harmonics
            if(king/10.ne.5) goto 500
 
c same treatment as cosine harmonics
            nmsg = 10
         else
c king=40-44  ki(i+*)=n1,m1,n2,m2  full range of harmonics
c king=45-49  ki(i+*)=m,n1,n2  resonant harmonics
            nmsg = 11
         endif
 
c see if any body set up for harmonics
         if(ktong.lt.0) goto 400
         do j = 1, (u_mxtes*(u_mxtes+1))/2-1
            tessfg(j) = .false.
         end do
         isv = Ki(i)/10
      else
         nmsg = 9
 
c see if any body set up for harmonics
         if(ktong.lt.0) goto 400
         n1 = Ki(i+1) - 1
         n2 = Ki(i+2) - 1
         do j = n1, n2
            Kount = Kount + 1
            Icntrl(Kount) = Ki(i)
            Iparp = Iparp + 1
            N     = N + 6
            if(kong.ne.Ncentr) then
               Itzone(ktong,j) = Kount
            else
               Iczone(j) = Kount
            endif
         end do
         i = i + 2
         goto 600
      endif
      do while( .true. )
         kong10 = Ki(i) - isv*10
         if(kong10.gt.4) then
            n1 = Ki(i+2)
            n2 = Ki(i+3)
            do j = n1, n2
               jdd = j*(j - 1)/2 + Ki(i+1) - 1
               tessfg(jdd) = .true.
            end do
            i = i + 3
         else
            n1 = (Ki(i+1)*(Ki(i+1)-1))/2 + Ki(i+2) - 1
            n2 = (Ki(i+3)*(Ki(i+3)-1))/2 + Ki(i+4) - 1
            do j = n1, n2
               tessfg(j) = .true.
            end do
            i = i + 4
         endif
         if(i.gt.Numki) goto 610
         if(Ki(i+1)/10.ne.isv) then
 
            kii = 10*isv + 1
            do j = 1, (u_mxtes*(u_mxtes+1))/2-1
               if(tessfg(j)) then
                  Kount = Kount + 1
                  Icntrl(Kount) = kii
                  Iparp = Iparp + 1
                  N     = N + 6
                  if(king.ge.50) then
 
c sine harmonics
                     if(kong.ne.Ncentr) then
                        Itsin(ktong,j) = Kount
                     else
                        Icsin(j) = Kount
                     endif
 
c cosine harmonics
                  else if(kong.ne.Ncentr) then
                     Itcos(ktong,j) = Kount
                  else
                     Iccos(j) = Kount
                  endif
               endif
            end do
            goto 600
         else
            i = i + 1
         endif
      end do
c*  start=2700
c
c no target planet set up for harmonic partial, abort here
  400 npmsg(4) = npmsg(nmsg)
      call SUICID(npmsg,16)
c*  start=2800
c usual parameter
  500 Kount = Kount + 1
      Icntrl(Kount) = Ki(i)
      Iparp = Iparp + 1

c need elliptic package for indirect term in partial w.r.t. mass of a
c target body
      j = Ki(i)
      if(j.gt.0 .and. j.le.30) then
         kt=Npltg(j)
         if(kt.gt.0) then
            Ltrg(kt)=1
         endif
      endif

c non-zero velocity initial condition for partial with respect to
c mass of integrated body or mass of central body (for elliptic
c orbital element initial conditions of equations of motion)
      if(j.le.0 .or. j.gt.30) then
         N = N + 6
      else if(Massfl(j).eq.0) then
         N = N + 6
 
c icnd=2 was never treated before
      else if(Icnd.lt.0) then
         N = N + 6
      else
         if(Massfl(j).lt.1) then
            quick = -0.5_10/Massp1
         else if(Massfl(j).eq.1) then
            quick = 0.5_10/Massp1
         else
            quick = 0.5_10/Mass(j)
         endif
         if(Icnd.eq.2) quick =quick*Goose/(Cond(4)*SQRT(Cond(1))+Goose)
         N = N + 3
         do j = 1, 3
            N     = N + 1
            V0(N) = quick*X0(j+3)
         end do
      endif
c*  start=2900
c end of loop i=8,numki
  600 if(i.lt.Numki) goto 200
      if(i.eq.Numki) goto 700
  610 call SUICID(
     .' INDEX FOR KI PARTIAL VECTOR EXCEEDS MAXIMUM, STOP IN SBSET ',15)
 
  700 if(Kount.gt.999) call SUICID(
     .        'PARTIAL COUNT GREATER THAN 999,STOP IN SBSET', 11)
c zero out remainder of Icntrl vector
      if(Kount.lt.999) then
         koumt = Kount + 1
         do i = koumt, 999
            Icntrl(i) = 0
         end do
      endif
c*  start=3000
c set up flags for grav.grad matrices  (dadx#)
c add input control here
      do i = 1,16
         lggfgs(i) = kggfgs(i)
      end do
      Npstep = -1
      Mpstep = Kkp(50)
      Pntflg = .false.
c
c numerical integration setup
      L1    = 0
      M     = Kp(90)
      T1    = Jdp1
      T1    = T1 + Epsp(1)
      T2    = Jdp2
      T2    = T2 + Epsp(2)
      xsign = SIGN(1._10, T2-T1)
      Nsign = xsign
      xh    = Kp(91)
      kk    = Kp(91)
      if(xh.le.0) then
         Hc = 2._10**kk*xsign
      else
         Hc = xh*xsign
      endif
      xn  = Kp(92)
      kkk = Kp(92)
      if(xn.le.0) then
         Hmn = 2._10**kkk
      else
         Hmn = xn
      endif
      if(Intp.le.0) then
         Hmx = (2._10**Intp)*xsign
      else
         Hmx = Intp*Nsign
      endif
      Epsi = Epsp(3)
      Epsz = Epsp(3)
c
c initial time setup
      T0   = Jdp0
      fcp0 = 0._10
      if(Intp1.ne.0) then
         fcp0 = Intp1
         fcp0 = fcp0*2._10**Intp2
         T0   = T0 + fcp0
      endif
      Tlpt0 = T0
      Pint5 = 4._10*Hmx
      Tstop = T2 + Hmx
      if(Kp(88).gt.0) Tstop = Tstop + Pint5
c variable nordsieck can stop right on the button
c
c extra setup for adams-moulton integration
      Npredt = Kp(89)
      Ncoret = Kp(89)
      Neq    = N
      Nh     = Kp(87)
      Inthmx = Intp
      Itype  = Kp(88)
c
c extra setup for nordsieck integration
      Iptr3 = 3
      if(Kp(88).gt.1) then
         Intm = 0
         Ints = 0
      else
         Intm = Kkp(88)
         Ints = Kkp(87)
         if(Intm.lt.0) Iptr3 = 2
         if(Kp(90).eq.6) then
            do i = 9,16
               lggfgs(i) = .false.
            end do
         endif
      endif
 
      if(Ict(42).gt.0) then
         if(Ki(1).le.0) call SUICID(
     .       'MUST INTEGRATE STATE & ICS, STOP SBSET  ', 10)
         do i = 2, 7
            if(Ki(i).lt.0) call SUICID(
     .       'MUST INTEGRATE PARTIALS WRT ICS, STOP SBSET ', 11)
         end do
      endif
c*  start=3200
c extra setup for pseudo corrector
c change initial conditions if this is checkpoint restart
c
      Both  = 0._10
      Iboth = -1
      if(Jdstrt.gt.0) then
         if(Iparp.ne.Iparst) call SUICID(
     .       'IPARP AND IPARST DO NOT AGREE, STOP IN SBSET', 11)
         if(Kp(100).ge.0) then
            T0 = (Jdstrt - Jdp0) + (Frstrt - fcp0)
            call ELIPT(1, t0)
         endif
         T0 = Jdstrt + Frstrt
         T1 = T0
         do i = 1, 6
            X0(i) = Rstrt(i,1)
            if(Kp(100).lt.0) then
               V0(i) = Rstrt(i,1)
            else
               V0(i) = Rstrt(i,1) - Ylpt(i)
            endif
         end do
         j1 = 1
         l  = 6
         if(Ki(1).ne.0) then
            do j = 1, 6
               if(Ki(j+1).ge.0) then
                  j1 = j1 + 1
                  do i = 1, 6
                     Dx0(i,j) = Rstrt(i,j1)
                     if(Ki(1).gt.0) then
                        l = l + 1
                        if(Kp(100).ge.0) then
                           V0(l) = Rstrt(i,j1) - Dylpt(i,j)
                           goto 810
                        endif
                     endif
                     V0(l) = Rstrt(i,j1)
  810             end do
               endif
            end do
         endif
         do while( j1.lt.Iparst )
            j1 = j1 + 1
            do i = 1, 6
               l     = l + 1
               V0(l) = Rstrt(i,j1)
            end do
         end do
c*  start=3400
c
c alter numerical integration setup for integrating forward
c and backward from initial epoch
      else if(ABS(T0-T1).ge.MIN(1._10,ABS(Hmx))) then
         Iboth = 0
         Nsign = -Nsign
         Hc    = -Hc
         Hmx   = -Hmx
         Pint5 = 4._10*Hmx
         hmxx  = Hmx
         if(ABS(Hmx).lt.1._10) hmxx = SIGN(1._10, Hmx)
         Tstop = T1 + 1E2_10*hmxx
 
         if(Jct(56).gt.0) call SUICID(
     .'CANNOT INTEGRATE BACKWARDS WITH FILTER, STOP SBSET  ', 13)
      endif
 
c setup kcnt for planet central body
      Kcnt = 0
      if(Ncentr.gt.0) then
         do i = 1,u_mxpl
            if(Nqlnt(i).eq.Ncentr) then
               Kcnt = i
               goto 900
            endif
         end do
      endif
c
c setup for time varying gravitational factor
  900 if(Kp(62).lt.1) then
      endif
c
c setup for violation of principle of equivalence
      Jpoev = 0
c metric parameters beta or gamma not equal to one give
c principle of equivalence violation (unless ratio of
c gravitational to total energy is zero), as do non-zero
c deltas or deltap
      eta    = 4._10*Betapr - Gamapr - 3._10
      Dltbod = deltap + eta*Con1(10)
      if(Dltbod.ne.0._10) Jpoev = 1
      Dltsun = Deltas + eta*Sunpoe
      if(Dltsun.ne.0._10) then
         Jpoev = 1
      else if(Deltas.ne.0._10 .or. deltap.ne.0._10) then
         Jpoev = 1
      else
 
c check for integration of eqns of motion and partials
         do i = 8, Numki
            if(Ki(i).eq.0) goto 1000
            if(Ki(i).eq.-20 .or. Ki(i).eq.40) then
               Jpoev = 1
               goto 1000
            else if(Ki(i).eq.43 .or. Ki(i).eq.44) then
               Jpoev = 1
               goto 1000
            endif
         end do
      endif
 1000 if(Jpoev.gt.0 .and. Ncentr.gt.0) call SUICID(
     . 'PRINCIPLE OF EQUIVALENCE VIOLATION NOT IMPLEMENTED FOR SATELLITE
     .S. STOP IN SBSET', 20)
c
c general relativity setup
      Rlfact = 0._10
      Index1 = 3
      Index2 = 3
      Jndex1 = 3
c these indices control velocity calculations for various vectors
c index1, index2 for motion and partials for probe wrt central body
c jndex1 for probe and central body wrt sun
      if(Kp(61).ge.0) then
         Index1 = 6
         if(Ncentr.gt.0) Jndex1 = 6
         Rlfact = Relfct
         if(Kp(97).gt.0) Rlfact = Con(24)
         if(Kp(61).gt.0) Index2 = 6
         Alph4  = 4._10*(Gauss*Aultsc/86400._10)**2
         Cvel2  = (86400._10/Aultsc)**2
         Alph16 = 4._10*Alph4
         B2g2   = 2._10*(Betapm + Gamapm)
         B2g21  = B2g2 + 1._10
         B2m1   = 2._10*Betapm - 1._10
         G21    = 2._10*Gamapm + 1._10
         G22    = G21 + 1._10
         G11    = Gamapm + 1._10
         G7     = (4._10*Gamapm + 3._10)/2._10
         A44    = Alph4/4._10
 
         do i = 1,10
            Xmass(i) = Mass(i)
            Kgr(i)   = Kp(30+i)
         end do
         if(Nplnt.le.9) Kgr(Nplnt)=-1
         Xmass(0) = 1._10
         Kgr(0) = 1
         Xmass(11) = 0._10
         if(Nplnt.le.30) Xmass(11)=Mass(Nplnt)
         Xmsb=Xmass(11)
         Kgr(11) = 1
         if(Kgr(10).ge.0) then
            Xmass(10)=Mass(3)*Mass(10)
            Xmass(3) =Mass(3)-Xmass(10)
            if(Nplnt.eq.3) Xmass(11)=Xmass(3)
         endif
      endif
c
c setup for second harmonic of sun
      if(Kp(63).ge.0) call SJ2SET
c
c setup for second harmonic of integrated planet
      if(Nplnt.le.30 .and. Kp(82).ge.0) then
         if(Nplnt.ne.3) then
            do i = 1, 4
               if(Nplnt.eq.Nplhar(i)) then
                  Phar2 = Pzhar(i,1)
                  goto 1020
               endif
            end do
         else
            Phar2 = Ezhar(1)
         endif
 1020    Phar2 = Phar2*(Con(1)/Ltvel/Aultsc)**2
      endif
c*  start=3600
c
c setup for perturbing asteroid or satellite orbits
c (elliptic orbits assumed for bodies between 11 and 30)
      Numast = 0
      do j = 11, 30
         if(Kp(j+30).ge.0) then
            if(j.ne.Nplnt) then
               Numast = Numast + 1
               if(Numast.gt.12) call SUICID(
     .        ' MORE THAN 12 PERTURBING ASTEROIDS, STOP IN SBSET   ',13)
               Kpast(Numast) = Kp(j+30)
               do k = 1, Numpln
                  if(Nqlnt(k).eq.j) then
                     Kast = k
                     kastr(Numast) = Kast
                     goto 1030
                  endif
               end do
               call SUICID(
     .            ' PERTURBING ASTEROID NOT INPUT, STOP IN SBSET   ',12)
 1030          if(Jcnd(Kast).ne.0) call SUICID(
     .' WRONG TYPE OF PERTURBING ASTEROID I.C., STOP IN SBSET  ',14)
 
c cross-reference with group of tape ast/sat bodies
               if(Nast.gt.0) then
                  do k = 1, Nast
                     if(Np2(k).eq.j) then
                        Nas(k) = Numast
                        goto 1050
                     endif
                  end do
               endif
 
c not on the tape, must get position in fmipt
               if(Kp(j+30).gt.1) then
                  write(Iout, 1040) j, Kp(j+30)
 1040             format(5X, 'BODY', I3,
     .              ' NOT FOUND ON ASTEROID TAPE DESPITE REQUEST, KP =',
     .              I4, ' OVERRULED')
                  Kp(j+30) = 1
                  Kpast(Numast) = 1
               endif
 1050          Nplast(Numast) = Nqlnt(Kast)
               Ncnast(Numast) = Npcent(Kast)
               if(Jdpl0(Kast).gt.0) then
                  Tlpast(Numast) = Jdpl0(Kast)
               else
 
c integration time variable is julian date + 0.5
                  Tlpast(Numast) = Con11(1,Kast+4) + 0.5_10
               endif
               Astmab(Numast) = Mass(j)
               Astmac(Numast) = 1._10
               if(Ncnast(Numast).gt.0) then
                  nc  = Ncnast(Numast)
                  di4 = Mass1(nc)
                  if(nc.eq.Ncentr) di4 = di4 + Mass(nc)*Mass(j)
                  gauss9(Numast) = Gauss*SQRT(di4)
                  Astmas(Numast) = Mass(j)*Mass(nc)
                  Astmac(Numast) = Mass(nc)
               else
                  gauss9(Numast) = Gauss*SQRT(1._10 + Mass(j))
                  Astmas(Numast) = Mass(j)
               endif
               call JNITL(gauss9(Numast), Pcond(1,Kast),
     .                    Ast999(1,Numast), 0, dum)
            endif
         endif
      end do
c*  start=4000
c
c setup for low-thrust forces
      Velflg = .false.
      Posflg = .false.
 
c setup ifleng to control whether sbeng is called by plnorb
      Ifleng = 0
      if(Nplnt.gt.30) then
         if(Kp(81).gt.0 .or. Kp(82).gt.0 .or. Kp(83).gt.0) Ifleng = 1
         Kreflt = Kp(81)/10
         Kdirct = Kp(81) - 10*Kreflt
         if(Kp(83).gt.0 .or. Kdirct.gt.0) then
            if(lcentr.le.0) then
 
c set-up for solprb
               Ltf11 = Con(3)*Con(1)
               Ltf21 = Con(4)*Con(1)
               Ltf31 = Con(5)*Con(1)
               Ltf12 = Con(3)*Con(2)
               Ltf22 = Con(4)*Con(2)
               Ltf32 = Con(5)*Con(2)
c*  start=4100
c
c setup atmospheric drag constants
            else
 
c set-up for plnorb
               deltar    = delta*Convd
               alphar    = alpha*Convd
               Canvec(1) = COS(deltar)*COS(alphar)
               Canvec(2) = COS(deltar)*SIN(alphar)
               Canvec(3) = SIN(deltar)
               if(Kp(83).gt.0) then
 
c setup gas leak quantities
                  Tx    = ta + tb - tc
                  Ty    = -ta + tb - tc
                  Tz    = tb + tc
                  betar = beta*Convd
                  Sbeta = SIN(betar)
                  Cbeta = COS(betar)
               endif
               if(Kdirct.gt.0) then
 
c setup direct radiation pressure quantities
                  beta2  = Con1(5)*Convd
                  Cbeta2 = COS(beta2)
                  Sbeta2 = SIN(beta2)
 
c setup planet radius for shadow model in plnorb
                  Ler  = -1
                  Lpsx = Lps
                  if(Ncentr.eq.3) Pcrad  = Econd(7)
                  if(Ncentr.eq.10) Pcrad = Mcond(7)
 
c setup reflected radiation quantities
                  if(Kcnt.gt.0) Pcrad = Pcond(7,Kcnt)
               endif
            endif
         endif
         if(Kcnt.le.0) then
c
c setup shadow and atmospheric drag for earth satellites
            if(Ncentr.eq.3) then
               Index1 = 6
c
c setup shadow for lunar orbiters
            else if(Ncentr.ne.10) then
            endif
         else if(Kp(82).gt.0) then
            Posflg   = .true.
            Velflg   = .true.
            Omegc    = Twopi/Pcond(13,Kcnt)/86400._10
            Index1   = 6
            Cdrg1    = Pcond(14,Kcnt)/360._10*Twopi
            Cdrg2    = Pcond(15,Kcnt)/360._10*Twopi
            Omega(1) = Omegc*COS(Cdrg1)*COS(Cdrg2)
            Omega(2) = Omegc*SIN(Cdrg1)*COS(Cdrg2)
            Omega(3) = Omegc*SIN(Cdrg2)
            Sh       = Pcond(27,Kcnt)
            Rhoz     = Pcond(28,Kcnt)
            Drgeps   = Con11(9,Kcnt+4)
c cdrg1=con11(10,kcnt+4)    drag coefficints into
c cdrg2=con11(11,kcnt+4)    probe con1 vector
            Cdrg1  = Con1(10)
            Cdrg2  = Con1(11)
            Rhzkm  = Con11(12,Kcnt+4) + Pcond(7,Kcnt)
            it1    = 100*Ncentr + 27
            it2    = it1 + 1
            Atmeps = 1E-18_10
            do i = 1, Kount
               if(Icntrl(i).gt.it2) goto 1100
               if(Icntrl(i).ge.it1) Atmeps = 0._10
 
c if drag partials are integrated always compute drag force
            end do
         endif
      endif
c*  start=4200
c
c setup for target planet quantities
 1100 do i = 1, i_mxtrg
         Itgast(i) = 0
         Masstc(i) = 1._10
         Ktrg(i)   = 0
         Ntrg(i)   = 0
      end do
      Numtar = 0
      do i = 1, i_mxtrg
         Ntrg(i) = Kp(i)
         if(Ntrg(i).le.0) goto 1300
         Numtar   = Numtar + 1
         nt       = Ntrg(i)
         Masst(i) = Mass(nt)
         if(Numast.gt.0) then
            do kt = 1, Numast
               if(nt.eq.Nplast(kt)) then
                  Itgast(i) = kt
                  Masstc(i) = Astmac(kt)
                  Masst(i)  = Astmas(kt)
               endif
            end do
         endif
         if(nt.eq.10) then
            Masstc(i) = Mass(3)
            Masst(i)  = Mass1(10)
         endif
         do j = 1,u_mxpl
            if(Ntrg(i).eq.Nqlnt(j)) then
               Ktrg(i) = j
               goto 1150
            endif
         end do
         if(Ntrg(i).eq.3) then
            Ktrg(i) = -3
         else if(Ntrg(i).eq.10) then
            Ktrg(i) = -2
         else
            call SUICID(' INCORRECT TARGET OR NO TARGET SET UP FOR '//
     .                  'PROBE IN SBSET, STOP IN SBSET ', 18)
         endif
 1150    if(Ltrg(i).ne.0) then
 
c must integrate some i.c. partial
            if(Ktrg(i).ne.-2) then
 
c nt=ntrg(i)
               do kt = 1, 9
                  if(nt.eq.Nplbd(kt)) then
                     call IMITL(Gauss, Mass(nt), Betabd(1,kt), 1, i)
 
c masst(i)= mass(nt)
                     T0mpt(i) = Jdbd0(kt)
                     goto 1200
                  endif
               end do
c perturbing planet not found, look for perturbing elliptic
c orbit
               if(Numast.gt.0) then
                  do kt = 1, Numast
                     if(nt.eq.Nplast(kt)) then
 
c itgast(i)= kt
                        Kast = kastr(kt)
 
c maybe see if it is on ast/sat tape for init. cond.
                        call IMITL(gauss9(kt),0._10,Pcond(1,Kast),1,i)
c masstc(i)= astmac(kt)
c masst(i)= astmas(kt)
                        T0mpt(i) = Tlpast(kt)
                        goto 1200
                     endif
                  end do
               endif
               call SUICID(
     .' TARGET PLANET FOR PARTIALS IS NOT PERTURBING PLANET, STOP IN SBS
     .ET ', 17)
            endif
            call LUNSET(Gauss,Mass(3),Jdbd0(10),0._10,-1)
c masstc(i)= mass(3)
c masst(i)= mass(3) * mass(10)
            T0mpt(i) = Jdbd0(10)
         endif
 1200 end do
c setup for extra print on kout
 1300 if(Kout.gt.0 .and. Intxpr.gt.0) then
         call ZFILL(Relacc,zxprcm)
         Xprnam(1)=Name
         if(Numtar.gt.0) Xprnam(2)=Aplnt(Ktrg(1))
         if(Kcnt.gt.0) then
            Xprnam(3)=Aplnt(Kcnt)
         else if(Ncentr.eq.3) then
            Xprnam(3)=Aplnt(-3)
         else if(Ncentr.eq.10 .or. (Nplnt.eq.3.and.Kp(40).ge.0)) then
            Xprnam(3)=Aplnt(-2)
         else
            Xprnam(3)='  ----  '
         endif
         write(Kout,1310)
 1310    format(' ACCELERATIONS AND OTHER INTERMEDIATE QUANTITIES')
      endif
c
c*  start=4400
c setup for distributed asteroidal perturbing force
      Nbelt = 0
      if(Kp(80).ge.0) then
         Nbelt = Nbelt + 1
         call ASBSET(Nbelt, 50, 49, 48, 47)
      endif
      if(Kp(76).ge.0) then
         Nbelt = Nbelt + 1
c---------------- temporary ? ? ? ? ?
c use prmter(46) for mass of 2nd belt, 73-75 for rad,inc,nod
         call ASBSET(Nbelt, 46, 73, 74, 75)
      endif
      if(Kp(75).ge.0) then
         Nbelt = Nbelt + 1
c---------------- temporary ? ? ? ? ?
c use prmter(45) for mass of 3rd belt, 76-78 for rad,inc,nod
         call ASBSET(Nbelt, 45, 76, 77, 78)
      endif
      if(Nbelt.gt.0) then
         call LEGNDR(0._10,1._10,32,0,Zz(2),Leg45,0._10,0._10)
         Zz(1) = 0._10
         cos45 = SQRT(0.5_10)
         call LEGNDR(cos45,cos45,32,0,Leg45(2),Leg145(2),0._10,0._10)
      endif
c
c set up limited asteroids as perturbing bodies
c the masses are assumed to be negligible in determining the orbits,
c and so the scale may be squeezed to get the mean motion right
      do j = 1, 3
         Sumsml(j) = 0._10
         do iden = 1, 5
            Dsmsml(j,iden) = 0._10
         end do
      end do
      if(Kp(30).ge.0 .and. Numsml.gt.0) then
         if(lcentr.le.0) then
            do i = 1, Numsml
               den = 2.0
               if(Denpts(i).gt.0) den = Prmter(Denpts(i))
               Smlvol(i)=Scond(7,i)**3*(4._10/3._10)*Pi*1E15_10/2E33_10
               Smlmas(i)=Smlvol(i)*den
               call JNITL(Gauss,Scond(1,i),Elptsm(1,i),0,dum)
            end do
         endif
      endif
c
c Set up cometary nongravitational forces
      if(Nplnt.le.30 .and. lcentr.le.0 .and. Kp(83).ge.0) then
         Nga1 = Con(14)
         Nga2 = Con(15)
         Ngm =  2.15_10
         Ngn =  5.093_10
         Ngk =  4.6142_10
         Ngalph=0.1113E-8_10
         Ngr0 = 2.808_10
      endif
c
c set up selectors for PBCOR
      do l=1,10
         Kpb(l)=Kp(l+30).ge.0
      end do
      do j=1,Numtar
         nt=Ntrg(j)
         if(nt.gt.0 .and. nt.le.10 .and. Kp(nt+30).lt.0) call SUICID(
     .    'TARGET PLANET NOT INCLUDED IN MOTION, STOP SBSET',12)
      end do
      if(Ncentr.gt.0 .and. Ncentr.le.10) Kpb(Ncentr)=.false.
      if(Nplnt.gt.0 .and. Nplnt.le.10) Kpb(Nplnt)=.false.
      if(Nplnt.eq.3) Kpb(10)=.false.
c
c*  start=5000
c
c set up central body harmonics
      Nczone = 0
      Nctess = 0
      Nczonp = 1
      Nctesp = 0
      Ntopc  = 0
      Ntopct = 0
      Ntpctp = 0
      Jczone = 0
      Jccos  = 0
      Jcsin  = 0
      Icrot  = 0
      if(Jct(78).gt.0) call LEGSET

      if(Ncentr.gt.0) then
         Jczone = Ncentr*100 + 31
         Jccos  = Ncentr*100 + 41
         Jcsin  = Ncentr*100 + 51
c
c earth is central body
         if(Ncentr.eq.3) then
            Nczone = Nezone
            Nctess = Netess
            kkcnt  = -3
c
c moon is central body
         else if(Ncentr.eq.10) then
            Nczone = Nmzone
            Nctess = Nmtess
            Nctesp = 1
            kkcnt  = -2
         else
c
c planet is central body
            do i = 1, Nmphar
               if(Ncentr.eq.Nplhar(i)) then
                  khar = i
                  if(Kcnt.le.0) call SUICID(
     .   ' CENTRAL BODY IS NOT INPUT PLANET, STOP IN SBSET', 12)
 
c*  start=5200
                  kkcnt = Kcnt
c set-up for partials w/r/t pole location and rotation rate.
c partials w/r/t con's 6,7,8,9 of central body
c icrot=0        done above
                  nrckm = 100*Ncentr + 12
                  nrckp = nrckm + 3
                  do k = 1, Kount
                    if(Icntrl(k).le.nrckp .and. Icntrl(k).ge.nrckm) then
                       Icrot = Icrot + 1
                       Icrotp(Icrot,1) = k
                       Icrotp(Icrot,2) = Icntrl(k) - nrckm + 1
                    endif
c icrotp( ,1)=i tells which partial is for rotation quan.
c icrotp( ,2)=j con(j+5) partial, j=1,4
                 end do
                  Nczone = Npzone(khar)
                  Nctess = Nptess(khar)
                  goto 1350
               endif
            end do
            Nczonp = 0
            goto 1400
         endif
 1350    Crad   = Pcond(7,kkcnt)/Aultsc/Ltvel
         Nczon1 = Nczone - 1
         if(Nczon1.gt.0) then
            do i = 1, Nczon1
               if(Ncentr.eq.3) then
                  har = Ezhar(i)
               else if(Ncentr.eq.10) then
                  har = Mzhar(i)
               else
                  har = Pzhar(khar,i)
               endif
               Czhar(i) = har
            end do
         endif
         Nctes1 = Nctess - 1
         if(Nctes1.gt.0) then
            k = 0
            do j = 2, Nctess
               do i = 1, j
                  k = k + 1
                  if(Ncentr.eq.3) then
                     har  = Echar(k)
                     hars = Eshar(k)
                  else if(Ncentr.eq.10) then
                     har  = Mchar(k)
                     hars = Mshar(k)
                  else
                     har  = Pchar(khar,k)
                     hars = Pshar(khar,k)
                  endif
                  if(Jct(78).gt.0) then
                     Cchar(k) = har
                     Cshar(k) = hars
                  else
                     Cchar(k) = har*LEGSCL(j, i)
                     Cshar(k) = hars*LEGSCL(j, i)
                  endif
               end do
            end do
         endif
c*  start=5300
c earth, moon rotation are special
         if(kkcnt.ge.1) then
            Alphc(1) = (90._10 + Pcond(15,Kcnt))*Convd
            Alphc(2) = 0._10
            Salphc   = SIN(Alphc(1))
            Calphc   = COS(Alphc(1))
            Deltc(1) = (90._10 - Pcond(14,Kcnt))*Convd
            Deltc(2) = 0._10
            Sdeltc   = SIN(Deltc(1))
            Cdeltc   = COS(Deltc(1))
            if(Nctess.gt.0) then
               Omegc = 0._10
               if(Pcond(13,Kcnt).ne.0._10) Omegc = Twopi/Pcond(13,Kcnt)
               Psic0  = Pcond(12,Kcnt)*Convd
               Epochc = Con11(1,Kcnt+4) + 0.5_10
            endif
         endif
c        add .5 to julian date of epoch because integration
c        variable s is julian date + .5
c*  start=5400
c
c        no.of central body harmonics included in partials integ.
         Ntopc  = max0(Nczon1, Nctes1)
         Ntopct = min0(3, Nczon1)
         if(Numtar.le.0) Ntopct = 0
         if(Kp(86).ge.0) then
            ktop   = Kp(86)
            ktopa  = ktop/100
            Nctesp = min0(ktopa-1, Nctes1)
            ktop   = ktop - ktopa*100 - 1
            Nczonp = min0(ktop, Nczon1)
            Ntpctp = min0(Ntopct, Nczonp)
         endif
      endif
c
c*  start=6000
c
c set up target body harmonics
 1400 dhsi = Kkp(85)
      dhsi = ABS(dhsi)/100._10
      do ll = 1, i_mxtrg
         do i = 2, 5
            Tchar(ll,i) = 0._10
            Tshar(ll,i) = 0._10
         end do
         Tzhar(ll,2) = 0._10
         Hsitb(ll)   = 0._10
         Ntzone(ll)  = 0
         Nttess(ll)  = 0
         Ntzonp(ll)  = 0
         Nttesp(ll)  = 0
         Jtzone(ll)  = 0
         Jtcos(ll)   = 0
         Jtsin(ll)   = 0
         Ntopt(ll)   = 0
         if(Ntrg(ll).gt.0) then
            Jtzone(ll) = Ntrg(ll)*100 + 31
            Jtcos(ll)  = Ntrg(ll)*100 + 41
            Jtsin(ll)  = Ntrg(ll)*100 + 51
c
c earth is target body
            if(Ntrg(ll).eq.3) then
               sma   = 1._10
               npz   = Nezone
               npt   = Netess
               kktrg = -3
            else if(Ntrg(ll).ne.10) then
c*  start=6300
c
c planet is target body
               do i = 1, Nmphar
                  if(Ntrg(ll).eq.Nplhar(i)) then
                     Kthar(ll) = i
                     kktrg     = Ktrg(ll)
                     sma = Pcond(1,kktrg)
                     mum = Kthar(ll)
                     npz = Npzone(mum)
                     npt = Nptess(mum)
                     goto 1420
                  endif
               end do
               goto 1500
c
c moon is the target body
            else
c harmonic effect of moon not implemented for EMbary integration
               if(Nplnt.eq.3) goto 1500
               sma   = 1._10
               npz   = Nmzone
               npt   = Nmtess
               kktrg = -2
            endif
 1420       Trad(ll)   = Pcond(7,kktrg)/Aultsc/Ltvel
            Ntzone(ll) = min0(npz, 5)
            Ntzon1(ll) = Ntzone(ll) - 1
            ntz = Ntzon1(ll)
            if(Ntzon1(ll).gt.0) then
               do i = 1, ntz
                  if(kktrg.eq.-3) then
                     har = Ezhar(i)
                  else if(kktrg.eq.-2) then
                     har = Mzhar(i)
                  else
                     har = Pzhar(mum,i)
                  endif
                  Tzhar(ll,i) = har
               end do
            endif
            Nttess(ll) = min0(npt, 3)
            Nttes2(ll) = (Nttess(ll)*(Nttess(ll)+1))/2 - 1
            Nttes1(ll) = Nttess(ll) - 1
            ntt = Nttess(ll)
            if(Nttes1(ll).gt.0) then
               k = 0
               do j = 2, ntt
                  do i = 1, j
                     k = k + 1
                     if(kktrg.eq.-3) then
                        har  = Echar(k)
                        hars = Eshar(k)
                     else if(kktrg.eq.-2) then
                        har  = Mchar(k)
                        hars = Mshar(k)
                     else
                        har  = Pchar(mum,k)
                        hars = Pshar(mum,k)
                     endif
                     if(Jct(78).gt.0) then
                        Tchar(ll,k) = har
                        Tshar(ll,k) = hars
                     else
                        Tchar(ll,k) = har*LEGSCL(j, i)
                        Tshar(ll,k) = hars*LEGSCL(j, i)
                     endif
                  end do
               end do
            endif
 
c earth, moon rotation are special
            if(kktrg.ge.1) then
               Alpht(ll,1) = (90._10 + Pcond(15,kktrg))*Convd
               Alpht(ll,2) = 0._10
               Salpht(ll)  = SIN(Alpht(ll,1))
               Calpht(ll)  = COS(Alpht(ll,1))
               Deltt(ll,1) = (90._10 - Pcond(14,kktrg))*Convd
               Deltt(ll,2) = 0._10
               Sdeltt(ll)  = SIN(Deltt(ll,1))
               Cdeltt(ll)  = COS(Deltt(ll,1))
               if(Nttess(ll).gt.0) then
                  Omegt(ll) = 0._10
                  if(Pcond(13,kktrg).ne.0._10) Omegt(ll) =
     .                        Twopi/Pcond(13,kktrg)
                  Psit0(ll)  = Pcond(12,kktrg)*Convd
                  Epocht(ll) = Con11(1,kktrg+4) + 0.5_10
               endif
            endif
c        add .5 to julian date of epoch because integration
c        variable s is julian date + .5
c*  start=6500
c
c        no.of target body harmonics included in partials integration
            Ntopt(ll) = max0(Ntzon1(ll), Nttes1(ll))
 
c find target body range of harmonic effect
            nt = Ntrg(ll)
            if(nt.le.30) then
              if(Ntzone(ll).gt.2 .or. Nttess(ll).ge.2) then
                 hsib = sma*sma*Masst(ll)*Trad(ll)*Trad(ll)
                 if(hsib.gt.0._10) then
                   hsia =hsib*MAX(ABS(Tchar(ll,2)),ABS(Tshar(ll,2)))
                   hsia =(hsia*hsftr)**(1.0/4.0)
                   hsib =hsib*Trad(ll)
     .                    *MAX(ABS(Tchar(ll,3)), ABS(Tchar(ll,4)),
     .                    ABS(Tchar(ll,5)), ABS(Tshar(ll,3)),
     .                    ABS(Tshar(ll,4)), ABS(Tshar(ll,5)),
     .                    ABS(Tzhar(ll,2)))
                   hsib =(hsib*hsftr)**(1.0/5.0)
                   Hsitb(ll) = MAX(hsia, hsib) + dhsi
                 endif
              endif
            endif
            if(Kp(85).ge.0) then
               ktop  = Kp(85)
               ktopa = ktop/100
               Nttesp(ll) = min0(ktopa-1, Nttes1(ll))
               ktop = ktop - ktopa*100 - 1
               Ntzonp(ll) = min0(ktop, Ntzon1(ll))
            endif
         endif
 
 1500 end do
c
c setup central body mascons
      call MSCSET(Ncentr)

c
c les-8/9 thrust and attitude history setup
      Les    = Kkp(91)
      Nrecls = 0
      Jendls = 0
c
c les-8/9 station keeping setup
      Azold = -1E10_10
      do i = 1, 6
         Itrns(i) = 0
      end do
      Irout  = 0
      Ainc   = Con1(2)*Convd
      Aomeg  = Con1(3)*Convd
      Sslong = Con1(4) + Con1(9)
c east longitude sought is corrected for discrepency between
c longitude derived from sun transits and actual longitude of
c center of les-8/9 figure eight ground track
      i = Sslong/360._10
      if(Sslong.lt.0._10) i = i - 1
      Sslong = Sslong - i*360._10
      Damp   = Con1(5)
      Coast  = Con1(6)
      Satdr  = Con1(7)
      Thlev  = Con1(8)*1E3_10
c
c difference between ephemeris and universal time
      i = Jdp0 - 2443145
      j = i/365
      if(i.lt.0) j = j - 1
      etutc = (48.184_10 + j)
c etutc=48.184 in 1977
c
c fire1 =(39-sslong/15)/24
c kkp(71-73)=utc hour,minute,second of station keeping firing time
      Fire1 = (Kkp(71)*3600._10+Kkp(72)*60._10+Kkp(73)+etutc)/864.E2_10
      fire0 = Fire1 - fcp0
      i     = fire0/Hmx
      Fire1 = fcp0 + i*Hmx
      i     = Fire1
      if(Fire1.lt.0._10) i = i - 1
      Fire1 = Fire1 - i
      Fire2 = Fire1 + 0.5_10
      if(Fire2.gt.1._10) Fire2 = Fire2 - 1._10
      if(Fire2.le.Fire1) then
         fire0 = Fire1
         Fire1 = Fire2
         Fire2 = fire0
      endif
      Lthrst = 0
      Kthrst = 0
      Jthrst = 0
      T0sav  = T0
c initialize central body tape reading
c kkp(84)=0  ignore central body tape during integration
c kkp(84)=1  read central body tape for positions, velocities only
c kkp(84)=2  read central body tape for partials as well
      lcplp = 0
 
c somehow suppress printout duplication
      call PLPRD1(lcplp)
c*  start=9000
c
c print out data and control constants for satellite
c probe integration
      call PLNHED
 
      if(Ncentr.gt.0) then
         ntaaa = Nctesp + 1
         nzaaa = Nczonp + 1
         if(Nctess.eq.0) ntaaa = 0
         if(Nczone.eq.0) nzaaa = 0
         call PAGCHK(60,3,0)
         write(Iout,1550) Ncentr,Kkp(1),Nczone,Nctess,nzaaa,ntaaa
 1550    format('- NCENTR =',I3,'  LCENTR =',I3,'   NCZONE   =',I3,
     .        '   NCTESS   =',I3,'   FOR PARTIALS ONLY:   NCZONE   =',
     .        I3,'   NCTESS   =',I3)
         if(Numtar.gt.0) then
            nzaaa = Ntopct + 1
            ntaaa = Ntpctp + 1
            if(nzaaa.le.1) nzaaa = 0
            if(ntaaa.le.1) ntaaa = 0
            call PAGCHK(60,1,0)
            write(Iout,1560) nzaaa,ntaaa
 1560       format(
     .           ' TARGET BODY EFFECT ON CENTER HARMONICS,  NCZONE-T =',
     .           I3,'   FOR PARTIALS ONLY:   NCZONE-T =',I3)
         endif
         call PAGCHK(60,1,0)
         write(Iout,1600) Ncentr,Mass1(Ncentr),Goose
 1600    format('  MASS1(',I2,')=',1PD22.15,5X,'SQRT(GM)=',
     .          D22.15)
      endif
      do i = 1,Numtar
         if(Ntrg(i).gt.0) then
            nzaaa = Ntzonp(i) + 1
            ntaaa = Nttesp(i) + 1
            if(Nttess(i).eq.0) ntaaa = 0
            if(Ntzone(i).eq.0) nzaaa = 0
            call PAGCHK(60,1,0)
            write(Iout,1700) i,Ntrg(i),i,Ntzone(i),i,Nttess(i),
     .                        i,Hsitb(i),i,nzaaa,i,ntaaa
         endif
      end do
 1700 format('  NTRG(',I1,')=',I3,'   NTZONE(',I1,')=',I3,
     .       '   NTTESS(',I1,')=',I3,'  HSITB(',I1,')=',F10.6,
     .       '   FOR PARTIALS ONLY:   NTZONE(',I1,')=',I3,
     .       '   NTTESS(',I1,')=',I3)
      if(lcentr.gt.0 .and. lcentr.ne.3 .and. lcentr.ne.10) then
         if(Nplnt.gt.30) then
 
c write out for low thrust forces in plnorb
            i = 2
            if(Kp(81).gt.0) i = 1
            write(Iout,1720) mesage(i)
 1720 format(10X,'FOR LOW THRUST FORCES:  RADIATION PRESSURE IS',A4)
            i = 2
            if(Kp(82).ne.0) i = 1
            write(Iout,1740) mesage(i)
 1740       format(44X,'AIR DRAG IS',A4)
            i = 2
            if(Kp(83) .gt. 0) i = 1
            call PAGCHK(60,1,0)
            write(Iout,1760) mesage(i)
 1760       format(42X,'GAS LEAKS ARE',A4)
         endif
      endif
c*  start=9900
c ensure end of page here
      if(Line.lt.60) Line = 60
      return
      end
