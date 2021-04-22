      subroutine RADPSR(nvel,norm)
 
      implicit none
c
c           subr radpsr - j.f.chandler - 1984 apr
c           compute pulsar pulse phase observables
c     in this routine 'tmdly' is used for the delay from the arrival
c     (free-space) at the solar-system barycenter to the observation
c     and may be either pos. or neg.  however, the eventually computed
c     observable is pulse phase.

c arguments
      integer*4 nvel,norm

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 dutrec
      equivalence (dutrec,Angdum(10))
      include 'difnct.inc'
      include 'ennips.inc'
      include 'eqnphs.inc'
c     xpp - unit vector in alpha direction from pulsar
c     xrho- unit vector in delta direction from pulsar
c     drp(1-3) - vector proper motion (rad/d)
c     drp(4-6) - partial of unit vector to pulsar w.r.t. ddelta
c     drp(7-9) - partial of unit vector to pulsar w.r.t. dalpha
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'ltrapx.inc'
      real*10 tmdly
      equivalence (Deriv(2,1),tmdly)
      include 'mnsprt.inc'
      include 'number.inc'
      include 'obscrd.inc'
      real*10 Ctrecf,Phobs,Freqis
      equivalence (Ctrecf,Dstf),(Phobs,Dstf(6)),(Freqis,Dstf(7))
      include 'kobequiv.inc'
      real*4 acctim
      equivalence (acctim, Estf)
      include 'param.inc'
      include 'prpgat.inc'
      include 'rtrdvl.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'
c
c local
      real*10 ca,ddc2,dinv,disp,dpsrdc,dpsrra,dra2,dum,
     . dxst(3),fbc,frprt,psrrelp,sa,sbcvel,sclf,stmdly,t,tdprlx,
     . tfct,tfct1,tfct2,tmdlyb,tmdly3,tp,xedus
      real*10 vcentr(3),
     . vc1950(3)/-0.882706_10,-13.488415_10, 5.855212_10/,
     . vc2000(3)/-0.760127_10,-13.498004_10, 5.850308_10/
      integer*4 i,j,jdbc,jdprt,loop1,nn,np,np1
c external functions
      real*10 ANCOR,DOT
c
c
c-----------initialization------------------------------------------
c*  start=1000
c
      if(Nice.ge.0) call SUICID(
     .      'CAN''T PROCESS PULSE PHASE RATE OBSERVABLE, STOP RADPSR  '
     .      , 14)
      if(Nplnt0.ge.0) call SUICID(
     .  'CAN''T PROCESS PULSE PHASE OBSERVABLE FOR PLANET, STOP RADPSR'
     .  , 15)
      if(nvel.le.0) call SUICID(
     .                       'MUST HAVE VELOCITIES FOR PULSAR PHASE   '
     .                       , 10)
 
      if(Nk1.lt.0) then
         dinv = Psrprm(1)*Convds/Aultsc
         D    = 1E30_10
         if(Psrprm(1).gt.0._10) D = 1._10/dinv
         sclf   = Convds/365.25_10
         dpsrra = Psrprm(2)*sclf
         dpsrdc = Psrprm(3)*sclf
         if(Dydphi(3,1).eq.0) then
            ca = -Dydphi(1,1)/Yspcd(3,1)
            sa = -Dydphi(2,1)/Yspcd(3,1)
         else
            ca = Yspcd(1,1)/Dydphi(3,1)
            sa = Yspcd(2,1)/Dydphi(3,1)
         endif
         Xpp(1) = -sa
         Xpp(2) = ca
         Xpp(3) = 0._10
         do i = 1, 3
            Xrho(i) = Dydphi(i,1)
            Drp(i)  = dpsrra*Xpp(i) + dpsrdc*Xrho(i)
            if(Jct(13).eq.0) then
               vcentr(i)=vc1950(i)/Ltvel
            else
               vcentr(i)=vc2000(i)/Ltvel
            endif
         end do
c ignore sun's peculiar velocity for now
c        sbcvel = DOT(vcentr,Yspcd)
        sbcvel = 0._10
      endif
c
c*  start=1200
c
c determine star coordinates in standard reference frame
c note: xspcd is not a unit vector here
      t = (Jd - Jdps0) + Ctrecf
      do i = 1, 3
         Xspcd(i,1) = Yspcd(i,1) + t*Drp(i)
      end do
c
c*  start=2200
c
      Nit(20) = Nit(20) + 1
 
c set up xsitep, etc. (do both, just in case)
      do j = 1, 2
         do i = 1, 3
            Xsitep(i,j) = -Xspcd(i,1)
 
c ignore pulsar velocity here
            Xsitep(i + 3,j) = Xsite(i + 3,1)
         end do
         call UVECTR(3,Xsitep(1,j),Rsitp(j),Xsitp0(1,j),dum)
      end do
c
c test for dummy observation below the horizon
      if(Idumob.eq.1) then
         nn = 1
         call HORIZN(Xsite,Xsitep,nn,Jd)
         if(Jd.le.0) return
      endif
c
c*  start=2500
c compute delay w.r.t. solar-system barycenter
      Tmdly1 = DOT(Xemlsc,Xsitp0)
      do i = 1, 3
         Xessbc(i,1)   = Xemlsc(i,1) - Xslcns(i,1)*Aultsc
         Xessbc(i+3,1) = Xemlsc(i+3,1)
      end do
      xedus = DOT(Xessbc,Xsitp0)
 
c correct for parallax
      do i = 1, 3
         dxst(i) = Xessbc(i,1) - xedus*Xsitp0(i,1)
      end do
      Tmdly2 = DOT(Xessbc,dxst)
      tdprlx = Tmdly2*dinv*0.5_10
      tmdly  = xedus + tdprlx
c
c compute interstellar transmission frequency and dispersion
      if(Nvlcns.gt.0) then
         do i = 4, 6
            Xessbc(i,1) = Xessbc(i,1) - Xslcns(i,1)*Aultvl
         end do
      endif
      freqis = Freq*(1._10 - sbcvel + DOT(Xessbc(4,1),Xsitp0))
      Freqi2 = 1._10/freqis**2
      Freqi2t= Freqi2*t
      disp   = (Psrprm(4)+t*Psrprm(9))*Freqi2
      tmdly  = tmdly + disp
c
c correction to delay for general relativity
      Rs1m = SQRT(DOT(Xemlsc,Xemlsc))
      Rs2m = Rs1m
      Raddum(1) = -2._10*LOG((Rs1m - Tmdly1)/Aultsc)
c add planet contributions
c
c
      tmdly = tmdly + Raddum(1)*Reltrm(1)
      if(MOD(Jct(6)/2048,2).eq.1) then
         if(Line.gt.56) call OBSPAG
         write(Iout,100) Jd,Dstf(1),Raddum(1)*Reltrm(1),disp
  100    format(1x,'PSRREL: JD.F=',i7,f13.12,' DT=',1pd22.15,' DISP=',
     .    d22.15)
         Line = Line + 1
      endif
c
c iterate position of pulsar companion, if any
c extra delay is conceived as pulsar "topography" passed in dstf(8)
c to be subtracted from calculated delay
c NOTE: if the orbit is not defined, but the delay is given in dstf(8),
c that pre-computed delay is applied
      if(Psrprm(10).gt.0._10 .or. Nmpex.gt.0) then
         tp = (Jd-Jdps0) + (Ctrecf - tmdly/Secday)
         np1=1
         if(Psrprm(10).gt.0._10) np1=0
         tmdlyb = 0._10
         tmdly3 = acctim+1.
         loop1  = 0
         do while(ABS(tmdly3-tmdlyb).gt.acctim)
c
            tmdly3 = tmdlyb
            Frect  = tp - tmdly3/Secday
            jdprt  = Jdps0+Frect
            frprt  = MOD(Frect,1._10)
            if(frprt.lt.0._10) frprt=frprt+1._10
            tmdlyb = 0._10
            do np=np1,Nmpex
               call JLIPTP(Frect,Elptn(1,np),7,Ynpt(1,np),Rynpt(np),
     .          Rynpt2(np),Rynpt3(np),Dynpt(1,1,np))
               call TRPLST(jdprt,frprt,0,np,'JLIPT(PSR)',Ynpt(1,np))
               tmdlyb = tmdlyb + Ynpt(3,np)*Aultsc*Elptn(40,np)
            end do
            loop1  = loop1 + 1
            if(loop1.gt.10)
     .       call SUICID('MORE THAN 10 BINARY PULSAR ITERATIONS   ',10)
         end do
         if(np1.eq.0 .and. Psrprm(15).gt.0._10) then
            Psrrel=-2._10*LOG((Rynpt(0)-Ynpt(3,0))/Elptn(8,0))
            psrrelp=Psrrel*Psrprm(16)
            tmdlyb = tmdlyb + psrrelp*Reltrm(1)
            Raddum(1) = Raddum(1)+psrrelp
            if(MOD(Jct(6)/2048,2).eq.1) then
               if(Line.gt.56) call OBSPAG
               write(Iout,110) psrrel*Reltrm(1),-tmdlyb
  110          format(1x,'PSRREL(BIN): DT=',1pe22.15,' TORB=',e22.15)
               Line = Line + 1
            endif
         endif
         Dstf(8)=-tmdlyb
         Nit(20) = Nit(20) - 1 + loop1
c set up for partials
         Opv3=1._10
         do np=np1,Nmpex
            Opv3=Opv3+Ynpt(6,np)*Aultvl
         end do
      endif
      tmdly=tmdly-Dstf(8)
c
c copy pulsar position into xsbsun
      do i = 1, 3
         Xsbsun(i) = -D*Xsitp0(i,1)
      end do
c
c propagation corrections
      call PROPCO(1,1)
      tmdly = tmdly + Sumcor(1)
c
c antenna offset
      tmdly = tmdly + ANCOR(1)
c
c biases
      if(Nrbias.gt.0) tmdly = tmdly + Rbsx(1)
      if(Neqnox.gt.0) tmdly = tmdly + Pnox +
     .                            dutrec*(Pequat + Plat*dutrec/2._10)
c
c calculate star coordinate partials
      call SPOTCV(0,Nplnt0,1)
c
c convert coordinate time delay to pulse phase
      stmdly = tmdly
      call PSRPHS(tmdly)
c
c extra print
      if(Ict(2).le.0) then
         if(Line.gt.57) call OBSPAG
         call TIMINC(Jd,ctrecf,jdbc,fbc,-stmdly/Secday - .5_10)
         jdbc = mod(jdbc,10000)
         write(Iout,50) jdbc,fbc,stmdly,tdprlx,freqis
   50    format(49x,'TBRY=', i4, f14.13, 22x, 1pd16.9, d10.3, d13.6)
         Line = Line + 1
      endif
 
      if(.not. (Ict(1).le.0 .or. (Idumob.eq.1 .and. Ict(3).le.0
     .    ))) then
c
c setup for partials
         do i = 1, 10
            Raddum(i) = Raddum(i)*Dpphdt
         end do
         Dptdsp = -Dpphdt*(1._10 - xedus*dinv)*Aultsc
         tfct1  = t*Rsitp(1)
 
c scale back to yr/arcsec (and compensate for au)
         tfct  = tfct1*sclf/Aultsc
         tfct2 = tfct1*tfct
         dra2  = dpsrra*tfct2
         ddc2  = dpsrdc*tfct2
 
         do i = 1, 3
            Xsitp0(i,1) = Xsitp0(i,1) + dxst(i)*dinv
            Drp(i + 3)   = Xpp(i)*tfct - Xspcd(i,1)*dra2
            Drp(i + 6)   = Xrho(i)*tfct - Xspcd(i,1)*ddc2
         end do
      endif
c
c*  start=9990
      return
      end
