      subroutine PSRPAR(kick)
 
      implicit none
 
c psrpar - j.f.chandler - 1984 jul
c compute partials w.r.t. pulsar parameters

c parameter
      integer*4 kick

c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'comdateq.inc'
      include 'empcnd.inc'
      include 'ennips.inc'
      include 'eqnphs.inc'
      include 'funcon.inc'
      include 'ltrapx.inc'
      integer*2 kind
      equivalence (kind,Numpar)
      include 'mtrapx.inc'
      include 'number.inc'
      include 'obscrd.inc'
      real*10 ctrecf,freqis
      equivalence (ctrecf,Dstf(1)), (freqis,Dstf(7))
      include 'param.inc'
      include 'prtpr9.inc'
      include 'psrcom.inc'
      include 'rtrdvl.inc'
      include 'yvect.inc'

c external functions
      real*10 DOT
c local
      real*10 cond(9),goose,mp3mc,pspfct(10:16),psrrx,sinc,trm,trmr,
     . val,combo(3),adtfos(3),dxdp(3)
      integer*4 iflag,i,iswtch,l,ll,m,mm,mt
c index for pulsar planet with parameters in PSRCN
      integer*4 n10/0/
c
c partials for exo-planet parameters
      mm=1
      do ll=1,Nmpex
c
c see if exo-planet is on obs.lib tape
         mt = 0
         if(Iabs1.gt.0 .and. mm.le.Mmpex .and.
     .    Mplex(mm).eq.Nplex(ll)) mt = 1
c
c start of loop for exo-planet ic partials
         l = 0
         call PCOPS(7, 'PEX ', mt)
         do while(l.le.6 .or. iflag.le.0)
            if(l.lt.6) then
               iflag = 0
               call PCOPY(l, 6, iflag, 1, Lpex(1,ll), Mpex(1,mm))
               if(iflag.le.0) then
                  if(kick.eq.1) then
                     if(l.eq.1) then
                        trm = Ynpt(3,ll)/Pcond(1,Klanex(ll))
                     else
                        trm = Dynpt(3,l,ll)
                     endif
                     Deriv(kind,1) = Dpphdt*trm*Aultsc/Opv3
                  else
c assume pulsar astrometry is too coarse to depend on companion orbit
                     Deriv(kind,1) = 0._10
                     Deriv(kind,2) = 0._10
                  endif
               endif
            else
               iflag=1
               call PCOPY(l, 30, iflag, 1, Lpex(1,ll), Mpex(1,mm))
               if(iflag.le.0) then
                  if(kick.eq.1) then
                     if(Lpex(l,ll).eq.3) then
                        trm = -Dynpt(3,6,ll)*Frect*Twopi/
     .                   Pcond(9,Klanex(ll))**2
                     else
                        call SUICID('UNDEFINED EXO-PLANET PARAMETER  ',
     .                   8)
                     endif
                     Deriv(kind,1) = Dpphdt*trm*Aultsc/Opv3
                  else
c assume pulsar astrometry is too coarse to depend on companion orbit
                     Deriv(kind,1) = 0._10
                     Deriv(kind,2) = 0._10
                  endif
               endif
            endif
         end do
      end do
      do i=1,3
         combo(i)=Ynpt(i,n10)-Elptn(36,n10)*Dynpt(i,1,n10)
         adtfos(i)=Elptn(37,n10)*Elptn(7,n10)*Ynpt(i,n10)/Elptn(8,n10)
         dxdp(i)=0._10
      end do
      psrrx=-2._10*Psrprm(16)*Reltrm(1)/(Rynpt(n10)-Ynpt(3,n10))
      m = 1
      l = 0
      call PCOPS(m, 'PSRX', iabs1)
      do while( .true. )
c
c copy/merge partials from input obslib tape
         iflag  = 1
         iswtch = 0
         call PCOPY(l, u_nmpsr, iflag, iswtch, Lpsrx, Mpsrx)
         if(iflag.gt.0) goto 999
         ll = Lpsrx(l)
 
c clear partial derivatives, in case
         Deriv(kind, 1) = 0._10
         Deriv(kind, 2) = 0._10
c
c ntype=0: model has p, p-dot, p-dot-dot
         if((ll.ge.10 .and. Psrprm(10).le.0._10) .or. ll.gt.18) call
     .    SUICID('CANNOT CALCULATE PARTIAL W.R.T. UNDEFINED PULSAR PARAM
     .ETER, STOP IN PSRPAR  ', 19)
c
         if(kick.eq.1) then
            trmr=0._10
            if(ll.eq.1) then
c parallax
               Deriv(kind,1) = 0.5_10*Tmdly2*Dpphdt*Convds/Aultsc
 
            else if(ll.eq.2) then
c right ascension motion
               Deriv(kind,1) = Dptdsp*DOT(Drp(4), Xessbc)

            else if(ll.eq.3) then
c declination motion
               Deriv(kind,1) = Dptdsp*DOT(Drp(7), Xessbc)
 
            else if(ll.eq.4) then
c dispersion
               Deriv(kind,1) = Dpphdt*Freqi2
            else if(ll.eq.5) then
 
c initial phase
               Deriv(kind,1) = 1._10
            else if(ll.eq.6) then
 
c pulse period
               trm = Psf0*Tpsr(1) + 2.0*(Psf1*Tpsr(2) + 
     .          (Psf2+Psf1**2*Psp0)*Tpsr(3))
 
            else if(ll.eq.7) then
c pulse period rate
               trm = Psf0*Tpsr(2) + 4.0*Psf1*Tpsr(3)
 
            else if(ll.eq.8) then
c pulse period acceleration
               trm = Psf0*Tpsr(3)

            else if(ll.eq.9) then
c dispersion time derivative
               Deriv(kind,1) = Dpphdt*Freqi2t
 
            else if(ll.eq.10) then
c pulsar companion semimajor axis (assumed non-secular)
               trm = Ynpt(3,n10)/Psrprm(10)
               if(Psmass.gt.0._10) trm = trm-adtfos(3)*Psmass*pspfct(10)
 
            else if(ll.eq.11) then
c pulsar companion eccentricity
               trm = Dynpt(3,2,n10)
               trmr=1._10
               do i=1,3
                  dxdp(i)=Dynpt(i,2,n10)
               end do
 
            else if(ll.eq.12) then
c pulsar companion periapse
               trm = Dynpt(3,5,n10)*Convd
               trmr=Convd
               do i=1,3
                  dxdp(i)=Dynpt(i,5,n10)
               end do
 
            else if(ll.eq.13) then
c pulsar companion initial mean anomaly
               trm = Dynpt(3,6,n10)*Convd
               trmr=Convd
               do i=1,3
                  dxdp(i)=Dynpt(i,6,n10)
               end do

            else if(ll.eq.14) then
c pulsar companion period
               if(Psmass.eq.0._10) then
                  trm = -Dynpt(3,6,n10)*Elptn(7,n10)*Twopi/Psrprm(14)**2
               else
                  trm=(combo(3)-adtfos(3)*Psmass/mp3mc)*pspfct(14)
                  trmr=pspfct(14)
                  do i=1,3
                     dxdp(i)=combo(i)-adtfos(i)*Psmass/mp3mc
                  end do
               endif

            else if(ll.eq.15) then
c pulsar companion inclination if specified
               trm = (-Ynpt(3,n10)*mp3mc+adtfos(3)*Psmass)*pspfct(15) +
     .          Dynpt(3,3,n10)*Convd
               trmr=Convd
               do i=1,3
                  dxdp(i)=Dynpt(i,3,n10)
               end do

            else if(ll.eq.16) then
c pulsar companion mass
               trm = -adtfos(3)/mp3mc

            else if(ll.eq.17) then
c pulsar companion signature time derivative
               trm = Elptn(7,n10)*Ynpt(3,n10)/Elptn(8,n10)/
     .          (Elptn(40,n10)*Elptn(32,n10)*Aultsc)

            else if(ll.eq.18) then
c pulsar companion periapse time derivative
               trm = Elptn(7,n10)*Dynpt(3,5,n10)*Convd

            endif
            if(ll.ge.6 .and. ll.le.8) then
               Deriv(kind,1) = -trm*Psf
            else if(ll.eq.16) then
               Deriv(kind,1) = Dpphdt*(trm*Aultsc*Elptn(40,n10) +
     .          Psrrel*Reltrm(1))
            else if(ll.ge.10 .and. ll.le.18) then
               Deriv(kind,1) = Dpphdt*(trm*Aultsc*Elptn(40,n10)
     .          + trmr*psrrx*(DOT(dxdp,Ynpt(1,n10))/Rynpt(n10)-dxdp(3)))
     .          /Opv3
            endif

         else if (kick.eq.2) then
 
            if(ll.eq.2) then
c right ascension motion
               Deriv(kind,1) = Stime/365.25_10/15._10/Cdeltc

            else if(ll.eq.3) then
c declination motion
               Deriv(kind,2) = Stime/365.25_10
            endif
         endif
c
c
c end of loop
c*  start=1500
         if(l.ge.u_nmpsr) goto 999
      end do
 
  999 return
c----------------------------------------------------------------------
      entry PSRBIN
c set up for pulsar companion with elliptic orbit

      if(Psrprm(14).le.0._10) call SUICID(
     . 'ZERO PULSAR ORBITAL PERIOD  ',7)
      cond(1)=Psrprm(10)/Aultsc
      cond(2)=Psrprm(11)
      cond(3)=90._10
      cond(4)=0._10
      cond(5)=Psrprm(12)
      cond(6)=Psrprm(13)
      goose=SQRT(cond(1)**3*(Twopi/Psrprm(14))**2)
      cond(7)=Psrprm(17)/Aultsc
      cond(8)=Psrprm(18)
      cond(9)=1._10
      Psmass=0._10
      sinc=1._10
      val=0._10
      mp3mc=1._10
      do i=10,16
         pspfct(i)=0._10
      end do
      pspfct(14)=-2._10/3._10/Psrprm(14)
c if companion mass and inclination are given, then we can calculate
c the actual orbit, aside from the node orientation
      if(Psrprm(15).gt.0._10 .and. Psrprm(16).gt.0._10) then
         cond(3)= Psrprm(15)
         sinc=SIN(Convd*cond(3))
         val=(Gauss*Psrprm(14)/Twopi)**2 * (sinc/cond(1))**3
         Psmass=1._10
         do ll=1,30
            Psmass=((Psrprm(16)+Psmass)**2/val)**(1._10/3._10)
         end do
         cond(9)=Psmass/(Psmass+Psrprm(16))
         cond(1)=cond(1)/cond(9)/sinc
         goose=Gauss*SQRT(Psmass+Psrprm(16))
         cond(7)=cond(7)/cond(9)/sinc
         mp3mc=Psmass+3._10*Psrprm(16)
         pspfct(10)=1._10/Psrprm(10)/mp3mc
c         pspfct(14)=-2._10*Psrprm(16)/Psrprm(14)/mp3mc
         pspfct(15)=SQRT(1._10-sinc*sinc)/sinc/mp3mc*Convd
         pspfct(16)=-1._10/mp3mc
      endif
      call JNITL(goose,cond,Elptn(1,n10),3,Dynpt(1,1,n10))
      return

      end
