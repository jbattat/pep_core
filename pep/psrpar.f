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
      real*10 cond(6), trm
      integer*4 iflag,iswtch,l,ll,m,mm,mt
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
      m = 1
      l = 0
      call PCOPS(m, 'PSRX', iabs1)
      do while( .true. )
c
c copy/merge partials from input obslib tape
         iflag  = 1
         iswtch = 0
         call PCOPY(l, 16, iflag, iswtch, Lpsrx, Mpsrx)
         if(iflag.gt.0) goto 999
         ll = Lpsrx(l)
 
c clear partial derivatives, in case
         Deriv(kind, 1) = 0._10
         Deriv(kind, 2) = 0._10
c
c ntype=0: model has p, p-dot, p-dot-dot
         if((ll.gt.8 .and. Psrprm(10).le.0._10) .or. ll.gt.14) call
     .    SUICID('CANNOT CALCULATE PARTIAL W.R.T. UNDEFINED PULSAR PARAM
     .ETER, STOP IN PSRPAR  ', 19)
c
         if(kick.eq.1) then
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
               Deriv(kind,1) = Dpphdt/freqis**2
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
 
            else if(ll.eq.10) then
c pulsar companion semimajor axis (explicitly non-secular)
               trm = Ynpt(3,n10)/Psrprm(10)
 
            else if(ll.eq.11) then
c pulsar companion eccentricity
               trm = Dynpt(3,2,n10)
 
            else if(ll.eq.12) then
c pulsar companion periapse
               trm = Dynpt(3,5,n10)*Convd
 
            else if(ll.eq.13) then
c pulsar companion initial mean anomaly
               trm = Dynpt(3,6,n10)*Convd
 
            else if(ll.eq.14) then
c pulsar companion period
               trm = -Dynpt(3,6,n10)*Frect*Twopi/Psrprm(14)**2

            endif
            if(ll.ge.6 .and. ll.le.8) then
               Deriv(kind,1) = -trm*Psf
            else if(ll.ge.10 .and. ll.le.14) then
               Deriv(kind,1) = Dpphdt*trm*Aultsc/Opv3
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
         if(l.ge.16) goto 999
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
      call JNITL(SQRT(cond(1)**3*(Twopi/Psrprm(14))**2),
     . cond, Elptn(1,n10), 1, Dynpt(1,1,n10))
      return

      end
