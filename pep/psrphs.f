      subroutine PSRPHS(tmdly)
 
      implicit none
 
c
c           subr psrphs - j.f.chandler - 1984 jun
c           compute pulsar pulse phase observables
c     in this routine the delay from the arrival (free-space) at
c     the solar-system barycenter to the observation is converted to
c     pulse phase.
c
c arguments
      real*10 tmdly

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdateq.inc'
      include 'eqnphs.inc'
      include 'obscrd.inc'
      real*10 ctrecf, freqis, phobs
      equivalence (ctrecf,Dstf),(freqis,Dstf(7))
      real*10 phi8
      equivalence (phi8,Dstf(2)),(phobs,Dstf(6))
      include 'psrcom.inc'
 
c local variables
      real*10 phi4
      real*10 third, one
      real*10 absdph, dph
c
c set up for series
      if(Nk1.lt.0) then
         Psp0 = Plsper
 
c form period from base and adjustable portion
         Psp0 = Psp0 + Psrprm(6)
 
c pulse frequency
         Psf0  = 1.0/Psp0
         Psf   = Psf0
         Psrat = -Psrprm(7)*Psf0
 
c pulse phase acceleration, etc.
         Psf1  = Psf0*Psrat
         Psf2  = 2.0*Psf1*Psrat - Psrprm(8)*Psf0**2
         one   = 1.0
         third = 3.0
         third = one/third
      endif
c
c set up extended precision time
      Tpsr(1) = Jd - Jdps0
      Tpsr(1) = (Tpsr(1) + ctrecf)*Secday - tmdly
      Tpsr(2) = Tpsr(1)**2*0.5
      Tpsr(3) = Tpsr(2)*Tpsr(1)*third
c
c compute cumulative pulse phase, and save in Dstf vector
      phi4 = Psrprm(5) + Psf0*Tpsr(1) + Psf1*Tpsr(2) + Psf2*Tpsr(3)
      phi8 = phi4
 
c save instantaneous pulse period
      Dstf(9) = Psp0 + Psrprm(7)*Tpsr(1) + Psrprm(8)*Tpsr(2)
c
c extract phase modulo cycle
      tmdly = MOD(phi4,one)
      if(tmdly.lt.0._10) tmdly = tmdly + 1._10
c
c correct for cycle slip
      if(Idumob.ne.1) then
         dph    = Result(1) - tmdly
         absdph = ABS(dph)
         if(ABS(absdph-1._10).le.0.4_10) tmdly = tmdly +
     .       SIGN(1._10,dph)
         if(phobs.ne.0._10) tmdly = phi4 - phobs
      endif
c
c set up for partials
      Dpphdt = -(Psf0 + Psf1*Tpsr(1) + Psf2*Tpsr(2))
 
      return
      end
