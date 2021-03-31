      subroutine MONROT(itype,psi,theta,phi,dpsi,dtheta,dphi)
 
      implicit none

c
c r.w.king   dec 1973   subroutine monrot
c calculate rotation matrix for lunar rotation integration and
c perform transformation between reference and m.o.d. ecliptic systems
c
c parameters
      integer itype
      real*10 psi, theta, phi, dpsi, dtheta, dphi
c
c           itype=-1  determine euler angles from m.o.d. ecliptic
c                      libration angles
c           itype=0  determine q matrix between equatorial and
c                      selenodetic coordinate systems
c           itype=1  determine m.o.d. ecliptic libration angles from
c                      euler angles
c           itpye=2  same as 1, but never print results to kout
c           psi, etc. are the euler angles and angle rates

c arguments for entry DMONRT
      real*10 drtdy(3,3,3),dwdy(3,6),spsit,cpsit,sthetat,cthetat,
     . sphit,cphit,dpsyt,dthetat,dphit
c these can be either for core or exterior

c array dimensions
      include 'globdefs.inc'

c common
      include 'cassin.inc'
      include 'empcnd.inc'
      real*10 meqinc
      equivalence (Mrcond(11),meqinc)
      include 'funcon.inc'
      include 'inodta.inc'
      include 'mnrtlb.inc'
      real*10 q(3,3)
      equivalence (Mrotlb(1,1),q(1,1))
      include 'precmn.inc'

c local variables
      real*10 sig
      real*10 u(3,3)
      real*10 psix, thetax, phix, dpsix, dthetx, dphix
      real*10 spsix, cpsix, sthetx, cthetx, sphix, cphix
      integer   i, j, line
 
      data line/60/
      if(itype.lt.0) then
c
c determine angles for reference epoch from m.o.d.
         theta = ACOS(Mrotlb(3,3))
         psi   = ATAN2(Mrotlb(3,1),-Mrotlb(3,2))
         phi   = ATAN2(Mrotlb(1,3),Mrotlb(2,3))
      else if(itype.ne.0) then
c
c determine u from q (=m) and a
         do i = 1, 3
            do j = 1, 3
               u(i,j) = Mrotlb(i,1)*Aa(j,1) + Mrotlb(i,2)*Aa(j,2)
     .                   + Mrotlb(i,3)*Aa(j,3)
            end do
         end do
 
         psix = ATAN2(-u(3,1),u(3,2))
         sig  = psix - Asc
         sig  = MOD(sig,Twopi)
         if(ABS(sig).gt.1._10) then
            if(sig.le.0) then
               sig = sig + Twopi
            else
               sig = sig - Twopi
            endif
         endif
         Isig   = meqinc*sig
         thetax = ACOS(u(3,3))
         Rhoc   = thetax - meqinc
         phix   = ATAN2(-u(1,3),-u(2,3))
         spsix  = SIN(psix)
         cpsix  = COS(psix)
         sthetx = SIN(thetax)
         cthetx = COS(thetax)
         sphix  = SIN(phix)
         cphix  = COS(phix)
         Tauc   = phix - Long - Pi + psix
         Tauc   = MOD(Tauc,Twopi)
         if(ABS(Tauc).gt.1._10) then
            if(Tauc.le.0) then
               Tauc = Tauc + Twopi
            else
               Tauc = Tauc - Twopi
            endif
         endif
         dpsix  = -(W1*sphix + W2*cphix + cthetx*(E1p*spsix-E2p*cpsix)
     .            + E3p*sthetx)/sthetx
         dthetx = -W1*cphix + W2*sphix + E1p*cpsix + E2p*spsix
         dphix  = W3 + sthetx*(E1p*spsix - E2p*cpsix)
     .            - cthetx*(E3p + dpsix)
         Drhoc  = dthetx
         Disig  = (dpsix - Ascd)*meqinc
         Dtauc  = dphix - Longd + dpsix
         if(Kout.gt.0 .and. itype.eq.1) then
            if(line.ge.58) then
               write(Kout,10)
   10          format('1',
     .' EXTENDED PRINTOUT FOR ROTATION INTEGRATION - CLASSICAL LIBRATION
     .ANGLES AND RATES')
               write(Kout,20)
   20          format(4x,'JED', 12x, ' TAU ', 15x, ' RHO ', 14x,
     .                'ISIG', 14x, 'DTAU/DT', 12x, ' DRHO/DT', 11x,
     .                'DISIG/DT')
               line = 2
            endif
            write(Kout,40) Tpr,Tauc,Rhoc,Isig,Dtauc,Drhoc,Disig
   40       format(f12.3,6F20.16)
            line = line + 1
         endif
         return
      endif
      Spsi   = SIN(psi)
      Cpsi   = COS(psi)
      Stheta = SIN(theta)
      Ctheta = COS(theta)
      Sphi   = SIN(phi)
      Cphi   = COS(phi)
      if(itype.eq.0) then
c
c determine q matrix
         q(1,1) = Cphi*Cpsi - Sphi*Ctheta*Spsi
         q(1,2) = Cphi*Spsi + Sphi*Ctheta*Cpsi
         q(1,3) = Sphi*Stheta
         q(2,1) = -Sphi*Cpsi - Cphi*Ctheta*Spsi
         q(2,2) = -Sphi*Spsi + Cphi*Ctheta*Cpsi
         q(2,3) = Cphi*Stheta
         q(3,1) = Stheta*Spsi
         q(3,2) = -Stheta*Cpsi
         q(3,3) = Ctheta
      else
         dpsi   = (W1*Sphi + W2*Cphi)/Stheta
         dtheta = W1*Cphi - W2*Sphi
         dphi   = W3 - dpsi*Ctheta
      endif
      return
c------------------------------------------------------------
c MNCROT functions as MONROT itype=0 for the lunar core
      entry MNCROT(psi,theta,phi,dpsi,dtheta,dphi)

      Spsic   = SIN(psi)
      Cpsic   = COS(psi)
      Sthetac = SIN(theta)
      Cthetac = COS(theta)
      Sphic   = SIN(phi)
      Cphic   = COS(phi)
c
c determine rotation matrix
      Mcrtlb(1,1) = Cphic*Cpsic - Sphic*Cthetac*Spsic
      Mcrtlb(1,2) = Cphic*Spsic + Sphic*Cthetac*Cpsic
      Mcrtlb(1,3) = Sphic*Sthetac
      Mcrtlb(2,1) = -Sphic*Cpsic - Cphic*Cthetac*Spsic
      Mcrtlb(2,2) = -Sphic*Spsic + Cphic*Cthetac*Cpsic
      Mcrtlb(2,3) = Cphic*Sthetac
      Mcrtlb(3,1) = Sthetac*Spsic
      Mcrtlb(3,2) = -Sthetac*Cpsic
      Mcrtlb(3,3) = Cthetac
      return
c------------------------------------------------------------
      entry DMONRT(drtdy,dwdy,spsit,cpsit,sthetat,cthetat,sphit,cphit,
     . dpsyt,dthetat,dphit)
c calculate partials of rotation matrix wrt euler angles
c Drtdy is transpose for convenience: Drtdy(i,j,k)=dMrotlb(j,i)/dy(k)
      drtdy(1,1,1) = -(spsit*cphit + cpsit*cthetat*sphit)
      drtdy(2,1,1) =  (cpsit*cphit - spsit*cthetat*sphit)
      drtdy(3,1,1) = 0._10
      drtdy(1,2,1) =  (spsit*sphit - cpsit*cthetat*cphit)
      drtdy(2,2,1) = -(cpsit*sphit + spsit*cthetat*cphit)
      drtdy(3,2,1) = 0._10
      drtdy(1,3,1) = cpsit*sthetat
      drtdy(2,3,1) = spsit*sthetat
      drtdy(3,3,1) = 0._10
      drtdy(1,1,2) =  spsit*sthetat*sphit
      drtdy(2,1,2) = -cpsit*sthetat*sphit
      drtdy(3,1,2) = cthetat*sphit
      drtdy(1,2,2) =  spsit*sthetat*cphit
      drtdy(2,2,2) = -cpsit*sthetat*cphit
      drtdy(3,2,2) = cthetat*cphit
      drtdy(1,3,2) = spsit*cthetat
      drtdy(2,3,2) = -cpsit*cthetat
      drtdy(3,3,2) = -sthetat
      drtdy(1,1,3) = -(cpsit*sphit + spsit*cthetat*cphit)
      drtdy(2,1,3) =  (cpsit*cthetat*cphit - spsit*sphit)
      drtdy(3,1,3) = sthetat*cphit
      drtdy(1,2,3) =  (spsit*cthetat*sphit - cpsit*cphit)
      drtdy(2,2,3) = -(spsit*cphit + cpsit*cthetat*sphit)
      drtdy(3,2,3) = -sthetat*sphit
      drtdy(1,3,3) = 0._10
      drtdy(2,3,3) = 0._10
      drtdy(3,3,3) = 0._10

c dependence of angular velocity on euler angles and rates
c order is psi, theta, phi, psidot, thetadot, phidot
      dwdy(1,1) = 0._10
      dwdy(2,1) = 0._10
      dwdy(3,1) = 0._10
      dwdy(1,2) = dpsyt*sphit*cthetat
      dwdy(2,2) = dpsyt*cphit*cthetat
      dwdy(3,2) = -dpsyt*sthetat
      dwdy(1,3) = -dthetat*sphit+dpsyt*cphit*sthetat
      dwdy(2,3) = -dthetat*cphit-dpsyt*sphit*sthetat
      dwdy(3,3) = 0._10
      dwdy(1,4) = sphit*sthetat
      dwdy(2,4) = cphit*sthetat
      dwdy(3,4) = cthetat
      dwdy(1,5) = cphit
      dwdy(2,5) = -sphit
      dwdy(3,5) = 0._10
      dwdy(1,6) = 0._10
      dwdy(2,6) = 0._10
      dwdy(3,6) = 1._10
      return
      end
