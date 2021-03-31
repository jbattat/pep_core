      subroutine MONLIB(mnrt1,ilib,inert)
 
      implicit none

c
c ash/freed   oct 1968    subroutine monlib
c koziel formulation of moon rotation-libration
c modified sept 1972 by r.w.king to calculate angles for lunar
c rotation integration

c parameters
      integer*4 ilib, mnrt1, inert
c           mnrt1= 0 lunorb or rotmat is calling routine
c                = 1 mrtfn is calling routine
c           ilib= 0 input angles are tau,rho,isig and derivatives
c               = 1 input angles are psi,theta,phi (mod) and derivatives
c           inert = 0 ecliptic rates are mod, noninertial
c           inert = 1 ecliptic rates are inertial w.r.t. epoch
c           moon rotation-libration matrix (output of this routine)

c array dimensions
      include 'globdefs.inc'

c commons
      include 'angles.inc'
      include 'empcnd.inc'
      real*10 meqinc
      equivalence (Mrcond(11),meqinc)
      include 'funcon.inc'
      include 'mnrtlb.inc'
      include 'precmn.inc'
      include 'prtcod.inc'
      real*4 rho,tau
      equivalence (Librat(1,1),tau),(Librat(1,2),rho)

      real*4 sigma,tausig
      real*10 u(3,3),tsigd,sigd
      integer   i, j
 
      if(ilib.le.0) then
         tausig = Librat(1,3)/meqinc
         sigma  = -tausig + tau
         Psy    = Asc + sigma
         Theta  = meqinc + rho
         if(mnrt1.ge.1) then
            Phi = Pi + Long - Asc + tausig
         else
            Phi = Pi + Anomx + tausig + Per
         endif
      endif
      Cpsi   = COS(Psy)
      Spsi   = SIN(Psy)
      Cphi   = COS(Phi)
      Sphi   = SIN(Phi)
      Ctheta = COS(Theta)
      Stheta = SIN(Theta)
      u(1,1) = Cpsi*Cphi - Spsi*Sphi*Ctheta
      u(1,2) = Spsi*Cphi + Cpsi*Sphi*Ctheta
      u(1,3) = -Sphi*Stheta
      u(2,1) = -Cpsi*Sphi - Spsi*Cphi*Ctheta
      u(2,2) = -Spsi*Sphi + Cpsi*Cphi*Ctheta
      u(2,3) = -Cphi*Stheta
      u(3,1) = -Spsi*Stheta
      u(3,2) = Cpsi*Stheta
      u(3,3) = Ctheta
      if(mnrt1.gt.0) then
         if(ilib.le.0) then
            tsigd  = Librat(2,3)/meqinc
            sigd   = -tsigd + Librat(2,1)
            Dpsy   = Ascd + sigd
            Dtheta = Librat(2,2)
            Dphi   = Longd - Ascd + tsigd
         endif
         if(inert.ge.1) then
            E1p = 0.0
            E2p = 0.0
            E3p = 0.0
         endif
         W1 = -Dtheta*Cphi - Dpsy*Sphi*Stheta +
     .        E1p*(Cphi*Cpsi - Sphi*Ctheta*Spsi)
     .        + E2p*(Cphi*Spsi + Sphi*Ctheta*Cpsi) - E3p*Sphi*Stheta
         W2 = Dtheta*Sphi - Dpsy*Cphi*Stheta -
     .        E1p*(Sphi*Cpsi + Cphi*Ctheta*Spsi)
     .        + E2p*(Cphi*Ctheta*Cpsi - Sphi*Spsi) - E3p*Cphi*Stheta
         W3 = Dpsy*Ctheta + Dphi - E1p*Stheta*Spsi + E2p*Stheta*Cpsi +
     .        E3p*Ctheta
      endif
      do i = 1, 3
         do j = 1, 3
            Mrotlb(i,j)=u(i,1)*Aa(1,j)+u(i,2)*Aa(2,j)+u(i,3)*Aa(3,j)
         end do
      end do
      return
      end
