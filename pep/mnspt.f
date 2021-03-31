      subroutine MNSPT(jda,fract,mnspt1,nvel,n,nlibp2)
 
      implicit none

c
c ash/freed   oct 1968   subroutine mnspt
c calculation of moon spot coordinates relative to center of moon
c modified nov 1972 by r.king to calculate partials w.r.t. lunar
c libration parameters
c
c parameters
      real*10 fract
      integer*4 jda,mnspt1,nvel,n
      integer*2 nlibp2
c
c           n=1 spot(nspot) on nplnt0
c           n=2 spot(nspot2) on nplnt2
c
c     nlibpr=-2 libration and partials determined from eckhardt's series
c           =-1 libration determined from series, no partials calculated
c           = 0 libration taken from moon tape, no partials calculated
c           = 1 libration taken from moon rotation tape, no partials
c               calculated
c           = 2 libration and partials taken from moon rotation tape
c
c     for nlibpr=1,2, kmr(100) (in /rotdta/) determines type of angles:
c           =-1 moon rotation angles are euler angles wrt coordinate
c               system (pep inegration)
c           = 0 moon rotation angles are cassini libration angles wrt
c               osculating initial conditions (jpl numerical integration
c           = 1 moon rotation angles are cassini libration angles wrt
c               mean elements (eckhardt analytic theory)

c array dimensions
      include 'globdefs.inc'

c commons
      include 'comdateq.inc'
      include 'coord.inc'
      include 'empcnd.inc'
      real*10 alpha, beta, gamma, meqinc
      equivalence (Mrcond(8),alpha),(Mrcond(9),beta),
     .            (Mrcond(10),gamma),(Mrcond(11),meqinc)
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'hoccom.inc'
      include 'mnsprt.inc'
      real*10 dxdisg(6,2),dxdrho(6,2),dxdtau(6,2)
      equivalence (Dxdpsi(1,1),dxdisg(1,1)),(Dxdth(1,1),dxdrho(1,1)),
     .            (Dxdphi(1,1),dxdtau(1,1))
      real*10 dxdbet(6,2),dxdgam(6,2)
      equivalence (Dxdpsi(1,1),dxdbet(1,1)),(Dxdphi(1,1),dxdgam(1,1))
      real*4    librat(2,3),tau,rho,dtau,drho
      equivalence (librat(1,1),tau),(librat(1,2),rho),
     .            (librat(2,1),dtau),(librat(2,2),drho)
      real*10 isig
      real*10 du(3,3),dudpsi(3,3),dudth(3,3),dudphi(3,3),
     .          d2udps(3,3),d2udth(3,3),d2udph(3,3)
      real*10 dmrtlb(3,3)
      real*10 dmdpsi(3,3),dmdth(3,3),dmdphi(3,3)
      real*10 d2psps(3,3),d2thth(3,3),d2phph(3,3),d2thps(3,3),
     .          d2phps(3,3),d2phth(3,3)
      equivalence (dudpsi(1,1),d2psps(1,1)),(dudphi(1,1),d2phph(1,1)),
     1   (du(1,1),d2thth(1,1)),(dmdpsi(1,1),d2thps(1,1)),
     2   (dmdth(1,1),d2phps(1,1))
      real*4    dlibrt(6,2),taubet,rhobet,isgbet,ibet,dtaubt,
     .          drhobt, disgbt, taugam, rhogam, isggam, igam, dtaugm,
     .          drhogm, disggm
      equivalence (dlibrt(1,1),taubet),(dlibrt(2,1),rhobet),
     .            (dlibrt(3,1),isgbet),(dlibrt(4,1),dtaubt),
     .            (dlibrt(5,1),drhobt),(dlibrt(6,1),disgbt),
     .            (dlibrt(1,2),taugam),(dlibrt(2,2),rhogam),
     .            (dlibrt(3,2),isggam),(dlibrt(4,2),dtaugm),
     .            (dlibrt(5,2),drhogm),(dlibrt(6,2),disggm)
      include 'param.inc'
      include 'rotdtamr.inc'
c
c local
      real*10 anomx,asc,clony,cphi,cpsi,cth,ctvcor,danom,dasc,dper,
     . dphi,dphidt,dpsi,dpsidt,dsig,dtausg,dth,dthdt,
     . dxdbt1,frt,per,per1,phi,psi,qq,qv,slony,sigma,sphi,spsi,sth,
     . tausig,theta,tt,x(6),xms(3),ysptmp(3)
      integer*4 i,ii,int,int1,j,jdt,jnt,nlibpr
      real*10 DOT
 
      nlibpr = nlibp2
c
c
c determination of time
c remove any ctvary correction
      ctvcor = Ctvary*(jda - Prm97 - 0.5_10 + fract)**2
      call TIMINC(jda,fract,jdt,frt,-ctvcor)
      int = jdt - 2415020
      tt  = int
      qq  = frt - 0.5_10
      tt  = tt + qq
 
c if(mnspt1.gt.0) goto 40
      mnspt1 = 1
 
      if(nlibpr.eq.0) then
c
c interpolation of libration from nbody or moon tape
         if(meqinc.eq.0._10) meqinc = 0.0268587_10
         call MNLIB(jdt,frt,librat,nvel)
c
c series for libration
      else if(nlibpr.ge.0) then
c
c libration read from tape
         call ROTCRD(jdt,frt,x,nvel,10,1)
         if(jdt.le.0) then
c
c
c not found on rotation tape
            jda = 0
            return
         else if(Kmr(100).lt.0) then
 
c get rot from euler angles
            psi   = x(1)
            theta = x(2)
            phi   = x(3)
            dpsi  = x(4)
            dth   = x(5)
            dphi  = x(6)
            cpsi  = COS(psi)
            spsi  = SIN(psi)
            cth   = COS(theta)
            sth   = SIN(theta)
            cphi  = COS(phi)
            sphi  = SIN(phi)
 
            Rot(1,1) = cphi*cpsi - sphi*cth*spsi
            Rot(1,2) = cphi*spsi + sphi*cth*cpsi
            Rot(1,3) = sphi*sth
            Rot(2,1) = -sphi*cpsi - cphi*cth*spsi
            Rot(2,2) = -sphi*spsi + cphi*cth*cpsi
            Rot(2,3) = cphi*sth
            Rot(3,1) = sth*spsi
            Rot(3,2) = -sth*cpsi
            Rot(3,3) = cth
            goto 300
         else
            librat(1,1) = x(1)
            librat(1,2) = x(2)
            librat(1,3) = meqinc*x(1) - x(3)
            if(nvel.gt.0) then
               librat(2,1) = x(4)
               librat(2,2) = x(5)
               librat(2,3) = meqinc*x(4) - x(6)
            endif
         endif
      else
         call DLIBRA(jdt,frt,librat,meqinc,nvel,beta,gamma)
      endif
c
c     the a matrix is no longer calculated in mnspt, but rather in
c     preces in order to make use of the precession matrix before it
c     is modified by the ad hoc small rotation angles.  this means that
c     any adjustment of these angles to account for errors in the
c     earth's orientation will not affect the orientation of the moon.
c
c        determination of ascending node
      per1 = 7.1995354166666667E-1_10 +
     .       tt*(-1.4709422833333333E-4_10 + tt*(4.325E-15_10+
     .       tt*1.3888888888888889E-22_10))
      int1 = 2
  100 jnt  = per1
      if(per1.lt.0) then
         jnt = jnt - 1
      else if(per1.eq.0) then
         goto 200
      endif
      qv   = jnt
      per1 = (per1 - qv)*Twopi
  200 if(int1.eq.1) then
         per = per1
c
c
c calculate psi, theta, sines, and cosines
         tausig = librat(1,3)/meqinc
         isig   = -librat(1,3) + meqinc*tau
         sigma  = -tausig + tau
         psi    = asc + sigma
         theta  = meqinc + rho
         if(Kmr(100).eq.1) theta = rho
         cpsi = COS(psi)
         spsi = SIN(psi)
         cth  = COS(theta)
         sth  = SIN(theta)
c
c determine mean anomaly
         int  = 1306*int
         int1 = int/36000
         per1 = int - int1*36000
         per1 = ((per1-6.3895392E3_10) + qq*1306._10 +
     .    tt*(0.49924465_10+tt*(6.889E-10_10+tt*2.99E-17_10)))/3.6E4_10
         int1 = 3
      else if(int1.eq.2) then
         asc = per1
c
c determination of perigee
         per1 = 2.08739669444444444444E-1_10 +
     .    tt*(4.56550006944444444444E-4_10 -
     .    tt*(2.58222222222222222222E-14_10+
     .    tt*8.61111111111111111111E-22_10))
         int1 = 1
      else
         anomx = per1
c
c calculate u matrix
         phi     = Pi + anomx + tausig + per
         cphi    = COS(phi)
         sphi    = SIN(phi)
         U(1,1) = cpsi*cphi - spsi*sphi*cth
         U(1,2) = spsi*cphi + cpsi*sphi*cth
         U(1,3) = -sphi*sth
         U(2,1) = -cpsi*sphi - spsi*cphi*cth
         U(2,2) = -spsi*sphi + cpsi*cphi*cth
         U(2,3) = -cphi*sth
         U(3,1) = -spsi*sth
         U(3,2) = cpsi*sth
         U(3,3) = cth
c
c calculate moon rotation-libration matrix
         do i = 1, 3
            do j = 1, 3
               Rot(i,j) = U(i,1)*A(1,j) + U(i,2)*A(2,j) + U(i,3)
     .                     *A(3,j)
            end do
         end do
         goto 300
      endif
      goto 100
c
c compute spot coordinates in standard coordinate system
  300 do i=1,3
         ysptmp(i)=Yspcd(i,n)
      end do
c apply ad hoc libration correction from external source
      call HOCRED(jda,fract)
      if(Ihoc.eq.0) then
         call PRODCT(Rot,Xm,xms,3,3,1)
         clony=Hocy/SQRT(xms(1)**2+xms(2)**2)
         slony=-xms(2)*clony
         clony=-xms(1)*clony
         ysptmp(1)=Yspcd(1,n)+Hocx*Yspcd(2,n)+clony*Yspcd(3,n)
         ysptmp(2)=-Hocx*Yspcd(1,n)+Yspcd(2,n)-slony*Yspcd(3,n)
         ysptmp(3)=-clony*Yspcd(1,n)+slony*Yspcd(2,n)+Yspcd(3,n)
      endif
      call PRODCT(Rot,ysptmp,Xspcd(1,n),-3,3,1)
      call TRPLST(jdt,frt,0,0,"MNSPT",Xspcd(1,n))

      return
      entry MNSPTV(nvel,n,nlibp2)
      nlibpr = nlibp2
c
c determine partial derivatives of spot position
c with respect to local radius
      if(Lspot(1,n).gt.0 .or. Jct(11).gt.0) then
         do i = 1, 3
            Dxspcd(i,1,n) = Xspcd(i,n)/Rspot(n)
         end do
      endif
c
c with respect to spot longitude
      if(Lspot(2,n).gt.0 .or. Jct(11).gt.0) then
         do i = 1, 3
            Dxspcd(i,2,n) = -Rot(1,i)*Yspcd(2,n) + Rot(2,i)*Yspcd(1,n)
         end do
      endif
c
c with respect to spot latitude
      if(Lspot(3,n).gt.0 .or. Jct(11).gt.0) then
         do i = 1, 3
            Dxspcd(i,3,n) = Rot(1,i)*Dydphi(1,n) + Rot(2,i)
     .                        *Dydphi(2,n) + Rot(3,i)*Dydphi(3,n)
         end do
      endif
      if(Kmr(100).ge.0) then
 
         if(nvel.le.0) then
            if(IABS(nlibpr).le.1) return
         endif
c
c partial of u w.r.t. psi
         dudpsi(1,1) = -U(1,2)
         dudpsi(1,2) = U(1,1)
         dudpsi(1,3) = 0.0_10
         dudpsi(2,1) = -U(2,2)
         dudpsi(2,2) = U(2,1)
         dudpsi(2,3) = 0.0_10
         dudpsi(3,1) = -U(3,2)
         dudpsi(3,2) = U(3,1)
         dudpsi(3,3) = 0.0_10
 
c partial of u w.r.t. theta
         dudth(1,1) = spsi*sth*sphi
         dudth(1,2) = -cpsi*sth*sphi
         dudth(1,3) = -cth*sphi
         dudth(2,1) = spsi*sth*cphi
         dudth(2,2) = -cpsi*sth*cphi
         dudth(2,3) = -cth*cphi
         dudth(3,1) = -spsi*cth
         dudth(3,2) = cpsi*cth
         dudth(3,3) = -sth
 
c partial of u w.r.t. phi
         dudphi(1,1) = U(2,1)
         dudphi(1,2) = U(2,2)
         dudphi(1,3) = U(2,3)
         dudphi(2,1) = -U(1,1)
         dudphi(2,2) = -U(1,2)
         dudphi(2,3) = -U(1,3)
         dudphi(3,1) = 0.0_10
         dudphi(3,2) = 0.0_10
         dudphi(3,3) = 0.0_10
c
c determine partial derivatives of rot w.r.t.
c psi, theta, phi
         if(iabs(nlibpr).le.1) goto 400
         call PRODCT(dudpsi,A,dmdpsi,3,3,3)
         call PRODCT(dudth,A,dmdth,3,3,3)
         call PRODCT(dudphi,A,dmdphi,3,3,3)
      else
 
c euler angles are used - find partials of m
         dmdpsi(1,1) = -Rot(1,2)
         dmdpsi(1,2) = Rot(1,1)
         dmdpsi(1,3) = 0.0_10
         dmdpsi(2,1) = -Rot(2,2)
         dmdpsi(2,2) = Rot(2,1)
         dmdpsi(2,3) = 0.0_10
         dmdpsi(3,1) = -Rot(3,2)
         dmdpsi(3,2) = Rot(3,1)
         dmdpsi(3,3) = 0.0_10
 
         dmdth(1,1) = sphi*sth*spsi
         dmdth(1,2) = -sphi*sth*cpsi
         dmdth(1,3) = sphi*cth
         dmdth(2,1) = cphi*sth*spsi
         dmdth(2,2) = -cphi*sth*cpsi
         dmdth(2,3) = cphi*cth
         dmdth(3,1) = cth*spsi
         dmdth(3,2) = -cth*cpsi
         dmdth(3,3) = -sth
 
         dmdphi(1,1) = Rot(2,1)
         dmdphi(1,2) = Rot(2,2)
         dmdphi(1,3) = Rot(2,3)
         dmdphi(2,1) = -Rot(1,1)
         dmdphi(2,2) = -Rot(1,2)
         dmdphi(2,3) = -Rot(1,3)
         dmdphi(3,1) = 0.0_10
         dmdphi(3,2) = 0.0_10
         dmdphi(3,3) = 0.0_10
         if(nvel.gt.0) then
 
c calculate spot velocities
            do i = 1, 3
               do j = 1, 3
                  dmrtlb(i,j) = (dmdpsi(i,j)*dpsi + dmdth(i,j)*dth +
     .                           dmdphi(i,j)*dphi)/8.64E4_10
               end do
            end do
            call PRODCT(dmrtlb,Yspcd(1,n),Xspcd(4,n),-3,3,1)
         endif
      endif
c
c determine partial derivative of spot position
c w.r.t. psi,theta,phi
      do j = 1, 3
 
c cpartl expects units of partials affecting motion to be au
         Dxdpsi(j,n) = (dmdpsi(1,j)*Yspcd(1,n) + dmdpsi(2,j)*Yspcd(2,n)
     .                  + dmdpsi(3,j)*Yspcd(3,n))/Aultsc
         Dxdth(j,n)  = (dmdth(1,j)*Yspcd(1,n) + dmdth(2,j)*Yspcd(2,n)
     .                  + dmdth(3,j)*Yspcd(3,n))/Aultsc
         Dxdphi(j,n) = (dmdphi(1,j)*Yspcd(1,n) + dmdphi(2,j)*Yspcd(2,n)
     .                  + dmdphi(3,j)*Yspcd(3,n))/Aultsc
      end do
      if(Kmr(100).lt.0) goto 500
c
c determine partial derivatives of spot position
c w.r.t. isig, rho,tau
      do i = 1, 3
         Dxdi(i,n) = Dxdth(i,n)
     .                + isig/meqinc**2*(Dxdphi(i,n) - Dxdpsi(i,n))
         if(Kmr(100).le.0 .and. nlibpr.ge.0) Dxdi(i,n) = 0.0_10
         dxdisg(i,n) = (Dxdpsi(i,n) - Dxdphi(i,n))/meqinc
      end do
c     dxdrho(i,n)=dxdth(i,n)           by virtue of
c     dxdtau(i,n)= dxdphi(i,n)         equivalent storage
c
c           determine partial derivatives of spot position w.r.t. beta
c           and gamma if series calculation used
      if(nlibpr.lt.0) then
         call DLIBP1(nvel,dlibrt,ibet,igam)
         do i = 1, 3
            dxdbt1 = dxdisg(i,n)*isgbet + dxdrho(i,n)
     .               *rhobet + dxdtau(i,n)*taubet + Dxdi(i,n)*ibet
            dxdgam(i,n) = dxdisg(i,n)*isggam + dxdrho(i,n)
     .                     *rhogam + dxdtau(i,n)*taugam + Dxdi(i,n)
     .                     *igam
            dxdbet(i,n) = dxdbt1
         end do
      endif
c
c determine spot velocity
  400 if(nvel.le.0) return
      danom = 2.2802713493961401E-1_10 +
     .        tt*(2.404714643398E-13_10 + tt*1.5655603390389E-20_10)
      dasc  = -9.24220294234919E-4_10 +
     .        tt*(5.434955290710E-14_10 + tt*2.617993877991E-21_10)
      dper  = 2.868588295626071E-3_10 -
     .        tt*(3.2449161453178E-13_10 + tt*1.62315620435472E-20_10)
 
c time derivatives of angles
      dtausg = librat(2,3)/meqinc
      dsig   = -dtausg + dtau
      dpsidt = (dasc + dsig)/8.64E4_10
      dthdt  = drho/8.64E4_10
      dphidt = (danom + dper + dtausg)/8.64E4_10
 
c time derivative of u
      do i = 1, 3
         do j = 1, 3
            du(i,j) = dudpsi(i,j)*dpsidt + dudth(i,j)
     .                 *dthdt + dudphi(i,j)*dphidt
         end do
      end do
 
c time derivative of moon rotation,libration matrix
      call PRODCT(du,A,dmrtlb,3,3,3)
 
c velocity of moon spot relative to center of moon in standard coords.
      call PRODCT(dmrtlb,Yspcd(1,n),Xspcd(4,n),-3,3,1)
      call TRPLST(jdt,frt,1,0,"MNSPT",Xspcd(1,n))
c
c        determine partial derivatives of spot velocity
c        with respect to local radius
c          test if derivs. of spot vel. needed
c  * caution *  code below not tested
  500 if(Lspot(1,n).gt.0 .or. Jct(11).gt.0) then
         do i = 4, 6
            Dxspcd(i,1,n) = Xspcd(i,n)/Rspot(n)
         end do
      endif
      if(Lspot(2,n).gt.0 .or. Jct(11).gt.0) then
c
c with respect to spot longitude
         do i = 4, 6
            j = i - 3
            Dxspcd(i,2,n) = -dmrtlb(1,j)*Yspcd(2,n) + dmrtlb(2,j)
     .                        *Yspcd(1,n)
         end do
      endif
c
c with respect to spot latitude
      if(Lspot(3,n).gt.0 .or. Jct(11).gt.0)
     .    call PRODCT(dmrtlb,Dydphi(1,n),Dxspcd(4,3,n),-3,3,1)
      if(Kmr(100).ge.0) then
         if(iabs(nlibpr).gt.1) then
            if(Ncode.gt.1) then
               call SUICID(
     .' PROBLEM WITH INDICES OF DM WRT LIBRATION PARAMETERS-STOP IN MNSP
     .T  ', 17)
c
c determine partial derivatives of spot velocty
c w.r.t. psi, theta, phi
               d2psps(1,1) = -U(1,1)
               d2psps(1,2) = -U(1,2)
               d2psps(1,3) = 0.0_10
               d2psps(2,1) = -U(2,1)
               d2psps(2,2) = -U(2,2)
               d2psps(2,3) = 0.0_10
               d2psps(3,1) = -U(3,1)
               d2psps(3,2) = -U(3,2)
               d2psps(3,3) = 0.0_10
               d2thth(1,1) = spsi*cth*sphi
               d2thth(1,2) = -cpsi*cth*sphi
               d2thth(1,3) = sth*sphi
               d2thth(2,1) = spsi*cth*cphi
               d2thth(2,2) = -cpsi*cth*cphi
               d2thth(2,3) = sth*cphi
               d2thth(3,1) = spsi*sth
               d2thth(3,2) = -cpsi*sth
               d2thth(3,3) = -cth
               d2phph(1,1) = -U(1,1)
               d2phph(1,2) = -U(1,2)
               d2phph(1,3) = -U(1,3)
               d2phph(2,1) = -U(2,1)
               d2phph(2,2) = -U(2,2)
               d2phph(2,3) = -U(2,3)
               d2phph(3,1) = 0.0_10
               d2phph(3,2) = 0.0_10
               d2phph(3,3) = 0.0_10
               d2thps(1,1) = cpsi*sth*sphi
               d2thps(1,2) = spsi*sth*sphi
               d2thps(1,3) = 0.0_10
               d2thps(2,1) = cpsi*sth*cphi
               d2thps(2,2) = spsi*sth*cphi
               d2thps(2,3) = 0.0_10
               d2thps(3,1) = -cpsi*cth
               d2thps(3,2) = -spsi*cth
               d2thps(3,3) = 0.0_10
               d2phps(1,1) = -U(2,2)
               d2phps(1,2) = U(2,1)
               d2phps(1,3) = 0.0_10
               d2phps(2,1) = U(1,2)
               d2phps(2,2) = -U(1,1)
               d2phps(2,3) = 0.0_10
               d2phps(3,1) = 0.0_10
               d2phps(3,2) = 0.0_10
               d2phps(3,3) = 0.0_10
               d2phth(1,1) = dudth(2,1)
               d2phth(1,2) = dudth(2,2)
               d2phth(1,3) = dudth(2,3)
               d2phth(2,1) = -dudth(1,1)
               d2phth(2,2) = -dudth(1,2)
               d2phth(2,3) = -dudth(1,3)
               d2phth(3,1) = 0.0_10
               d2phth(3,2) = 0.0_10
               d2phth(3,3) = 0.0_10
 
c following do's could be replaced by calls to prodct
               do i = 1, 3
                  do j = 1, 3
                     d2udps(i,j) = d2psps(i,j)*dpsidt + d2thps(i,j)
     .                              *dthdt + d2phps(i,j)*dphidt
                     d2udth(i,j) = d2thps(2,j)*dpsidt + d2thth(i,j)
     .                              *dthdt + d2phth(i,j)*dphidt
                     d2udph(i,j) = d2phps(2,j)*dpsidt + d2phth(i,j)
     .                              *dthdt + d2phph(i,j)*dphidt
                  end do
               end do
               call PRODCT(d2udps,A,dmdpsi,3,3,3)
               call PRODCT(d2udth,A,dmdth,3,3,3)
               call PRODCT(d2udph,A,dmdphi,3,3,3)
               do j = 4, 6
                  Dxdpsi(j,n)= DOT(dmdpsi(1,j-3),Yspcd(1,n))/Aultvl
                  Dxdth(j,n) = DOT(dmdth(1,j-3),Yspcd(1,n))/Aultvl
                  Dxdphi(j,n)= DOT(dmdphi(1,j-3),Yspcd(1,n))/Aultvl
               end do
c
c determine partial derivatives of spot velocity
c w.r.t. isig, rho, tau
               do i = 4, 6
                  Dxdi(i,n)= Dxdth(i,n) + isig/meqinc**2*(Dxdphi(i,n) -
     .             Dxdpsi(i,n) + Dxdphi(ii,n) - Dxdpsi(ii,n))
                  dxdisg(i,n) = (Dxdpsi(i,n) - Dxdphi(i,n))/meqinc
               end do
c     dxdrho(i,n)=dxdth(i,n)           by virtue of
c     dxdtau(i,n)= dxdphi(i,n)         equivalent storage
c
c           determine partial derivatives of spot velocity w.r.t. beta
c           and gamma if series calculation used
               if(nlibpr.eq.-2) then
                  do i = 4, 6
                     ii     = i - 3
                     dxdbt1 = dxdisg(i,n)*isgbet + dxdrho(i,n)
     .                        *rhobet + dxdtau(i,n)*taubet + Dxdi(i,n)
     .                        *ibet + disgbt*dxdisg(ii,n)
     .                        + drhobt*dxdrho(ii,n)
     .                        + dtaubt*dxdtau(ii,n)
                     dxdgam(i,n) = Dxdi(i,n)*igam + dxdisg(i,n)
     .                              *isggam + dxdrho(i,n)
     .                              *rhogam + dxdtau(i,n)
     .                              *taugam + disggm*dxdisg(ii,n)
     .                              + drhogm*dxdrho(ii,n)
     .                              + dtaugm*dxdtau(ii,n)
                     dxdbet(i,n) = dxdbt1
                  end do
               endif
            endif
         endif
      endif
      return
      end
