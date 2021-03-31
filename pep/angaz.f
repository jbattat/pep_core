      subroutine ANGAZ(kathy)
 
      implicit none

c j.f.chandler - 1979 january - subroutine angaz
c calculate azimuth/elevation observables
c derived from old subr. angle

c arguments
      integer*4 kathy
c          kathy=-1  azimuth,elevation observation (topocentric)
c          kathy= 0  meridian circle observation (geocentric referred
c                    to true equinox and equator of date)
c          kathy= 1  photographic observation (topocentric without
c                    elliptic aberration removed referred to the mean
c                    equinox and equator of ref. epoch) (astrometric)
c          kathy= 2  photographic observation (topocentric with
c                    all aberration removed referred to the mean
c                    equinox and equator of ref. epoch) (astrographic)

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdateq.inc'
      include 'coord.inc'
      include 'coordxeq.inc'
      include 'coordoeq.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'ltrapobs.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*4    acctim,accdst,accprc
      equivalence (acctim,Estf),(accdst,Estf(2)),
     .            (accprc,Estf(3))
      real*10 dist,fdev
      equivalence (dist,Dstf(4)),(fdev,Dstf(10))
      integer*2 itime,ntime
      equivalence (Kob(17),itime),(Kob(18),ntime)
      include 'jdfrequv.inc'
      include 'param.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'

c external functions
      real*10 DOT

c local
      real*10 qxn,qxn2,qxo1,qxo2,rltsc,rsq,rsqh,xstp(3)
      integer*4 i,iter,j,ksu,lsw,lswm
      real*4    elref,dref,eltrue
c
c*  start=1000
c begin iteration to get observed body over meridian and/or
c to make light time correction (planetary aberration)
      lswm = -1
      lsw  = 1
      iter = 0
      do while( .true. )
         iter = iter + 1
         if(Kst1.gt.1) call TIMINC(Jd,fr(1),Jdy,fr(3),-dist/Secday)
         if(Klan.eq.u_mxpl+1) then
 
c moon or lunar orbiter
            fr(2) = fr(3)
            Jdx   = Jdy
            call MNTRP(1,Jdx,fr(2),0,lswm,1)
            if(Jdx.le.0) then
               Jd = 0
               return
            else if(Mnau.ne.1._10) then
               do i = 1,3
                  Xm(i,1) = Xm(i,1)*Mnau
               end do
            endif
 
         else if(Klan.gt.0 .and. Klan.le.u_mxpl) then
c
c read planet tape, interpolate to get planet position
c (planet velocity needed for partial derivatives)
            call PLTRP(1,Jdy,fr(3),0,lsw)
            if(Jdy.le.0) then
c
c*  start=9000
               Jd = 0
               return
            else
 
c get planet relative to earth
               do j = 1,3
                  x(j,4) = Xp(j) + x(j,2)
               end do
            endif
         endif
c
c*  start=1500
c determine if probe or satellite is observed body
         if(Klanb.gt.0) then
c
c determine satellite position
            call SBTRP(1,Jdy,fr(3),0,lsw)
            if(Jdy.le.0) then
               Jd = 0
               return
c
c determine probe or satellite relative to earth
            else if(Ncp0.ne.3) then
               ksu = 1
               if(Ncp0.ne.10) then
                  ksu = 4
                  if(Ncp0.eq.0) ksu = 2
               endif
               do j = 1,3
                  x(j,1) = Xsb(j)*Cmfct + x(j,ksu)
               end do
            else
               do j = 1,3
                  x(j,1) = Xsb(j)
               end do
            endif
         endif
c
c topocentric correction
         do j = 1,3
            x(j,Kst) = x(j,Kst) - Xsitau(j)
         end do
c
c see if there is to be further iteration
         r = SQRT(DOT(x(1,Kst),x(1,Kst)))
 
c check for sun - no distance check needed
         if(Kst1.ne.1) then
            rltsc = r*Aultsc
            if(ABS(dist-rltsc).gt.accdst) then
c
c set up quantities for next iteration
               dist = rltsc
               lsw  = -1
               if(iter.lt.20) goto 100
               call SUICID(
     .    'MORE THAN 20 ITERATIONS NEEDED FOR LIGHT TIME, STOP IN ANGAZ'
     .    ,15)
            endif
         endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c*  start=2000
c           termination calculations
         Nit(20) = Nit(20) + iter
c
c test for dummy observation below horizon
         if(Idumob.eq.1) then
            do i = 1,3
               xstp(i) = -x(i,Kst)
            end do
            call HORIZN(Xsitau,xstp,1,Jd)
            if(Jd.le.0) return
         endif
         goto 200
  100 end do
 
  200 dist = rltsc
 
c aberration correction
      do j = 1,3
         x(j,Kst) = x(j,Kst) + r*Ves(j)
         xo(j)     = x(j,Kst)
      end do
c
c is this needed?
c need planet velocity for partials
      if(Ict(3).ge.0 .or. Idumob.ne.1) then
         if(Ict(1).gt.0) call ANGVP
      endif
c
c*  start=2500
c calculate theoretical value of elevation
      rsq  = DOT(xo,xo)
      r    = SQRT(rsq)
      do i=1,3
         Xsitp0(i,1) = -xo(i)/r
      end do
      qxn  = DOT(Sitnrm,xo)
      rsqh = rsq - qxn**2
      qxn1 = qxn/r
      Deriv(2,2) = ASIN(qxn1)/Convd
 
c set up vector for partials
      qxn2 = qxn/rsq
      do i=1,3
         xnp(i) = Sitnrm(i,1) - qxn2*xo(i)
      end do
      rr1 = -SQRT(rsqh)*Convd
      if(Neqnox.gt.0) Deriv(2,2) = Deriv(2,2) + Plat
c
c add refraction correction to theoretical value of elevation
c (radio frequency refraction table from millstone)
      if(fdev.le.(1._10-900E-10_10)) then
         eltrue = Deriv(2,2)
         if(eltrue.lt.0.) eltrue = 0.
         dref = 0.
         do i = 1,3
            elref = eltrue + dref
            call DELL(elref,dref)
         end do
         Deriv(2,2) = Deriv(2,2) + dref
      endif
c (only one call to dell would be needed to remove refraction
c from observed value of elevation)
c
c calculate theoretical value of azimuth
      if(Nice.le.0) then
 
c pole unit vector is nutpr(3,.)
         do i = 1,3
            xmerid(i) = Nutpr(3,i) - Snrm(1)*Sitnrm(i,1)
         end do
         call CROSS(Sitnrm,xo,xop)
         qxo2 = DOT(xmerid,xop)
         qxo1 = DOT(xmerid,xo)
         Deriv(2,1) = ATAN2(qxo2,qxo1)/Convd
         if(Deriv(2,1).lt.0.0_10) Deriv(2,1) = Deriv(2,1)+360._10
         cake = rsqh*Convd
         if(Neqnox.gt.0) Deriv(2,1) = Deriv(2,1) + Pnox
      endif
 
c*  start=9990
      return
      end
