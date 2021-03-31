      subroutine CHNCNE(nplnt,ncentr,cond,jd0,fract0,lcnd,jtype,
     .                  incnd,icnd,iepoch)

      implicit none

c     m.e.ash    april 1973    subroutine chncne
c  Transform orbital elements to mean equatorial frame of epoch if
c  necessary.  Change initial orbital elements of type 'incnd' to type
c  'icnd' (possibly by way of Cartesian).
c

c parameters
      integer*2 nplnt,ncentr,jtype,icnd
      integer*4 jd0,lcnd,incnd,iepoch
      real*10 cond(6),fract0

c NPLNT    = planet number
c NCENTR   = central body number
c LCND = 0 nothing
c LCND = 1 input adbarv cond(4) is velocity and not difference between
c          velocity and circular orbit velocity  (units au & au/day)
c          change to difference velocity
c
c JTYPE specifies the coordinate frame of input elements - to be
c       converted to mean equinox and equator of reference epoch.
c JTYPE = 0 mean equinox and equator of input date
c       = 1 true equinox and equator of input date
c       = 2 mean equinox and ecliptic of input date
c       = 4 mean lunar plane of date (x axis along intersection of mean
c           lunar and ecliptic planes of date)
c       = 5 mean equinox and true equator of input date
c
c INCND=-1 input initial conditions are cartesian
c      = 0 input initial conditions are elliptic elements
c      = 1 input initial conditions are elliptic with angle sums
c      = 2 input initial conditions are adbarv
c ICND =   same logic as incnd: specifies output initial conditions
c      Both INCND and ICND are ignored for rotation i.c.'s.
c
c IEPOCH=1 Reference epoch is B1950.0
c       =2 Reference epoch is J2000.0

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'funcon.inc'
      include 'param.inc'
c
c local
      integer*4 i,icndx,icndt,jd9
      real*10 casc,cinc,cper,cphi,cpsi,cth,DOT,drot(3,3),dum,elpt(18),
     . fract9,gaussx,nutpr(3,3),r,r2,r3,rot(3,3),sasc,sinc,sper,sphi,
     . spsi,sth,tmprot(3,3),v,vn,vp(3),vpep(3),xasc,xinc,xmass,xper,
     . x0(6),ynorm(3),yper(3)
      equivalence (elpt(1),nutpr(1,1)),(elpt(10),x0(1)),(elpt(16),vp(1))
c
c set epoch to B1950 or J2000
      jd9 = 2433282
      fract9 = 0.923_10
      if(iepoch.eq.2) then
         jd9 = 2451545
         fract9 = 0.5_10
      endif
c
c ICNDX tracks the current form of orbital elements,
c starting from INCND and ending at ICND
      icndx = incnd

      if(nplnt.gt.0) then
         icndt = icnd
         if(icnd.eq.1) icndt = 0
         if(icndx.eq.1) then
            if(icndx.eq.icnd .and. jd0.eq.jd9 .and.
     .       ABS(fract9-fract0).lt.1.E-5_10) goto 100

c convert input elliptic to normal form
            cond(6) = cond(6) - cond(5)
            if(cond(6).lt.0._10) cond(6) = cond(6) + 360._10
            cond(5) = cond(5) - cond(4)
            if(cond(5).lt.0._10) cond(5) = cond(5) + 360._10
            icndx = 0
         endif
c
c calculate constants
         xmass = 0._10
         if(nplnt.le.30) xmass = Mass(nplnt)
         if(ncentr.gt.0) then
            xmass = Mass(ncentr)
            if(nplnt.ne.10) then
               if(ncentr.eq.3) xmass  = Mass(3)*(1._10 - Mass(10))
               if(ncentr.eq.10) xmass = Mass(3)*Mass(10)
            endif
         else
            xmass = 1._10 + xmass
         endif
         gaussx = Gauss*SQRT(xmass)
c
c fix velocity in adbarv initial conditions
         if(lcnd.ge.1) cond(4) = cond(4) - gaussx/SQRT(cond(1))
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c shall we change from true of date to mean of epoch
      if(jtype.gt.0 .or. jd0.ne.jd9 .or.
     . ABS(fract9-fract0).gt.1.E-5_10) then
         call PRCNT9(jd0,fract0,jtype,nutpr,jd9,fract9)
         if(nplnt.lt.0) then
c transform rotation i.c.'s (Euler angles)
c  psi   = cond(1)
c  theta = cond(2)
c  phi   = cond(3)
c  dpsi  = cond(4)
c  dth   = cond(5)
c  dphi  = cond(6)
            if(cond(1).eq.0._10 .and. cond(2).eq.0._10
     .       .and. cond(3).eq.0._10) return
            cpsi  = COS(cond(1))
            spsi  = SIN(cond(1))
            cth   = COS(cond(2))
            sth   = SIN(cond(2))
            cphi  = COS(cond(3))
            sphi  = SIN(cond(3))
c
c get body-frame-to-original-frame matrices
            rot(1,1) = cphi*cpsi - sphi*cth*spsi
            rot(1,2) = cphi*spsi + sphi*cth*cpsi
            rot(1,3) = sphi*sth
            rot(2,1) = -sphi*cpsi - cphi*cth*spsi
            rot(2,2) = -sphi*spsi + cphi*cth*cpsi
            rot(2,3) = cphi*sth
            rot(3,1) = sth*spsi
            rot(3,2) = -sth*cpsi
            rot(3,3) = cth
            do i=1,3
               drot(1,i)= sphi*rot(3,i)*cond(5) + rot(2,i)*cond(6)
               drot(2,i)= cphi*rot(3,i)*cond(5) - rot(1,i)*cond(6)
               drot(3,i)= -(sphi*rot(1,i)+cphi*rot(2,i))*cond(5)
            end do
            do i=1,3
               drot(i,1)= drot(i,1) - rot(i,2)*cond(4)
               drot(i,2)= drot(i,2) + rot(i,1)*cond(4)
            end do
c
c get corresponding matrices for new reference frame and extract new
c Euler angles and rates
            call PRODCT(rot,nutpr,tmprot, 3,-3,3)
            cth = tmprot(3,3)
            sth = SQRT(tmprot(1,3)**2+tmprot(2,3)**2)
            if(ABS(sth).lt.1.E-16_10) call SUICID(
     . 'DEGENERATE EULER ANGLES, STOP IN CHNCNE ',10)
            sphi= tmprot(1,3)/sth
            cphi= tmprot(2,3)/sth
            spsi= tmprot(3,1)/sth
            cpsi= -tmprot(3,2)/sth
            cond(1) = ATAN2(spsi,cpsi)
            cond(2) = ACOS(cth)
            cond(3) = ATAN2(sphi,cphi)

            call PRODCT(drot,nutpr,tmprot, 3,-3,3)
            cond(4) = (tmprot(3,1)*cpsi+tmprot(3,2)*spsi)/sth
            cond(5) = -tmprot(3,3)/sth
            cond(6) = (tmprot(1,3)*cphi-tmprot(2,3)*sphi)/sth

         else if(icndx.lt.0) then
c
c transform cartesian coordinates
            do i = 1, 6
               x0(i) = cond(i)
            end do
            call PRODCT(nutpr, x0, cond, -3, 3, 2)
         else if(icndx.eq.0) then
c
c transform euler angles
            cinc = cond(3)*Convd
            sinc = SIN(cinc)
            cinc = COS(cinc)
            casc = cond(4)*Convd
            sasc = SIN(casc)
            casc = COS(casc)
            cper = cond(5)*Convd
            sper = SIN(cper)
            cper = COS(cper)
            ynorm(1) = sasc*sinc
            ynorm(2) = -casc*sinc
            ynorm(3) = cinc
            yper(1)  = casc*cper - sasc*sper*cinc
            yper(2)  = sasc*cper + casc*sper*cinc
            yper(3)  = sper*sinc
            call PRODCT(nutpr, ynorm, x0, -3, 3, 1)
            cinc = x0(3)
            sinc = SQRT(1._10 - cinc**2)
            casc = -x0(2)/sinc
            sasc = x0(1)/sinc
            xinc = ACOS(x0(3))/Convd
            xasc = ATAN2(x0(1), -x0(2))/Convd
            if(xasc.lt.0._10) xasc = xasc + 360._10
            call PRODCT(nutpr, yper, x0, -3, 3, 1)
            cper = casc*x0(1) + sasc*x0(2)
            sper = cinc*(-sasc*x0(1) + casc*x0(2)) + sinc*x0(3)
            xper = ATAN2(sper, cper)/Convd
            if(xper.lt.0._10) xper = xper + 360._10
            cond(3) = xinc
            cond(4) = xasc
            cond(5) = xper
         else
c
c transform adbarv to cartessian before applying precession-nt
            call RADVAF(cond, x0, dum, gaussx, 0)
            icndx = -1
            call PRODCT(nutpr, x0, cond, -3, 3, 2)
         endif
      endif
      if(nplnt.lt.0) return
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c shall we change from orbital elements to cartesian coord.
      if(icndx.ne.icndt) then
         if(icndx.eq.0) then
            call JNITL(gaussx,cond,elpt,0,dum)
            call JLIPT(0._10,elpt,1,cond,r,r2,r3,dum)
            icndx = -1
         else if(icndx.eq.2) then

c convert adbarv to cartesian
            call RADVAF(cond,x0,dum,gaussx,0)
            do i = 1, 6
               cond(i) = x0(i)
            end do
            icndx = -1
         endif
c now elements are cartesian, possibly need another type
c (code merged from subroutine chncnd)
         if(icndx.ne.icndt) then
            do i = 1, 6
               x0(i) = cond(i)
            end do
            if(icndt.ne.0) then

c convert to adbarv
               cond(1) = SQRT(DOT(x0,x0))
               cond(2) = ATAN2(x0(2), x0(1))/Convd
               if(cond(2).lt.0._10) cond(2) = cond(2) + 360._10
               cond(3) = ASIN(x0(3)/cond(1))/Convd
               v = SQRT(DOT(x0(4),x0(4)))
               cond(4) = v - gaussx/SQRT(cond(1))
               do i = 1, 3
                  ynorm(i) = x0(i)/cond(1)
               end do
               yper(1) = -ynorm(3)*ynorm(1)
               yper(2) = -ynorm(3)*ynorm(2)
               yper(3) = 1._10 - ynorm(3)**2
               vn = DOT(x0(4), ynorm)
               do i = 1, 3
                  vp(i) = x0(i + 3) - vn*ynorm(i)
               end do
               call CROSS(vp, yper, vpep)
               cond(5) = ATAN2(DOT(vpep,ynorm), DOT(vp,yper))/Convd
               if(cond(5).lt.0._10) cond(5) = cond(5) + 360._10
               cond(6) = ACOS(DOT(ynorm,x0(4))/v)/Convd
               icndx   = 2
            else

c convert to elliptic
               call CHNCNC(gaussx, x0, cond)
               icndx = 0
            endif
         endif
      endif

      if(icndx.eq.0 .and. icnd.eq.1) then
         cond(5) = cond(4) + cond(5)
         if(cond(5).ge.360._10) cond(5) = cond(5) - 360._10
         cond(6) = cond(5) + cond(6)
         if(cond(6).ge.360._10) cond(6) = cond(6) - 360._10
         icndx = 1
      endif
  100 if(icndx.ne.icnd) call SUICID(
     .   'LOGIC ERROR IN INITIAL CONDITION CONVERSION, STOP IN CHNCNE '
     .   , 15)
      return
      end
