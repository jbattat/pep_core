      subroutine CHNCNC(goose,x,cond)
 
      implicit none
 
 
c*** start of declarations inserted by spag
 
c*** end of declarations inserted by spag
 
 
c
c ash / amuchastegui - august 1970 - subroutine chncnc
c initial set-up for conic section motion from position and velocity
c revised sept 1971 for use in pep
c optimized 1977 jan  j.f.chandler
c
c arguments
      real*10 goose,x(6),cond(6)

c goose       = input value of square root of gravitational constant
c               times mass of central body. units are l**3/2*(t**-1)
c x(1-6)      = input cartesian coordinates
c     cond(1) = a    = semi-major axis
c     cond(2) = e    = eccentricity
c     cond(3) = inc  = inclination
c     cond(4) = asc  = right ascension of ascending node
c     cond(5) = per  = argument of perigee
c     cond(6) = anom = initial mean anomaly
c
c common
      include 'funcon.inc'

c local
      real*10 rcv(3),inc,mu,y(2),dy(2)
      real*10 a,a2rt,absa,anom,asc,casc,cecc,cinc,cper,cpsi,
     . e,e2,ecc,g,g2,p,p2rt,per,psr,quan3,quan4,quan5,r,r2,
     . rdot,sasc,secc,sinc,sper,spsi,stheta,tpsi,v,v2

c external functions
      real*10 DOT

c the inverse hyperbolic sine function is defined in some compilers, but
c not in others
      real*10 XASINH,xdu
      XASINH(xdu) = LOG(xdu+SQRT(1._10+xdu**2))
 
c
      mu     = goose**2
      call CROSS(x(1), x(4), rcv)
      g2     = DOT(rcv(1), rcv(1))
      g      = SQRT(g2)
      p      = g2/mu
      p2rt   = g/goose
      quan5  = g
      r2     = DOT(x(1), x(1))
      r      = SQRT(r2)
      v2     = DOT(x(4), x(4))
      v      = SQRT(v2)
      stheta = g/r/v
c
c straight line motion
      if(stheta.gt.1.E-14_10) goto 100
      e = 1E75_10
      cond(1) = r
      cond(2) = e
      cond(3) = ATAN2(x(2), x(1))/Convd
      if(cond(3).lt.0._10) cond(3) = cond(3) + 360._10
      cond(4) = ASIN(x(3)/cond(1))/Convd
      cond(5) = v
      cond(6) = 0._10
      call SUICID('STRAIGHT LINE MOTION NOT ALLOWED, STOP IN CHNCNC',12)
      return
 
c
c calculate inclination and ascending node
  100 cinc = rcv(3)/g
      sinc = 1._10 - cinc**2
      if(sinc.le.0._10) cinc = SIGN(1._10, cinc)
      inc  = ACOS(cinc)
      casc = -rcv(2)
      sasc = rcv(1)
      if(casc.ne.0._10 .and. sasc.ne.0._10) then
         asc = ATAN2(sasc, casc)
         if(asc.lt.0._10) asc = asc + Twopi
      else
         asc = 0._10
      endif
c
c test a
      a = 2._10/r - v2/mu
      if(ABS(a).lt.1.E-14_10) then
c
c parabolic motion
         a     = 1E75_10
         e     = 1._10
         psr   = p/r
         spsi  = SQRT((2._10-psr)*psr)
         cpsi  = psr - 1._10
         y(1)  = r*cpsi
         y(2)  = r*spsi
         rdot  = goose*spsi/p2rt
         dy(1) = -rdot
         dy(2) = rdot*spsi + quan5/r*cpsi
         tpsi  = spsi/psr
         anom  = ((tpsi**3)/3._10 + tpsi)*0.5_10
      else
         e2   = 1._10 - p*a
         a    = 1._10/a
         absa = ABS(a)
         a2rt = SQRT(absa)
c
c elliptic or hyperbolic motion
         if(e2.le.1.E-16_10) then
 
c circular motion
            e     = 0._10
            e2    = 0._10
            y(1)  = r
            y(2)  = 0._10
            dy(1) = 0._10
            dy(2) = v
            anom  = 0._10
         endif
 
c ellipse or hyperbola
         quan3 = p2rt*a2rt
         quan4 = a2rt*goose
         if(e2.gt.0._10) then
c bug in SQRT function clobbers padding
            e     = SQRT(e2)+0._10
            secc  = DOT(x(1),x(4))/(quan4*e)
            cecc  = (1._10 - r/a)/e
            y(1)  = a*(cecc - e)
            y(2)  = quan3*secc
            dy(1) = -quan4*secc/r
            dy(2) = quan5*cecc/r
            if(a.lt.0._10) then
c
c hyperbolic motion
               ecc  = XASINH(secc)
               anom = (e*secc - ecc)
            else
               ecc  = ATAN2(secc, cecc)
               anom = (ecc - e*secc)
               if(anom.lt.0._10) anom = anom + Twopi
            endif
         endif
      endif
c
c determines argument of perigee
      if(sinc.gt.0._10) then
         sper = x(3)*dy(2) - x(6)*y(2)
         cper = x(6)*y(1) - x(3)*dy(1)
      else
         sper = x(1)*dy(1) - x(4)*y(1)
         cper = x(1)*dy(2) - x(4)*y(2)
      endif
      per = ATAN2(sper, cper)
      if(per.lt.0._10) per = per + Twopi
c
c
c move orbital elements into output vector
      cond(1) = a
      cond(2) = e
      cond(3) = inc/Convd
      cond(4) = asc/Convd
      cond(5) = per/Convd
      cond(6) = anom/Convd
 
      return
      end
