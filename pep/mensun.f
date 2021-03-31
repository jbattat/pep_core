      subroutine MENSUN(jd, fract, xsun, k)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      a, botm, capt, chck, cnvdr, d, dele, delg, diff, ecanm,
     .          ecc, ecdot, edot, epdot, eps, etry, g, gtry, h
      integer   i, ij, j, jd, k, kj, m, n
      real      r, ratio, rdot, siv, t, top, v, vadot, xlam, xldot, xn,
     .          xovr, yovr, zovr
 
c*** end of declarations inserted by spag
 
 
c
c     m. slade  feb. 1970           subroutine mensun
c     computation of approximate  solar position in rectangular
c     coordinates (equator and equinox of ref.epoch) for use in phase
c     corrections and  star occultations by moon
c
c
      real*10 fract
      real*10 delta
      real*4    xdate(3)
      real*4    vdate(3)
      real*10 xsun(6)
      real*10 gg, xl, dt, dd, xkj
 
      include 'nutprc.inc'
 
      data cnvdr/1.745329E-2/
      data a/1.0000002/
c
c     if k=-1, velocity only
c     if k=0, position only
c     if k=1, position and velocity
c     distances in a.u.
c     velocity in  a.u./ephemeris days
c
 
      eps = Moblq
 
      dt = (jd - 2415020) + (fract - 0.5_10)
      t  = dt
      dd = dt*1E-4_10
      d  = dd
c
c see p. 98 explanatory supplement
c
      ecc   = 0.01675104-1.1444E-5*d-9.4E-9*d**2
      xl    = 279.696678+0.98564733543*dt+2.267E-5_10*dd**2
      gg    = 358.475845+0.985600267*dt-1.12E-5_10*dd**2-7E-8_10*dd**3
      kj    = gg/360._10
      xkj   = kj
      gg    = gg - 360._10*xkj
      ij    = xl/360._10
      xkj   = ij
      xl    = xl - 360._10*xkj
      delta = xl - gg
      diff  = delta*cnvdr
      g     = gg
      g     = g*cnvdr
      xlam  = xl*cnvdr
 
c iterate to solve kepler's equation
      n    = 1
      etry = g + ecc*SIN(g)
      do while( .true. )
         gtry = etry - ecc*SIN(etry)
         delg = g - gtry
         h    = abs(g) - 0.01
         if(h.lt.0.) ratio = delg
         if(h.ge.0.) ratio = delg/g
         ratio = abs(ratio)
         if(ratio.ge.1.E-6) then
            if(n.gt.100) call SUICID(
     .' ITERATION FOR ECC ANOM  HAS EXCEEDED 100 PASSES-STOP IN MENSUN '
     ., 16)
            dele = delg/(1. - ecc*COS(etry))
            n    = n + 1
            etry = etry + dele
            go to 100
         endif
         ecanm = etry
         r     = a*(1. - ecc*COS(ecanm))
         top   = SQRT(1. - ecc**2)*SIN(ecanm)
         botm  = COS(ecanm) - ecc
         chck  = abs(botm) - 0.01
         if(chck.lt.0.) then
            siv = SQRT(1. - ecc**2)*SIN(ecanm)/(1. - ecc*COS(ecanm))
            v   = ASIN(siv)
         else
            v = ATAN2(top, botm)
         endif
         xlam     = diff + v
         xovr     = COS(xlam)
         yovr     = SIN(xlam)*COS(eps)
         zovr     = SIN(xlam)*SIN(eps)
         xsun(1)  = 0.0
         xsun(2)  = 0.0
         xsun(3)  = 0.0
         xdate(1) = xovr*r
         xdate(2) = yovr*r
         xdate(3) = zovr*r
         if(k.ge.0) then
            do i = 1, 3
               do j = 1, 3
                  xsun(i) = Prec(j, i)*xdate(j) + xsun(i)
                  end do
               end do
         endif
         if(k.ne.0) then
            xsun(4) = 0.0
            xsun(5) = 0.0
            xsun(6) = 0.0
            capt    = d/3.6525
            xn    = 3.141593/(365.25964 + 3.04E-6*capt)
            xn    = xn*2.
            ecdot = -1.1444E-9 - 18.8E-13*d
            epdot = -3.5626E-7 - 2.46E-11*d + 3.09E-12*d**2
            epdot = epdot*cnvdr
            edot  = xn/(1. - ecc*COS(ecanm))
            rdot  = -a*ecdot*COS(ecanm) + a*ecc*SIN(ecanm)*edot
            vadot = xn*(1. + ecc*COS(v))**2/SQRT((1.-ecc**2)**3)
            xldot = 4.70684E-5 + 6.78E-9*d + 2.10E-11*d**2
            xldot = xldot*cnvdr + vadot
            vdate(1) = -r*SIN(xlam)*xldot + xovr*rdot
            vdate(2) = r*(COS(xlam)*COS(eps)*xldot - SIN(xlam)*SIN(eps)
     .                 *epdot) + yovr*rdot
            vdate(3) = r*(COS(xlam)*SIN(eps)*xldot + SIN(xlam)*COS(eps)
     .                 *epdot) + zovr*rdot
            do i = 4, 6
               do j = 4, 6
                  m = i - 3
                  n = j - 3
c
c xsun(i)=prec(n,m)*vdate(n) +dprec(n,m)*xdate(n)   +  xsun(i)
c if needed to include effects  of precession
c
                  xsun(i) = Prec(n, m)*vdate(n) + xsun(i)
                  end do
               end do
         endif
         return
  100    end do
      end
