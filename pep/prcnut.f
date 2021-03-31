      subroutine PRCNUT(jd, fract)
 
      implicit none

c
c m.e.ash   nov 1966     subroutine prcnut
c the nutation-precession matrix is calculated
c
 
c parameters
      integer*4 jd
      real*10 fract
c jd    =julian day number
c fract = fraction of day from midnight

c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'argfun.inc'
      include 'comdateq.inc'
      include 'comnut.inc'
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'funcon.inc'
      include 'leon.inc'
      include 'nutprc.inc'
      real*4    dpsi(2), deps(2)
      equivalence (dpsi, Nutat), (deps, Nutat(3))
c
c           kind  =0 nutation-precession determined
c           kind  =1 nutation-precession and its derivative determined
c           nutpr =nutation times precession matrices
c           dnutpr=derivative of nutation times precession matrices
c           prec  =precession matrix
c           dprec =derivative of precession matrix
c           moblq,dmoblq=mean obliquity of ecliptic and its derivative
c           dpsi  =nutation in longitude and its derivative (radians)
c           deps  =nutation in obliquity and its derivative (radians)
c           (unit of time for derivatives is the second)
c
c
c   coefficients for tidal corrections to nutation
      real*10 arg(7), scof(7)/-2.085E-8_10, -19.926E-8_10, 9.357E-8_10,
     . 21.380E-8_10, -.8242E-8_10, 1.745E-8_10, -.1454E-8_10/,
     . ccof(7)/1.066E-8_10, 9.793E-8_10, 2.715E-8_10, -4.169E-8_10,
     . 0.388E-8_10, 0.1939E-8_10, .097E-8_10/
c
c coefficients for everett interpolation
      real*10 evcf(9)/1.5333333333333333333_10, -0.3_10,
     . 3.3333333333333333333E-2_10, -5.8333333333333333333E-1_10,
     . 3.3333333333333333333E-1_10, -4.1666666666666666667E-2_10,
     . 0.05_10, -3.3333333333333333333E-2_10,
     . 8.3333333333333333333E-3_10/
      integer*2 kcof(4)/9, 10, 33, 34/, ntab, nr
c
c internal to this subroutine
      real*10 corj2,dnutm(3,3),f1,f2,fract1,ndp(3,3),npsig(3),
     . nut(3,3),oblq,p(4),sig(3),t,tjd,wobarg,xnp,xoz,y(3,4),ynp
      integer i,ir,j,jct28,jd1,k,k1
      real*4 xws,yws
      logical iprnt

c external functions
      real*10 ARGMNT,SXY06

c arithmetic statement function
      real*10 DBL10
      real*4 x4
      DBL10(x4)=x4

      iprnt = mod(Jct(6)/128,2) .ne. 0
c
c adjust time quantities
      i = fract
      if(fract.lt.0._10) i = i - 1
      fract1 = fract - i
      jd1    = jd + i
      t    = jd1
      p(1) = fract1 - 0.5_10
      t    = t + p(1)
c
c determine precession matrix and mean obliquity
      call PRECES(t)
c
c determine nutation p vector and tabular indices
      ntab = 2*Inter*(jd1 - Jder(2))
      if(p(1).lt.0._10) then
         p(1) = fract1
      else
         ntab = ntab + Inter
      endif
      p(1) = p(1)/Dnut
      if(p(1).lt.0._10) then
         p(1) = p(1) + 1._10
         ntab = ntab - 1
      endif
      p(3) = 1._10 - p(1)
      p(2) = p(1)**2
      p(4) = p(3)**2
c
c determine nutation y vectors
      do k = 1, 4
         nr = ntab + kcof(k)
         f1 = DBL10(Psidx(nr+1)) + DBL10(Psidx(nr-1))
         f2 = DBL10(Psidx(nr+2)) + DBL10(Psidx(nr-2))
         do i = 1, 3
            j = i*3 - 2
            y(i, k) = evcf(j)*DBL10(Psidx(nr))
     .                + (evcf(j+1)*f1 + evcf(j+2)*f2)
         end do
      end do
c
c everett fourth difference interpolation for nutation
      do k = 1, 3, 2
         k1 = k + 1
         Nutat(k) = p(3)*(y(1,k) + p(4)*(y(2,k)+p(4)*y(3,k))) + p(1)
     .              *(y(1,k1) + p(2)*(y(2,k1)+p(2)*y(3,k1)))
         if(Kindnp.gt.0) Nutat(k1)
     .       = (y(1,k1) + p(2)*(y(2,k1)*3._10+p(2)*y(3,k1)*5._10)
     .       - (y(1,k)+p(4)*(y(2,k)*3._10+p(4)*y(3,k)*5._10)))/Dnut
      end do
      if(Jct(21).eq.2) then
c correct nutations assumed to be IAU 2000 model to be consistent with
c IAU 2006 precession
         corj2 = -2.7774e-6_10*((jd-2451545)+(fract-0.5_10))/36525._10
         dpsi(1) = dpsi(1) + dpsi(1)*(corj2+0.4697e-6_10)
         deps(1) = deps(1) + corj2*deps(1)
      endif
c
c           add corrections to nutation theory to account for the
c           earth's elasticity and fluid core
c             if mod(jct28,2).ne.0 apply melchior's corrections
c                 to woolards's rigid-body theory
c             if jct(29).gt.0 apply corrections input as con(3-18)
c                 note:  these two options are not mutually exclusive
      jct28 = Jct(28)
      if(mod(jct28,2).ne.0 .or. Jct(29).gt.0) then
         tjd = jd + fract - .5_10
 
c get angles in brown lunar theory (units are revolutions)
         call FUNARG(tjd)
         arg(1) = ARGMNT(2.*F + 2.*Ascm)
         arg(2) = ARGMNT(2.*F - 2.*D + 2.*Ascm)
         arg(3) = ARGMNT(Lp)
         arg(4) = ARGMNT(Ascm)
         arg(5) = ARGMNT(Lp + 2.*F - 2.*D + 2.*Ascm)
         arg(6) = ARGMNT(L)
         arg(7) = ARGMNT(L + 2.*F + 2.*Ascm)
         do i = 1, 7
            Carg(i) = COS(arg(i))
            Sarg(i) = SIN(arg(i))
         end do
         if(mod(jct28,2).ne.0) then
            do i = 1, 7
               deps(1) = deps(1) + ccof(i)*Carg(i)
               dpsi(1) = dpsi(1) + scof(i)*Sarg(i)
            end do
            if(Jct(29).le.0) go to 100
         endif
         do i = 9, 15, 2
            ir = i/2 - 3
            deps(1) = deps(1) + Ercond(i)*Convds*Carg(ir)
         end do
         do i = 10, 16, 2
            ir = i/2 - 4
            dpsi(1) = dpsi(1) + Ercond(i)*Convds*Sarg(ir)
         end do
 
c out-of-phase terms
         do i = 17, 23, 2
            ir = i/2 - 7
            deps(1) = deps(1) + Ercond(i)*Convds*Sarg(ir)
         end do
         do i = 18, 24, 2
            ir = i/2 - 8
            dpsi(1) = dpsi(1) + Ercond(i)*Convds*Carg(ir)
         end do
      endif
 
c add free core nutation
  100 if(Jct(29).ne.0 .and. Ercom(2).ne.0._10) then
         wobarg = MOD(t, Ercom(2))/Ercom(2)*Twopi
 
c ercom(2) is the period of the free mode, in days
         Cfnut   = COS(wobarg)
         Sfnut   = SIN(wobarg)
         deps(1) = deps(1) + (Ercond(7)*Cfnut - Ercond(8)*Sfnut)*Convds
         dpsi(1) = dpsi(1) + (Ercond(7)*Sfnut + Ercond(8)*Cfnut)
     .             /SIN(Moblq)*Convds
      endif
c
c determine pc,ps
      deps(1) = deps(1)*Nutn/Nutn0
      dpsi(1) = dpsi(1)*Nutn/Nutn0
      oblq    = Moblq
      Cobliq  = COS(oblq)
      Sobliq  = SIN(oblq)
      Pc(1)   = Cobliq*Nutat(1)
      Ps(1)   = Sobliq*Nutat(1)
      if(Kindnp.gt.0) then
         deps(2) = deps(2)/8.64E4*Nutn/Nutn0
         dpsi(2) = dpsi(2)/8.64E4*Nutn/Nutn0
         Ps(2)   = Dmoblq + deps(2)
         Pc(2)   = Cobliq*dpsi(2) - Ps(1)*Ps(2)
         Ps(2)   = Sobliq*dpsi(2) + Pc(1)*Ps(2)
      endif
c
c calculate nutation matrix
      nut(1,1) = 1._10 - 0.5_10*dpsi(1)**2
      nut(1,2) = -Pc(1)
      nut(1,3) = -Ps(1)
      nut(2,1) = Pc(1) - deps(1)*Ps(1)
      nut(2,2) = 1._10 - 0.5_10*(deps(1)**2 + Pc(1)**2)
      nut(2,3) = -deps(1) - 0.5_10*Pc(1)*Ps(1)
      nut(3,1) = Ps(1) + deps(1)*Pc(1)
      nut(3,2) = deps(1) - 0.5_10*Pc(1)*Ps(1)
      nut(3,3) = 1._10 - 0.5_10*(deps(1)**2 + Ps(1)**2)
c
c calculate partials of nutation matrix w.r.t. deps and dpsi
      Dnutde(1,1) = 0._10
      Dnutde(1,2) = 0._10
      Dnutde(1,3) = 0._10
      Dnutde(2,1) = -Ps(1)
      Dnutde(2,2) = -deps(1)
      Dnutde(2,3) = -1._10
      Dnutde(3,1) = Pc(1)
      Dnutde(3,2) = 1._10
      Dnutde(3,3) = -deps(1)
c
      Dnutdp(1,1) = -dpsi(1)
      Dnutdp(1,2) = -Cobliq
      Dnutdp(1,3) = -Sobliq
      Dnutdp(2,1) = Cobliq - deps(1)*Sobliq
      Dnutdp(2,2) = -Pc(1)*Cobliq
      Dnutdp(2,3) = -Ps(1)*Cobliq
      Dnutdp(3,1) = Sobliq + deps(1)*Cobliq
      Dnutdp(3,2) = -Pc(1)*Sobliq
      Dnutdp(3,3) = -Ps(1)*Sobliq
c
c print nutation matrix
      if(mod(Jct(6)/512,2).ne.0) then
         if(Line.gt.54) call OBSPAG
         write(Iout,150) ((nut(i,j),j=1,3), i = 1, 3)
  150    format(' PRCNUT:     NUT=', (t18,3F20.16))
         Line = Line + 3
      endif
c
c calculate product of nutation and precession matrices
      call PRODCT(nut, Prec, Nutpr, 3, 3, 3)
c
c calculate derivative of nutprc
      if(Kindnp.gt.0) then
         dnutm(1, 1) = 0._10
         dnutm(1, 2) = -Pc(2)
         dnutm(1, 3) = -Ps(2)
         dnutm(2, 1) = Pc(2)
         dnutm(2, 2) = 0._10
         dnutm(2, 3) = -deps(2)
         dnutm(3, 1) = Ps(2)
         dnutm(3, 2) = deps(2)
         dnutm(3, 3) = 0._10
         call PRODCT(dnutm, Prec, Dnutpr, 3, 3, 3)
         call PRODCT(nut, Dprec, ndp, 3, 3, 3)
         do i = 1, 3
            do j = 1, 3
               Dnutpr(i, j) = Dnutpr(i, j) + ndp(i, j)
            end do
         end do
c
c print time derivative of precession/nutation matrix
         if(mod(Jct(6)/512,2).ne.0) then
            if(Line.gt.54) call OBSPAG
            write(Iout,155) ((Dnutpr(i,j),j=1,3), i = 1, 3)
  155       format(' PRCNUT:  DNUTPR=', (t18,3F20.16))
            Line = Line + 3
         endif
      endif
c
c determine earth wobble
      if(Iwob.gt.0) then
         call WOBBLE(jd, fract, xws, yws)
         Xwob = xws*Convds
         Ywob = yws*Convds
         if(iprnt) then
            if(Line.gt.56) call OBSPAG
            write(Iout, 160) Xwob, Ywob, xws, yws
  160       format(' PRCNUT: XWOB,YWOB =', 1p, 2E16.7, ' (RAD)   =',
     .             2E16.7, ' (ARCSEC)')
            Line = Line + 1
         endif
      endif

c Determine correction for Greenwich sidereal time.
c For old models, use just the projection of nutation in longitude to
c correct mean to true sidereal time.
c For 2000 and later models, we use GST-GMST = -(ERA-GST)-(GMST-ERA)
c where ERA-GST= Equation of origins, and GMST-ERA was saved earlier
c when the mean sidereal time at midnight was computed.   See
c Wallace & Capitaine A&A 2006 and Capitaine, Wallace, & Chapront 2003
      if(Jct(21).le.0) then
         Dgst = Pc(1)
      else
         xnp = Nutpr(3,1)
         ynp = Nutpr(3,2)
         xoz = xnp/(1._10+Nutpr(3,3))
         sig(1) = 1._10-xnp*xoz
         sig(2) = -ynp*xoz
         sig(3) = -xnp
         call PRODCT(Nutpr,sig,npsig,3,3,1)
         Dgst = -SXY06(jd,fract)+xnp*ynp*0.5_10 - Dera
         if(npsig(1).ne.0._10 .or. npsig(2).ne.0._10)
     .    Dgst=Dgst+ATAN2(npsig(2),npsig(1))
      endif

      if(iprnt) then
         if(Line.gt.56) call OBSPAG
         write(Iout, 200) Nutat(1), Nutat(3)
  200    format(' PRCNUT: PSI,EPS (RADIANS) =', 1p, 2E18.7)
         Line = Line + 1
      endif
      return
      end
