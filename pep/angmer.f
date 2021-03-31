      subroutine ANGMER(kathy)
 
      implicit none

c j.f.chandler  -  1979 february  -  subroutine angmer
c calculation of theoretical value of right ascension and/or
c declination for optical observations of sun,moon,planets,stars

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
      real*10 dist
      equivalence (dist,Dstf(4))
      integer*2 itime,ntime
      equivalence (Kob(17),itime),(Kob(18),ntime)
c nice=-1 observation of right ascension only
c nice= 0 observation of right ascension and declination
c nice= 1 observation of declination only
      include 'jdfrequv.inc'
      include 'param.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'

c external functions
      real*10 DOT

c local
      real*10 vcentr(3)/3*0._10/
      real*10 alphat,deltat,eist,hrang,rr2,rsq
      integer   i,iter,j,ksu,lsw,lswem,mm
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=1000
c
c begin iteration to get observed body over meridian
      hrang = 0._10
      lswem = -1
      lsw   = 1
      iter  = 0
  100 iter  = iter + 1
      if(iter.gt.1) call TIMINC(Jds,Fract,Jd,fr(1),Ctut)
      Jdx   = Jd
      fr(2) = fr(1)
      if(Kst1.gt.1) call TIMINC(Jd,fr(1),Jdy,fr(3),-dist/Secday)
c
c determine nutation-precession
      if(ABS(hrang).gt.accprc) then
         call PRCNUT(Jd,fr(1))
         Sidtm = Sidtm0 + Dgst - Longr(1)
         hrang = 0._10
      endif
c
c find type of central body
      if(Klan.eq.u_mxpl+1) then
 
c moon or lunar orbiter
         fr(2) = fr(3)
         Jdx   = Jdy
         call MNTRP(1,Jdx,fr(2),0,lswem,1)
         if(Jdx.le.0) then
            Jd = 0
            return
         else if(Mnau.ne.1._10) then
            do i = 1,3
               Xm(i,1) = Xm(i,1)*Mnau
            end do
         endif
c
c need earth position,velocity for trans-lunar object
      else if(Kst1.ne.2) then
 
c already done once in angctl
         if(iter.ne.1) then
            call MNTRP(1,Jdx,fr(2),0,lswem,1)
            if(Jdx.le.0) then
c
c*  start=9000
               Jd = 0
               return
            else
               call MNTRP(1,Jdx,fr(2),-1,0,1)
c
c read earth-moon barycenter tape, perform interpolation
c to determine sun relative to earth (vel. vice-versa)
               call EMTRP(1,Jd,fr(1),0,lswem,1)
               if(Jd.le.0) then
                  Jd = 0
                  return
               else
                  call EMTRP(1,Jd,fr(1),-1,0,1)
                  do j = 1,3
                     x(j,2)     = Mnfct*Xm(j,1) - Xem(j,1)
                     x(j + 3,2) = -Mnfct*Xm(j + 3,1) + Xem(j + 3,1)
                     Ves(j)      = (x(j+3,2) - vcentr(j))*Aultvl
                  end do
               endif
            endif
         endif
 
         if(Klan.gt.0 .and. Klan.le.u_mxpl) then
c
c read planet tape, interpolate to get planet position
            call PLTRP(1,Jdy,fr(3),0,lsw)
            if(Jdy.le.0) then
               Jd = 0
               return
            else
 
c get planet relative to earth
               do j = 1,3
                  x(j,4) = Xp(j) + x(j,2)
               end do
            endif
         endif
      endif
c
c*  start=1500
c determine if probe or satellite is observed body
      if(Klanb.gt.0) then
c
c determine satellite position, velocity
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
c*  start=1800
c aberration correction
      r = SQRT(DOT(x(1,Kst),x(1,Kst)))
 
c no correction for cis-lunar object
      if(Kst1.ne.2) then
         do j = 1,3
            x(j,Kst) = x(j,Kst) + r*Ves(j)
         end do
      endif
c
c determine vector pointing to observed body referred to true
c equinox and equator of date
      call CORCHN(xo,x(1,Kst))
c
c determine theoretical value of right ascension
      alphat = ATAN2(xo(2),xo(1))
c
c determine time correction for meridian circle observation
      hrang = Sidtm + Sidvel*Utrec - alphat
      do while( hrang.lt.-Pi )
         hrang = hrang + Twopi
      end do
      do while( hrang.ge.Pi )
         hrang = hrang - Twopi
      end do
      Utrec = Utrec - hrang/Sidvel
      Fract = Utrec/Secday
      Jd    = Jds
c
c see if there is to be further iteration
      if(Kst1.ne.1) then
         r = r*Aultsc
         if(ABS(dist-r).gt.accdst) goto 200
      endif
      if(ABS(hrang).le.acctim) goto 300
c
c set up quantities for next iteration
  200 dist = r
      lsw  = -1
      if(iter.lt.100) goto 100
      call SUICID(
     .    ' MORE THAN 100 ITERATIONS NEEDED TO GET BODY OVER MERIDIAN  '
     .    ,15)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c*  start=2000
c           termination calculations
  300 Nit(20) = Nit(20) + iter
 
      dist = r
      rsq  = DOT(xo,xo)
      r    = SQRT(rsq)
      do i=1,3
         Xsitp0(i,1) = -xo(i)/r
      end do
 
      mm   = Utrec
      eist = Utrec - mm
      Ihr  = mm/3600
      mm   = mm - Ihr*3600
      Imin = mm/60
      Sec  = mm - Imin*60 + eist
 
      rr2 = xo(1)**2 + xo(2)**2
      Angdum(10) = SQRT(rr2)
      rr1 = rr2*Convhs
 
c need velocity for angle rates
      if(Jct(39).ge.0) then
         if(Ict(3).lt.0 .and. Idumob.eq.1) goto 400
         if(Ict(1).le.0) goto 400
      endif
      fake = xo(3)/rsq
      cake = Angdum(10)*Convds
 
c need planet velocity for partials
      call ANGVP
 
c get velocity in system of date
      do i = 1,3
         xop(i) = Xsitep(i+3,1)
      end do
      call CORCHN(Xsitep(4,1),xop(1))
      Angdum(8) = rr2*Sidvel + (xo(1)*Xsitep(5,1) - xo(2)*Xsitep(4,1))
c*  start=2500
c calculate observables
c
c calculate angle rates
  400 if(Jct(39).ge.0) then
c
c correct for star catalog errors
         deltat = ATAN2(xo(3),Angdum(10))
         call SKYCOR(alphat,deltat)
c
c calculate hour angle instead of right ascension
         if(Jct(39).le.1) then
c
c get right ascension between 0 and 360 degrees
            if(alphat.lt.0._10) alphat = Twopi + alphat
         else
            alphat = Twopi - (alphat - (Sidtm+Sidvel*Utrec))
            mm     = alphat/Twopi
            if(alphat.lt.0._10) mm = mm - 1
            alphat = alphat - mm*Twopi
         endif
         Deriv(2,1) = alphat/Convhs
c
c calculate declination
         Deriv(2,2) = deltat/Convds
c equator-equinox and phase corrections on right ascension
c and declination (calculate following for observation
c library tape even if neqnox=0)
         Salph = SIN(alphat)
         Calph = COS(alphat)
         Tdelt = TAN(deltat)
 
c is this a meridian circle observation with limb correction
         if(ntime.gt.0) call LMBCOR(alphat,deltat)
         if(Neqnox.gt.0 .or. Nphase.gt.0) call EQPCOR(kathy)
      else
         Deriv(2,1) = (xo(2)*Xsitep(4,1)-xo(1)*Xsitep(5,1))/rr1
         Angdum(9)  = -DOT(xo,Xsitep(4,1))
         Deriv(2,2) = (-Xsitep(6,1) - xo(3)/rsq*Angdum(9))/Angdum(10)
     .    /Convds
      endif
 
c*  start=9990
      return
      end
