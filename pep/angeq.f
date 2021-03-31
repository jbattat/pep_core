      subroutine ANGEQ(kathy)
 
      implicit none

c j.f.chandler - 1979 february - subroutine angeq
c compute photographic observables
c derived from subr. angle

c arguments
      integer*4 kathy
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
      include 'mnsprt.inc'
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
      include 'yvect.inc'

c external functions
      real*10 DOT

c local
      real*10 alphat,deltat,rltsc,rr2,rsq,xstp(3)
      integer   i,iter,j,ksu,lsw,lswm,mm
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c*  start=1000
c           begin iteration to make light time correction
      lsw  = 1
      lswm = -1
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
         else if(Nplnt0.eq.-4) then
            call SPOTCD(Jd,fr(1),0,0,Nplnt0,0._10,0._10,1)
            do j=1,3
               x(j,3)=Xspcd(j,1)
            end do
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
 
c topocentric correction
         if(Nplnt0.ne.-4) then
            do j = 1,3
               x(j,Kst) = x(j,Kst) - Xsitau(j)
            end do
         endif
c
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
     .        'MORE THAN 20 LIGHT TIME ITERATIONS NEEDED, STOP IN ANGEQ'
     .        ,14)
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
c for non-meridian circle observation
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
      if(Kst1.ne.2) then
         if(kathy.lt.2) then
            do j = 1,3
               x(j,Kst) = x(j,Kst) + r*Ves(j)
            end do
         endif
      endif
c determine vector pointing to observed body referred to mean
c equinox and equator of reference epoch or true of date
      if(Jct(39).gt.0) then
         call CORCHN(xo,x(1,Kst))
      else
c jct(39).gt.0 means topocentric right ascension-declination
c referred to the true equinox and equator of date for dummy
c observation prediction
         do i = 1,3
            xo(i) = x(i,Kst)
         end do
      endif
      rsq = DOT(xo,xo)
      r   = SQRT(rsq)
      rr2 = xo(1)**2 + xo(2)**2
      rr1 = rr2*Convhs
      Angdum(10) = SQRT(rr2)
      do i=1,3
         Xsitp0(i,1) = -xo(i)/r
      end do
c
c determine theoretical value of right ascension
      alphat = ATAN2(xo(2),xo(1))
c
c*  start=2200
      if(Jct(39).ge.0) then
         if(Ict(3).lt.0 .and. Idumob.eq.1) goto 300
         if(Ict(1).le.0) goto 300
      endif
      fake = xo(3)/rsq
      cake = Angdum(10)*Convds
c
c need planet velocity for partials or angle rates
      call ANGVP
c
c*  start=2500
c calculate angle rates
      if(Jct(39).lt.0) then
         Deriv(2,1) = (xo(2)*Xsitep(4,1)-xo(1)*Xsitep(5,1))/rr1
         Angdum(9)  = -DOT(xo,Xsitep(4,1))
         Deriv(2,2) = -(Xsitep(6,1)+xo(3)/rsq*Angdum(9))/Angdum(10)
     .    /Convds
         return
      endif
c
c correct for star catalog errors
  300 deltat = ATAN2(xo(3),Angdum(10))
      call SKYCOR(alphat,deltat)
c
c calculate hour angle instead of right ascension
      if(Jct(39).gt.1) then
         alphat = Twopi - (alphat - (Sidtm+Sidvel*Utrec))
         mm     = alphat/Twopi
         if(alphat.lt.0._10) mm = mm - 1
         alphat = alphat - mm*Twopi
c
c get right ascension between 0 and 360 degrees
      else if(alphat.lt.0) then
         alphat = Twopi + alphat
      endif
      Deriv(2,1) = alphat/Convhs
c
c calculate declination
      Deriv(2,2) = deltat/Convds
c
c equator-equinox and phase corrections on right ascension
c and declination (calculate following for observation
c library tape even if neqnox=0)
      Salph = SIN(alphat)
      Calph = COS(alphat)
      Tdelt = TAN(deltat)
      if(Neqnox.gt.0 .or. Nphase.gt.0) call EQPCOR(kathy)

      call SPOTCV(0,Nplnt0,1)

c*  start=9990
      return
      end
