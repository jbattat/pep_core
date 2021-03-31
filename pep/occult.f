      subroutine OCCULT(ncodet)
 
      implicit none

c     ash/forni  june 1968   subroutine occult
c     revisions:
c     slade 1971,  ng 1976,  chandler 1977
c  type of observation is specified by ncodf and further by ib
c arguments
      integer*4 ncodet
c     ncodet= ncodf - 8
c     ncodet=   -1   occultation
c     ncodet=    0   transit
c     ncodet.ge. 1   photographs as black dot crosses sun
c
c       ib       event                    logical center
c     ib=1 moon occults star                    observer
c     ib=2 planet occults star                  observer
c     ib=3 moon occults planet                  observer
c     ib=4 planet transits sun                  sun
c     ib=5 (not used)
c     ib=6 planet occults planet (mutual event) distant body
c     ib=7 (reserved for black dot photos)
c     ib=8 planet occults planet (mid-time)     observer
c     ib=9 planet eclipses planet (mid-time)    sun
c
c  all coordinates are to be expressed in au, au/day
c  possible exception - moon
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 xe(6)
      equivalence (Xem, xe)
      real*10 xo(3),rr1,r,xz(3),rzsq,rz
      equivalence (Angdum, xo), (Angdum(4), rr1), (Angdum(5), r),
     .            (Angdum(6), xz), (Angdum(9), rzsq),
     .            (Angdum(10), rz)
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'ltrapobs.inc'
      include 'mnsprt.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 xra,xmotra,xdec,xmotdc,prlx,rvel
      character*8 numst
      equivalence (Save(28), xra), (Save(29), xmotra), (Save(30), xdec),
     .            (Save(31), xmotdc), (Save(32), prlx),
     .            (Save(33), rvel), (Save(34), numst)
      real*4 acctim, accdst, accprc
      equivalence (acctim, Estf), (accdst, Estf(2)),
     .            (accprc, Estf(3))
      real*10 dist,distz
      equivalence (dist, Dstf(4))
      equivalence (distz, Dstf(5))
      include 'jdfrequv.inc'
      integer*2 itime
      equivalence (itime, Kob(17))
c nice=-1 observation of crossing of first planet limb
c nice= 0 observation of crossing of first & second planet limb
c nice= 1 observation of crossing of second planet limb
      include 'param.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'trnocc.inc'
      real*10 mrad0, mrad
      equivalence (Rhos, mrad0), (Rhom, mrad)
      include 'watstf.inc'
      include 'yvect.inc'
      real*10 xema(3), xmna(3), xpla(3), xsba(3), xsca(3)
      equivalence (Accm, xmna), (Accem, xema),
     .            (Accp, xpla)
c
c observed body coordinates always placed in xsbs
c note: these coordinates are eventually expressed relative to
c the 'logical center' of the event, which depends on ib (q.v.)
c
c    xo   = vector pointing from earth to observed body
c    rr1  = xo(1)**2+xo(2)**2+xo(3)**2
c    r    = SQRT(rr1)
c    xz = vector pointing to second (distant) body
c    rz = SQRT(xz(1)**2+xz(2)**2+xz(3)**2 )
c
c           local variables
      real*10 cdec,cortim,cs1,ctht,dec,dfdt3r,dtprds,dum,f,fk,fl,fm,
     . hrang,onebd,ra,raab,rehrm,rhsm,rmx,s1,sdec,tf1950,usomu,vex,vmx
      integer*4 i,iso,iter,j,jmn,jpl,lsw,lsw0,mgo,mnspt1,
     . ndi,nrtv,nvs
      real*10 dist0(10)/300._10, 150._10, 0._10, 300._10, 2200._10,
     .       4250._10, 9090._10, 14530._10, 19220._10, 1.28_10/
      real*10 rmp(3)
      real*10 usshat(3)
      real*10 muhat(3), usohat(3), ushat(3, 2)
      equivalence (ushat(1,1), usshat), (ushat(1,2), usohat)
      real*10 rhomn(3)
      integer*2 nlibp2/-1/
 
c unit vector derivatives
      real*10 dreh(3, 2)

c external functions
      real*10 A1UT1,CTUTF,DOT

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c this would be unnecessary if limb correction were programmed
      if(ncodet.lt.0) mrad = mrad0
 
c set up (jd,fr) indices
      jmn = 1
      if(Nplnt0.eq.10) jmn = 2
      jpl  = jmn + 1
      lsw  = -1
      mgo  = 1
      rhsm = Rhos + Rhom
 
c set up dummy spot on moon for rotation routine
      Yspcd(1, 1) = 0._10
      Yspcd(2, 1) = 0._10
      Yspcd(3, 1) = 0._10
      mnspt1 = -1
      nvs    = 1
 
c need site acceleration for mutual occultation
      if(Ibtrn.eq.8) nvs = 2
c*  start=100
c
c first guess at time of transit or occultation
  100 if(Nice.lt.0) then
         if(mgo.eq.1) goto 200
         if(mgo.eq.2) goto 1000
      else if(Nice.eq.0) then
         goto 200
      endif
      if(mgo.eq.1) goto 1000
  200 Utrec = Deriv(2, mgo)
      lsw0  = mgo - 2
      hrang = 1E19_10
      if(dist.le.0._10) dist = Tguess
      if(Nps1.gt.0) then
         if(distz.le.0._10) then
            distz = dist
            ndi   = Nps1
            if(Ncs1.gt.0) ndi = Ncs1
            if(ndi.gt.0 .and. ndi.le.10) then
               if(ndi.ne.Nplnt(Klan)) distz = dist0(ndi)
            endif
         endif
      endif
c
c determine time quantities
      Fract = Utrec/Secday
      if(itime.ge.0) then
         Ctat = 32.15_10
         if(itime.gt.1) goto 300
         if(Jds.ge.2435490) then
            Atuts = A1UT1(Jds, Fract)
            goto 300
         endif
      endif
      Ctat  = CTUTF(Jds, Fract)
      Atuts = 0._10
  300 Ctut  = (Ctat + Atuts)/Secday
c
c determine sidereal time quantities
      call SIDTIM(Jds,Ctat+Atuts,Sidtm0,Sidvel,Dera)
      call PLATMO(Jds)
 
c read moon tape for precession-nutation values
      Jd = Jds
      call MNREED(Jd)
      if(Jd.le.0) goto 1200
      iter = 0
      if(Ibtrn.le.2) then
 
c set up star direction unit vector
c N.B.  This assumes star positions are given for mean epoch 1950,
c       which is not necessarily true.  More code will be needed if
c       lunar occultations are seriously used.
         tf1950 = Jds + Fract - 0.5_10 - 2433282.423_10
         onebd  = prlx
         onebd  = onebd*(1._10 - onebd*rvel*tf1950)
         ra     = xra + xmotra*tf1950
         dec    = xdec + xmotdc*tf1950
 
c correct for aberration due to local standard of rest
         sdec = SIN(dec)
         cdec = COS(dec)
         raab = 168.75_10*Convd + ra
         dec  = dec + (0.341_10*COS(raab)*sdec + 0.029_10*cdec)*Convds
         ra   = ra + 0.341_10*SIN(raab)*Convds/cdec
         usshat(2) = COS(dec)
         usshat(1) = COS(ra)*usshat(2)
         usshat(2) = SIN(ra)*usshat(2)
         usshat(3) = SIN(dec)
      endif
c*  start=1000
c
c start of iteration to determine time of transit or occultati
  400 iter = iter + 1
      call TIMINC(Jds, Fract, jde(1), fr(1), Ctut)
 
c get precession-nutation matrix
      if(ABS(hrang).ge.accprc) then
         call PRCNUT(jde(1), fr(1))
         Sidtm = Sidtm0 + Dgst
      endif
 
c get site position, velocity, and maybe acceleration
      call SITCOR(Sidtm + Sidvel*Utrec, 1, nvs, 0)
      do i = 1, 3
         j = i + 3
         Xsitau(i) = Xsite(i, 1)/Aultsc
         Xsitau(j) = Xsite(j, 1)/Aultvl
      end do
      call TIMINC(jde(1), fr(1), jde(2), fr(2), -dist/Secday)
      if(Nps1.gt.0) call TIMINC(jde(1), fr(1), jde(3), fr(3),
     .                              -distz/Secday)
 
c get moon, embary position and velocity
      call ETRP(1, jde(1), fr(1), 0, lsw, 1, 1)
      if(jde(1).le.0) goto 1200
      call ETRP(1, jde(1), fr(1), -1, 0, 1, 1)
c
c need observer acceleration for mutual occultation
      if(Ibtrn.eq.8) then
         call MNTRPA(xmna)
         call EMTRPA(xema)
         do i = 1, 3
            xema(i) = xema(i) + Xsite(i + 6, 1)/Auacc - xmna(i)*Mnfct
         end do
      endif
c
c get site coordinates
      do i = 1, 6
         xe(i) = Xemlsc(i, 1) + Xsitau(i)
      end do
      if(Ibtrn.ne.1) then
         call PLTRP(1, jde(jpl), fr(jpl), 0, lsw)
         if(jde(jpl).le.0) goto 1200
         call PLTRP(1, jde(jpl), fr(jpl), -1, 0)
 
c also get accelerations for mutual events
         if(Ibtrn.gt.7) call PLTRPA(xpla)
         if(Klanb.gt.0) then
 
c get coordinates for satellite observed body
            call SBTRP(1, jde(jpl), fr(jpl), 0, lsw)
            if(jde(jpl).le.0) goto 1200
            call SBTRP(1, jde(jpl), fr(jpl), -1, 0)
            do i = 1, 6
               Xsbsun(i) = Xsb(i) + Xp(i)
            end do
 
c also get accelerations for mutual events
            if(Ibtrn.gt.7) then
               call SBTRPA(xsba)
               do i = 1, 3
                  xsba(i) = xpla(i) + xsba(i)
               end do
            endif
         else
 
c save planet stuff as xsbs
            do i = 1, 6
               Xsbsun(i) = Xp(i)
            end do
            do i = 1, 3
               xsba(i) = xpla(i)
            end do
         endif
         if(Nps1.gt.0) then
c get coordinates of distant planet at retarded time
c
            do i = 1, 3
               Xsc(i, 1)     = 0._10
               Xsc(i + 3, 1) = 0._10
               xsca(i)       = 0._10
            end do
 
c see if distant body same as central planet
            if(Klans1.gt.0) then
               call SCTRP(1, jde(3), fr(3), 0, lsw, 1)
               if(jde(3).le.0) goto 1200
               call SCTRP(1, jde(3), fr(3), -1, 0, 1)
 
c also get accelerations for mutual events
               if(Ibtrn.gt.7) call SCTRPA(xsca)
 
c see if distant body is a satellite
               if(Ncs1.le.0) then
 
c if not satellite, then coordinates are already relative to sun
                  do i = 1, 6
                     Xscsun(i, 1) = Xsc(i, 1)
                  end do
                  goto 500
               endif
            endif
 
c if so, get its center and coordinates relative to sun
            call PLTRP(1, jde(3), fr(3), 0, 0)
            if(jde(3).le.0) goto 1200
            call PLTRP(1, jde(3), fr(3), -1, 0)
            do i = 1, 6
               Xscsun(i, 1) = Xsc(i, 1) + Xp(i)
            end do
            if(Ibtrn.gt.7) then
 
c add on planets acceleration
               call PLTRPA(xpla)
               do i = 1, 3
                  xsca(i) = xsca(i) + xpla(i)
               end do
            endif
         endif
      endif
c
c*  start=2000
  500 lsw = lsw0
      if(Ibtrn.le.2) then
c moon or planet occults star
c correct star direction for parallax effect
         dtprds = DOT(usshat, xe)
         do i = 1, 3
            usohat(i) = usshat(i) + onebd*(usshat(i)*dtprds - xe(i))
         end do
      endif
      if(Nplnt0.ne.10) then
c correct velocity for retarded time
c set up observer coordinates
         do i = 1, 6
            Xsitau(i) = xe(i)
         end do
      else
 
c convert moon to au, au/day (rel. to site)
         call MNTRP(1, jde(2), fr(2), 0, 0, 1)
         if(jde(2).le.0) goto 1200
         call MNTRP(1, jde(2), fr(2), -1, 0, 1)
         do j = 1, 6
            Xsbsun(j) = Xm(j, 1)*Mnau
         end do
      endif
 
      nrtv = 9
      do i = 1, 3
         Xsitep(i, 1) = Xsitau(i) - Xsbsun(i)
         xo(i) = -Xsitep(i, 1)
      end do
      call UVECTR(3, Xsitep(1,1), Rsitp(1), Xsitp0(1,1), dum)
      r   = Rsitp(1)
      rr1 = r**2
      if(Nps1.gt.0) then
         nrtv = 11
         do i = 1, 3
            Xsitep(i, 2) = Xscsun(i, 1) - Xsbsun(i)
            xz(i) = Xscsun(i, 1) - xe(i)
         end do
         call UVECTR(3, Xsitep(1,2), Rsitp(2), Xsitp0(1,2), dum)
         if(Ibtrn.ne.9) then
            rzsq = DOT(xz, xz)
            rz   = SQRT(rzsq)
         endif
      endif
      call VLRTRD(Xsitau(4), Xsbsun(4), Xscsun(4,1), Aultvl, nrtv, 0)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=2500
c
c           determine correction to time of transit or occultation
c        dispatch to proper observable section
      if(Ibtrn.eq.2) then
      else if(Ibtrn.eq.4) then
         goto 600
      else if(Ibtrn.eq.5) then
         goto 1100
      else if(Ibtrn.eq.6) then
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=3000
c           ncodet=    0   transit
c           vector   r     sub e =  xe(1-3)
c           vector   r dot sub e =  xe(4-6)
c           vector   r     sub m = xps(1-3)
c           vector   r dot sub m = xps(4-6)
c        transit of one planet across another (mutual occultation)
c           substitute coordinates wrt occulted body
         do i = 1, 6
            xe(i)     = xe(i) - Xscsun(i, 1)
            Xsbsun(i) = Xsbsun(i) - Xscsun(i, 1)
         end do
         goto 600
      else if(Ibtrn.eq.7) then
         goto 1100
      else if(Ibtrn.eq.8) then
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=4000
c        mid-time of mutual occultation or eclipse
c           set up arrays for closest approach computation
c        mid-time occultation
c           phase correction to position of occulted body
         call PHSCOR(1, xz, xo, xz, Rhom, Rhos)
         do i = 1, 6
            Xsbsun(i)    = Xsbsun(i) - xe(i)
            Xscsun(i, 1) = Xscsun(i, 1) - xe(i)
         end do
         do i = 1, 3
            xsba(i) = xsba(i) - xema(i)
            xsca(i) = xsca(i) - xema(i)
         end do
         goto 700
      else if(Ibtrn.eq.9) then
c mid-time eclipse
c compute extra delay to eclipsing body
         rz = r + Rsitp(2)
 
c phase correction to position of eclipsed body
         call PHSCOR(0, xo, Xscsun, Xsbsun, Rhos, Rhom)
         goto 700
      else
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c           ncodet=   -1   occultation
c
c        moon occults star or planet
c        set up selenographic coordinates frame
         call MNSPT(jde(jmn), fr(jmn), mnspt1, 0, 1, nlibp2)
         if(Nps1.gt.0) then
c moon occults planet
c set up quantities like stellar occultation
            call UVECTR(2, xz, rz, usohat, dum)
            onebd = 1._10/rz
         endif
      endif
      do i = 1, 6
         Xsbsun(i) = Xsbsun(i) - Xsitau(i)
      end do
 
c muhat is unit vector from observer to moon
      do i = 1, 3
         muhat(i) = -Xsitp0(i, 1)
      end do
      usomu = DOT(usohat, muhat)
      s1    = SQRT(1._10 - usomu**2)
      cs1   = usomu/s1
 
c calculate unit vector center-to-limb
      do i = 1, 3
         rhomn(i) = usohat(i)*cs1 - muhat(i)/s1
      end do
c        rotate rhomn to selenographic frame and read watt's limb corect
c           note: some special treatment needed here for planet case
c     call prodct(rot,rhomn,rhomn0,3,3,1)
c           compute correction to limb
c     mrad=watlim(rhomn0)*mrad0
c        compute correction to occultation time
      cortim = (r*s1 - mrad)/DOT(rhomn, Xsbsun(4))
      goto 800
 
c transit across sun (or mutual occ. with coordinates moved)
  600 call UVECTR(3, xe, Re(mgo), Reh(1,mgo), dum)
      rmx = DOT(Reh(1,mgo), Xsbsun)
 
c get perpendicular component of planet position
      do i = 1, 3
         rmp(i) = Xsbsun(i) - rmx*Reh(i, mgo)
      end do
      rehrm = SQRT(DOT(rmp,rmp))
      fk    = rehrm*Rhos
      fl    = fk/Re(mgo) - rmx
      fm    = (rehrm - rhsm)*rehrm
      do i = 1, 3
         Xpr(i, mgo)  = fl*rmp(i) + fm*Reh(i, mgo)
         Vect(i, mgo) = Re(mgo)*rmp(i) + fk*Reh(i, mgo)
      end do
      dfdt3r = DOT(Xpr(1,mgo), xe(4)) + DOT(Vect(1,mgo), Xsbsun(4))
 
c cortim = -f/dfdt3
      cortim = rehrm*(Re(mgo)*(rhsm-rehrm) - Rhos*rmx)/dfdt3r
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c ncodet.ge. 1   photographs as black dot crosses sun
c not implemented
      goto 800
 
c compute correction to median time
  700 call UVECTR(7, Xsbsun, Re(1), Reh(1,1), dreh(1,1))
      call UVECTR(7, Xscsun, Re(2), Reh(1,2), dreh(1,2))
      f    = DOT(Reh(1,1), dreh(1,2)) + DOT(Reh(1,2), dreh(1,1))
      ctht = DOT(Reh(1,1), Reh(1,2))
      vex  = DOT(Xsbsun(4), Reh(1,1))
      vmx  = DOT(Xscsun(4,1), Reh(1,2))
      do i = 1, 3
         Xpr(i, 1)  = (Reh(i,2) - ctht*Reh(i,1))/Re(1)
         Xpr(i, 2)  = (Reh(i,1) - ctht*Reh(i,2))/Re(2)
         Vect(i, 1) = (dreh(i,2) - ctht*dreh(i,1) - f*Reh(i,1)
     .                - vex*Xpr(i,1))/Re(1)
         Vect(i, 2) = (dreh(i,1) - ctht*dreh(i,2) - f*Reh(i,2)
     .                - vmx*Xpr(i,2))/Re(2)
      end do
      Dfdt3(1) = DOT(Vect(1,1), Xsbsun(4)) + DOT(Xpr(1,1), xsba)
     .           + DOT(Vect(1,2), Xscsun(4,1)) + DOT(Xpr(1,2), xsca)
      cortim   = -f/Dfdt3(1)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=5000
c apply correction and test for convergence
  800 hrang = cortim*Secday
      Utrec = Utrec + hrang
      Fract = Utrec/Secday
c
c see if there is to be further iteration
      r = r*Aultsc
      if(ABS(dist-r).le.accdst .and. ABS(hrang).le.acctim) goto 900
c
c set up quantities for next iteration
      dist = r
      if(Nps1.gt.0) distz = rz*Aultsc
      if(iter.lt.100) goto 400
      call SUICID(
     .' MORE THAN 100 ITERATIONS NEEDED TO GET TRANSIT OR OCCULTATION TI
     .ME, STOP IN OCCULT ', 21)
c
c loop on mgo
  900 Deriv(2, mgo) = Utrec
      Nit(20) = Nit(20) + iter
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=6000
c calculate quantities for partial
      if(Ict(1).gt.0) then
         if(Ibtrn.eq.3) then
 
c moon occults a planet
            iso    = 2
            dtprds = r*usomu
         else if(Ibtrn.eq.4 .or. Ibtrn.eq.6) then
 
c transits
            Rmxc(mgo) = rmx
            do i = 1, 3
               Vect(i, mgo) = Vect(i, mgo)/rehrm
               Xpr(i, mgo)  = Xpr(i, mgo)/rehrm
            end do
            Rmrx(mgo)  = rehrm
            Dfdt3(mgo) = dfdt3r/rehrm/Secday
            goto 1000
         else if(Ibtrn.eq.5 .or. Ibtrn.eq.7) then
            goto 1200
         else if(Ibtrn.eq.8 .or. Ibtrn.eq.9) then
c
c set up for partials for mid-time observable
            Dfdt3(1) = Dfdt3(1)/Secday
            if(Ict(1).gt.0) call PHSCRS
            if(Nice.ge.0) then
c separation in terms of sum of radii
c get sine of angle using difference of unit vectors
               Rmxc(1) = SQRT(DOT(Xpr,Xpr))*Re(1)
               Rmrx(1) = Rhom/Re(1) + Rhos/Re(2)
               if(Rmrx(1).le.0._10) Rmrx(1) = 1._10
               Rmrx(2)     = ASIN(Rmxc(1))/Rmrx(1)
               Deriv(2, 2) = Rmrx(2)
               if(Ict(1).gt.0) then
 
c quantities only for separation partials
                  Ftct(1) = Rhom/Re(1)**2
                  Ftct(2) = Rhos/Re(2)**2
                  Rmxc(2) = (Ftct(1)*DOT(Reh(1,1),Xsbsun(4)) + Ftct(2)
     .                      *DOT(Reh(1,2),Xscsun(4,1)))/Secday
               endif
            endif
            goto 1100
         else
c occultations
c
c occultation of a star
            iso    = 1
            dtprds = DOT(xo, usshat)
         endif
         do i = 1, 3
            Xpr(i, mgo)  = xo(i) - dtprds*ushat(i, iso)
            Vect(i, mgo) = rhomn(i)
            Reh(i, mgo)  = ushat(i, iso)
         end do
         Rmrx(mgo) = cs1
         Re(mgo)   = onebd
         dtprds    = DOT(Xpr(1,mgo), xe(4))
         if(Ibtrn.eq.3) dtprds = dtprds -
     .                               DOT(Xpr(1,mgo), Xscsun(4,1))
 
c 2 - rho dot
         Dfdt3(mgo) = (-DOT(rhomn,Xsbsun) + cs1*onebd*dtprds)/Secday
      endif
 1000 if(mgo.lt.2 .and. Ibtrn.lt.8) then
         mgo  = 2
         rhsm = Rhos - Rhom
         goto 100
      endif
c
c*  start=9000
c termination calculations
c
 1100 return
 1200 Jd = 0
      return
      end
