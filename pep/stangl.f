      subroutine STANGL(xsun)
 
      implicit none
c
c m.e. ash  june 1970   subroutine stangl
c calculate theoretical value of look angle
c observation of one satellite by another (both with same
c central body)
c
c interpolation through sbtrp,sctrp added 1977 oct - j.f.chandler

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 eyaw(3),epitch(3),eroll(3),r3yaw,r2roll,r2r3,r1r,r1rr,
     . ryaw,rpitch,rr1,fake,cake
      equivalence (Raddum,epitch),(Raddum(4),eroll),
     .            (Raddum(7),eyaw)
      equivalence (Angdum,r3yaw),(Angdum(2),r2roll),
     .            (Angdum(3),r2r3),(Angdum(4),r1r),
     .            (Angdum(5),r1rr)
      equivalence (Angdum(6),rpitch),(Angdum(7),ryaw)
      equivalence (Angdum(4),rr1),(Angdum(6),fake),(Angdum(7),cake)
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'ltrapobs.inc'
      include 'mnsprt.inc'
      include 'number.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      real*4    acctim
      equivalence (acctim,Estf)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'
c
c external functions
      real*10 A1WWV,CTATF,DOT

c local
      real*10 alphat,ctuts,deltat,fractr,r1ptch
      real*10 xsun(13),xmon(6),tmdls1,tmdls2,frcts,xsnmn(3)
      integer   i,iter1,j,lsw,mnspt1
      integer*4 itstvl/0/
c
c determine time quantities
      lsw   = 1
      Utrec = Ihr*3600 + Imin*60
      Utrec = Utrec + Sec
      Fract = Utrec/Secday
      Atuts = A1WWV(Jd,Fract)
      Ctat  = CTATF(Jd,Fract + Atuts/Secday,1,1)
      Ututs = 0.
      ctuts = Ctat + Atuts
      Ctrec = Utrec + ctuts
      do while( .true. )
         if(Ctrec.lt.0) then
            Jd    = Jd - 1
            Ctrec = Ctrec + Secday
         else if(Ctrec.eq.0) then
            goto 100
         else
            do while( Ctrec.ge.Secday )
               Jd    = Jd + 1
               Ctrec = Ctrec - Secday
            end do
            goto 100
         endif
      end do
  100 ctrecf = Ctrec/Secday
      Tmdly1 = 0._10
      if(Dstf(4).gt.0._10) Tmdly1 = Dstf(4)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           calculate position of observing satellite at
c           receive time
c
      call SCTRP(1,Jd,ctrecf,0,lsw,1)
      if(Jd.le.0) goto 9900
 
      do i = 1,3
         Xscsun(i,1) = Xsc(i,1)*Aultsc
      end do
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           iteration to calculate position of observed satellite
c           at reflection time
 
      if(Nplnt0.ge.0) then
c
c           determine reflection time
         iter1  = 0
         mnspt1 = 0
         tmdls1 = Aultsc
      else
c see if star is observed
         if(Nspot.le.0) call SUICID(' NO STAR, STOP IN STANGL',6)
         mnspt1 = 0
         call SPOTCD(Jd,ctrecf,mnspt1,0,Nplnt0,ctuts,Atuts,1)
         do i = 1,3
            Xsbsun(i)   = Xspcd(i,1)
            Xsitep(i,1) = -Xspcd(i,1)
         end do
         Rsitp(1) = 1._10
 
c no test for dummy star observations occulted by central body
         goto 700
      endif
  200 call TIMINC(Jd,ctrecf,Jdx,fractr,-Tmdly1/Secday)
c
c see if central body is observed
      if(Ncs1.ne.Nplnt0) then
c
c see if sun is observed
         if(Nplnt0.eq.0) then
            if(Ncs1.ne.3) call SUICID(' CANNOT PROCESS SATELLITE BASED O
     .BSERVATIONS OF SUN FOR CENTRAL BODY NOT EARTH, STOP STANGL ', 23)
c
c (must be used if sun sensor accuracy is better than .005 dg)
c true sun interpolation
            if(Ict(26).eq.0) then
c
c mean sun calculations
               call PRECES(Jdx + (fractr-0.5_10))
               call MENSUN(Jdx, fractr, xsun, itstvl)
               call MENMON(Jdx, fractr, xmon, itstvl)
               do i = 1, 3
                  Xsbsun(i) = (xsun(i) + Mass(10)*xmon(i))*Aultsc
               end do
            else
               call ETRP(1, Jdx, fractr, 0, 0, 1, 2)
               do i = 1, 3
                  Xsbsun(i) = -Xemlsc(i, 1)
               end do
            endif
            goto 400
c
c see if probe is observed
         else if(Klanb.gt.0) then
c
c read probe peripheral data set
            call SBTRP(1, Jdx, fractr, 0, lsw)
            if(Jdx.le.0) then
               Jd = 0
               return
            else
               do i = 1, 3
                  Xsbsun(i) = Xsb(i)*Aultsc
               end do
               goto 400
            endif
c see if planet or moon is observed
         else if(Mnplnt.ne.0) then
            call PLTRP(1,Jdx,fractr,0,lsw)
            if(Jdx.le.0) goto 9900
            do i=1,3
              Xsb(i)=Xp(i)
              Xsbsun(i)=Xsb(i)*Aultsc
           end do
            goto 400
         endif
 
c central body is observed
      else if(Nspot.gt.0) then
         call SPOTCD(Jdx, fractr, mnspt1, 0, Nplnt0, ctuts, Atuts, 1)
         if(Jdx.gt.0) goto 9900
         do i = 1, 3
            Xsbsun(i) = Xspcd(i, 1)
         end do
         goto 400
      else
         do i = 1, 3
            Xsbsun(i) = 0._10
         end do
         goto 400
      endif
c
c true moon interpolation
      call MNTRP(1, Jdx, fractr, 0, 0, 1)
      do i = 1, 3
         Xsbsun(i) = Xm(i, 1)*Mnltsc
      end do
c
c sun for special moon-sun mode
      if(Jct(50).gt.0) then
         frcts = (Ctrec - tmdls1)/Secday
         call ETRP(1, Jdx, frcts, 0, 0, 1, 2)
         do i = 1, 3
            xsun(i) = Xscsun(i, 1) - Xemlsc(i, 1)
         end do
         xsun(7) = xsun(1)**2 + xsun(2)**2 + xsun(3)**2
         xsun(8) = SQRT(xsun(7))
      endif
c
c determine reflection to receive time delay, decide
c whether to re-iterate for reflection time
  400 do i = 1, 3
         Xsitep(i, 1) = Xscsun(i, 1) - Xsbsun(i)
      end do
      Rsitp(1) = SQRT(DOT(Xsitep(1,1),Xsitep(1,1)))
      Tmdly2   = Rsitp(1)
      iter1    = iter1 + 1
c
c decision on sun light time in special sun-moon mode
      if(Mnplnt.eq.0 .and. Jct(50).gt.0) then
         tmdls2 = tmdls1
         tmdls1 = xsun(8)
         if(ABS(tmdls2-tmdls1).gt.acctim) goto 500
      endif
c
c decision on observed body light time
      if(ABS(Tmdly2-Tmdly1).le.acctim) goto 600
  500 if(iter1.gt.1000) call SUICID(
     .'MORE THAN 1000 REFLECTION ITERATIONS, STOP IN STANGL', 13)
      lsw    = -1
      Tmdly1 = Tmdly2
      goto 200
  600 Nit(20) = Nit(20) + iter1
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c test for dummy observation occulted by central body
      if(Idumob.eq.1) then
         if(Ncs1.eq.Nplnt0 .and. Nspot.le.0) goto 700
         call CNTOCC(Xsb, Xsitep, 1, Jd)
         if(Jd.le.0) return
      endif
c
c calculate spot partials
  700 if(Nspot.gt.0) call SPOTCV(0, Nplnt0, 1)
c
c calculate velocites
c
c calculate velocity of observing satellite at receive
      call SCTRP(1, Jd, ctrecf, 1, 0, 1)
      do i = 4, 6
         Xscsun(i, 1) = Xsc(i, 1)*Aultvl
      end do
c see if velocity of observed body is to be calculated
c (itstvl is always zero)
      if(itstvl.gt.0) then
c
c see if star is observed
         if(Nplnt0.lt.0) then
            do i = 4, 6
               Xsbsun(i) = 0._10
            end do
c
c see if central body observed
         else if(Ncs1.eq.Nplnt0) then
            do i = 4, 6
               Xsbsun(i) = 0._10
 
c logic to go here for observing spot on central body
            end do
c
c see if sun is observed
         else if(Nplnt0.eq.0) then
            do i = 4, 6
c this is not valid if subroutine mensun was not called
c but itstvl is always zero anyway
               Xsbsun(i) = (xsun(i) + Mass(10)*xmon(i))*Aultvl
            end do
c
c see if planet or moon is observed
         else if(Klanb.le.0) then
c           planet is observed
            call PLTRP(1,Jdx,fractr,1,0)
            do i=4,6
               Xsb(i)=Xp(i)
               Xsbsun(i)=Xsb(i)*Aultvl
            end do
         else
c
c calculate velocity of observed satellite at reflection time
            call SBTRP(1, Jdx, fractr, 1, 0)
            do i = 4, 6
               Xsbsun(i) = Xsb(i)*Aultvl
            end do
         endif
c
c calculate relative velocities
         do i = 4, 6
            Xsitep(i, 1) = Xsbsun(i) - Xscsun(i, 1)
         end do
      endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c decide if look angles are in satellite reference system
      if(Ncodf.gt.4) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           calculation of theoretical values of look angles
c           relative to stars
         if(Ncodf.gt.5) then
c
c astrographic observations
            do i = 1, 3
               Ves(i) = Xscsun(i + 3, 1)
            end do
         else
c
c astrometric observation (elliptic aberration not removed
c from comparison star catalogue)
            do i = 1, 3
               Ves(i) = Xscsun(i + 3, 1) + Velipt(i)
            end do
         endif
c
c     ves should be incremented by velocity of earth relative to sun
c     for sun observations and relative to center of mass of solar
c     system for star observations if look angles are being calculated
c     for comparison with observations made in satellite based
c     inertial refernce frame, rather than for comparison with
c     observations against the star background with reduction using
c     a comparison star catalogue as is implicitely assumed here
c
c           aberration correction
         do i = 1, 3
            Xsitep(i, 1) = Xsitep(i, 1) - Rsitp(1)*Ves(i)
c negative sign in aberration correction because xsitep is
c minus the observed body relative to the observing body
         end do
c
c calculate theoretical value of right ascension and
c declination referred to the mean equinox and equator
c of the reference epoch
         rr1    = Xsitep(1, 1)**2 + Xsitep(2, 1)**2
         alphat = ATAN2(-Xsitep(2,1), -Xsitep(1,1))
         if(alphat.lt.0._10) alphat = alphat + Twopi
         Deriv(2, 1) = alphat/Convd
         deltat = ATAN2(-Xsitep(3,1), SQRT(rr1))
         Deriv(2, 2) = deltat/Convd
 
         if(Mnplnt.eq.0 .and. Jct(50).gt.0) call SUICID(
     .'CANNOT CALCULATE SUN RA,DEC IN SPECIAL MOON-SUN MODE, STOP IN STA
     .NGL', 17)
c
c calculate partial derivative quantities
         Rsitp(1) = SQRT(rr1 + Xsitep(3,1)**2)
         fake     = -Xsitep(3, 1)/Rsitp(1)
         cake     = Rsitp(1)*SQRT(1._10 - fake**2)*Convd/Aultsc
         fake     = fake/Rsitp(1)
c
c constant biases in right ascension and declination
         rr1 = rr1*Convd/Aultsc
      else
c
c           correct for aberration due to movement of observing satellit
c     for sun observation must increment xscsun(4-6) with velocity of
c        earth relative to sun  $ no, because of light time correction
c     for star observation must increment xscsun(4-6) with velocity of
c        earth relative to center of mass of solar system
c     however, this accuracy (20 arc sec) not needed in earth oriented
c        frame look angle observations since the attitude control of
c        the satellite is not that accurate
         do i = 1, 3
            Xsitep(i, 1) = Xsitep(i, 1) - Rsitp(1)*Xscsun(i + 3, 1)
c negative sign in aberration correction because xsitep is
c minus the observed body relative to the observing body
         end do
         Rsitp(2) = Rsitp(1)*Ltvel
c
c calculate unit vector
         Rsitp(1) = SQRT(Xsitep(1,1)**2 + Xsitep(2,1)**2 + Xsitep(3,1)
     .              **2)
         do j = 1, 3
            Xsitp0(j, 1) = Xsitep(j, 1)/Rsitp(1)
         end do
c
c calculate unit yaw, pitch and roll vectors
         ryaw = SQRT(DOT(Xscsun(1,1),Xscsun(1,1)))
         do i = 1, 3
            eyaw(i) = -Xscsun(i, 1)/ryaw
         end do
         call CROSS(Xscsun(1,1), Xscsun(4,1), epitch)
         rpitch = SQRT(DOT(epitch,epitch))
         do i = 1, 3
            epitch(i) = epitch(i)/rpitch
         end do
         call CROSS(eyaw, epitch, eroll)
c
c sun angles in special sun-moon mode
         if(Mnplnt.eq.0 .and. Jct(50).gt.0) then
            do i = 1, 3
               xsun(i) = xsun(i) - xsun(8)*Xscsun(i + 3, 1)
            end do
            xsun(7)  = xsun(1)**2 + xsun(2)**2 + xsun(3)**2
            xsun(8)  = SQRT(xsun(7))
            xsun(9)  = ATAN2(-DOT(xsun,eroll), -DOT(xsun,eyaw))/Convd
            xsun(10) = ASIN(-DOT(xsun,epitch)/xsun(8))/Convd
c
c fraction of moon disk illuminated by sun as
c seen from satellite
            do i = 1, 3
               xsnmn(i) = Xsitep(i, 1) - xsun(i)
            end do
            xsun(11) = DOT(xsnmn, Xsitep)
     .                 /(SQRT(DOT(xsnmn,xsnmn))*Rsitp(1))
            xsun(11) = (Pi - ACOS(xsun(11)))/Pi
            xsun(12) = xsun(11)*180._10
            xsun(13) = 3476._10/Rsitp(2)/Convd
         endif
c
c calculate components of look vector in pitch, roll,
c yaw coordinate system
         r1ptch = -DOT(Xsitep(1,1), epitch)
         r2roll = -DOT(Xsitep(1,1), eroll)
         r3yaw  = -DOT(Xsitep(1,1), eyaw)
c
c calculation of theoretical values of look angles
c in satellite reference system
         if(Nice.le.0) then
            Deriv(2, 1) = ATAN2(r2roll, r3yaw)/Convd
            r2r3 = (r2roll**2 + r3yaw**2)*Convd/Aultsc
         endif
 
c always calculate roll angle for error in pitch
         r1r = r1ptch/Rsitp(1)
         Deriv(2, 2) = ASIN(r1r)/Convd
         r1rr = Rsitp(1)*SQRT(1._10 - r1r**2)*Convd/Aultsc
      endif
c
c constant biases in pitch and roll
      if(Neqnox.gt.0) then
         if(Nice.le.0) Deriv(2, 1) = Deriv(2, 1) + Pnox
         if(Nice.ge.0) Deriv(2, 2) = Deriv(2, 2) + Plat
      endif
      return
 
c missing ephemeris data
c return with error condition
 9900 Jd = 0
      return
      end
