      subroutine LESTHP(ncall)
 
      implicit none
c
c m.e.ash   aug 1975   subroutine lesthp
c calculate les-8/9 thrust and radiation pressure accelerations
c
c arguments
      integer*4 ncall
c        ncall=-2 setup once per iteration of a given step for partials
c        ncall=-1 setup once per iteration of a given step for motion
c        ncall= 0 setup once per step
c        ncall= k evaluate acceleration for equation k.gt.0
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'funcon.inc'
      include 'incon.inc'
      integer*4 lstrt1
      equivalence (lstrt1,L1)
      include 'intstf.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'sbstuf.inc'
      include 'sscon.inc'
      include 'sbthng.inc'

c local
      real*10 az,az1,az2,azold1,cpt1,crl2,cyw3,dtrns,el,flevel,
     .          frold,frold1,radfor(3),rlam,spt1,srl2,syw3
      integer   i,i8or9,if1,ii,itherf,itrns6,itrnst,itstep,
     .          jdold,jdold1,ntrnst
      real*10 eyaw(3),epitch(3),eroll(3),rpitch
      real*10 pp,bcorp(3),satmas,factlb
      character*4 les8/'LES8'/
      real*4    ther(4),thrst,er1pt,er2rl,er3yw
      equivalence (ther,thrst),(ther(2),er1pt),(ther(3),er2rl),
     .            (ther(4),er3yw)
      real*10 offset(2),of1pt,of2rl
      equivalence (offset,of1pt),(offset(2),of2rl)
 
      real*10 lesrdp(3),lesths(3),frdp,fths
c for les-8 and les-9 earth satellites in subroutine lesthp
c     con1(1) = weight of satellite in pounds (if zero, default launch
c               values used in in computing radiation pressure and
c               thrust accelerations)
c     con(17) = solvable parameter multiplying les-8/9 thrust
c               (if zero, assumed one)
c     con(18) = solvable parameter multiplying les-8/9 radiation
c               pressure (if zero, assumed one)
c
      equivalence (ntrnst,Itrns(5)),(itrns6,Itrns(6))
c     lstrt1=0 during starting proceedure
c
c kthrst = 0 steady state integration
c kthrst = 1 starting proceedure needed during thrust initiation
c            or termination
c
c jthrst = 0 no thrust at start, thrust at end of thrust  interval
c            as these points are approached during integration
c jthrst = 1 vice-versa during starting proceedure for thrust
c            initiation and termination
c
c t0sav  =   saved value of initial epoch to be restored after end
c            of starting proceedure after thrust initiation or
c            termination
c
c external functions
      real*10 DOT

      if(ncall.lt.0) then
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c set-up once per iteration of a given step for motion
         if(ncall.ge.-1) then
            itstep = itstep + 1
c
c     xlong is satellite longitude measured easterly in radians
c     if(nonce.gt.0) go to 311
c     nonce =1
c     if(fract.gt.1.0E-9_10) go to 311
c     xlong=atan2(cslng(1),cclng(1))       /   1.7453292519943296E-2_10
c     write(6,305)   jd,fract,xlong
c 305 format(i8,1p, 2d22.14)
c 311 continue
c
c           calculate unit yaw, pitch and roll vectors
            do i = 1, 3
               eyaw(i) = -Sbcor(i)/Rsb
            end do
            call CROSS(Sbcor(1),Sbcor(4),epitch)
            rpitch = SQRT(DOT(epitch,epitch))
            do i = 1, 3
               epitch(i) = epitch(i)/rpitch
            end do
            call CROSS(eyaw,epitch,eroll)
c
c les-8 is les-9 turned upside down
            if(Name(1:4).eq.les8) then
               i8or9 = 8
               if(Con1(1).le.0._10) satmas = 1092.52_10/32.15_10
               do i = 1, 3
                  epitch(i) = -epitch(i)
                  eroll(i)  = -eroll(i)
               end do
            endif
c
c correct for pitch and roll offsets and for
c attitude errors in pitch, roll and yaw
            if(itherf.gt.0) then
               call LESROT(cpt1,spt1,eroll,eyaw)
               call LESROT(crl2,srl2,eyaw,epitch)
               call LESROT(cyw3,syw3,epitch,eroll)
            endif
c see les-8/9 drawing d-60069, sheet 6 of 15
c
c is radiation pressure to be calculated
            if(Kp(81).gt.0) then
c
c           determine if satellite in shadow (lambda in /sbthng/)
c           sunlight  lambda=1
c           penumbra  lambda between 0 and 1
c           umbra     lambda=0
               call SHADOW
               rlam = Lambda/Rb2
            endif
c
c           calculate sun azimuth and elevation in satellite frame
c           azimuth measured from yaw axis positive to roll
c           elevation is angle above yaw-roll plane
c     bcor = coordinates of satellite relative to sun
            el = -ASIN(DOT(Bcor,epitch)/Rb)/Convd
            pp = -DOT(Bcor,epitch)
            do i = 1, 3
               bcorp(i) = -Bcor(i) - pp*epitch(i)
            end do
            az = ATAN2(DOT(bcorp,eroll),DOT(bcorp,eyaw))/Convd
c
c * * * * * * * * * *
c
c are station keeping computations to be done
            if(Kkp(94).gt.0) then
c
c no computations before station keeping enabling time,
c of which con1(10) is julian ephemeris date (half day less
c than jd,fract)
               if(Jd + Fract.ge.Con1(10) + 0.5_10) then
c
c no calculations if in midst of starting proceedure
                  if(lstrt1.ne.0) then
c
c no calculations if this is first iteration of
c predictor-corrector numerical integration technique
                     if(itstep.le.1) then
                        if(Kp(88).ne.3) goto 20
                     endif
c
c skip if there is no previous azimuth value
                     if(Azold.gt.-1.E9_10) then
c
c was the step size halved (integration must be forward in
c time, because station keeping thruster decisions depend
c on rememberance of things past)
                        if(jdold.ge.Jd) then
                           if(frold.ge.Fract) then
                              jdold = jdold1
                              frold = frold1
                              Azold = azold1
c itrnst=0 no transit last step
c itrnst=1 there was a transit last step (to be cancelled if
c step redone with half the step size)
                              if(itrnst.gt.0) ntrnst = ntrnst - 1
                           endif
                        endif
                        itrnst = 0
c
c calculate sun transit time (az=90 or -90 deg)
                        if(az.le.170._10) then
                           if(az.ge.-170._10) then
                              az1 = Azold - 90._10
                              az2 = az - 90._10
                              if((az1*az2).gt.0.) then
                                 az1 = Azold + 90._10
                                 az2 = az + 90._10
                                 if((az1*az2).gt.0.) goto 10
                                 ntrnst = ntrnst + 1
                                 Itrns(ntrnst) = -1
                              else
                                 ntrnst = ntrnst + 1
                                 Itrns(ntrnst) = 1
                              endif
                              dtrns = ABS(az2)
                              dtrns = dtrns/(dtrns + ABS(az1))
                              Frtrns(ntrnst) = Fract -
     .                           dtrns*((Jd-jdold) + (Fract-frold))
                              i = Frtrns(ntrnst)
                              if(Frtrns(ntrnst).lt.0._10) i = i - 1
                              Jdtrns(ntrnst) = i + Jd
                              Frtrns(ntrnst) = Frtrns(ntrnst) - i
                              itrnst = 1
                           endif
                        endif
                     else
                        itrnst = 0
                        jdold  = 0
                        frold  = 0
                     endif
c
c update the previous azimuth saved value
   10                jdold1 = jdold
                     frold1 = frold
                     azold1 = Azold
                     Azold  = az
                     jdold  = Jd
                     frold  = Fract
c
c determine thrusting strategy on the step after the fourth
c sun transit occured (in case step was halved)
                     if((ntrnst.ge.4) .and. (itrnst.le.0))
     .                  call STATON(i8or9,Jd,Fract,Hmx)
                  endif
               endif
            endif
c
c * * * * * * * * * *
c
c           calculate solar radiation pressure force on les-8/9
c     radfor(1) = force in x direction (z2 roll)   (micropounds)
c     radfor(2) = force in y direction (z3 yaw)    (micropounds)
c     radfor(3) = force in z direction (z1 pitch)  (micropounds)
c
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c        set-up once per a given step for partials
   20       if(Kp(81).gt.0) call FLES(az,el,rlam,i8or9,radfor)
         endif
      else if(ncall.eq.0) then
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c        set-up once per step
c     nonce =0
         itstep = 0
         i8or9  = 9
         frdp   = 1._10
         fths   = 1._10
         if(Con(17).ne.0._10) fths = Con(17)
         if(Con(18).ne.0._10) frdp = Con(18)
         satmas = 949.09_10/32.15_10
         if(Con1(1).gt.0._10) satmas = Con1(1)/32.15_10
         factlb = .3048E-9_10/(Aultsc*Ltvel)*8.64E4_10**2
c factor for conversion from force in micropounds to acceleration
c in au/day**2 is factlb/satmas
c
c read thrust and attitude history
         call LESTIN(Jd,Fract,ther,offset,itherf)
c itherf = 0 no thrust or attitude output from subroutine lestin
c itherf = 1 there is the following thrust and attitude output
c thrst  = signed thrust in postive roll direction (micropounds)
c          if thrust positive, thrusters on -z2 roll face are firing
c          with dv in +z2 roll direction
c              (west face on les9, dv with orbital velocity)
c              (east face on les8, dv against orbital velocity)
c          if thrust negative, thrusters on +z2 roll face are firing
c          with dv in -z2 roll direction
c              (east face on les9, dv against orbital velocity)
c              (west face on les8, dv with orbital velocity)
c note: les8 is turned upside down from les9
c er1pt  = attitude error about z1 pitch axis (radians)
c er2rl  = attitude error about z2 roll axis (radians)
c er3yw  = attitude error about z3 yaw axis (radians)
c of1pt  = offset angle about z1 pitch axis (radians)
c                (between -180 and +180 deg)
c of2rl  = offset angle about z2 roll axis (radians)
c                (between -11.25 and +11.25 deg)
c
c           calculate sines and cosines of attitude angles
         if(itherf.gt.0) then
            cpt1 = er1pt + of1pt
            spt1 = SIN(cpt1)
            cpt1 = COS(cpt1)
            crl2 = er2rl + of2rl
            srl2 = SIN(crl2)
            crl2 = COS(crl2)
            cyw3 = er3yw
            syw3 = SIN(cyw3)
 
            cyw3 = COS(cyw3)
         endif
      else
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c compute acceleration for motion
         if(Kkk.gt.0) then
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c compute accelerations for partials
            if(Icntrl(Kkk).eq.-17) Fn(1) = Fn(1) + lesths(Kk)
            if(Icntrl(Kkk).eq.-18) Fn(1) = Fn(1) + lesrdp(Kk)
            return
         else
c
c radiation pressure
            lesrdp(Kk) = 0.0_10
            if(Kp(81).gt.0) then
               lesrdp(Kk) = (radfor(1)*eroll(Kk) + radfor(2)*eyaw(Kk)
     .                      + radfor(3)*epitch(Kk))*factlb/satmas
               Fn(1) = Fn(1) + lesrdp(Kk)*frdp
            endif
c
c thrusting
            lesths(Kk) = 0.0_10
            if(itherf.gt.0) then
               lesths(Kk) = thrst*eroll(Kk)*factlb/satmas
               Fn(1) = Fn(1) + lesths(Kk)*fths
               if(itstep.le.1) then
                  if(Kp(88).lt.3) goto 40
               endif
               if((Kkp(92).gt.0) .and. (Kk.eq.3)) write(6,30)
     .           Jd, Fract, lesths, (Sbcor(i),i = 1,3),thrst,
     .           factlb, satmas
   30          format(i8,f17.14,1p,9D12.5)
            endif
c
c see if there is station keeping thrusting
   40       if(Kkp(94).gt.0) then
               if(itrns6.gt.0) then
                  flevel = 1.0_10
                  do i = 1, itrns6
                     ii  = i
                     if1 = 0
                     if(Jd.lt.Jdfire(i,1)) return
                     if(Jd.eq.Jdfire(i,1)) then
                        if(Fract.lt.Frfire(i,1)) return
                        if(Fract.eq.Frfire(i,1)) goto 60
                     endif
                     if1 = 1
                     if(Jd.lt.Jdfire(i,2)) goto 80
                     if(Jd.eq.Jdfire(i,2)) then
                        if(Fract.lt.Frfire(i,2)) goto 80
                        if(Fract.eq.Frfire(i,2)) goto 60
                     endif
                  end do
               endif
            endif
            return
c
c special treatment of start and end of thrust in constant
c step size integration methods
   60       if(Kp(88).gt.1) then
               flevel = 0.5_10
               goto 100
c
c variable step size method needs starting proceedure
c at start and end of thrust
            else if(if1.gt.0) then
               if(Jthrst.gt.0) return
               Kthrst = 1
            else if(Jthrst.le.0) then
               Kthrst = 1
               return
            endif
c
c station keeping thrusting
   80       Fn(1) = Fn(1) + Vthrst(ii)*eroll(Kk)*factlb/satmas*flevel
            if(lstrt1.eq.0) return
            if(itstep.le.1) then
               if(Kp(88).lt.3) return
            endif
         endif
  100    if((Kk.eq.3) .and. (Kkp(94).gt.1)) write(6,150) Jd,
     .     Fract, ii, Vthrst(ii),eroll,factlb,satmas,flevel
 
  150    format(i8,f17.14,i3,1p,7D13.6)
      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      return
      end
