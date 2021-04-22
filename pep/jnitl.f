      subroutine JNITL(goose,cond,setp,kind,dy)
      implicit none

c Subroutine JNITL   J.F.Chandler   1977 November 30
c Set up for elliptic or hyperbolic orbit calculations
c
c arguments to JNITL, JLIPT
      real*10 goose,cond(9),setp(50),y(6),dy(6,6),t,ry,ry2,ry3
      integer*4 kind,nv
c        goose - the square root of (g*(mass of body+mass of center))
c        cond - array of initial conditions + optional extras
c        cond(7)=time derivative of cond(1)
c        cond(8)=time derivative of cond(5)
c        cond(9)=mass factor for scaling output
c        Usage of array 'setp' --- storage for saved quantities from
c        setting up.  It is effectively equivalenced to the following:
c      b(3,2),tsv,a,e,anom0,motpi,mu2,secc,cecc,quan2,quan3,quan4,quan5,
c      1      7   8 9  10    11    12  13   14   15    16    17    18
c
c      quan12(3),quan13(3),quan11,motion,mu,sasc,casc,spci,cpci
c         19        22      25     26    27  28   29   30   31
c
c      sinc,cinc,sper,cper,a0,a1,per0,per1,pscale
c       32   33   34   35  36 37  38   39   40
c
c        Note: all variables up thru quan5 (setp(18)) are used in all
c        elliptic calculations, but those from there on are needed only
c        if partials with respect to the initial conditions are needed.
c
c        kind - indicates type of setup,
c                0 => only want coordinates
c                1 => also want partials w.r.t. initial conditions
c                2 => also need sinc,cinc,sper,cper in extended setp
c                3 => also set up for moving periapse and semimajor axis
c        y - array for output coordinates
c        ry,ry2,ry3 - output radial distance, square, and cube
c        dy - output array of partials
c
c        input orbital elements
c     equivalence (cond(1),ap),(cond(2),ep),(cond(3),incp),
c    1 (cond(4),ascp),(cond(5),perp),(cond(6),anomp)
c
      include 'funcon.inc'

c local variables
      real*10 anom,anoms,casc,cf,cinc,cpci,cper,qq1,qq2,qq2n,qq4,qq5,
     . qq6,qq7,qq8,qqh,qql,qqs,quan1,sasc,sinc,spci,sper,th
      integer   i,int,iyb,j,nitr,nv1,type
      real*10 motion, ybar(6), absa
      real*4    anom4

      setp(7)  = 1E60_10
      setp(8)  = cond(1)
      absa     = ABS(setp(8))
      setp(9)  = cond(2)
      setp(10) = cond(6)/360._10
      quan1    = 1._10 - setp(9)**2
      setp(15) = SQRT(ABS(quan1))
      setp(16) = absa*setp(15)
      setp(17) = SQRT(absa)
      motion   = goose/absa/setp(17)
      setp(11) = motion/Twopi
      setp(12) = goose**2
      setp(17) = setp(17)*goose
      setp(18) = setp(17)*setp(15)
      setp(17) = -setp(17)
      th   = cond(3)*Convd
      cinc = COS(th)
      sinc = SIN(th)
      th   = cond(4)*Convd
      casc = COS(th)
      sasc = SIN(th)
      th   = cond(5)*Convd
      sper = SIN(th)
      cper = COS(th)
      spci = sper*cinc
      cpci = cper*cinc
      setp(1) = casc*cper - sasc*spci
      setp(4) = -casc*sper - sasc*cpci
      setp(2) = sasc*cper + casc*spci
      setp(5) = casc*cpci - sasc*sper
      setp(3) = sper*sinc
      setp(6) = cper*sinc

c Decide if partials setup needed
      if(kind.eq.0) return

c Fill rest of setp array for partials
      setp(25) = setp(9)/quan1
      setp(26) = motion
      setp(27) = goose
      setp(28) = sasc
      setp(29) = casc
      setp(30) = spci
      setp(31) = cpci
      do i = 1, 3
         setp(i+18) = setp(i)*setp(8)
         setp(i+21) = setp(i+3)*setp(25)
      end do
      dy(3,4) = 0._10
      dy(6,4) = 0._10
      if(kind.le.1) return
      setp(32) = sinc
      setp(33) = cinc
      setp(34) = sper
      setp(35) = cper
      if(kind.le.2) return
      setp(36)=setp(8)
      setp(37)=cond(7)
      setp(38)=th
      setp(39)=cond(8)*Convd
      setp(40)=cond(9)
      return

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry JLIPTP(t,setp,nv,y,ry,ry2,ry3,dy)
c same as JLIPT, but with moving periapse and semimajor axis
      type=1
      goto 40
c
      entry JLIPT(t,setp,nv,y,ry,ry2,ry3,dy)
c Compute elliptic/hyperbolic quantities for time t from initial epoch.
c Corresponding setp array must already have been set up.
c nv - operation indicator.
c      1 do position and velocity
c     >1 do partials as well (up to dy(1-6,nv-1) )
c      0 do just position
c     -1 do just velocity, in y(1-3)
c     -2 do just acceleration, in y(1-3)

      type=0
   40 continue
c        Take short cut if same epoch as last time
      if(t.ne.setp(7)) then
         setp(7) = t

c Get mean anomaly
         anom= setp(11)*t + setp(10)
         if(setp(8).gt.0._10) then
c elliptic orbits
            int= anom
            if(anom.lt.0._10) int = int - 1
            anoms= Twopi*(anom - int)
            anom4= anoms
            if(anoms.lt.Pi) then
c range is 0 to pi
               qqh= Pi
               qql= 0._10
            else if(setp(9).le.0.9_10) then
c range is pi to 2 pi (this branch is kept for compatibility of results)
               qqh= Twopi
               qql= Pi
            else
c range is -pi to 0
               anoms  = anoms - Twopi
               qqh = 0._10
               qql = -Pi
            endif

c Solve Kepler's eqn for eccentric anomaly
            nitr = 0
            if(setp(9).gt.0.7_10 .and. ABS(anoms).lt.0.3_10) then
               qq1 = anoms/(ABS(1._10-setp(9)) + ABS(anoms))
            else
               qq1 = anoms + setp(9)
     .          *(sin(anom4)+setp(9)*sin(2.*anom4)*0.5_10)
            endif

   50       nitr = nitr+1
            qq2n = anoms - qq1 + setp(9)*SIN(qq1)
            qq2 = qq2n/(1._10 - setp(9)*COS(qq1))
c Damp oscillations if starting value is very bad
            if(ABS(qq2).gt.Pi) qq2 = SIGN(1._10,qq2)
            qqs = qq1
            qq1 = qq1+qq2
            if(qq2.lt.0._10) then

c Step is negative
               qqh = qqs
               if(qq1.lt.qql) qq1 = qql
            else

c Step is positive
               qql = qqs
               if(qq1.gt.qqh) qq1 = qqh
            endif
            if(ABS(qq2).gt.5.E-15_10) then
               if(nitr.le.9) goto 50
               if(ABS(qq2n).gt.5.E-14_10) then
                  if(nitr.lt.20) goto 50
                  call SUICID(
     .             'CAN''T SOLVE KEPLER EQUATION, STOP JNITL ',10)
               endif
            endif
c Save results of solution as sine,cosine
            setp(13) = SIN(qq1)
            setp(14) = COS(qq1)
         else

c hyperbolic orbits
            anoms=anom*Twopi
            qq1=anoms/(setp(9)-1._10)
            if(setp(9).gt.1.E5_10) then
            else if(ABS(qq1).lt.0.8) then
               qq1=qq1 - setp(9)*qq1**3/6._10/(setp(9)-1._10) +
     .          (9._10*setp(9)**2+setp(9))*qq1**5/120._10
     .          /(setp(9)-1._10)
            else
               qq1=SIGN(LOG(2._10*ABS(qq1)),qq1)
            endif
            nitr=0
            qq2=1._10
   60       nitr = nitr+1
            if(ABS(qq2).lt.0.2_10) then
               qq2n = anoms + qq1 - setp(9)*SINH(qq1)
               qq2 = -qq2n/(1._10 - setp(9)*COSH(qq1))
            else
               qq2n=2._10*ABS(anoms+qq1)/setp(9)+EXP(-ABS(qq1))
               qq2=SIGN(LOG(qq2n),qq1)-qq1
            endif
            qq1=qq1+qq2
            if(ABS(qq2).gt.5.E-15_10) then
               if(nitr.le.9) goto 60
               if(ABS(qq2n).gt.5.E-14_10) then
                  if(nitr.lt.30) goto 60
                  call SUICID(
     .             'CAN''T SOLVE HYP. KEPLER EQUATION, STOP JNITL',11)
               endif
            endif
c Save results of solution as hyperbolic sine,cosine
            setp(13) = SINH(qq1)
            setp(14) = COSH(qq1)
         endif

c Get radial distance including change of periapse and semimajor axis
         if(type.eq.1) then
            setp(8)=setp(36)+t*setp(37)
            setp(16) = setp(8)*setp(15)
            th=setp(38)+t*setp(39)
            sper=SIN(th)
            cper=COS(th)
            sasc=setp(28)
            casc=setp(29)
            sinc=setp(32)
            cinc=setp(33)
            spci=sper*cinc
            cpci=cper*cinc
            setp(30)=spci
            setp(31)=cpci
            setp(34)=sper
            setp(35)=cper
            setp(1) = casc*cper - sasc*spci
            setp(4) = -casc*sper - sasc*cpci
            setp(2) = sasc*cper + casc*spci
            setp(5) = casc*cpci - sasc*sper
            setp(3) = sper*sinc
            setp(6) = cper*sinc
            do i = 1, 3
               setp(i+18) = setp(i)*setp(8)
               setp(i+21) = setp(i+3)*setp(25)
            end do
         endif
         ry = setp(8)*(1._10 - setp(9)*setp(14))
         ry2= ry*ry
         ry3= ry2*ry
      endif

c Compute vector in orbit plane coordinates
  100 ybar(1) = setp(8)*(setp(14) - setp(9))
      ybar(2) = setp(16)*setp(13)
      ybar(3) = setp(17)*setp(13)/ry
      ybar(4) = setp(18)*setp(14)/ry
      if(nv.gt.0) then
c Compute position and velocity
c Get position and velocity
         call PRODCT(setp,ybar(1),y(1),3,2,1)
         call PRODCT(setp,ybar(3),y(4),3,2,1)
         if(nv.gt.1) then
            cf  = -setp(12)/ry3
            qq1 = t*1.5_10
            qq2 = cf*qq1
            qq4 = setp(13)/setp(26)
            if(setp(8).lt.0._10) qq4=-qq4
            qq5 = setp(8)*setp(14)/ry
            qq6 = cf*qq4
            qq7 = -ybar(4)*setp(25)
            do i=1,3
               j = i+3

c Partials w.r.t. a
               dy(i,1) = (y(i) - qq1*y(j))/setp(8)
               dy(j,1) = (-0.5_10*y(j) - qq2*y(i))/setp(8)

c Partials w.r.t. e
               dy(i,2) = qq4*y(j) - setp(i+18) - setp(i+21)*ybar(2)
               dy(j,2) = qq5*y(j) + qq6*y(i) + qq7*setp(j)
            end do
            if(nv.gt.3) then

c Partials w.r.t. inc (radians)
               dy(1,3) = setp(28)*y(3)
               dy(2,3) = -setp(29)*y(3)
               dy(3,3) = setp(30)*ybar(1) + setp(31)*ybar(2)
               dy(4,3) = setp(28)*y(6)
               dy(5,3) = -setp(29)*y(6)
               dy(6,3) = setp(30)*ybar(3) + setp(31)*ybar(4)
               if(nv.gt.4) then

c Partials w.r.t. asc (radians)
                  dy(1,4) = -y(2)
                  dy(2,4) = y(1)
                  dy(4,4) = -y(5)
                  dy(5,4) = y(4)
                  if(nv.gt.5) then
                     qq8 = cf/setp(26)
                     do i=1,3
                        j = i+3

c Partials w.r.t. per (radians measured from node)
                        dy(i,5) = setp(j)*ybar(1) - setp(i)*ybar(2)
                        dy(j,5) = setp(j)*ybar(3) - setp(i)*ybar(4)

c Partials w.r.t. anom0 (radians from periapse)
                        dy(i,6) = y(j)/setp(26)
                        dy(j,6) = y(i)*qq8
                     end do
                  endif
               endif
            endif
         endif
         return
      else
c Compute only position or velocity or acceleration
c leave resulting vector in y(1-3)
         nv1 = iabs(nv)+1
         if(nv1.eq.1) then
         else if(nv1.eq.2) then

c Get velocity
            iyb = 3
            goto 150
         else if(nv1.eq.3) then
c Get acceleration
c Use radial force law
            cf = -setp(12)/ry3
            ybar(5) = ybar(1)*cf
            ybar(6) = ybar(2)*cf
            iyb     = 5
            goto 150
         else
            call SUICID('BAD NV, STOP IN JLIPT   ', 6)
         endif

c Get position
         iyb = 1

c Rotate vector to standard frame
  150    call PRODCT(setp,ybar(iyb),y,3,2,1)
         return
      endif
      end
