      subroutine LMBCOR(alphat,deltat)
 
      implicit none

c
c     make limb corr.  to theoretical value of meridian circle obs.
c     specialized for corrections to  r. a., dec. of center observed
c     when the limb is on the meridian. the declination  code is then
c     without value and if the corresponding r.a. code is not available
c     as in the greenwich series which have r.a. and dec. separate,
c      then trial and error  for dec. correction is the only solution.
c
c parameters
      real*10 alphat,deltat

c array dimensions
      include 'globdefs.inc'

c commons
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 x(6,3)
      equivalence (Xm(1,1),x(1,1))
c           x(.,1) =coordinates of moon relative to earth (retarded time
c                      if moon the observed body, in which case there is
c                      no velocity)
c           x(.,2) =coordinates of sun relative to earth (velocity
c                      vice-versa)
c           x(.,3) =coordinates of planet relative to earth at retarded
c                      time (no velocity)
      include 'empcnd.inc'
      real*10 mnr0
      equivalence (mnr0,Mcond(7))
      include 'funcon.inc'
      include 'ltrapx.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      include 'kobequiv.inc'
      real*10 dist
      equivalence (dist,Dstf(4))
      include 'param.inc'
      include 'sitcrd.inc'
      character*4 sitf4
      equivalence (Sitf,sitf4)
 
c locals
      integer   mnrad, nlimb
      real*10 test1,test2,chck1,chck2,sgn
      character*2 north/'N '/, south/'S '/, upper/'U '/,lower /'L '/,
     1      east/'E '/, west/'W '/, one/'1 '/, two/'2 '/, center/'C '/,
     2      top /'T '/, botm /'B '/, blank/'  '/
      character*4 rname /'RADC'/
      character*4 tcon /'GTOK'/, tname
      real*10 rphipr/0.90011_10/
 
c rphipr is in radians
      real*10 rrho/0.997950734_10/
 
c rrho is  radius  to radcliffe site in units of earth radii
      real*10 conspr/1.659251067E-2_10/
      real*10 condis/1.282221187_10/
c condis is in light-seconds
c conspr is the constant of sin(parallax) at mean distance condis
c sin(prlx)= 3422.451 sec of arc in radians
      real*10 tphipr/0.619424_10/
      real*10 trho/0.998876_10/
      character*4 mern /'MUSN'/, mrser /'M561'/
      real*10 mrho/0.99870_10/
      real*10 mphipr/0.676024_10/
      character*4 gname /'GREN'/, grnser /'1312'/
      real*10 gphipr/0.8845846_10/
      real*10 grho/0.997994_10/
      real*10 tfract
      real*10 z(6,2)
      real*10 caps
      real*10 dradst
      real*10 ddcdst
      real*10 sinpdf
      real*10 dellim
 
      nlimb = 0
      tname = tcon
      sgn  = 1._10
      if(Series.eq.grnser) trho   = grho
      if(Series.eq.grnser) tphipr = gphipr
      if(Series.eq.grnser) tname  = gname
      if(Series.eq.mrser) trho    = mrho
      if(Series.eq.mrser) tphipr  = mphipr
      if(Series.eq.mrser) tname   = mern
      tfract = Sidvel*8.640E4_10/Twopi
      call CORCHN(z(1,1),x(1,1))
      call CORCHN(z(4,1),x(4,1))
      mnrad  = mnr0/Ltvel
      caps   = mnrad/dist
      caps   = ASIN(caps)
      dradst = (z(1,1)*z(5,1) - z(2,1)*z(4,1))/(z(1,1)**2 + z(2,1)**2)
      dradst = dradst/(Convhs*tfract)
      dradst = dradst/4.320E4_10
c
c make limb correction to right ascension
      if(Nice.le.0) then
         if(Limb(1).ne.center) then
            if((Limb(1).eq.one) .or. (Limb(1).eq.west)) then
               sgn = 1._10
            else if(Limb(1).eq.blank) then
               sgn = 1._10
            else
               if((Limb(1).ne.two) .and. (Limb(1).ne.east))
     .            goto 100
               sgn = -1._10
            endif
            nlimb = 1
            if(sitf4.eq.rname) then
               Deriv(2,1) = Deriv(2,1)
     .                       - caps*sgn/(Convhs*(1.0_10-dradst)
     .                       *COS(deltat))
            else if(sitf4.ne.tname) then
 
               Deriv(2,1) = Deriv(2,1)
     .                       - dradst*caps*sgn/(Convhs*(1.0_10-dradst)
     .                       *COS(deltat))
               if(Limb(1).eq.blank) then
                  test1 = Deriv(2,1)
                  test2 = Deriv(2,1)
     .                    + 2.0_10*dradst*caps*sgn/(Convhs*(1.0_10-
     .                    dradst)*COS(deltat))
                  chck1 = ABS(Result(1) - test1)
                  chck2 = ABS(Result(1) - test2)
                  if(chck1.lt.chck2) then
                     Deriv(2,1) = test1
                  else
                     Deriv(2,1) = test2
                  endif
               endif
            else
               Deriv(2,1) = Deriv(2,1)
     .                       - caps*sgn/(Convhs*(1.0_10-dradst)
     .                       *COS(deltat))
               test1 = Deriv(2,1)
               test2 = Deriv(2,1)
     .                 + 2._10*caps*sgn/(Convhs*(1.0_10-dradst)
     .                 *COS(deltat))
               chck1 = ABS(Result(1) - test1)
               chck2 = ABS(Result(1) - test2)
               if(chck1.lt.chck2) then
                  Deriv(2,1) = test1
               else
                  Deriv(2,1) = test2
               endif
            endif
         endif
      endif
 
c make limb correction to declination
  100 if(Nice.ge.0) then
         if(Limb(2).ne.center) then
            if((Limb(2).eq.upper) .or. (Limb(2).eq.north) .or.
     .         (Limb(2).eq.top)) then
               sgn = -1._10
            else if(Limb(2).eq.blank) then
               sgn = -1._10
            else
               if((Limb(2).ne.lower) .and. (Limb(2).ne.south)
     .            .and. (Limb(2).ne.botm)) goto 200
               sgn = 1._10
            endif
            nlimb  = 1
            ddcdst = (z(1,1)**2*z(6,1) + z(2,1)**2*z(6,1) - z(1,1)
     .               *z(3,1)*z(4,1) - z(2,1)*z(3,1)*z(5,1))
     .               /(z(1,1)**2 + z(2,1)**2)
            ddcdst = ddcdst/SQRT(z(1,1)**2 + z(2,1)**2 - z(3,1)**2)
            ddcdst = ddcdst*15.0_10/(Convhs*tfract)
            ddcdst = ddcdst/4.320E4_10
            if(sitf4.eq.rname) then
               dellim = Result(2)*4.848136811E-6_10
               sinpdf = rrho*conspr*condis/dist*SIN(rphipr - dellim)
               Deriv(2,2) = Deriv(2,2)
     .                       - ASIN(sinpdf + sgn*SIN(caps))
     .                       /4.848136811E-6_10
            endif
            if(sitf4.eq.tname) then
               if(Series.ne.mrser) then
                  dellim = Result(2)*4.848136811E-6_10
                  sinpdf = trho*conspr*condis/dist*SIN(tphipr - dellim)
                  test1  = Deriv(2,2) - ASIN(sinpdf + sgn*SIN(caps))
     .                     /4.848136811E-6_10
                  test2  = Deriv(2,2) - ASIN(sinpdf - sgn*SIN(caps))
     .                     /4.848136811E-6_10
                  chck1  = ABS(Result(2) - test1)
                  chck2  = ABS(Result(2) - test2)
                  if(chck1.lt.chck2) then
                     Deriv(2,2) = test1
c
c go to 601
c if time is not at  transit of  1st or 2nd limb
c
                  else
                     Deriv(2,2) = test2
                  endif
               endif
            endif
            if(Limb(2).ne.top) then
               if(Limb(2).ne.botm) then
                  test1 = Deriv(2,2)
     .                    + ddcdst*caps*sgn/(Convhs*(1.0_10-dradst)
     .                    *COS(deltat))
                  test2 = Deriv(2,2)
     .                    - ddcdst*caps*sgn/(Convhs*(1.0_10-dradst)
     .                    *COS(deltat))
                  chck1 = ABS(Result(2) - test1)
                  chck2 = ABS(Result(2) - test2)
                  if(chck1.lt.chck2) then
                     Deriv(2,2) = test1
                  else
                     Deriv(2,2) = test2
                  endif
               endif
            endif
         endif
      endif
 
  200 if(nlimb.eq.0) write(6,300)
  300 format(
     .' CENTER OR UNRECOGNIZABLE CODES FOUND-ASSUME OBSERVATION IS OF CE
     .NTER AT TRANSIT OF CENTER')
      return
      end
