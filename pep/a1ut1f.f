      real*10 function A1UT1F(jd,fract)
 
      implicit none

c
c r.w.king    mar 1973    real*10 function a1ut1f
c
 
c calculate corrections to usno or bih a1-ut1 from analytic models

c parameters
      integer*4 jd
      real*10 fract
 
c array dimensions
      include 'globdefs.inc'

c commons
      include 'comdat.inc'
      real*10 jdut0
      equivalence (Comcon(25),jdut0)
c jdut0 is epoch for seasonal & polynomial terms of ut,
c input as 17 jan for given year.  also used in wobblf.
      include 'dtparm.inc'
      include 'empcnd.inc'
      real*10 r0,r1,r2
      equivalence (r0,Ercond(17)),(r1,Ercond(18)),(r2,Ercond(19))
      real*10 c0(2,2),s0(2,2)
      equivalence (c0(1,1),ercond(20)),(s0(1,1),ercond(24))
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'leon.inc'
c
c local
      real*10 xnum(2)/1._10,2._10/,td,a1,a2,h,tf,tfy
      integer   i,nn
 
      A1UT1F = 0.0_10
c
c determine a1-ut1 from piecewise linear function
      if(Numdt.gt.0) then
         if(Jddt0.ge.0) then
            I1 = 0
            I2 = 0
            if(jd.ge.Jddt(1)) then
               if(jd.lt.Jddt(Numdt)) then
                  do i = 1,Numdt
                     if(jd.lt.Jddt(i)) then
                        I2 = i
                        goto 10
                     endif
                  end do
                  call SUICID(' JD NOT IN DT TABLE, STOP IN A1UT1F ',9)
   10             I1 = I2 - 1
                  a2 = Dt(I2)
                  a1 = Dt(I1)
                  h  = Jddt(I2) - Jddt(I1)
                  Ff = ((jd-Jddt(I1)) + fract)/h
                  if(Jct(12).lt.0) Ff = 0._10
                  A1UT1F = A1UT1F + (a2 - a1)*Ff + a1
               endif
            endif
         endif
      endif
c
c
c determine if erotat con(11-21) are ut1 or nutation terms
      if(Jct(29).le.0) then
c
c determine time from epoch
         if(jdut0.gt.0._10) then
            td    = (jd - jdut0) + fract
            Ttutf = td/36525._10
c
c evaluate seasonal terms
c similar calc. for tcos & tsin done in wobblf.
            if(Ict(34).gt.0) then
               tfy = Twopi*MOD(td/365.2421988_10 + .04654_10, 1._10)
               nn  = Ict(34)
               do i = 1, nn
                  tf = tfy*xnum(i)
                  Tcos(i) = COS(tf)
                  Tsin(i) = SIN(tf)
                  A1UT1F  = A1UT1F + (c0(1,i) + c0(2,i)*Ttutf)*Tcos(i)
     .                      + (s0(1,i) + s0(2,i)*Ttutf)*Tsin(i)
               end do
            endif
c
c evaluate polynomial terms
            A1UT1F = A1UT1F + r0 + Ttutf*(r1 + Ttutf*r2)
         endif
      endif
      return
      end
