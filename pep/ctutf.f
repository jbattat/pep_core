      real*10 function CTUTF(jd,fract)
 
      implicit none

c
c m.e.ash    nov 1969    real*10 function ctutf
c determine coordinate time minus ut1 universal time
c
c parameters
      integer*4 jd
      real*10 fract

c array dimensions
      include 'globdefs.inc'

c commons
      include 'dtparm.inc'
      include 'empcnd.inc'
      real*10 c0(2,2),s0(2,2)
      equivalence (c0(1,1),Ercond(20)),(s0(1,1),Ercond(24))
      real*10 r0,r1,r2
      equivalence (r0,Ercond(17)),(r1,Ercond(18)),(r2,Ercond(19))
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'leon.inc'

c external functions
      real*10 ETUTF, UT2UT1
c
c local
      real*10 a1,a2,h,td,tf,tfy
      real*10 xnum(2)/1._10,2._10/
      integer   i,nn
c
c determine ct-ut1 from table and standard seasonal effects
      if(Numdt.gt.0) then
c
c determine et-ut2
         do i = 1, Numdt
            if(jd.ge.Jddt(i)) then
               I2 = i
               goto 50
            endif
         end do
         call SUICID(' JD NOT IN DT TABLE, STOP IN CTUTF  ', 9)
   50    I1    = I2 - 1
         a2    = Dt(I2)
         a1    = Dt(I1)
         h     = Jddt(I2) - Jddt(I1)
         Ff    = ((jd-Jddt(I1)) + fract)/h
         CTUTF = (a2 - a1)*Ff + a1
c
c evaluate seasonal terms
         if(Jct(29).le.0) then
 
c erotat con(11-21) are nutation terms if jct(29).gt.0
            if(Ict(34).gt.0) then
 
c 2435490 = 17.0 jan 1956 = jddt0
               td    = (jd - 2435490) + fract
               tfy   = Twopi*MOD(td/365.2421988_10 + .04654_10, 1._10)
               Ttutf = td/36525._10
               nn    = Ict(34)
               do i = 1, nn
                  tf = tfy*xnum(i)
                  Tcos(i) = COS(tf)
                  Tsin(i) = SIN(tf)
                  CTUTF   = CTUTF + (c0(1,i) + c0(2,i)*Ttutf)*Tcos(i)
     .                      + (s0(1,i) + s0(2,i)*Ttutf)*Tsin(i)
               end do
            endif
         endif
      else
         CTUTF = ETUTF(jd,fract) + UT2UT1(jd,fract)
      endif
c
c add polynomial terms
      if(Jct(34).gt.0) then
         td    = (jd - 2435490) + fract
         Ttutf = td/36525._10
         CTUTF = CTUTF + r0 + Ttutf*(r1 + Ttutf*r2)
      endif
 
      return
      end
