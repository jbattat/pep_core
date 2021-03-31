      subroutine ECLPRC(jd, fract, mnrt1)
 
      implicit none

c
c        ash/freed   oct 1968   subroutine eclprc
c     Determine mean lunar orbit orientation angles and AA trans-
c     formation matrix from ecliptic of date to equator of reference
c     epoch because mean lunar orbit package was not called.
c     Modified Aug 1972 by R.King to calculate orbit quantities for
c     lunar rotation integration.

c arguments
      integer*4 jd,mnrt1
      real*10 fract
c
c MNRT1= 0 LUNORB or ROTMAT is calling routine, determine
c          anomx, asc, per, and aa matrix
c        1 MORFN is calling routine, determine long, asc,
c          their time derivatives, and aa matrix
      include 'aacoff.inc'
      include 'funcon.inc'
      include 'orbstf.inc'
      include 'precmn.inc'

c local
      real*10 qq, ppp
      real*10 cappi, spi
      integer i, int, int1, j

c
c determination of time
      int = jd - 2415020
      qq  = fract - 0.5_10
      Tt  = int + qq
      if(mnrt1.le.0) then
c
c determination of mean anomaly
         int  = 1306*int
         int1 = int/36000
         ppp  = int - int1*36000
         ppp  = ((ppp-6.3895392E3_10) + qq*1.306E3_10 +
     .    Tt*(0.49924465_10+Tt*(6.889E-10_10+Tt*2.99E-17_10)))/3.6E4_10
         int  = ppp
         Anomx= (ppp - int)*Twopi
      else
c
c determine mean longitude and time derivatives
         int  = 1317*int
         int1 = int/36000
         ppp  = int - int1*36000
         ppp  = ((ppp+2.70434164E4_10) + qq*1.317E3_10 +
     .    Tt*(0.63965268_10+Tt*(-8.5E-11_10+Tt*3.9E-18_10)))/3.6E4_10
         int  = ppp
         Long = (ppp - int)*Twopi
         Longd  = (13.1763965268_10 - Tt*(1.7E-12_10+
     .    Tt*1.17E-19_10))*Convd
         Longdd = (1.7E-12_10 + Tt*2.34E-19_10)*Convd
c
c determine angular velocity of mod ecliptic coord system
         cappi= 3.036017687_10 + 4.361929415E-7_10*Tt
         spi  = (2.283957251E-6_10 + 9.291432628E-14_10*Tt)/365.25_10
         E1p  = spi*COS(cappi)
         E2p  = spi*SIN(cappi)
         E3p  = -(2.436499028E-4_10 + 2.946711490E-12_10*Tt)/365.25_10
      endif
c
c determination of ascending node
      ppp = 7.19953541666666666667E-1_10 +
     . Tt*(-1.47094228333333333333E-4_10 +
     . Tt*(4.325E-15_10+Tt*1.38888888888888888889E-22_10))
      int = ppp
      Asc = (ppp - int)*Twopi
      Ascd  = (-1.4709422833333333333E-4_10 +
     .        Tt*(8.65E-15_10+Tt*4.16666666666666667E-22_10))*Twopi
      Ascdd = (8.65E-15_10 + Tt*8.3333333333333333333E-22_10)*Twopi
c
c determination of perigee
      ppp = 2.08739669444444444444E-1_10 +
     . Tt*(4.56550006944444444444E-4_10 -
     . Tt*(2.58222222222222222222E-14_10+
     . Tt*8.61111111111111111111E-22_10))
      int = ppp
      Per = (ppp - int)*Twopi
c
c determination of aa matrix
      Tpr = Tt - 18262.423_10
      if(Kepoch.eq.2) Tpr = Tt - 36524.5_10
      do j = 1, 3
         do i = 1, 3
            Aa(i,j) = Tpr*(Aa1(1,i,j,Kepoch)
     .       +Tpr*(Aa1(2,i,j,Kepoch)+Tpr*(Aa1(3,i,j,Kepoch)
     .       +Tpr*(Aa1(4,i,j,Kepoch)+Tpr*(Aa1(5,i,j,Kepoch))))))
            end do
         end do
      Aa(1,1) = Aa(1,1) + 1._10
      Aa(2,2) = Aa(2,2) + Cob0(Kepoch)
      Aa(2,3) = Aa(2,3) + Sob0(Kepoch)
      Aa(3,2) = Aa(3,2) - Sob0(Kepoch)
      Aa(3,3) = Aa(3,3) + Cob0(Kepoch)
      return
      end
