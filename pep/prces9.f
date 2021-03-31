      subroutine PRCES9(tt, prec, moblq, tref)

      implicit none


c*** start of declarations inserted by spag
      integer   i

c*** end of declarations inserted by spag


c
c m.ash      sept 1972    subroutine prces9
c calculate precession matrix and mean obliquity
c
      real*10 prec(3,3), moblq, tt, tref
      real*10 z(3), cz(3), sz(3), q1, q2, tcntry, tf, t, DOTN,
     .       zcof(3,3), ocof(3), ttrm(3), oblq1
      include 'funcon.inc'

      if(tref.gt.2433283._10) then
         tcntry = 36525._10
         tf     = (tref - 2451545._10)/tcntry
         oblq1  = (84381.448_10 - 46.8150_10*tf - 0.00059_10*tf**2 +
     .          0.001813_10*tf**3)*Convds
         ocof(1) = -46.8150_10 - 0.00117_10*tf + 0.005439*tf**2
         ocof(2) = -0.00059_10 + 0.005439_10*tf
         ocof(3) = 0.001813_10
c
         zcof(1, 1) = 2306.2181_10 + 1.39656_10*tf - 0.000139_10*tf**2
         zcof(1, 2) = 0.30188_10 - 0.000345_10*tf
         zcof(1, 3) = 0.017998_10
         zcof(2, 1) = zcof(1,1)
         zcof(2, 2) = 1.09468_10 + 0.000066_10*tf
         zcof(2, 3) = 0.018203_10
         zcof(3, 1) = 2004.3109_10 - 0.85330_10*tf - 0.000217_10*tf**2
         zcof(3, 2) = -0.42665_10 - 0.000217_10*tf
         zcof(3, 3) = -0.041833_10
      else

c use old expressions (tropical centuries, etc.)
         tcntry = 36524.21988_10
         oblq1 = 84404.8363494530058_10*Convds
         zcof(1,1) = 2304.948_10
         zcof(1,2) = 0.302_10
         zcof(1,3) = 0.0179_10
         zcof(2,1) = zcof(1,1)
         zcof(2,2) = 1.093_10
         zcof(2,3) = 0.0192_10
         zcof(3,1) = 2004.255_10
         zcof(3,2) = -0.426_10
         zcof(3,3) = -0.0416_10
         ocof(1) = -46.8485418486769944_10
         ocof(2) = -3.18487539448290358E-3_10
         ocof(3) = 1.80988402545302038E-3_10
      endif
c
c compute precession
      t = (tt-tref)/tcntry
      ttrm(1) = Convds*t
      ttrm(2) = ttrm(1)*t
      ttrm(3) = ttrm(2)*t
      call PRODCT(zcof, ttrm, z, 3, 3, 1)
      moblq = oblq1 + DOTN(ocof, ttrm, 3)
c
      do i = 1, 3
         cz(i) = COS(z(i))
         sz(i) = SIN(z(i))
      end do
c
c prec(j,k)= precession matrix  j,k=1,2,3
      q1 = cz(2)*cz(3)
      q2 = sz(2)*cz(3)
      prec(1,1) = cz(1)*q1 - sz(1)*sz(2)
      prec(1,2) = -sz(1)*q1 - cz(1)*sz(2)
      prec(1,3) = -cz(2)*sz(3)
      prec(2,1) = cz(1)*q2 + sz(1)*cz(2)
      prec(2,2) = cz(1)*cz(2) - sz(1)*q2
      prec(2,3) = -sz(2)*sz(3)
      prec(3,1) = cz(1)*sz(3)
      prec(3,2) = -sz(1)*sz(3)
      prec(3,3) = cz(3)

      return
      end
