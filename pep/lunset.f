      subroutine LUNSET(gauss, mass, jd, fract, nn)

      implicit none


c*** start of declarations inserted by spag
      integer   i, int, j, jd, k, nn

c*** end of declarations inserted by spag


c
c m.e.ash   sept 1965    subroutine lunset
c setup of initial quantities for browns mean lunar orbit
c
      real*10 gauss, mass, fract
c GAUSS = square root of gravitational constant times sun mass
c MASS  =(mass of earth+moon)/(mass of sun)
c
c common
      include 'aacoff.inc'
      include 'orbstf.inc'
      real*10 pq1,pq2,pr1,pr2,qr1,qr2,pqr,dk
      equivalence (Ww,pq1), (Ww(2),pq2), (Ww(3),pr1), (Ww(4),pr2),
     .            (Ww(5),qr1), (Ww(6),qr2), (Ww(7),pqr), (Ww(7),dk)
      include 'precmn.inc'
c
c calculation of miscellaneous quantities used in lunorb
      E22   = E/2._10
      E3    = 3._10*E
      Ert2  = 1._10 - E*E
      Ert   = SQRT(Ert2)
      Erta  = A*Ert
      Erte  = E/Ert
      Ertae = A*Erte
      Ert23 = 3._10*Ert2
      Mot3  = 3._10*gauss*SQRT(mass/A)/A
      Mot32 = Mot3/2._10
c
c calculation of coefficients daa1(k,i,j) and ddaa1(k,i,j) for
c the derivatives of the aa(i,j) matrix
      do i = 1, 3
         do j = 1, 3
            dk = 0._10
            do k = 1, 5
               dk = dk + 1._10
               Daa1(k,i,j) = dk*Aa1(k,i,j,Kepoch)
               end do
            dk = 0._10
            do k = 1, 4
               dk = dk + 1._10
               Ddaa1(k,i,j) = dk*Daa1(k+1,i,j)
               end do
            end do
         end do
c
c calculation of initial position,velotity and initial time
      Luncon = 1
      call LUNORB(jd,fract,-1)
      Luncon = 2
      T0 = Tt
      if(nn.eq.0) return

c Calculation of 3x3 matrix DINTL = the matrix of partial derivatives of
c orientation angles for mean lunar orbit relative to mean equinox and
c ecliptic at initial time with respect to orientation angles for mean
c lunar orbit relative to mean equinox and equator of reference epoch at
c initial time.  Orientation angles are inclination, longitude of
c ascending node, argument of perigee.
c
c           partial derivatives w.r.t. inclination on equator
      pq1 = Sasc*Sinc
      pq2 = -Casc*Sinc
      qr1 = Aa(1,3)*pq1 + Aa(2,3)*pq2 + Aa(3,3)*Cinc
      qr2 = SQRT(1._10 - qr1*qr1)
      pr1 = Aa(1,1)*pq1 + Aa(2,1)*pq2 + Aa(3,1)*Cinc
      pr2 = Aa(1,2)*pq1 + Aa(2,2)*pq2 + Aa(3,2)*Cinc
      pq1 = pr1/qr2
      pq2 = pr2/qr2
      pqr = qr1/qr2
      Prec(1,1) = pq1*C(3,1)
      Prec(1,2) = pq1*C(3,2)
      Prec(1,3) = pq1*qr1
      Prec(2,1) = pq2*C(3,1)
      Prec(2,2) = pq2*C(3,2)
      Prec(2,3) = pq2*qr1
      Prec(3,1) = pqr*C(3,1)
      Prec(3,2) = pqr*C(3,2)
      Prec(3,3) = -qr2

      do int = 1, 3
c
c similar calculation for all partial derivatives
         pqr = 0._10
         pq1 = 0._10
         pq2 = 0._10
         qr1 = 0._10
         qr2 = 0._10
         do i = 1, 3
            pqr = pqr - Aa(3,i)*Prec(i,3)
            pq1 = pq1 + Aa(1,i)*Prec(i,3)
            pq2 = pq2 + Aa(2,i)*Prec(i,3)
            qr1 = qr1 + Aa(3,i)*Prec(i,1)
            qr2 = qr2 + Aa(3,i)*Prec(i,2)
            end do
         Dintl(1,int) = pqr/Sinc
         Dintl(2,int) = (Casc*pq1 + Sasc*pq2)/Sinc
         Dintl(3,int) = (Cper*qr1 - Sper*qr2)/Sinc
         if(int.eq.1) then
c
c partial derivatives w.r.t. ascending node on equator
            do i = 1, 2
               Prec(1,i) = -C(2,i)
               Prec(2,i) = C(1,i)
               Prec(3,i) = 0._10
               end do
            Prec(1,3) = -pr2
            Prec(2,3) = pr1
            Prec(3,3) = 0._10
         else if(int.eq.2) then
c
c partial derivatives w.r.t. argument of perigee from equator
            do i = 1, 3
               Prec(i,1) = C(i,2)
               Prec(i,2) = -C(i,1)
               Prec(i,3) = 0._10
               end do
         endif
         end do
c
c calculation of initial values of partial derivatives of
c position,velocity
      if(nn.gt.0) call LUNBRO

      return
      end
