      subroutine YCOFF(y)
 
      implicit none
c
c m.e.ash  feb 1966   determination of interpolation y vectors
c update 2014 oct - j.f.chandler - add entry for 14-point interpolation
c
c argument array to fill with y vectors
c positions 6-7 are ignored for 10-point interpolation, but the
c dimension is set to 7 for compatibility with the 14-point version
      real*10 y(7,3)

c common      
      include 'tabval.inc'

c local
      real*10 f(6)
      integer*4 i,k,n,nr
 
c interpolation coefficients for 10-point interpolator
      real*10 evcf(5,5)/
     . 1.78730158730158730159_10, -0.496031746031746031746_10,
     . 0.120634920634920634921_10, -1.98412698412698412698E-2_10,
     . 1.58730158730158730159E-3_10,
     . -0.935956790123456790123_10, 0.605709876543209876543_10,
     . -0.163271604938271604938_10, 2.77998236331569664903E-2_10,
     . -2.25970017636684303351E-3_10,
     . 0.158217592592592592593_10, -0.117129629629629629630_10,
     . 4.60648148148148148148E-2_10, -8.79629629629629629630E-3_10,
     . 7.52314814814814814815E-4_10,
     . -9.75529100529100529101E-3_10, 7.60582010582010582011E-3_10,
     . -3.50529100529100529101E-3_10, 8.59788359788359788360E-4_10,
     . -8.26719576719576719577E-5_10,
     . 1.92901234567901234568E-4_10, -1.54320987654320987654E-4_10,
     . 7.71604938271604938272E-5_10, -2.20458553791887125220E-5_10,
     . 2.75573192239858906526E-6_10/

c interpolation coefficients for 14-point interpolation
      real*10 evc7(7,7)/
     . 1.95513375513375513376_10, -6.37723387723387723377E-1_10,
     . 2.05133755133755133755E-1_10, -5.43900543900543900544E-2_10,
     . 1.06893106893106893107E-2_10, -1.35975135975135975136E-3_10,
     . 8.32500832500832500833E-5_10,
     . -1.18373456790123456790_10, 8.14922839506172839506E-1_10,
     . -2.88089726631393298071E-1_10, 7.88745590828924162289E-2_10,
     . -1.57341269841269841270E-2_10, 2.01789722623055956386E-3_10,
     . -1.24158249158249158249E-4_10,
     . 2.50447530864197530864E-1_10, -1.95034722222222222222E-1_10,
     . 9.26008597883597883598E-2_10, -2.78829548794826572601E-2_10,
     . 5.80798059964726631402E-3_10, -7.62511022927689594342E-4_10,
     . 4.75823045267489711923E-5_10,
     . -2.28127755731922398596E-2_10, 1.86441798941798941799E-2_10,
     . -1.01159474206349206349E-2_10, 3.58428865373309817757E-3_10,
     . -8.10322971781305114630E-4_10, 1.11331569664902998236E-4_10,
     . -7.14193856554967666064E-6_10,
     . 9.85725308641975308642E-4_10, -8.25617283950617283951E-4_10,
     . 4.81219686948853615521E-4_10, -1.90145502645502645503E-4_10,
     . 4.83630952380952380952E-5_10, -7.16490299823633156961E-6_10,
     . 4.82253086419753086420E-7_10,
     . -1.98162177328843995506E-5_10, 1.68350168350168350168E-5_10,
     . -1.02400493025493025493E-5_10, 4.34236545347656458775E-6_10,
     . -1.21502725669392336055E-6_10, 2.00416867083533750197E-7_10,
     . -1.46137298915076692856E-8_10,
     . 1.48385565052231718902E-7_10, -1.27187627187627187627E-7_10,
     . 7.94922669922669922670E-8_10, -3.53298964410075521190E-8_10,
     . 1.05989689323022656354E-8_10, -1.92708526041859375196E-9_10,
     . 1.60590438368216145993E-10_10/
c
c     if subroutine angle is the calling program, we always have
c     ntab1=2 and ntab2=3. the same is true if subroutine partl is the
c     calling program and partl was called by subroutine optic or trnsit
c     if subroutine deldop is the calling program, we usually have
c     ntab1=2 and ntab2=3, but the values ntab1=1 and ntab2=3 or ntab1=1
c     and ntab2=2 are possible.  if subroutine partl is the calling
c     program and partl was called by subroutine radar, then we usually
c     have ntab1=2 and ntab2=3, but ntab1=1 and ntab2=3 is also possible
c     only the values tabvl(i) for i=ntab1 to ntab2+8 are used in the
c     computations in subroutine ycoff.  the dimension of y in the
c     calling program is sometimes (7,3,3) and the calling statement
c     looks like;
c           call YCOFF(y(1,1,m))
      do k = Ntab1, Ntab2
         nr = k + 4
         do n=1,4
            f(n) = Tabvl(nr+n) + Tabvl(nr-n)
         end do
         do i = 1,5
            y(i,k) = evcf(1,i)*Tabvl(nr) +
     .       (evcf(2,i)*f(1) +
     .       (evcf(3,i)*f(2) +
     .       (evcf(4,i)*f(3) +
     .       (evcf(5,i)*f(4)))))
         end do
      end do
      return
c
c fill all 7 positions in the y vector, for 14-point interpolation
      entry YCOF14(y)
      do k = Ntab1,Ntab2
         nr = k+6
         do n=1,6
            f(n) = Tabvl(nr+n) + Tabvl(nr-n)
         end do
c note: the f vector could, in principle, be dotted with evc7(2:7,i),
c but the magnitude of the coefficients drops off with increasing
c index.  therefore, roundoff errors are minimized by adding the
c terms in reverse order
         do i = 1,7
            y(i,k) = evc7(1,i)*Tabvl(nr) +
     .       (evc7(2,i)*f(1) +
     .       (evc7(3,i)*f(2) +
     .       (evc7(4,i)*f(3) +
     .       (evc7(5,i)*f(4) +
     .       (evc7(6,i)*f(5) +
     .       (evc7(7,i)*f(6)))))))
         end do
      end do
      return
      end
