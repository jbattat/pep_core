      subroutine YPRTCD(x,y,n)
 
      implicit none

c
c m.e.ash   may 1966   subroutine yprtcd
c computation of everett interpolation y-vectors
c
c arguments
      real*10 x(6,84),y(5,9,41)
      integer*4 n
 
c common
      include 'prtpin.inc'

c local
      integer   i, j, k, l, m, nr
      real*10 f1, f2, f3, f4
 
c interpolation coefficients
      real*10 evcf(25)/
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
      real*10 fact(2:9)/3._10,5._10,7._10,9._10,2._10,4._10,6._10,8._10/
c
c determination of position y-vectors
      do k = 1,n
         nr = k+4
         do l = 1,3
            f1 = x(l,nr+1)+x(l,nr-1)
            f2 = x(l,nr+2)+x(l,nr-2)
            f3 = x(l,nr+3)+x(l,nr-3)
            f4 = x(l,nr+4)+x(l,nr-4)
            do i = 1,5
               j = (i-1)*5+1
               y(i,l,k) = evcf(j)*x(l,nr)
     .                     +(evcf(j+1)*f1+(evcf(j+2)
     .                      *f2+(evcf(j+3)*f3+evcf(j+4)*f4)))
            end do
         end do
      end do
c
c determination of velocity,acceleration y-vectors (if rel.)
      if((Ler.ge.0) .or. (Nptrp.ne.Lps)) then
         if(Ler.le.0) return
      endif
c future option to calculate velocity y-vectors from
c coordinates if they exist rather than by numerical
c differention
      do l = 1,3
         j = l+3
         m = j+3
         do k = 1,n
            y(1,j,k) = y(1,l,k)
            do i = 2,5
               y(i,j,k) = y(i,l,k)*fact(i)
               y(i,m,k) = y(i,j,k)*fact(i+4)
            end do
         end do
      end do
      return
      end
