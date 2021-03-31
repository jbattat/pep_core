      subroutine YPLPCD(x,y,n,lvl,k2,kr)
 
      implicit none

c  j.f.chandler - 1976 sep
c   set up y vectors  for everett interpolation for planet parials

c arguments
      real*10 y(5,6,2,1),x(6,24,1)
      integer*4 n,lvl,k2
      integer*2 kr(n)
c  x  - input array of tabular pts
c  yv - array of everett y-vectors for 6tuples (two tabular pts)
c  n  - number of 6tuples in the arrays
c  lvl- number of coordinates within each 6tuple
c  k2 - number of output tabular pts (1 or 2)
c  kr - pointers from y-vector into input array
 
c local
      real*10 f1, f2, f3, f4
      integer   i, j, k, l, m, mx, nr
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

      do k = 1, k2
         nr = k+4
         do mx = 1, n
            m = kr(mx)
            if(m.gt.0) then
               do l = 1, lvl
                  f1 = x(l,nr+1,m)+x(l,nr-1,m)
                  f2 = x(l,nr+2,m)+x(l,nr-2,m)
                  f3 = x(l,nr+3,m)+x(l,nr-3,m)
                  f4 = x(l,nr+4,m)+x(l,nr-4,m)
                  do i = 1, 5
                     j = (i-1)*5+1
                     y(i,l,k,mx) = evcf(j)*x(l,nr,m)
     .                +(evcf(j+1)*f1+(evcf(j+2)
     .                *f2+(evcf(j+3)*f3+evcf(j+4)*f4)))
                  end do
               end do
            endif
         end do
      end do
      return
      end
