      subroutine YNUTLB(nlb,ynlb,k2,kiss)
 
      implicit none
c
c m.e.ash   aug 1968   subroutine ynutlb
c calculate y-vectors for nutation and libration fourth difference
c interpolation
c
c arguments
      real*4 nlb(120)
      real*10 ynlb(3,2,2)
      integer k2,kiss

c local
      real*10 f1, f2
      real*10 evcf(9)/1.5333333333333333333_10,-0.3_10,
     . 3.3333333333333333333E-2_10, -0.58333333333333333333_10,
     . 0.33333333333333333333_10, -4.1666666666666666667E-2_10,
     . 0.05_10, -3.3333333333333333333E-2_10,
     . 8.3333333333333333333E-3_10/
      real*10 fact(3)/0._10,3._10,5._10/
      integer   i, j, k, nr
 
      do k = 1,k2
         nr = k + 4
         f1 = nlb(nr+1)
         f1 = f1 + nlb(nr-1)
         f2 = nlb(nr+2)
         f2 = f2 + nlb(nr-2)
         do i = 1,3
            j = 3*i - 2
            ynlb(i,1,k) = evcf(j)*nlb(nr) + evcf(j + 1)
     .                      *f1 + evcf(j + 2)*f2
         end do
      end do
 
      if(kiss.gt.1) then
         do k = 1,k2
            ynlb(1,2,k) = ynlb(1,1,k)
            do i = 2,3
               ynlb(i,2,k) = ynlb(i,1,k)*fact(i)
            end do
         end do
      endif
 
      return
      end
