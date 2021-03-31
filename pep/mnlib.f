      subroutine MNLIB(jd,fract,librt,kind)
 
      implicit none
c
c m.e.ash    oct 1968    subroutine mnlib
c determine physical libration of moon by everette fourth
c difference interpolation at half day tabular interval
c
c arguments
      integer*4 jd,kind
      real*4    librt(2,3)
      real*10 fract

c commons
      include 'comnut.inc'
      real*10 dlib
      integer*4 jdmr(3),intmr
      equivalence (Dnut,dlib),(Jder,jdmr),(Inter,intmr)
 
c everett fourth difference interpolation coefficients
      real*10 evcf(9)/1.5333333333333333333_10,-0.30_10,
     . 3.3333333333333333333E-2_10, -.58333333333333333333_10,
     . .33333333333333333333_10, -4.1666666666666666667E-2_10,
     . 5.000000000000000E-2_10, -3.3333333333333333333E-2_10, 
     . 8.3333333333333333333E-3_10/
      real*10 p(4),f1,f2,y(3,2,3), fract1
      integer   i, j, jd1, k, l, nr, ntab
c
c adjust time quantities
      jd1    = jd
      fract1 = fract
      do while( .true. )
         if(fract1.lt.0) then
            jd1    = jd1 - 1
            fract1 = fract1 + 1.0_10
         else if(fract1.eq.0) then
            goto 100
         else
            do while( fract1.ge.1.0_10 )
               jd1    = jd1 + 1
               fract1 = fract1 - 1.0_10
            end do
            goto 100
         endif
      end do
c
c determine p vector and tabular indices
  100 p(1) = fract1 - 0.5_10
      ntab = 2*intmr*(jd1 - jdmr(2)) + 8
      if(p(1).lt.0) then
         p(1) = fract1
      else
         ntab = ntab + intmr
      endif
      p(1) = p(1)/dlib
      if(p(1).lt.0) then
         p(1) = p(1) + 1.0_10
         ntab = ntab - 1
      endif
      p(3) = 1.0_10 - p(1)
      p(2) = p(1)**2
      p(4) = p(3)**2
c
c determine libration y-vectors
      do k = 1, 2
         nr = ntab + k
         do l = 1, 3
            f1 = Librat(nr+1,l)
            f1 = f1 + Librat(nr-1,l)
            f2 = Librat(nr+2,l)
            f2 = f2 + Librat(nr-2,l)
            do i = 1, 3
               j = i*3 - 2
               y(i,k,l) = evcf(j)*Librat(nr,l)
     .                      + (evcf(j+1)*f1 + evcf(j+2)*f2)
            end do
         end do
      end do
c
c everett fourth difference interpolation for libration
      do k = 1, 3
         librt(1,k) = p(3)*(y(1,1,k) + p(4)*(y(2,1,k)+p(4)*y(3,1,k)))
     .                 + p(1)*(y(1,2,k) + p(2)*(y(2,2,k)+p(2)*y(3,2,k)))
         if(kind.gt.0) librt(2,k)
     .       = (y(1,2,k)+p(2)*(y(2,2,k)*3.0_10+p(2)*y(3,2,k)*5.0_10)
     .       - (y(1,1,k)+p(4)*(y(2,1,k)*3.0_10+p(4)*y(3,1,k)*5.0_10)))
     .       /dlib
      end do
 
      return
      end
