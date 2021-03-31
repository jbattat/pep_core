      subroutine SCALE(npz, npt, leg, gleg)
 
      implicit none

c
c r. j. cappallo     august 1970    subroutine scale
c
c scale calculates and multiplies by the Legendre normilization coef
c leg and gleg are dimensioned in calling program (pshape)
c formulas for normilization constants is from same source as
c legendre polynomial formulas as used in subroutine legndr

c arguments
      integer*4 npz, npt
      real*10 leg(1), gleg(1)

c common
      include 'funcon.inc'

c local
      logical*1 switch/.false./
      real*10 fact(0:30), a(20), b(120)
      real*10 ds2n1
      integer*4 m, n, npt1
 
      if(.not.switch) then
         switch = .true.
 
         fact(0) = 1._10
         do n = 1, 30
            fact(n) = fact(n-1)*n
         end do
         do n = 1, 20
            ds2n1 = (2*n+1)/Twopi
            a(n) = SQRT(ds2n1)
         end do
         do n = 1, 15
            ds2n1 = (2*n+1)/Twopi
            do m = 1, n
               b(n*(n-1)/2+m) = SQRT(ds2n1*fact(n-m)/fact(n+m))
            end do
         end do
      endif
      npt1 = npt*(npt+1)/2
      do n = 1, npz
         leg(n) = leg(n)*a(n)
      end do
      do n = 1, npt1
         gleg(n) = gleg(n)*b(n)
      end do
      return
      end
