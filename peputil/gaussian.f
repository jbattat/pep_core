      double precision function gaussian(idummy)
      implicit none
      integer idummy
      real*8 v1,v2,s, random

c Return gaussian distributed random value with 0 mean and unit variance

      s=99d0
      do while (s.ge.1d0 .or. s.le.0d0)
         v1=2d0*random(idummy)-1d0
         v2=2d0*random(idummy)-1d0
         s=v1**2+v2**2
      end do
      gaussian=v1*dsqrt(-2d0*dlog(s)/s)
      return
      end
