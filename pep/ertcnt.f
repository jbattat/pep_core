      subroutine ERTCNT(jd, fract, x, nvel, ntype)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, jd, jo, ntype, nvel
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   dec 1971    subroutine ertcnt
c
      real*10 fract, x(6), y(6)
c     x      = output coordinates of earth or earth-moon barycenter
c              relative to center of mass of solar system
c              determined from n-body tape by calling solcnt & ertsnt
c     output units are au and au/day
c     input quantities are
c        jd     = julian day number
c        fract  = coordinate time fraction of day
c        nvel   = 0 x(1-3) position output
c        nvel   = 1 x(4-6) velocity output
c        ntype  = 0 earth-moon barycenter relative to center of mass
c                   of solar system
c        ntype  = 1 earth relative to center of mass of solar system
c
      call ERTSNT(jd, fract, x, nvel, ntype)
      call SOTRP(jd, fract, y, nvel)
      jo = 3*nvel
      do i = 1, 3
         x(i + jo) = x(i + jo) - y(i + jo)
      end do
 
      return
      end
