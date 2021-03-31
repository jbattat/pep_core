      subroutine MONORB(ncall)
 
      implicit none

c
c ash / amuchastegui - october 1969 - subroutine monorb
c computation of forces acting on artificial lunar orbiter
c
c parameter
      integer   ncall
c        ncall=-2 setup once per iteration of a given step for partials
c        ncall=-1 setup once per iteration of a given step for motion
c        ncall= 0 setup once per step
c        ncall= k evaluate acceleration for equation k.gt.0

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'empcnd.inc'
      include 'sbthng.inc'
 
      if(ncall.lt.0) then
c
c set-up once per iteration of a given step for motion
         if(ncall.lt.-1) then
c
c set-up once per a given step for partials
         endif
      else if(ncall.ne.0) then
c
c compute acceleration for motion
         if(Kkk.gt.0) then
         endif
c
c set-up once per step
      endif
c
c compute accelerations for partials
c
      return
      end
