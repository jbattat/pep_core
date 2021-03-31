      real*10 function ARGMNT(angle)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   int
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash  sept 1967   function subroutine argmnt
c get angle in revolutions between -twopi and twopi radians
c
      include 'funcon.inc'
      real*10 angle, tint
c angle =  input angle in revolutions
c
      int    = angle
      tint   = int
      ARGMNT = (angle - tint)*Twopi
      return
      end
