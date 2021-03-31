      real*10 function FUNCOF(t, a, intb, b, c, d)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   int, intb
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash  sept 1967   function subroutine funcof
c third order polynomials for fundamental arguments of brown lunar
c theory are evaluated
c
      real*10 t, a, b, c, d, tint, tt, ab, bb, poly
c
c reduce number of revolutions
      int  = t
      tint = int
      tt   = t - tint
      int  = int*intb
      ab   = int - ((int/10000)*10000)
      bb   = intb
c
c evaluate polynomial
      poly = a + (ab + bb*tt)*1.0E-4_10 + t*(b + t*(c+t*d))
c
c get between -1 and 1 revolutions
      int    = poly
      tint   = int
      FUNCOF = poly - tint
      return
      end
