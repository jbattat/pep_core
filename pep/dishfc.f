      subroutine DISHFC(alph, amu, anu, fz, fy)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      cos2a, cos3a, cosa, sin2a, sin3a, sina
 
c*** end of declarations inserted by spag
 
 
      real*4    alph, amu, anu, fz, fy
c
c        r.b. goldstein  dec 1978
c        remove in-line code from dish and make into subroutine
c        calculates equation 7 of tm391-240
c
      include 'dish1.inc'
 
      cosa  = cos(alph)
      sina  = sin(alph)
      cos2a = cosa*cosa - sina*sina
      sin2a = 2*sina*cosa
      cos3a = cosa*cos2a - sina*sin2a
      sin3a = sina*cos2a + cosa*sin2a
 
      fy = A(5)*amu*sina + A(3)*anu*sin2a + A(4)*amu*sin3a
 
      fz = A(1)*anu + (1. + A(2)*amu)*cosa + A(3)*anu*cos2a + A(4)
     .     *amu*cos3a
 
      return
      end
