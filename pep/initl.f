      subroutine INITL(gauss, mass, cond, kind)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 Elpte, g2
      integer   kind
 
c*** end of declarations inserted by spag
 
 
c m.e.ash   april 1965   initial elliptic orbit computations
c optimized 1977 j.f.chandler
c code moved to jnitl/jlipt  1977 dec  j.f.chandler
      real*10 cond(6), gauss, mass
      include 'ellips.inc'
 
c elliptic orbit quantities (encke labeled common)
      common /ENCKE / Elpte(31)
 
      g2 = mass
      if( kind .gt. 0 ) g2 = 1._10 + mass
 
c note: the 2 parameter to jnitl is only to match routine change
      call JNITL(gauss*SQRT(g2), cond, Elpte, 2, Dylpt)
      return
      end
