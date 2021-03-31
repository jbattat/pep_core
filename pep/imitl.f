      subroutine IMITL(gauss, mass, cond, kind, jx)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 g2
      integer   jx, kind
 
c*** end of declarations inserted by spag
 
 
c ash/becker  may 1968   subroutine imitl
c code moved to jnitl/jlipt  1977 dec  j.f.chandler
c initial eliptic orbit computations for target bodies
      real*10 cond(6), gauss, mass
c
c common
      include 'emcke.inc'
      include 'emmips.inc'
 
      g2 = mass
      if( kind .gt. 0 ) g2 = 1._10 + mass
      call JNITL(gauss*SQRT(g2), cond, Elptg(1,jx), 1, Dympt(1,1,jx))
      return
      end
