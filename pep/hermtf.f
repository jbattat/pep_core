      real*10 function HERMTF(yh, ntype, mtype, psb)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   mtype, ntype
 
c*** end of declarations inserted by spag
 
 
c
c m.ash,friedman,smith - sept. 1969 - hermtf
c evaluate hermite interpolation polynomial
c
      real*10 yh(12), psb
c ntype=1 position interp.    ntype=2 velocity interp.
c mtype=0 coordinate itself   mtype=1 partial of corrdinate
c
c fifth degree polynomial
      HERMTF = yh(1)
     .         + psb*(yh(2) + psb*(yh(3)+psb*(yh(4)+psb*(yh(5)+psb*yh(6)
     .         ))))
      return
      end
