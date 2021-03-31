      real*10 function DTRPF(p,y,d)
 
      implicit none

c function dtrpf, d2trpf - j.f.chandler  1977 june
c perform numerical differentiation from position everett
c y vectors.  returns the 1st derivative (DTRPF) or 2nd
c derivative (D2TRPF)
c
c updated 2014 oct to include entry points for 14-point interpolation
c DTRP14 and D2TRP14 - same calling sequence

c arguments
      real*10 p(4),y(7,2),d
c p - p,q for desired epoch
c y - y vectors
c d - signed time interval between tabular points

c local
      real*10 D2TRPF,DTRP14,D2TRP14,s1,s2
      integer*4 kk
 
c coefficients for differentiation
c note: the fact array in COMCON is too short and cannot be used
      real*10 fact(2:7)/3._10,5._10,7._10,9._10,11._10,13._10/
      real*10 ff(2:7)/6._10,20._10,42._10,72._10,110._10,156._10/

      s1 = 0._10
      s2 = 0._10
      do kk = 5,2,-1
         s1 = s1*p(2) + y(kk,2)*fact(kk)
         s2 = s2*p(4) + y(kk,1)*fact(kk)
      end do
      DTRPF = (s1*p(2) + y(1,2) - s2*p(4) - y(1,1))/d
      return

c
c entry d2trpf
c
      entry D2TRPF(p,y,d)
 
      s1 = 0._10
      s2 = 0._10
      do kk = 5,2,-1
         s1 = s1*p(2) + y(kk,2)*ff(kk)
         s2 = s2*p(4) + y(kk,1)*ff(kk)
      end do
 
      D2TRPF = (s1*p(1) + s2*p(3))/d**2
 
      return

c-----------------------------------------------
c 14-point interpolation functions

      entry DTRP14(p,y,d)

      s1 = 0._10
      s2 = 0._10
      do kk = 7,2,-1
         s1 = s1*p(2) + y(kk,2)*fact(kk)
         s2 = s2*p(4) + y(kk,1)*fact(kk)
      end do
      DTRP14 = (s1*p(2) + y(1,2) - s2*p(4) - y(1,1))/d
      return

c
c entry d2trp14
c
      entry D2TRP14(p,y,d)
 
      s1 = 0._10
      s2 = 0._10
      do kk = 7,2,-1
         s1 = s1*p(2) + y(kk,2)*ff(kk)
         s2 = s2*p(4) + y(kk,1)*ff(kk)
      end do
 
      D2TRP14 = (s1*p(1) + s2*p(3))/d**2
 
      return
      end
