      subroutine ROTOPS
 
      implicit none
 
      real*10 cang, sang, xin(3), xout(3)
c        r. goldstein   july 1977
c        program to do rotations for phi,psi, or inc
c           cang=cosine of phi,psi, or inc
c           sang=sine of same
c           xin=body fixed rectangular coords, of spot
c           xout=inertial coords (some frame), of spot
c        for reference see "rotation of mars" by reasenberg & king
c
c
      entry ROT3(cang, sang, xin, xout)
 
      xout(3) = xin(3)
      xout(1) = xin(1)*cang - xin(2)*sang
      xout(2) = xin(2)*cang + xin(1)*sang
      return
 
      entry ROT1(cang, sang, xin, xout)
 
      xout(1) = xin(1)
      xout(2) = cang*xin(2) + sang*xin(3)
      xout(3) = cang*xin(3) - sang*xin(2)
 
      return
      end
