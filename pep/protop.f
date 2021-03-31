      subroutine PROTOP
 
      implicit none
 
      real*10 cang, sang, xin(3), xout(3)
c
c        r. goldstein   may 1977
c        program to apply partials of rotation operations for phi,psi,
c           or inc
c        for definition of input parameters see rotops
c
      entry PROT3(cang, sang, xin, xout)
 
      xout(3) = 0._10
      xout(2) = xin(1)*cang - xin(2)*sang
      xout(1) = -xin(1)*sang - xin(2)*cang
      return
c
c
      entry PROT1(cang, sang, xin, xout)
 
      xout(1) = 0._10
      xout(2) = xin(3)*cang - xin(2)*sang
      xout(3) = -xin(2)*cang - xin(3)*sang
c
c
      return
      end
