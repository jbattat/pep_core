      subroutine ASTERP(y, dist, n)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 dist, y
      integer   i, n
 
c*** end of declarations inserted by spag
 
 
c j.f.chandler - 1976 july
c everett interpolation for perturbing satellites
      dimension y(5, 9, 41), dist(2)
c           y   =input y vector
c           dist=input tabular interval and its square
c                             now not needed, since no vel, acc.
c
c     common
      include 'fmmips.inc'
      include 'prtpin.inc'
c
c determine position coordinates for body n
      do i = 1, 3
         Yast(i, n) = P(1)
     .                *(y(1,i,Mtab1) + P(2)*(y(2,i,Mtab1)+P(2)*(y(3,i,
     .                Mtab1)+P(2)*(y(4,i,Mtab1)+P(2)*y(5,i,Mtab1)))))
     .                + P(3)
     .                *(y(1,i,Ntab1) + P(4)*(y(2,i,Ntab1)+P(4)*(y(3,i,
     .                Ntab1)+P(4)*(y(4,i,Ntab1)+P(4)*y(5,i,Ntab1)))))
      end do
      return
      end
