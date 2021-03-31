      subroutine PRENUT
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   j
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   aug 1967   subroutine prenut
c calculation of the product of the nutation and precession matrices
c
 
      include 'nutces.inc'
      include 'precmn.inc'
      include 'prtcod.inc'
 
      Pc = Dpsi(1)*Cmoblq
      Ps = Dpsi(1)*Smoblq
      do j = 1, 3
         Nutprc(1, j) = Prec(1, j) - (Pc*Prec(2,j) + Ps*Prec(3,j))
         Nutprc(2, j) = Prec(2, j) + (Pc*Prec(1,j) - Deps(1)*Prec(3,j))
         Nutprc(3, j) = Prec(3, j) + (Ps*Prec(1,j) + Deps(1)*Prec(2,j))
      end do
      return
      end
