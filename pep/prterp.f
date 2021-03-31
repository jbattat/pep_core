      subroutine PRTERP(y,dist)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 dist, DOT, y
      integer   i, j
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   may 1966   subroutine prterp
c everett interpolation for perturbing planets
c
      dimension y(5,9,41),dist(2)
c y   =input y vector
c dist=input tabular interval and its square
c
c common
      include 'prtcod.inc'
      include 'prtpin.inc'
c
c determine position coordinates for body n
      do i = 1, 3
         Xpert(i,Nptrp) = P(1)
     .                     *(y(1,i,Mtab1) + P(2)*(y(2,i,Mtab1)+P(2)*(y
     .                     (3,i,Mtab1)+P(2)
     .                     *(y(4,i,Mtab1)+P(2)*y(5,i,Mtab1))))) + P(3)
     .                     *(y(1,i,Ntab1) + P(4)
     .                     *(y(2,i,Ntab1)+P(4)*(y(3,i,Ntab1)+P(4)
     .                     *(y(4,i,Ntab1)+P(4)*y(5,i,Ntab1)))))
      end do
c
c determine velocity,acceleration for body n (if there is rel)
      if((Ler.ge.0) .or. (Nptrp.ne.Lps)) then
         if(Ler.le.0) goto 100
      endif
      do i = 4, 6
         Xpert(i,Nptrp) = (y(1,i,Mtab1) + P(2)*(y(2,i,Mtab1)+P(2)*(y(3,
     .                     i,Mtab1)+P(2)
     .                     *(y(4,i,Mtab1)+P(2)*y(5,i,Mtab1))))
     .                     - y(1,i,Ntab1) - P(4)
     .                     *(y(2,i,Ntab1)+P(4)*(y(3,i,Ntab1)+P(4)
     .                     *(y(4,i,Ntab1)+P(4)*y(5,i,Ntab1)))))/dist(1)
         if(Ler.gt.1) then
            j = i + 3
            Xpert(j,Nptrp) = (P(1)*(y(2,j,Mtab1)+P(2)*(y(3,j,Mtab1)+P(2
     .                       )*(y(4,j,Mtab1)+P(2)*y(5,j,Mtab1))))
     .                        + P(3)
     .                        *(y(2,j,Ntab1)+P(4)*(y(3,j,Ntab1)+P(4)
     .                        *(y(4,j,Ntab1)+P(4)*y(5,j,Ntab1)))))
     .                        /dist(2)
         endif
      end do
c
c determination of further quantities for body n
  100 if(Lps.le.0) then
         if(Nptrp.ne.10) then
            Rpert2(Nptrp) = DOT(Xpert(1,Nptrp),Xpert(1,Nptrp))
            Rpert(Nptrp)  = SQRT(Rpert2(Nptrp))
            Rpert3(Nptrp) = Rpert2(Nptrp)*Rpert(Nptrp)
            do i = 1, 3
               Xpert3(i,Nptrp) = Xpert(i,Nptrp)/Rpert3(Nptrp)
            end do
         endif
      endif
      return
      end
