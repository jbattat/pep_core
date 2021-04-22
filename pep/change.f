      subroutine CHANGE(goose,mass,x,kind)
 
      implicit none
c
c m.e.ash  april 1965   subroutine change
c optimized 1977 j.f.chandler
c
c arguments
      real*10 goose,mass,x(6)
      integer kind

c eliptic orbit quantities are determined from the input
c position and velocity (x(1),...,x(6))

c common
      include 'ellips.inc'
      include 'funcon.inc'
 
c elliptic orbit quantities (encke labeled common)
      common /ENCKE/ B(3,2),Tsv,A,E,Anom0,Motpi,Mu2,Secc,Cecc,
     .        Quan2,Quan3,Quan4,Quan5,Quan12(3),Quan13(3),
     .        Quan11,Motion,Mu,Sasc,Casc,Spci,Cpci,Sinc,
     .        Cinc,Sper,Cper
      real*10 B,Tsv,A,E,Anom0,Motpi,Mu2,Secc,Cecc,Quan2,Quan3,
     . Quan4,Quan5,Quan12,Quan13,Quan11,Motion,Mu,Sasc,
     . Casc,Spci,Cpci,Sinc,Cinc,Sper,Cper

      real*10 DOT,e2,quan1,ybar(4),cf,qq4,qq5,qq6,qq7,qq8,vv,rv
      integer   i,j
 
      if(kind.gt.0) then
         Mu2 = 1.0_10 + mass
      else
         Mu2 = mass
      endif
      Mu    = goose*SQRT(Mu2)
      Mu2   = goose*goose*Mu2
      Rylpt2= DOT(x,x)
      Rylpt = SQRT(Rylpt2)
      Rylpt3= Rylpt2*Rylpt
      rv    = DOT(x,x(4))
      vv    = DOT(x(4),x(4))
      do i = 1,6
         Ylpt(i) = x(i)
      end do
      A     = Mu2/(2.0_10*Mu2/Rylpt - vv)
      Quan4 = SQRT(A)
      Motion = Mu/A/Quan4
      Sinc   = Quan4/Mu
      Quan4  = Quan4*Mu
      Secc   = rv/Quan4
      Cecc   = 1.0_10 - Rylpt/A
      e2     = Secc**2 + Cecc**2
      E = SQRT(e2)
      if(E.le.1.E-15_10) then
         E     = 0._10
         Anom0 = 0._10
         Secc  = 0._10
         Cecc  = 1._10
      else
         Anom0 = (ATAN2(Secc,Cecc) - Secc)/Twopi
         Secc  = Secc/E
         Cecc  = Cecc/E
      endif
      Motpi  = Motion/Twopi
      quan1  = 1.0_10 - e2
      Quan2  = SQRT(quan1)
      Quan3  = A*Quan2
      Quan5  = Quan2*Quan4
      Quan4  = -Quan4
      Quan11 = E/quan1
      Casc   = Cecc/Rylpt
      Sasc   = Secc/Rylpt
      Cper   = Sinc*(Cecc - E)
      Sper   = Sinc*Secc
      do i = 1,3
         j = i + 3
         B(i,1)   = Ylpt(i)*Casc - Ylpt(j)*Sper
         B(i,2)   = (Ylpt(i)*Sasc + Ylpt(j)*Cper)/Quan2
         Quan12(i) = B(i,1)*A
         Quan13(i) = B(i,2)*Quan11
      end do
      Sinc = SQRT(B(3,1)**2 + B(3,2)**2)
      Cinc = B(1,1)*B(2,2) - B(1,2)*B(2,1)
      if(Sinc.le.1.0E-15_10) then
         Casc = 1.0_10
         Sasc = 0.0_10
         Cper = B(1,1)
         Sper = -B(1,2)
         Sinc = 0.0_10
         Cinc = SIGN(1.0_10,Cinc)
      else
         Sper = B(3,1)/Sinc
         Cper = B(3,2)/Sinc
         Sasc = B(2,1)*Cper - B(2,2)*Sper
         Casc = B(1,1)*Cper - B(1,2)*Sper
      endif
      Spci = Sper*Cinc
      Cpci = Cper*Cinc
c compute partial derivatives
      ybar(1) = A*(Cecc - E)
      ybar(2) = Quan3*Secc
      ybar(3) = Quan4*Secc/Rylpt
      ybar(4) = Quan5*Cecc/Rylpt
      cf  = -Mu2/Rylpt3
      qq4 = Secc/Motion
      qq5 = A*Cecc/Rylpt
      qq6 = cf*qq4
      qq7 = -ybar(4)*Quan11
      qq8 = cf/Motion
      do i=1,3
         j = i+3

c Partials w.r.t. a
         Dylpt(i,1) = Ylpt(i)/A
         Dylpt(j,1) = -0.5_10*Ylpt(j)/A

c Partials w.r.t. e
         Dylpt(i,2) = qq4*Ylpt(j) - Quan12(i) - Quan13(i)*ybar(2)
         Dylpt(j,2) = qq5*Ylpt(j) + qq6*Ylpt(i) + qq7*B(i,2)
      end do

c Partials w.r.t. inc (radians)
      Dylpt(1,3) = Sasc*Ylpt(3)
      Dylpt(2,3) = -Casc*Ylpt(3)
      Dylpt(3,3) = Spci*ybar(1) + Cpci*ybar(2)
      Dylpt(4,3) = Sasc*Ylpt(6)
      Dylpt(5,3) = -Casc*Ylpt(6)
      Dylpt(6,3) = Spci*ybar(3) + Cpci*ybar(4)

c Partials w.r.t. asc (radians)
      Dylpt(1,4) = -Ylpt(2)
      Dylpt(2,4) = Ylpt(1)
      Dylpt(3,4) = 0._10
      Dylpt(4,4) = -Ylpt(5)
      Dylpt(5,4) = Ylpt(4)
      Dylpt(6,4) = 0._10
      do i=1,3
         j = i+3

c Partials w.r.t. per (radians measured from node)
         Dylpt(i,5) = B(i,2)*ybar(1) - B(i,1)*ybar(2)
         Dylpt(j,5) = B(i,2)*ybar(3) - B(i,1)*ybar(4)

c Partials w.r.t. anom0 (radians from periapse)
         Dylpt(i,6) = Ylpt(j)/Motion
         Dylpt(j,6) = Ylpt(i)*qq8
      end do
      return
      end
