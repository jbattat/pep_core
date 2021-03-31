      subroutine CHANGE(goose,mass,x,kind)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i,j
 
c*** end of declarations inserted by spag
 
 
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
     .        Cinc,Sper,Cper,Vv,Rv,Ry2
      real*10 B,Tsv,A,E,Anom0,Motpi,Mu2,Secc,Cecc,Quan2,Quan3,
     . Quan4,Quan5,Quan12,Quan13,Quan11,Motion,Mu,Sasc,
     . Casc,Spci,Cpci,Sinc,Cinc,Sper,Cper,Vv,Rv,Ry2

      real*10 DOT,e2,quan1
 
      if(kind.gt.0) then
         Mu2 = 1.0_10 + mass
      else
         Mu2 = mass
      endif
      Mu    = goose*SQRT(Mu2)
      Mu2   = goose*goose*Mu2
      Ry2   = DOT(x,x)
      Rylpt2= Ry2
      Rylpt = SQRT(Ry2)
      Rylpt3= Ry2*Rylpt
      Rv    = DOT(x,x(4))
      Vv    = DOT(x(4),x(4))
      do i = 1,6
         Ylpt(i) = x(i)
      end do
      A     = Mu2/(2.0_10*Mu2/Rylpt - Vv)
      Quan4 = SQRT(A)
      Motion = Mu/A/Quan4
      Sinc   = Quan4/Mu
      Quan4  = Quan4*Mu
      Secc   = Rv/Quan4
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
      return
      end
