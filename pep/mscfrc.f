      real*10 function MSCFRC(kk,gamat3,sbforc)

      implicit none

c m.ash  jan 1993   evaluate mascon force on satellite plus
c        central body force reduced by fractional mass of masscons

c if there are masscons, this function subroutine MSCFRC called from SBFN1 
c where central force is evaluated

      integer*4 kk
      real*10 gamat3,sbforc(3)

c input:
c    kk         = equation number (1, 2, or 3)
c    gamat3     = -(Gauss**2)*Mass(ncentr) for zero mass spacecraft
c                  plus time variation of gravitational constant Gauss
c    sbforc(kk) = kk coordinate of satellite / r**3
c
c output:
c    MSCFRC  =  function return with reduced central force + masscon force

      include 'mascon.inc'
      include 'mascn1.inc'

      integer*4 i,j
      real*10 sum, DOT

c*****************************
c input quantites to entry point MSCFR1
      real*10 sbcor(3),cntrot(3,3),MSCFR1
c sbcor  = satellite coordinates relative to central body
c cntrot = central body rotation matrix

c quantities setup at entry point MSCFR1
      real*10 zsbmsc(3,100),rsbms2,rsbms3(100)
c zsbmsc(1-3,j) = vector pointing from masscon j to satellite in inertial
c                 frame (masscon coordinates transformed from central 
c                 body frame for this computation)
c*****************************

c function MSCFRC not called from SBFN1 unless Nmmsc1 > 0

      sum  = 0._10
      do j=1,Nmmsc1
         sum  = sum + Masmsc(j)*zsbmsc(kk,j)/rsbms3(j)
         end do

      MSCFRC  = gamat3*(Fcnmsc*sbforc(kk) + sum)
      return

c***************************************************************************

c computations done only once for a given iteration of a given step
c called from SBFN to setup quantities for function call of MSCFRC from SBFN1

      entry MSCFR1(sbcor,cntrot)

      MSCFR1=0._10

      do j=1,Nmmsc1
         do i=1,3
            zsbmsc(i,j) = sbcor(i) - DOT(cntrot(1,i), Xmsc(1,j))
            end do
         rsbms2  = DOT(zsbmsc(1,j),zsbmsc(1,j))
         rsbms3(j) = rsbms2 * SQRT(rsbms2)
         end do

      return
      end
