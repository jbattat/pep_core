      subroutine PVBUS(fract,x,klbi,kick)
 
      implicit none
c
c      r.king              march 1979         subroutine pvbus
c     calculate antenna offset for pioneer-venus multiprobe bus space-
c     craft, using (adjustable) input values for the spin frequency,
c     amplitude, and phase, and (non-adjustable) input values
c     for the direction of the spin axis.
c

c parameters 
      real*10 fract,x(6)
      integer   kick
      integer*2 klbi

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      real*10 uhat(3),vhat(3)
      equivalence (Sbcom(4),uhat),(Sbcom(7),vhat)
      include 'coord.inc'
      real*10 tpvbus,spvbus,cpvbus
      equivalence (tpvbus,Angdum(5)),(spvbus,Angdum(6)),
     .            (cpvbus,Angdum(7))
      include 'empcnd.inc'
      include 'obscrd.inc'

c local variables
      real*10 phi,r,tpvm0
      integer   i
      integer*2 klb
 
      klb = klbi
      if(Nk1.eq.-1) then
         tpvm0 = (-Sbcom(1) - 2443852._10)*Secday
c this operation to reduce roundoff error assumes all bus
c tracking data were taken on the same day ( 9 dec 78 = 2443852 )
c
         r   = Pcond(8,klb)
         phi = Pcond(9,klb)
         if(kick.eq.4) then
            r   = 1.34_10*Pcond(8,klb)
            phi = (Pcond(9,klb) + 0.72082_10)
         endif
      endif
      tpvbus = fract*Secday - tpvm0
      spvbus = SIN(-Pcond(7,klb)*tpvbus + phi)
      cpvbus = COS(-Pcond(7,klb)*tpvbus + phi)
 
      do i = 1,3
         x(i) = x(i) + r*(spvbus*uhat(i) + cpvbus*vhat(i))
      end do
      return
 
      entry PVBUSV(x)
      do i = 1,3
         x(i + 3) = x(i + 3) - r*Pcond(7,klb)
     .              *(cpvbus*uhat(i) - spvbus*vhat(i))
      end do
 
      return
      end
