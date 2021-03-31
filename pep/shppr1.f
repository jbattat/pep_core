      subroutine SHPPR1(ll,dlyrad)
 
      implicit none

c ash/cappallo   subroutine shppr1     september 1970
c shppr1 takes partial wrt first couple of shape harmonic coeffs.
c
c arguments
      integer*2 ll
      real*10 dlyrad

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdateq.inc'
      include 'ltrapx.inc'
      integer*2 kind
      equivalence (kind,Numpar)
      include 'param.inc'
      include 'shpcom.inc'
      include 'shphar.inc'

c local variables 
      real*10 dphist
c
c only delay partials are computed
      Deriv(kind,2) = 0._10
      if(Ncode.le.2) then
         if(Nshp.le.0) then
c
c harmonic expansion for shape
c
c partial wrt  radius
            if(ll.le.1) then
               Deriv(kind,1) = dlyrad*Aultsc
               return
c
c partial wrt  j1
            else if(ll.le.2) then
               Deriv(kind,1) = Leg(1)*Aultsc*dlyrad
               return
c
c partial wrt  c11
            else if(ll.le.3) then
               Deriv(kind,1) = Gleg(1)*Cosmln(1)*Aultsc*dlyrad
               return
c
c partial wrt s11
            else if(ll.gt.4) then
c
c planet rotation rate
               call SUICID(
     .    'CANNOT CALCULATE PLANET ROTATION PARTIAL, STOP IN SHPPR1',14)
            else
               Deriv(kind,1) = Gleg(1)*Sinmln(1)*Aultsc*dlyrad
               return
            endif
         endif
c
c
c trigonometric series or altitude grid for shape
c partial wrt  radius
         do while( ll.gt.1 )
c
c partial wrt  degree of flattening
            if(ll.gt.2) then
               call SUICID(
     .    'CANNOT CALCULATE PLANET ROTATION PARTIAL, STOP IN SHPPR1',14)
            else
               dphist = 2.0_10*Tphis/(Pflat1*(1.0_10+Tphis**2))
               Deriv(kind,1) = dlyrad*Aultsc*Pradls/Quaf*(Quaf12*Pflat1
     .                          *Sphis2 - 2.0_10*Pflat3*Sphis2 +
     .                          Cphis*Sphis*((Pflat4-1.0_10)
     .                          -Quaf12*(Pflat2-1.0_10))*dphist)/Quaf2
               return
            endif
         end do
         Deriv(kind,1) = dlyrad*Quaf*Aultsc
      endif
c note: the foregoing is equivalent to  the following -
c (why isn't it used?)
c deriv(kind,1)=-dlyrad*aultsc*pradls*quaf*sphis2*pflat1/quaf2
c
      return
      end
