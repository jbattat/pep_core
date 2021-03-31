      real*10 function SKALE(index)
 
      implicit none
c
c d. white  april 1973  real*8 function skale
c
c         returns the appropriate scale factor for the apriori coeff.
c         of the parameter indicated by index:
c index .gt. 0 => index to prmter array (-30)
c index .lt. 0 => -index div 100, -index mod 100 are
c             indices to a tesseral harmonic coefficient
c index .eq.71 => site or spot radius or z comp
c index .eq.72 => site or spot long or lat
c index .gt. 80 => 1 of 6 factors for body ic's & params
c entry SKALEINIT initializes scale factors for prmter array
c
c         parameter
      integer*4 index, idummy

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'funcon.inc'
      include 'param.inc'
c
c locals and external function
      real*10 LEGSCL
      real*10 SKALEINIT
      real*10 ssk(86)
      integer*4 i, j
c
      if(index.eq.0 .or. index.gt.86) goto 100
c
      if(index.lt.0) then
c tesserals
         i = -index/100
         j = mod(-index, 100)
         SKALE = 1._10/LEGSCL(i, j)
      else
c solar system parameters
         SKALE = ssk(index)
      endif
      return
c
c anything else
  100 SKALE = 1._10
      return
 
 
c
c initialization
c
      entry SKALEINIT(idummy)
 
      SKALEINIT = 0._10
 
      do i = 1, 86
         ssk(i) = 1._10
         end do
      ssk(1)  = (8.64E4_10/(Gauss*Aultsc))**2
      ssk(2)  = 1._10/365.25_10
      ssk(3)  = (Aultsc*Ltvel/Sunrad)**2
      ssk(11) = ssk(1)/Aultsc
      ssk(12) = ssk(11)
      ssk(17) = 1._10/Convd
      ssk(18) = ssk(17)
      ssk(21) = Aultsc
      ssk(23) = ssk(11)
      ssk(71) = Aultsc*Ltvel
      ssk(72) = Aultsc/Convd
      ssk(81) = 1._10/Convd
      ssk(82) = ssk(1)
      ssk(83) = Aultsc*Ltvel
      ssk(84) = Aultsc
      ssk(85) = Aultsc/Convd
      ssk(86) = -Aultsc/Convd
 
      return
      end
