      real function EIONDL(theta, dmx, amx, freq, nflag)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 freq
      real      amx, degthe, del, dmx, theta, var
      integer   nflag
 
c*** end of declarations inserted by spag
 
 
c
c pfeiffer/friedman   july 1969    function eiondl
c
c
c          this function calculates the delay caused by the earth's
c     ionosphere for the one way transit of a signal.  the answer is
c     returned in sec.  the dummy variables in the order in which
c     they appear are;
c          1. theta - the zenith angle of the source in radians
c          2. dmx - the peak electron density in electrons /c.c.
c          3. amx - the altitude above the surface of the earth of
c     the peak density in km.  if it is negative, a value of 350 km.
c     is used
c          4. freq - the frequency of the source signal in hz
c          5. flag - if flag is less than 0.0, the calculation is for
c     the phase delay, otherwise the group delay is calculated
c          if theta is less than 0 or greater than pie, theta is
c     assumed to be 0.

c portable single-precision arithmetic statement function
      real*4 SNG10
      real*10 x
      SNG10(x)=x
c
      degthe = theta*57.2957795
      var    = 1.
      if(nflag.lt.0) var = -1.
      if(theta.lt.0.0 .or. degthe.gt.90.) then
         del = 97.718
      else if(degthe .lt. 35.) then
         del = 97.718 + 85.1*(1./cos(theta) - 1.)
      else if(degthe .lt. 55.) then
         del = 115.52 + 67.8*(1./cos(theta) - 1.2208)
      else if(degthe .lt. 60.) then
         del = 151.04 + 2.85*(degthe - 55.)
      else if(degthe .lt. 65.) then
         del = 165.28 + 3.43*(degthe - 60.)
      else if(degthe .lt. 80.) then
         del = 182.45 + 4.42*(degthe - 65.)
      else
         del = 248.77 + 2.6*(degthe - 80.)
      endif
      if(amx .lt. 0.0) then
         EIONDL = var*del*dmx*1.E3/(SNG10(freq)**2*2.997925)
         return
      else if(degthe .le. 46.) then
         EIONDL = (1. - 1.2E-6*(46.-degthe)*(350.-amx))
     .            *del*dmx*var*1.E3/(SNG10(freq)**2*2.997925)
         return
      else if(degthe .le. 80.) then
         EIONDL = (1. - 1.825E-5*(degthe-46.)*(350.-amx))
     .            *del*dmx*var*1.E3/(SNG10(freq)**2*2.997925)
         return
      else
         EIONDL = (1. + (degthe-80.)*(350.-amx)*1.8E-4 + .031)
     .            *del*dmx*var*1.E3/(SNG10(freq)**2*2.997925)
         return
      endif
      end
