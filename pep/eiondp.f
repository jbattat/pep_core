      real function EIONDP(theta,dmx,amx,freq,derthe,derdmx)
 
      implicit none

c
c pfeiffer/friedman   july 1969    function eiondp
c
c
c          this function calculates the doppler shift of a signal
c     caused by the earth's ionosphere.  the answer is returned in
c     hz/hz.  the dummy variables in the order in which they appear
c     are;
c          1. theta - the zenith angle of the source in radians
c          2. dmx - the peak electron density in electrons /c.c.
c          3. amx - the altitude above the surface of the earth of
c     the peak density in km.  if it is negative, a value of 350 km.
c     is used
c          4. freq - the frequency of the source signal in hz
c          5. derthe - the derivitive of theta in radians/sec
c          6. derdmx - the derivitive of dmx in electrons/sec
c          if theta is less than 0 or greater than pie, theta is
c     assumed to be 0.

c arguments
      real*10 freq
      real*4 theta,dmx,amx,derthe,derdmx

c local variables
      real*4 degthe,dop1,dop2,err,gthe,xma

c external functions
      real*4 EIONDL

c portable single-precision arithmetic statement function
      real*4 SNG10
      real*10 x
      SNG10(x)=x
c
      degthe = theta*57.2957795
      dop1   = EIONDL(theta, derdmx, amx, freq, +1)
      xma    = amx
      if(amx .le. 0.0) xma = 350.
      if((theta .le. 0.0) .or. (degthe .gt. 90.0) ) then
         EIONDP = dop1
         return
      else
         if(degthe .le. 45.) then
            dop2 = 523.*sin(theta)/cos(theta)**2 + 6.04*(degthe)
         else if(degthe .le. 55.) then
            gthe = 18.*(degthe - 45.)/57.2957795
            dop2 = 100.55 + 4.92*(degthe - 45.) - .70*sin(gthe)
            dop2 = 10.*dop2
         else if(degthe .le. 60.) then
            dop2 = 149.7 + 6.24*(degthe - 55.)
            dop2 = 10.*dop2
         else if(degthe .le. 72.) then
            dop2 = 180.9 + 6.98*(degthe - 60.)
            dop2 = 10.*dop2
         else if(degthe .le. 73.) then
            dop2 = 2720.
         else
            dop2 = 167.5*(90. - degthe)
         endif
         if(degthe .le. 45.) then
            err = 1. + .0025 + 1.178E-5*degthe*(350. - xma)
         else
            err = 1. + (6.5*sin(theta)/cos(theta) - .7)*(350. - xma)
     .            *1.E-4
         endif
         EIONDP = dop1 + err*dop2*derthe*dmx*1.E2/
     .                   (SNG10(freq)**2*2.997925)
         return
      endif
      end
