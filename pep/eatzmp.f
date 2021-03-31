      subroutine EATZMP(i,temp,press,rh,wetzi,dryzi)

c     To do:
c     .    - be aware that xc and cco2 are fixed values in this code
c     .      may want to look into making them adjustable...

c     computes the total zenith path delay due to the Earth's atmosphere
c     using the Mendes and Pavlis model of 2004
c     Geophysical Research Letters, Vol 31, L14602, 2004
c     which I also call Paper A

c     Paper A cites work by Ciddor, 1996
c     Ciddor, P. E (1996) Refractive index of air: 
c     New equations for the visible and near infrared, 
c     Applied Optics, 35, 1566-1573
c     which I also call Paper B

c     Paper B (Ciddor) cites the work of James C. Owens
c     J. C. Owens, "Optical refractive index of air: 
c     dependence on pressure, temperature and composition," 
c     Applied Optics, 6, 51-59, 1967

c     James Battat July 20, 2005
c     jbattat@cfa.harvard.edu

      implicit none

c external variables
      integer i
      real*4 temp,press,rh
      real*4 wetzi,dryzi
c temp, press, rh are local temperature (K), pressure (mb), and relative
c     humidity (fractional)
c dryzi is the hydrostatic     zenith delay (eqn 26), in sec
c wetzi is the non-hydrostatic zenith delay (eqn 38), in sec
c i is the site number that gets passed from eatzdl.f
      
c local variables
      real*4 xc,presPa,es
      real*10 c2lat,hgt
      real*10 k0,k1,k2,k3,cco2,sigsq
      real*10 fhlam,flath,fnhlam
      real*10 w0,w1,w2,w3

c     define some parameters for the Mendes and Pavlis model
c     required to compute the hydrostatic delay
      parameter (k0=238.0185_10, k1=19990.975_10)
      parameter (k2=57.362_10,   k3=579.55174_10)
c     k0 is k0  in Paper A and has units of microns**-2
c     k1 is k1* in Paper A and has units of microns**-2
c     k2 is k2  in Paper A and has units of microns**-2
c     k3 is k3* in Paper A and has units of microns**-2

c     Now define the parameters needed for the non-hydrostatic
c     component (for the standard phase and group refractivity of 
c     water vapor).
      parameter (w0=295.235,w1=2.6422,w2=-0.032380,w3=0.004028)
c     These parameters are taken from Paper A, page 2
c     w0 is a unitless parameter
c     w1 has units of um**-2
c     w2 has units of um**-4
c     w3 has units of um**-6

c array dimensions
      include 'globdefs.inc'

c common blocks
      include 'param.inc'
      include 'sitcrd.inc'
c sitcrd.inc --> Shgt(2), Freq
c param.inc  --> Ltvel

c     Convert pressure to pascals from mBar
      presPa = press*100.0

c     A formula to convert relative humidity and temperature into 
c     the water vapor pressure at the site (surface water vapor pressure)
c     Using the formalism of Marini and Murray, 1973
c     NASA Document X-591-73-351
      es  = 100.0*rh*6.11*10**(7.5*(temp-273.15)/(237.3+(temp-273.15)))

c     es     = water vapor pressure at laser site in Pascals

c     define observing site parameters 
      c2lat = Cnrm(i)*Cnrm(i)-Snrm(i)*Snrm(i)
      hgt   = Shgt(i)
c     Cnrm and Snrm are the cos and sin of the geodetic (as opposed to the 
c     .             geocentric) latitude.
c     c2lat = cos(2*phi) where phi is the geodetic latitude of the 
c     .       observing site 
c     hgt   = height of the tracking station above sea level (in km)

c     compute sigma^2 of the observation from the global Freq
      sigsq  = 1.0E-18*(Freq/Ltvel)**2
c Freq   = frequency of observation in Hz
c Ltvel  = velocity of light in km/sec
c sigsq  = square of the wavenumber 
c     the wavenumber is the inverse of vacuum wavelength in microns

c     *************************************
c     *** COMPUTE THE HYDROSTATIC DELAY ***
c     *************************************

c     compute the CO2 content function 
      xc   = 375.0
      cco2 = 1._10+0.534E-6_10*(xc-450.0)
c     xc is the carbon dioxide content in parts per million
c        and is typically set to 375 as per IAG recommendations
c     cco2 is a function of the carbon dioxide content

c     calculate the quantity flambda which is used in the group 
c     refractivity
      fhlam  = 0.01 * ( k1*(k0+sigsq)/(k0-sigsq)**2 +
     .                  k3*(k2+sigsq)/(k2-sigsq)**2) *
     .                  cco2

c     calculate the quantity f(phi,H) which depends on the
c     latitude and height of the site
      flath = 1.0_10-0.00266*c2lat-0.00028*hgt

c     Finally, calculate the zenith hydrostatic delay
      dryzi = 2.416579E-5_10*(fhlam/flath)*presPa
c     dryzi = the zenith hydrostatic delay (eqn 26 of Paper A)
c                 and has units of meters.


c     *****************************************
c     *** COMPUTE THE NON-HYDROSTATIC DELAY ***
c     *****************************************

c     Start with the dispersion formula for the non-hydrostatic 
c     atmospheric component
      fnhlam = 0.003101*(w0 + 3.0*w1*sigsq    + 
     .                        5.0*w2*sigsq**2 + 
     .                        7.0*w3*sigsq**3)
c     fnhlam is unitless and is defined in Paper A, eqn 32

c     Finally compute the zenith non-hydrostatic delay 
      wetzi = 1.0E-6_10 * (5.316*fnhlam-3.759*fhlam) *
     .          (es/flath)
c     wetzi is the zenith, non-hydrostatic delay 
c     .       it has units of meters as defined in
c     .       Paper A, equation 38


c     ***** CONVERT THE DELAYS FROM METERS INTO SECONDS *****
c     dryzi and wetzi are in meters here but eatzdl.f wants them
c     to be in seconds.
      wetzi = wetzi/(Ltvel*1e3_10)
      dryzi = dryzi/(Ltvel*1e3_10)

      return
      end
