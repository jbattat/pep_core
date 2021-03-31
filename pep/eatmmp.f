      real function EATMMP(i, wetz, dryz, zen)
c     the old input parameters were (active, elev, mapping)

c     This mapping function depends only on:
c     1. Temperature
c     2. Orthometric height of the ranging station
c     3. Latitude of the ranging station
c     4. (of course) Elevation Angle of the line of sight

c     This code is called from eatctl.f
c     depending on the value of Ict(23)
c     Ict(23) is binary coded
c                  Mapping    Zenith
c            0       EATMDL   EATZDL
c            1       EATMMP   EATZDL
c            2       EATMDL   EATZMP
c            3       EATMMP   EATZMP

c     Here we compute the mapping function 
c     using the Mendes,Prates,Pavlis^2,Langley model of 2002
c     Geophysical Research Letters, Vol 29, No. 10 1414, 2004
c     Improved Mapping Functions for Atmospheric Refraction 
c     Correction in SLR.
c     which I also call Paper X
c     In this paper, they present two models:  
c        FCULa is for data that has location + meteorology
c        FCULb is for data that has only location and day of year
c     Here we only compute FCULa
c     so we do require meteorology.

c     James Battat July 22, 2005
c     jbattat@cfa.harvard.edu

c     INPUTS:
c     i    = site number as determined in eatctl.f
c     wetz = zenith delay of the wet component of the atm. (sec)
c     dryz = zenith delay of the dry component of the atm. (sec)
c     zen  = zenith angle of line of sight  (radians)

      implicit none

c     external variables
      real*4   zen, elevRad
      real*4   wetz, dryz
c     zen is the zenith angle in radians

c     local variables
      integer  i, is
      logical  smlsav
      real*4   tempCel, hgtOrth, xc
      real*10 a1, a2, a3
      real*10 a10,a11,a12,a13
      real*10 a20,a21,a22,a23
      real*10 a30,a31,a32,a33
      real*10 cLat, numer, denom, mapping 

c     THESE PARAMETERS ARE USED FOR FCULa ONLY!
c     define some parameters used to compute the actual a_i
c     parameters in the continued fraction mapping function
      parameter(a10=12100.8E-7_10,a11=1729.5E-9_10)
      parameter(a12=319.1E-7_10,  a13=-1847.8E-11_10)
      parameter(a20=30496.5E-7_10,a21=234.6E-8_10)
      parameter(a22=-103.5E-6, a23=-185.6E-10_10)
      parameter(a30=6877.7E-5_10, a31=197.2E-7_10)
      parameter(a32=-345.8E-5_10, a33=106.0E-9_10)
c     The parameter values are given in Table 1 of Paper X.


c array dimensions
      include 'globdefs.inc'

c COMMON BLOCKS
      include 'obscrd.inc'
      include 'param.inc'
      include 'sitcrd.inc'
c     sitcrd.inc --> Cnrm(2), Snrm(2), Shgt(2), Freq
c     param.inc  --> Ltvel
c     obscrd.inc --> Save()

c     define the local temperature
      is = 0
      if (i .eq. 2) is = 3
      smlsav = (Numsav .lt. is + 43)
      tempCel = Save(is + 41)
      if( smlsav .or. tempCel .lt. 150. .or. tempCel .gt. 350. ) 
     . tempCel = 273.15
c     convert from Kelvin to celsius
c     tempCel = temperature in degrees Celsius
      tempCel = tempCel - 273.15

c     define observing site parameters 
      hgtOrth = Shgt(i)*1000.0
      elevRad = 1.57079632679-zen
      cLat    = Cnrm(i)
c     Shgt(i) = height of tracking station above reference ellipsoid (km)
c     hgtOrth = orthometric height of the station in meters
c     orthometric height is the height above mean sea level
c     because we do not know the distance between the ellipsoid and 
c     the geoid, we take the orthometric height to be the height above 
c     the ellipsoid (i.e. hgtOrth=Shgt(i)*1000)
c     elevRad = elevation angle in radians

c     **************************
c     *** COMPUTE a1, a2, a3 ***
c     **************************
      a1 = a10 + a11*tempCel + a12*cLat + a13*hgtOrth
      a2 = a20 + a21*tempCel + a22*cLat + a23*hgtOrth
      a3 = a30 + a31*tempCel + a32*cLat + a33*hgtOrth
c     a_i are computed from Equation 5 in Paper X

c     Finally, calculate the FCULa mapping function 
c     via the continued fraction method
      numer = 1+(a1/(1+(a2/(1+a3))))
      denom = SIN(elevRad)+(a1/(SIN(elevRad)+(a2/(SIN(elevRad)+a3))))
      mapping = numer/denom
c     mapping = m(epsilon) in Equation 4 of Paper X
c              and is a unitless scaling of the total zenith delay

      EATMMP = mapping*(wetz+dryz)

      return
      end
