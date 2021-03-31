      subroutine DELL(e, d)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      a, aa, an, b, bb, d, del, e
      integer   ix
 
c*** end of declarations inserted by spag
 
 
c
c a.freed   dec 1972   subroutine dell (millstone refraction table)
c
      dimension a(27), b(27)
 
      data an/306.123/
c     data pp,tc,tk,pw,an/984.,9.92,283.,7.79,306.123/
c     pp =absolute pressure in millibars
c     tc =temperature in degrees centigrade
c     tk =temperature in degrees kelvin   =tc+273.18
c     pw =water vapor partial pressure in millibars
c     an =surface refractivity = (refractive index-1.)*1.e6
c        = 77.6*(pp+4810.*pw/tk)/tk
c
c     a,b=refraction table
c
c        table is set up in form a+b*ns for a sequence of fixed elevatio
c        angles where ns = surface refractivity (n-units) = an based
c        on albany weather bureau data for february and august.
c
c
c        entries a (mdeg) and b (mdeg/n-unit) correspond to elevation
c        intervals of 0.5 deg between 0 deg and 6 deg elevation,
c        intervals of 2.0 deg between 6 deg and 20 deg elevation,
c        intervals of 10. deg between 20 deg and 90 deg elevation.
c        model fails if unrealistic refractivity is specified.
c
      data a/1210., 512.291, 268.304, 160., 95.9147, 57., 40.9796, 29.,
     .     20.5, 14.15, 10.2317, 7.1, 5.0, 1.26, 0.30535, .131, -.043,
     .     -.217, -.391, -.56519, -.456, -.346, -.23687, +1.2, 2.8, 3.5,
     .     4.8/
      data b/6.25, 3.4733, 2.3716, 1.78, 1.4094, 1.16, 0.98541, .87,
     .     .77, .68, .60963, .555, .511, .386, .308628, .256, .218, .19,
     .     .168, .150541, .094, .064, .0460739, .034, .026, .02, .016/
 
      if( e .gt. 20. ) then
         ix  = 18 + INT(e/10.)
         del = (amod(e,10.))*.1
      else if( e .gt. 6. ) then
         ix  = 13 + INT((e-6.)/2.)
         del = (amod(e,2.))*.5
      else
         ix  = 1 + INT(2*e)
         del = (amod(e,0.5))*2.0
      endif
      aa = a(ix) - (a(ix) - a(ix+1))*del
      bb = b(ix) - (b(ix) - b(ix+1))*del
      d  = .001*(-aa + bb*an)
      return
      end
