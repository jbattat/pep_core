      subroutine FLATC(flat,lat)
 
      implicit none
c
c        r.b. goldstein oct. 1978
c        code removed from harshp and put into this routine
c        which is now called from harshp and grdshp
c        performs flattening calculations for shape models
c           removed from subr. flats 1979, j.f.chandler

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdat.inc'
      real*10 pradls
      equivalence (Comcon(100),pradls)
      include 'funcon.inc'
      include 'shpcom.inc'

c local variables
      real*10 cphis2,flat,lat,quaf1
c
c        compute flattening at observation latitude given flattening
c        coeficient pflat and equatorial radius pradls (in light seconds
c        equations from m.ash 1972-5 (blue book) 19 apr. 1972 p. 250
c
      Tphis  = TAN(lat*Convd)/Pflat2
      cphis2 = 1._10/(1._10 + Tphis**2)
      Sphis2 = 1._10 - cphis2
      Cphis  = SQRT(cphis2)
      Sphis  = Tphis*Cphis
      quaf1  = cphis2 + Pflat4*Sphis2
      Quaf2  = cphis2 + Pflat2*Sphis2
      Quaf12 = quaf1/Quaf2
      Quaf   = SQRT(Quaf12)
      flat   = pradls*(Quaf - 1._10)
      return
      end
