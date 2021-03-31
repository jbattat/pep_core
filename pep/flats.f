      subroutine FLATS(pflat)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 pflat
 
c*** end of declarations inserted by spag
 
 
 
c * * * implicit real*8 * * *
c
c        r.b. goldstein oct. 1978
c        code removed from harshp and put into this routine
c           set up flattening calculations, called from plrd1
c
      include 'shpcom.inc'
 
      Pflat1 = 1._10 - pflat
      Pflat2 = Pflat1**2
      Pflat3 = Pflat1*Pflat2
      Pflat4 = Pflat2**2
      Quaf12 = 1._10
      Quaf   = 1._10
      return
      end
