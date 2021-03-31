      subroutine SKYCOR(ra,dec)
 
      implicit none

c           subr. skycor - j.f.chandler - 1983 jan

c arguments
      real*10 ra,dec

c     correct ra and dec for star catalog errors (assumed tied to
c     one independent variable, namely, ra).
c  ra  - input right ascension in radians (corrected on return)
c  dec - input declination in radians (corrected on return)
c     note: correction coefficients are in arc seconds

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'funcon.inc'
      include 'ltrapx.inc'
      include 'shpcom.inc'
      real*10 sncs(20,2),csky(20),ssky(20)
      equivalence (csky,Leg,sncs),(ssky,Leg1)
      include 'skymap.inc'
c
c local
      real*10 cor(4),csk,ssk,sum
      integer   i,ii,k,n,nt
 
      if(Nskyc.le.0) return
      nt = (Nskyc + 3)/4
 
c compute cos(n*ra), sin(n*ra) to desired order
      csk     = COS(ra)
      ssk     = SIN(ra)
      csky(1) = csk*Convds
      ssky(1) = ssk*Convds
      do i = 2,nt
         csky(i) = csk*csky(i - 1) - ssk*ssky(i - 1)
         ssky(i) = csk*ssky(i - 1) + ssk*csky(i - 1)
      end do
c
c sum series: cosine,sine for both ra and dec
      k = 1
      do n = 1,4
         sum = 0._10
         ii  = n
         do i = 1,nt
            sum = sum + Sky(ii)*sncs(i,k)
            ii  = ii + 4
         end do
         k = 3 - k
         cor(n) = sum
      end do
 
c apply correction to input position
      ra  = ra + cor(1) + cor(2)
      dec = dec + cor(3) + cor(4)
      return
      end
