      subroutine TIMINC(jd, fr, jdn, frn, t)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   ijd, jd, jdn
 
c*** end of declarations inserted by spag
 
 
c
c        r.goldstein    dec. 1975
c        adds t to jd.fr to give jdn.frn
c
c        modified july,1977 to allow no restrictions on t. t now can
c        range from -infinity  to  +infinity
c
c        future modifications will allow normalization of jd.fract
c        where fract is out of range
c
      real*10 fr, frn, t
 
      frn = fr + t
      if( frn .ge. 1._10 ) then
 
         ijd = frn
         frn = frn - ijd
         jdn = jd + ijd
         return
      else if( frn .lt. 0._10 ) then
 
         ijd = frn
         frn = 1._10 + frn - ijd
         jdn = jd + ijd - 1
         return
      else
         jdn = jd
         return
      endif
 
      end
