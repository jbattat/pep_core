      real*10 function UT2UT1(jd, fract)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   int, jd, ntab
      real      s, s2
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   oct 1968    function ut2ut1
c determine the difference ut2-ut1
c
      real*10 fract, t
      real*4 a(4,2)/-.017, .022, .006, -.007, -.012, .022, .007, -.006/
      include 'funcon.inc'
 
      t = jd
      t = t + fract - 0.5_10
      t = t - 2437665.5_10
 
c (negative before 1962 jan 1.0, postive after)
      ntab = 1
      if(t.ge.0) ntab = 2
      t   = t/365.2421988_10
      int = t
      t   = t - int
      s   = t*Twopi
      s2  = s*2.0E0
      UT2UT1 = a(1, ntab)*cos(s) + a(2, ntab)*sin(s) + a(3, ntab)
     .         *cos(s2) + a(4, ntab)*sin(s2)
 
      return
      end
