      real*10 function LEGSCL(n,ih)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, if, ih, m, n
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   aug 1971    real*10 function legscl
c calculate scale factors for legendre functions
c for normalization to four pi
c
      real*10 f, quan
 
      if(ih.gt.0) then
 
         if = n - ih
         f  = if + 1
         m  = 2*ih
         do i = 2, m
            f = f*(if+i)
            end do
         quan = 4*n+2
         LEGSCL = SQRT(quan/f)
      else
         quan = 2*n+1
         LEGSCL = SQRT(quan)
      endif
 
      return
      end
