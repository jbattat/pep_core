      real*10 function TERPF(p, y)
 
      implicit none
 
c
c m.e.ash  feb 1966   10-point everett interpolation
c revised 2014 oct to add entry point for 14-point interpolation
c
      real*10 p(4), y(7,2)
c     the dimensions of p and y in the calling program can be different
c     from the dimensions in this subroutine. for instance,
c           real*10 p(4,2),y(7,3,3).
c     in this case the calling statement might look like
c           x=terpf(p(1,m),y(1,1,n))
c     positions 6-7 are ignored here, since this is the 10-point
c     interpolator, but they are used for the 14-point interpolator
c     NOTE: if these unused positions are filled with zeroes, then
c     a call to TERP14 should give exactly the same result
c
      real*10 TERP14

c perform 10-point interpolation
      TERPF = p(3)
     .        *(y(1,1) + p(4)*(y(2,1)+p(4)*(y(3,1)+p(4)*(y(4,1)+p(4)*y
     .        (5,1))))) + p(1)
     .        *(y(1,2) + p(2)*(y(2,2)+p(2)*(y(3,2)+p(2)*(y(4,2)+p(2)
     .        *y(5,2)))))
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c perform 14-point everett interpolation
      entry TERP14(p,y)

      TERP14 = p(3)*(y(1,1)+p(4)*(y(2,1)+p(4)*(y(3,1)+p(4)*(y(4,1)+
     . p(4)*(y(5,1)+p(4)*(y(6,1)+p(4)*y(7,1)))))))
     . + p(1)*(y(1,2)+p(2)*(y(2,2)+p(2)*(y(3,2)+p(2)*(y(4,2)+
     . p(2)*(y(5,2)+p(2)*(y(6,2)+p(2)*y(7,2)))))))
      return
      end
