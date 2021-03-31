      subroutine PRDHAR(lhar, klnsiz, klnhar, mhar, length)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, klnsiz
 
c*** end of declarations inserted by spag
 
 
c
c     m.e.ash   june 1969    subroutine prdhar
c     pick out derivatives w.r.t. gravitational potential harmonic coef.
c     also for target body initial conditions and harmonics.
c     this subroutine differs from prdbdy only in lhar being doubly
c     dimensioned and in allowing gaps in the mhar vector compared to
c     the lhar vector for target bodies.
c
      integer*2 lhar(klnsiz, 100), mhar(100), length, klnhar
c
c common
      include 'wrkcompm.inc'
 
      if( length .gt. 0 ) then
         do i = 1, length
            if( klnhar .gt. 0 ) then
               if( lhar(klnhar,i) .gt. 0 ) then
                  Jfk = Jfk + 1
                  if( mhar(i) .gt. 0 ) then
                     Lbj = Lbj + 1
                     Ipts(Jfk) = Lbj
                  endif
                  go to 50
               endif
            endif
            if( mhar(i) .gt. 0 ) Lbj = Lbj + 1
   50    end do
      endif
 
      return
      end
