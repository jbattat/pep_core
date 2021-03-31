      subroutine LVTHAR(lharx,mharx,lhar,klnsiz,klnhar,n1,n2,
     .                  nsize,msize,maxsiz)
 
      implicit none

c           j.f.chandler - 1980 dec - subr. lvthar
c           based on m.e.ash routines lvthar,thrvct of 1969 june
c           calculate lharx vector from mharx,lhar vectors
c
c        parameters
      integer*4 klnsiz,klnhar,n1,n2
      integer*2 lharx(99),mharx(99),lhar(klnsiz,99),nsize,msize
c lharx - controls for gravitational harmonic partials, output
c mharx - ditto from input obslib tape
c lhar  - ditto from input stream
c klnsiz- dimension of lhar array
c klnhar- index into lhar array
c n1    - starting index for lhar copy
c n2    - stopping index for lhar copy
c nsize - returned index to last meaningfull item in lharx
c msize - size of mharx array
c maxsiz- limiting size of lharx array (apply to n2 for target bodies)

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdat.inc'

c local
      integer   i,maxsiz,n2t,nsizt

      if(n2.gt.0) then
 
c copy from lhar (input controls)
         n2t = n2
         if(n2t.gt.maxsiz) n2t = maxsiz
         nsizt = 0
         do i = n1,n2t
            if(lhar(klnhar,i).gt.0) then
               lharx(i) = 1
               nsizt    = i
            endif
         end do
         if(nsizt.gt.nsize) nsize = nsizt
      endif
      if(Iabs1.gt.0) then
c
c check for items in lhar but not mharx
         if(Lnotm.le.0 .and. nsize.gt.0) then
            do i = 1,nsize
               if(lharx(i).gt.0) then
                  if(i.gt.msize .or. mharx(i).le.0) then
                     Lnotm = 1
                     goto 50
                  endif
               endif
            end do
         endif
c
c copy from mharx (old obslib tape)
   50    if(msize.gt.0) then
            do i = 1,msize
               if(mharx(i).gt.0) lharx(i) = 1
            end do
            if(msize.gt.nsize) nsize = msize
         endif
      endif
 
      return
      end
