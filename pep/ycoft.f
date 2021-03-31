      subroutine YCOFT(yv,body,ndim,icall,index)
 
      implicit none

c        subr. ycoft - j.f.chandler - 1979 nov
c        derived from subr. sctrp - j.f.chandler, 1977 may
c arguments
      real*10 yv(7,3,index),body(6,ndim,1)
      integer*4 ndim,icall,index
c          yv      array for interpolation y-vectors
c                  note: most interpolators use only 1-5, but the moon
c                  uses 1-7
c          body    coordinates from integration tape
c          ndim    dimension of 'body' - max. # of partials + 1
c          icall   indicates which is the calling routine
c                   1 - emtrp,  2 - mntrp,  3 - pltrp,  4 - sbtrp,
c                   5 - sctrp,  6 -(prtrp), 7 -(ertrp), 8 - sotrp,
c                   9 - ertsnt, 10 - sotrp for mercury
c          index   dimensions of y to set up (usually 3)
c
c
c common
      include 'tabval.inc'
      include 'trpcom.inc'
c local variables
      integer*4 i,idir,j,nr,nr1,ntb2p8
c
c determine y vectors
      ntb2p8 = Ntab2 + 8
c moon used 14-point interpolator, others use 10
      if(icall.eq.2) ntb2p8=ntb2p8+4
      idir   = Idirb(icall)
      nr1    = Nbtrp(icall) + idir*Ntab1
      do j = 1, index
         nr = nr1
         do i = Ntab1, ntb2p8
            Tabvl(i) = body(j,1,nr)
            nr = nr + idir
         end do
         if(icall.eq.2) then
            call YCOF14(yv(1,1,j))
         else
            call YCOFF(yv(1,1,j))
         endif
      end do
      return
      end
