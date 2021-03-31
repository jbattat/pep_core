      subroutine WFILTR(b, buff, bthts, nparam, kode, index, npnp,
     .                  wtrans, has, nd)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, k, nd
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1974  subroutine wfiltr
c paul macneil january, 1978 modified for iteration
c
c write one epoch of filter data set
c
c parameters
      integer*4 nparam, kode, index
      integer*4 npnp
      real*10   b(1), bthts(nd, npnp), wtrans(nparam, npnp)
      real*10   has(npnp, npnp)
      real*10 buff(nparam)
c
c common
      include 'filtds.inc'
c
c id
      write(Filter) kode, Jd1fil, Jd2fil, index
c
c initial condition transformation matrix
      if( Iconof .gt. 0 ) then
         do i = 1, nparam
            write(Filter) (wtrans(i,j), j = 1, npnp)
         end do
      endif
c
c
      do i = 1, npnp
         write(Filter) (has(i,j), j = 1, npnp)
      end do
c
c gain matrix
      do i = 1, nd
         write(Filter) (bthts(i,j), j = 1, npnp)
      end do
c
c propagated info matrix
      do i = 1, nparam
         buff(i) = 0.0_10
      end do
      do i = 1, nparam
         k = i*(i - 1)/2
         do j = 1, i
            buff(j) = b(k + j)
         end do
         write(Filter) buff
      end do
      return
      end
