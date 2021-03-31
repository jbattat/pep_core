      subroutine RFILTR(b, buff, bthts, nparam, npnp, ipoch, wtrans,
     .                  has, nd)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, k, kode, nd
 
c*** end of declarations inserted by spag
 
 
c
c paul macneil january, 1978 modified for iteration
c
c read one epoch of filter data set
c
c parameters
      integer*4 ipoch
      integer*4 nparam
      integer*4 npnp
      real*10   b(1), bthts(nd, npnp), wtrans(nparam, npnp)
      real*10   has(npnp, npnp)
      real*10 buff(nparam)
c
c common
      include 'filtds.inc'
c
c local
      real*10 jd1, jd2
c
c id
      read(Filter) kode, jd1, jd2, i
c
c check epoch
      if( i .ne. ipoch ) call SUICID('EPOCH CHECK IN RFILTR   ', 6)
 
      if( Iconof .gt. 0 ) then
 
c i(k/k-1) matrix
         do i = 1, nparam
            read(Filter) (wtrans(i,j), j = 1, npnp)
         end do
      endif
c
c
      do i = 1, npnp
         read(Filter) (has(i,j), j = 1, npnp)
      end do
 
      do i = 1, nd
         read(Filter) (bthts(i,j), j = 1, npnp)
      end do
c
c info matrix
      do i = 1, nparam
         read(Filter) buff
         k = i*(i - 1)/2
         do j = 1, i
            b(k + j) = buff(j)
         end do
      end do
      return
      end
