      subroutine FRMWTR(b, wtrans, npnp, np, nptr)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, index, j, k, l, np, npnp
 
c*** end of declarations inserted by spag
 
 
c
c paul macneil january, 1978 subroutine frmwtr
c
      real*10   b(1), wtrans(np, npnp)
      integer*2 nptr(npnp)
c form initial condition offset transformation matrix
c wtrans (subset of i(k/k-1))
c
      do i = 1, np
         k = i
         do j = 1, npnp
            l = nptr(j)
            if( k .ge. l ) index = k*(k - 1)/2 + l
            if( l .gt. k ) index = l*(l - 1)/2 + k
            wtrans(i, j) = b(index)
         end do
      end do
      return
      end
