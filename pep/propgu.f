      subroutine PROPGU(bthts, u, temp, nptr, rptr, npnp, np, has, nd)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, id, in, j, jpt, nd
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1974  subroutine propgu
c
c         propagate rhs from gain matrix
c        perform part of equation 20 of memo
c        g=i(k/k-1)*s
c        u(k)=(e-g)*u(k-1)
c
c         parameters
      integer*4 npnp, np
      real*10   bthts(nd, npnp), temp(np), has(npnp, npnp)
      real*10 u(np)
      integer*2 nptr(npnp), rptr(np)
c
c zero temp
      do i = 1, np
         temp(i) = 0.0
      end do
 
      id = 0
      in = 0
      do i = 1, np
         if( rptr(i) .gt. 0 ) then
c
c process noise parameters
            in = in + 1
            do j = 1, npnp
               jpt = nptr(j)
               if( jpt .gt. 0 ) temp(i) = temp(i) + has(in, j)*u(jpt)
            end do
         else
 
c deterministic parameter
            id = id + 1
            do j = 1, npnp
               jpt = nptr(j)
               if( jpt .gt. 0 ) temp(i) = temp(i) + bthts(id, j)*u(jpt)
            end do
         endif
 
      end do
 
      do i = 1, np
         u(i) = u(i) - temp(i)
      end do
      return
      end
