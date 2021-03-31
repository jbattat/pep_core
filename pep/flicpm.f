      subroutine FLICPM(iconof, mfile, w, npnpx2, nepoch, itrwnd)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, iconof, ind, indsav, j, k, mfile, nepoch, npnp, 
     .          npnpx2
 
c*** end of declarations inserted by spag
 
 
c
c paul macneil january, 1978
c
c transfer initial conditions and offsets from direct access
c to permanent file
      real*10 w(npnpx2), prev(40)
      real*10 backdf(40)
      integer*2 itrwnd(99)
 
      rewind iconof
      npnp = npnpx2/2
      itrwnd(iconof) = 1
 
c skip epochs
      read(iconof)
      do i = 1, nepoch
         ind    = 2*i - 1
         indsav = ind
         read(mfile, rec = ind) w
         write(6, 50) w
   50    format(' FR', 4(5D15.5,/))
         if( i .ne. 1 ) then
            do j = 1, npnp
               w(j + npnp) = w(j) - prev(j) + w(j + npnp)
            end do
            write(iconof) w
            write(6, 60) w
   60       format(' FW', 4(5D15.5,/))
            ind = indsav
         endif
 
         do j = 1, npnp
            prev(j) = w(j)
         end do
 
      end do
c
c
c backward differences
      do k = 1, nepoch
         i   = nepoch - k + 1
         ind = 2*i - 1
         indsav = ind
         read(mfile, rec = ind) w
         do j = 1, npnp
            backdf(j) = w(j + npnp)
         end do
         if( i .ne. nepoch ) then
            do j = 1, npnp
               w(j + npnp) = w(j) - prev(j)
            end do
            write(iconof) w
         endif
 
         do j = 1, npnp
 
c include current oldsum in prev for backward differences
            prev(j) = w(j) + backdf(j)
         end do
 
      end do
 
      endfile iconof
      rewind iconof
      itrwnd(iconof) = 0
 
      return
      end
