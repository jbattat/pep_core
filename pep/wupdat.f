      subroutine WUPDAT(iconof, mfile, w, npnpx2, nepoch, itrwnd,
     .                  oldsum, pvoflg)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, iconof, ind, indsav, j, mfile, nepoch, npnp, npnpx2
 
c*** end of declarations inserted by spag
 
 
c transfer process noise parameter offsets (summed over
c previous iterations) to direct access (mfile)
c paul macneil  november, 1978
      real*10 w(npnpx2)
      real*10 oldsum(npnpx2)
      integer*2 itrwnd(99)
      logical*4 pvoflg
 
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
   50    format(' WWR', 4(5D15.5/))
         if(i.eq.1) then
 
            do j = 1, npnpx2
               w(j) = 0._10
               end do
            ind = indsav
            if(pvoflg) write(mfile, rec = ind) w
            if(pvoflg) write(6, 80) w
         else
            read(iconof) oldsum
            write(6, 60) oldsum
   60       format(' WOR', 4(5D15.5/))
            do j = 1, npnp
               if(pvoflg) w(j) = 0._10
               w(j + npnp) = oldsum(j + npnp)
               end do
            ind = indsav
            write(mfile, rec = ind) w
            write(6, 80) w
   80       format(' WWW', 4(5D15.5/))
         endif
 
         end do
c
c
      rewind iconof
      itrwnd(iconof) = 0
 
      return
      end
