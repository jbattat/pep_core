      subroutine WOTSNE(b, side, buff, mparam, kode, ipoch)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ipoch, irow, j, l
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1974  subroutine wotsne
c
c         write rhs and lhs in filter format sne
c         kode=1 => write forward only (last epoch)
c         kode=2 => write backward only (first epoch)
c         kode=3 => write forward & backward (intermediate epoch)
c
c         parameters
      integer*4 kode
      integer*4 mparam
      real*10   b(1)
      real*10 side(mparam), buff(mparam)
 
      include 'filtda.inc'
      include 'filtds.inc'
c
c id record
      write(Outsne) kode, Jd1fil, Jd2fil, ipoch
c
c         write out passed norm eqn
c         (may be last epoch forward, first or some intermediate epoch
c         backwards)
c
c         rhs
      write(Outsne) side
c
c lhs
      do i = 1, mparam
         buff(i) = 0._10
         end do
      do i = 1, mparam
         irow = i*(i - 1)/2
         do j = 1, i
            buff(j) = b(irow + j)
            end do
         write(Outsne) buff
         end do
c
c if doing both forward & backward, read in forward from
c direct access and write
      if(kode.le.2) return
c
c get index to direct access
      l = (ipoch - 1)*(mparam + 1) + 1
c
c rhs
      read(Mfile, rec = l) buff
      write(Outsne) buff
c
c lhs
      do i = 1, mparam
         read(Mfile, rec = l + i) buff
         write(Outsne) buff
         end do
      return
      end
