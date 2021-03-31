      subroutine PREPDB(buff, dim, nr, file)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i
 
c*** end of declarations inserted by spag
 
 
c
c d. white  november 1974  subroutine prepdb
c
c zero direct access
c
      integer*4 nr, dim, file
      real*10 buff(dim)
 
      do i = 1, nr
         write(file, rec = i) buff
      end do
      return
      end
