      subroutine LIBCND(jdb, bcond, l, jd, cond)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   jd, jdb
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   jan 1970    subroutine libcnd
c move adjusted initial condtions from observation library tape
c into correct storage locations
c
      real*10 bcond(30), cond(30)
      integer*2 l(30)
 
      integer*4 i, j
 
      if( bcond(1) .ne. 0.0_10 ) then
         do i = 1, 6
            if( l(i) .gt. 0 ) then
               do j = 1, 6
                  cond(j) = bcond(j)
               end do
               jd = jdb
               return
            endif
         end do
      endif
 
      return
      end
