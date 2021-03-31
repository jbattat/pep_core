      subroutine PRDPRM(lpr, mpr, length)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, length, mice, mml
 
c*** end of declarations inserted by spag
 
 
c
c ash/friedman    nov 1968    subroutine prdprm
c
      integer*2 lpr(length), mpr(length)
      include 'wrkcompm.inc'
 
      mice = 1
      do i = 1, length
         if( lpr(i) .le. 0 ) go to 200
         Jfk = Jfk + 1
c
c mice is counter for mpr vector
         do while( mice .le. length )
            if( mpr(mice) .le. 0 ) then
 
               mice = length + 1
               go to 100
            else
               mml = mpr(mice) - lpr(i)
               if( mml .gt. 0 ) go to 100
               Lbj  = Lbj + 1
               mice = mice + 1
               if( mml .ge. 0 ) then
                  Ipts(Jfk) = Lbj
                  go to 100
               endif
            endif
         end do
  100 end do
 
  200 if( mice .le. length ) then
         do j = mice, length
            if( mpr(j) .le. 0 ) return
            Lbj = Lbj + 1
         end do
      endif
 
      return
      end
