      subroutine FRMBDY(rst, lbdy, mbdy, nsize, ntop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, msize, nsize, ntop
 
c*** end of declarations inserted by spag
 
 
c
c m.e. ash february 1970 subroutine frmbdy
c body logic routine for forming normal equations from
c saved equations
c
      real*10 rst(1000)
      integer*2 lbdy(30), mbdy(30)
      include 'restor.inc'
 
      msize = 6
      if( nsize .lt. 0 ) msize = -nsize
 
      do i = 1, msize
         if( mbdy(i) .gt. 0 ) then
            Nsav = Nsav + 1
            if( lbdy(i) .gt. 0 ) then
               Nrst = Nrst + 1
               if( Nrst .le. ntop ) rst(Nrst) = Sav(Nsav)
            endif
         else if( lbdy(i) .gt. 0 ) then
            Nrst = Nrst + 1
         endif
      end do
 
      if( nsize .ge. 0 ) then
         msize = nsize - 6
         call FRMPRM(rst, lbdy(7), mbdy(7), msize, ntop)
      endif
 
      return
      end
