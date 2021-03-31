      subroutine FRMPRM(rst, lprm, mprm, nsize, ntop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   l, m, nsize, ntop
 
c*** end of declarations inserted by spag
 
 
c
c m.e. ash february 1970 subroutine frmprm
c logic routine for forming normal equations from
c saved equations
c
      real*10 rst(1000)
      integer*2 lprm(100), mprm(100)
      include 'restor.inc'
c ntop = positive move saved vector into restored vector
c ntop = negative find restored row corresponding to saved row(ntop)
c
      l = 1
      m = 0
      do while( .true. )
         m = m + 1
         if( m .gt. nsize ) go to 200
         if( mprm(m) .le. 0 ) go to 200
         Nsav = Nsav + 1
         do while( l .le. nsize )
            if( lprm(l) .le. 0 ) go to 100
            if( lprm(l) .gt. mprm(m) ) go to 100
            l    = l + 1
            Nrst = Nrst + 1
 
            if( lprm(l-1) .ge. mprm(m) ) then
               if( Nrst .le. ntop ) rst(Nrst) = Sav(Nsav)
               go to 100
            endif
         end do
  100 end do
 
  200 do while( l .le. nsize )
         if( lprm(l) .le. 0 ) return
         l    = l + 1
         Nrst = Nrst + 1
      end do
 
      return
      end
