      subroutine FRMHAR(rst, lhar, mhar, krst, ksav, kdim, n1, m1, ntop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, kdim, krst, ksav, m1, n0, n1, nm, ntop
 
c*** end of declarations inserted by spag
 
 
c
c m.e. ash february 1970 subroutine frmhar
c gravitational potential harmonic logic routine for forming
c normal equations from saved equations
c
      integer*2 lhar(kdim, 100), mhar(kdim, 100)
      real*10 rst(1000)
      include 'restor.inc'
c ntop = positive, move saved vector into restored vector
c ntop = negative, find restored row corresponding to saved row
c
      n0 = 1
 
      if( m1 .gt. 0 ) then
         if( n1 .gt. 0 ) then
 
            nm = min0(m1, n1)
            do i = 1, nm
               if( mhar(ksav,i) .gt. 0 ) then
                  Nsav = Nsav + 1
                  if( lhar(krst,i) .gt. 0 ) then
                     Nrst = Nrst + 1
                     if( Nrst .le. ntop ) rst(Nrst) = Sav(Nsav)
                  endif
               else if( lhar(krst,i) .gt. 0 ) then
                  Nrst = Nrst + 1
               endif
            end do
 
            n0 = nm + 1
            if( m1 .le. nm ) then
               if( n1 .gt. nm ) go to 100
               return
            endif
         endif
         do i = n0, m1
            if( mhar(ksav,i) .gt. 0 ) Nsav = Nsav + 1
         end do
         return
 
      else if( n1 .le. 0 ) then
         return
      endif
  100 do i = n0, n1
         if( lhar(krst,i) .gt. 0 ) Nrst = Nrst + 1
      end do
 
      return
      end
