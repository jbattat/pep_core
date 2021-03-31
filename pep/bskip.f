      subroutine BSKIP(mat, mparam)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   jmat, mat, mparam, nsav, nsav0
 
c*** end of declarations inserted by spag
 
 
c
c r.reasenberg  oct 1970  subroutine bskip
c
c
c        mat.gt.0  skip over rhs and matrix records
c        mat.lt.0  skip over matrix record only
c        read from fortran logical iabs(mat)
c
      jmat = iabs(mat)
      nsav = 0
      if( mat .lt. 0 ) nsav = 1
  100 do while( .true. )
         read(jmat, err=200) nsav0
         if( nsav0 .ne. -1 ) go to 300
      end do
  200 read(jmat) nsav0
  300 if( nsav .gt. 0 ) then
         if( nsav0 .lt. mparam ) go to 100
         return
      else
         nsav = 1
         go to 100
      endif
      end
