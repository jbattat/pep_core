      subroutine QREAD(nit, n1, a, num)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   n1, nit, num
 
c*** end of declarations inserted by spag
 
 
c
c r.reasenberg  oct 1970 subroutine qread
c quick read for subroutine nrmfrm
c
 
c use fortran 'h' opt=2
      real*10 a(num)
      read(nit) n1, a
      return
      end
