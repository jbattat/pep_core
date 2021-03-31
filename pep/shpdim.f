      subroutine SHPDIM(tl, del, itxt, ldim)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      dl, td, tst
      integer   idl
 
c*** end of declarations inserted by spag
 
 
      character*4 itxt
      real*4    tl(2), del
      integer*4 ldim
c
c        r.b. goldstein  june, 1978
c        calculate latitude and longitude dimensions of shape
c        grid from input parameters. alters input
c        parameters to create an integral grid if necessary
c        called from bodred
c
c
      td  = tl(2) - tl(1)
      dl  = td/del + .001
      idl = dl
      tst = td - float(idl)*del
      if( abs(tst) .ge. 1.E-3 ) then
 
c bad grid
         write(6, 50) itxt, del
   50    format('0', 20x, 'SHPDIM..CHANGING', a4, 'IN. OLD VALUE: ',
     .          f10.3)
         del = del + tst/float(idl)
         write(6, 100) del
  100    format(' ', 46x, 'NEW VALUE: ', f10.3)
      endif
 
      ldim = dl + 1.
      return
      end
