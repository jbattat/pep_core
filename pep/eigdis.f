      subroutine EIGDIS(a, n, messag, nmess, wr, wi)
 
      implicit none
 
 
c paul macneil  oct 13, 1978
c debug for kalman filter
 
      integer*4 n, nmess
      real*10 a(n, n), wr(n), wi(n)
      character*4 messag(nmess)
      real*10 fv1(50)
      integer*4 ierr, iv1(20)
 
      write(6, 100) messag
  100 format('0MATRIX ', 10A4)
c
c note that eispac destroys A !!!!!!!!!!!
c
c     call EISPAC(n, n, MATRIX('REAL',a), VALUES(wr,wi))
      call rg(n,n,a,wr,wi,0,fv1,iv1,fv1,ierr)
      write(6, 200)
  200 format(' EIGENVALUES (REAL AND IMAGINARY PARTS):')
      write(6, 300) wr
  300 format(' ', 8(d13.5,2x))
      write(6, 300) wi
      write(6, 400)
  400 format(' ')
 
      return
      end
