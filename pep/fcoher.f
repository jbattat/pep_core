      subroutine FCOHER(a, n, wr, wi, zp)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, imax, j, k, lp, n
 
c*** end of declarations inserted by spag
 
 
c
c     r. w. babcock      sept. 3, 1980
c
c     calculates coherence time for kalman filter
c
c     needs ddname eispaclb defined to eispac library
c
c     a is the product of 2 real symmetric matrices, so it has
c     only real eigenvalues (garbage will be printed if this is
c     not true)
c
      include 'filtim.inc'
      include 'inodta.inc'
c
c local and calling sequence
      real*10 a(n, n), wr(n), wi(n), zp(n, n), sr, si, sum
      real*10 fv1(20)
      integer*4 iv1(20), ierr
c
c note that eispac destroys the input matrix a
c
c     call EISPAC(n, n, MATRIX('REAL',a), VALUES(wr,wi), VECTOR(zp))
      call rg(n,n,a,wr,wi,1,zp,iv1,fv1,ierr)
c
c     sort the eigenvalues (and vectors)
c     results from eispac are in approximately descending order,
c     but eigenvalues of approximately the same order of magnitude
c     may be interchanged.
c
      do j = 2, n
         if(wr(j-1).lt.wr(j)) then
            sr = wr(j)
            si = wi(j)
            do i = 1, n
               a(i, 1) = zp(i, j)
               end do
            k = j
            do while( .true. )
               wr(k) = wr(k - 1)
               wi(k) = wi(k - 1)
               do i = 1, n
                  zp(i, k) = zp(i, k - 1)
                  end do
               k = k - 1
               if(k.le.1 .or. wr(k-1).ge.sr) then
                  wr(k) = sr
                  wi(k) = si
                  do i = 1, n
                     zp(i, k) = a(i, 1)
                     end do
                  go to 100
               endif
               end do
         endif
  100    end do
c
c normalize the eigenvectors
c
      do j = 1, n
         imax = 1
         sum  = 0._10
         do i = 1, n
            sum = sum + zp(i, j)**2
            if(ABS(zp(i,j)).gt.ABS(zp(imax,j))) imax = i
            end do
         sum = SIGN(SQRT(sum), zp(imax,j))
         do i = 1, n
            zp(i, j) = zp(i, j)/sum
            end do
         end do
c
c convert eigenvalues to time in days
c
      do i = 1, n
         a(i, 1) = delta/wr(i)
         end do
c
c print results
c
      lp = 1 + (n - 1)/8
      call PAGCHK(60, 2 + 2*lp, 0)
      write(Iout, 200) wr
  200 format('0REAL AND IMAGINARY PARTS OF EIGENVALUES OF S*A ',
     .       '(SHOULD BE PURE REAL)'/(1p,8D16.7))
      write(Iout, 300) wi
  300 format(1p, 8D16.7)
      call PAGCHK(60, 2 + lp, 0)
      write(Iout, 400) (a(i,1), i = 1, n)
  400 format('0CORRESPONDING COHERENCE TIMES IN DAYS'/(1p,8D16.7))
      call PAGCHK(60, 2, 0)
      write(Iout, 500)
  500 format('0EIGENVECTORS (READ DOWN)')
      do i = 1, n
         call PAGCHK(60, lp, 0)
         write(Iout, 300) (zp(i,j), j = 1, n)
         end do
      return
      end
