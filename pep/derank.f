      subroutine DERANK(b,side,nparm,elim,eigval,temp,iuse,ierinv)
      implicit none
c derank - j.f.chandler - 2007 nov
c Solve normal equations by eigenvalue decomposition and elimination of
c a specified number of the smallest eigenvalues.

c calling parameters
      integer*4 nparm,ierinv,elim
      real*10 b(nparm,nparm,2),side(nparm),eigval(nparm,2),temp(nparm)
      integer*4 iuse(nparm)

c nparm - number of parameters
c elim  - number of eigenvalues to eliminate
c b     - on input, a scaled lower-triangular-form nparm*nparm matrix
c         used for two square matrices during solution
c         returned as inverse matrix in lower-triangular form
c side  - on input, right-hand side of normal equations
c         returned as solution vector
c eigval- work array, filled with eigenvalues (real and imaginary parts)
c temp  - work array, filled on return with eigensolutions
c iuse  - work array
c ierinv- completion code

c local variables
      integer*4 n0,n,i,j,ij,elimx
      logical change

c overrule any request to eliminate all eigenvalues
      if(elim.lt.nparm) then
         elimx=MAX(elim,0)
      else
         elimx=nparm-1
         call SUICID(
     .    'REQUESTED ELIMINATION OF ALL EIGENVALUES. KEEPING ONE ANYHOW'
     .    ,-15)
      endif

c expand matrix to square form, starting with last row
      n0=(nparm*(nparm+1))/2
      do i=nparm,1,-1
         n=n0
         do j=nparm,1,-1
            if(j.le.i) then
               b(j,i,1)=b(n,1,1)
               n=n-1
            else
               b(j,i,1)=b(i,j,1)
            endif
         end do
         n0=n0-i
      end do

c decompose into eigenvectors
      call RG(nparm,nparm,b,eigval(1,1),eigval(1,2),1,b(1,1,2),iuse,
     . temp,ierinv)
      ierinv=ABS(ierinv)

c verify all real eigenvalues
      if(ierinv.eq.0) then
         do i=1,nparm
            if(eigval(i,2).ne.0._10) ierinv=999
         end do
      endif

c normalize eigenvectors and ensure that the largest component of each
c eigenvector is positive
      do i=1,nparm
         temp(1)=0._10
         temp(2)=0._10
         do j=1,nparm
            temp(1)=temp(1)+b(j,i,2)**2
            if(ABS(b(j,i,2)).gt.ABS(temp(2))) temp(2)=b(j,i,2)
         end do
         if(temp(1).gt.0._10) temp(1)=1._10/SQRT(temp(1))
         if(temp(2).lt.0._10) temp(1)=-temp(1)
         do j=1,nparm
            b(j,i,2)=b(j,i,2)*temp(1)
         end do
      end do

c find and eliminate the smallest eigenvalues
c array of eigenvalues is already aproximately in decreasing order
c use bubble sort to enforce strict ordering and keep eigenvectors
c lined up with corresponding eigenvalues
      n=nparm+1
      change=.true.
      do while (change .and. n.gt.2)
         n=n-1
         change=.false.
         do i=2,n
            if(ABS(eigval(i,1)).gt.ABS(eigval(i-1,1))) then
               do ij=1,2
                  temp(1)=eigval(i-1,ij)
                  eigval(i-1,ij)=eigval(i,ij)
                  eigval(i,ij)=temp(1)
               end do
               do j=1,nparm
                  temp(1)=b(j,i-1,2)
                  b(j,i-1,2)=b(j,i,2)
                  b(j,i,2)=temp(1)
               end do
               change=.true.
            endif
         end do
      end do

c form inverse of eigenvalues, i.e., variances of eigenparameters
      do i=1,nparm
         if(eigval(i,1).ne.0._10) eigval(i,1)=1._10/eigval(i,1)
      end do

c multiply by eigenvector matrix to get covariance of regular parameters
      ij=0
      do i=1,nparm
         do j=1,i
            ij=ij+1
            temp(1)=0._10
            do n=1,nparm-elimx
               temp(1)=temp(1)+b(i,n,2)*b(j,n,2)*eigval(n,1)
            end do
            b(ij,1,1)=temp(1)
         end do
      end do

c transform right-hand side to orthogonal form and multiply by variances
c this is the eigensolution vector
      call PRODCT(b(1,1,2),side,temp,-nparm,nparm,1)
      do i=1,nparm
         temp(i)=temp(i)*eigval(i,1)
      end do

c multiply by eigenvector matrix to obtain the conventional solution
c (but omit the adjustments corresponding to the smallest eigenvalues)
      call PRODCT(b(1,1,2),temp,side,nparm,nparm-elimx,1)

      return
      end
