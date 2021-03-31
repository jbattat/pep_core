      subroutine MULITR(b, f, roa, v, buf, n, row, col)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 col, col0, row
      integer   i, i0, ij, ik, j, j0, j1, jk, k, k1, l, n
 
c*** end of declarations inserted by spag
 
 
c
c     f. amuchastegui/m. ash - october 1968 - subroutine mulitr
c     multer multiplies two matrices  a and  b, the first  a
c     is on a disk, and the other b is a triangular symmetric
c     matrix in storage (the calculated inverse of a), to calculate
c     each column of error matrix  f = i - a*b
c     as each column of f is calculated, half row of b is cleaned
c     up by  b = b + b*f      the cleaned-up part of b is not used
c     in subsequentes calculations of columns of f
c
c
 
      include 'inodta.inc'
      real*10 roa(1), buf(1), f(1), v(1)
      real*10 b(1), sum(2)/2*0._10/
 
      do l = 1, n
         buf(l) = 0._10
      end do
      col = 0._10
 
      i0 = 0
      do i = 1, n
c
c determination of column i of error matrix
         do j = 1, n
            read(Ibuf) k, (roa(l), l = 1, n)
            sum(1) = 0._10
            sum(2) = 0._10
            k1     = 1
            ik     = i0
            do k = 1, n
               ik = ik + k1
               if(k.ge.i) k1 = k
 
c sum=sum+roa(k)*b(ik)
               call XLOAD8(roa(k))
               call XMUL8(b(ik))
               call XADD(sum)
               call XSTORE(sum)
            end do
            if(j.eq.i) call XSUB8(1._10)
            call STORND(f(j))
            f(j) = -f(j)
         end do
c
c evaluation of error matrix's row and column norm
         col0 = 0._10
         do l = 1, n
            buf(l) = buf(l) + ABS(f(l))
            col0   = col0 + ABS(f(l))
         end do
         col = MAX(col, col0)
         if(i.ge.n) then
            row = 0._10
            do l = 1, n
               row = MAX(row, buf(l))
            end do
         endif
c
c multiplication of all n rows of inverse matrix
c by column i of error matrix
         j1 = 1
         ij = i0
         j0 = 0
         do j = 1, n
            ij = ij + j1
            if(j.ge.i) j1 = j
            sum(1) = 0._10
            sum(2) = 0._10
            jk     = j0
            j0     = j0 + j
            k1     = 1
            do k = 1, n
               jk = jk + k1
               if(k.ge.j) k1 = k
 
c sum = sum+b(jk)*f(k)
               call XLOAD8(b(jk))
               call XMUL8(f(k))
               call XADD(sum)
               call XSTORE(sum)
               end do
            call XLOAD8(b(ij))
            call XADD(sum)
            call STORND(v(j))
         end do
c
c correction of lower diagonal half  of row i
c of inverse matrix
         do j = 1, i
            b(i0 + j) = v(j)
         end do
 
         rewind Ibuf
         i0 = i0 + i
      end do
 
      rewind Ibuf
      return
      end
