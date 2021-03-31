      subroutine MULIT2(h, hinv, f, v, buf, n, row, col)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, k, l
 
c*** end of declarations inserted by spag
 
 
c
c d. white  subroutine mulit2  september 1974
c
c modified version of mulitr to do cleanup on calculated
c inverse, where both matrices are full n*n and in core
c
c parameters
      integer*4 n
      real*10 h(n, n), hinv(n, n), f(n), v(n), buf(n)
      real*10 row, col
c
c local
      real*10 col0
      real*10 sum(2)
 
      do l = 1, n
         buf(l) = 0._10
         end do
c
c determination of column i of error matrix
      col = 0._10
      do i = 1, n
         do j = 1, n
            sum(1) = 0._10
            sum(2) = 0._10
            do k = 1, n
 
c sum = sum + hinv(i, k)*h(k, j)
               call XLOAD8(h(k,j))
               call XMUL8(hinv(i,k))
               call XADD(sum)
               call XSTORE(sum)
               end do
            call STORND(f(j))
            if(j.eq.i) then
               call XLOAD(sum)
               call XSUB8(1._10)
               call STORND(f(j))
            endif
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
         do j = 1, n
            sum(1) = 0._10
            sum(2) = 0._10
            do k = 1, n
 
c sum = sum + f(k)*hinv(k, j)
               call XLOAD8(hinv(k,j))
               call XMUL8(f(k))
               call XADD(sum)
               call XSTORE(sum)
               end do
            call XLOAD8(hinv(i,j))
            call XADD(sum)
            call STORND(v(j))
            end do
c
c correction of row i
c of inverse matrix
         do j = 1, n
            hinv(i, j) = v(j)
            end do
 
         end do
 
      return
      end
