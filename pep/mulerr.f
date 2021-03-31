      subroutine MULERR(b, f, roa, buf, n)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 col, row, row0
      integer   i, j, j0, jfk, jk, k, k1, l, linexp, n
 
c*** end of declarations inserted by spag
 
 
c
c f amuchastegui - august 1969 - subroutine mulerr
c evaluate and write out error matrix for calculated
c inverse of coefficient matrix of normal equations
c
      real*10 b(1), f(1), roa(1), buf(1)
      real*10 sum(2)
 
      include 'fcntrl.inc'
      include 'inodta.inc'
 
      character*8 words(2)/' IS     ',' (CONT) '/
 
      row = 0._10
 
      jfk    = 1
      Line   = 64
      linexp = 2 + (n - 1)/16
      do l = 1, n
         buf(l) = 0._10
      end do
c
c evaluation of error matrix
      do i = 1, n
         read(Ibuf) j, (roa(l), l = 1, n)
         j0 = 0
         do j = 1, n
            sum(1) = 0._10
            sum(2) = 0._10
            k1     = 1
            jk     = j0
            j0     = j0 + j
            do k = 1, n
               jk = jk + k1
               if(k.ge.j) k1 = k
 
c sum = sum + roa(k)*b(jk)
               call XLOAD8(roa(k))
               call XMUL8(b(jk))
               call XADD(sum)
               call XSTORE(sum)
            end do
            if(j.eq.i) call XSUB8(1._10)
            call STORND(f(j))
            f(j) = -f(j)
         end do
c
c evaluation of row and column norm
         row0 = 0._10
         do l = 1, n
            buf(l) = buf(l) + ABS(f(l))
            row0   = row0 + ABS(f(l))
         end do
         row = MAX(row, row0)
         if(i.ge.n) then
            col = 0._10
            do l = 1, n
               col = MAX(col, buf(l))
            end do
         endif
c
c write out row of error matrix
         if(Line + linexp.gt.60) then
            call NEWPG
            write(Iout, 20) n, words(jfk)
   20       format('0PRODUCT OF THE COEFFICIENT MATRIX OF ORDER',
     .             i4,
     .          ' AND ITS CALCULATED INVERSE MINUS THE IDENTITY MATRIX'
     .          , a8)
            jfk  = 2
            Line = 3
         endif
         write(Iout, 50) i, (f(k), k = 1, n)
   50    format('0', i3, (t5,1p,16E8.1))
         Line = Line + linexp
 
      end do
      rewind Ibuf
c
c printout row and column norm
      if(Line.gt.52) then
         call NEWPG
         write(Iout, 20) n, words(2)
         Line = Line + 2
      endif
      write(Iout, 100) row, col
  100 format('-ROW NORM=', 1pd11.4, '  COL NORM=', 1pd11.4)
      Line = Line + 3
      call TIMRIT(' CALCULATING & WRITING OUT ERROR MATRIX ', 10)
 
      return
      end
