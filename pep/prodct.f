 
      subroutine PRODCT(a, b, c, k, m, n)
 
c-----------------------------------------------------------------------
c  matrix multiplication
c
c  j.f.chandler  -  1977 dec 13
c  modified slightly 11/25/91 mam (sao)
c
c  a and b are the input matrices and c is the result: c = a * b
c
c  a is  k x m,  b is  m x n,  c is  k x n
c
c  If any of the dimensions are negative, then the corresponding array
c  is actually its own transpose: a(m,k) or b(n,m) or c(n,k).  The
c  absolute values of k, m, and n are used as the actual dimensions for
c  the arrays.
c
c  For example, suppose we have A(3,2), B(2,8), and C(3,8). A call
c  would look like
c
c     call PRODCT( A, B, C, 3, 2, 8 )
c
c  On the other hand, suppose A and C are dimensioned A(2,3) and C(8,3).
c  Then
c
c     call PRODCT( A, B, C, -3, 2, -8 )
c
c  would put transpose(A) * B into transpose(C).
c
c  If any dimensions are zero, the results may be incorrect, but there
c  are no other adverse effects.
c
c-----------------------------------------------------------------------
 
      implicit none
 
c*** start of declarations inserted by spag
      integer   ia, ia0, iak, iam, ib, ib0, ibm, ibn, ibx, ic, ic0, ick,
     .          icn, k, ka, kk, m, ma, mm, n
      integer   na, nn
      real*10 s
c*** end of declarations inserted by spag
 
      real*10 a(*), b(*), c(*)
 
      ka  = IABS(k)
      ma  = IABS(m)
      na  = IABS(n)
      iak = 1
      iam = ka
      if(k.le.0) then
         iak = ma
         iam = 1
      endif
      ia0 = 1 - iak
      ibn = ma
      ibm = 1
      if(m.le.0) then
         ibn = 1
         ibm = na
      endif
      ib0 = 1 - ibn
      ick = 1
      icn = ka
      if(n.le.0) then
         ick = na
         icn = 1
      endif
      ic0 = 1 - ick - icn
 
      do kk = 1, ka
         ia0 = ia0 + iak
         ic0 = ic0 + ick
         ic  = ic0
         ibx = ib0
         do nn = 1, na
            ibx = ibx + ibn
            ic  = ic + icn
            s   = 0._10
            ia  = ia0
            ib  = ibx
            do mm = 1, ma
               s  = s + a(ia)*b(ib)
               ia = ia + iam
               ib = ib + ibm
               end do
            c(ic) = s
            end do
         end do
 
      return
      end
