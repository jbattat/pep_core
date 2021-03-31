      subroutine MULER2(h,hinv,f,buf,n)
 
      implicit none
 
 
c*** start of declarations inserted by spag
 
c*** end of declarations inserted by spag
 
 
c
c d. white  subroutine muler2  september 1974
c
c         modified version of mulerr to print error matrix for inverse,
c         where both matrices are full n*n and in core
c
c
c         parameters
      integer*4 n
      real*10 h(n,n),hinv(n,n),f(n),buf(n)
c
c common
      include 'fcntrl.inc'
      integer*2 itmax
      equivalence (Ict(49),itmax)
      include 'inodta.inc'
c
c local
      real*10 sum(2),col,row,row0
      integer   i,j,jfk,k,l,lintst
      character*8 words(2)/' IS     ',' (CONT) '/
 
      row = 0._10
 
      jfk    = 1
      lintst = 59 - (n - 1)/16
      do l = 1,n
         buf(l) = 0._10
      end do
c
c evaluation of error matrix
      do i = 1,n
         do j = 1,n
            sum(1) = 0._10
            sum(2) = 0._10
            do k = 1,n
 
c sum = sum + hinv(i,k)*h(k,j)
               call XLOAD8(h(k,j))
               call XMUL8(hinv(i,k))
               call XADD(sum)
               call XSTORE(sum)
            end do
            call STORND(f(j))
            if(j.eq.i) then
               call XLOAD(sum)
               call XSUB8(1._10)
               call XSTORE(sum)
               call STORND(f(j))
            endif
            f(j) = -f(j)
         end do
c
c evaluation of row and column norm
         row0 = 0._10
         do l = 1,n
            buf(l) = buf(l) + ABS(f(l))
            row0   = row0 + ABS(f(l))
         end do
         row = MAX(row,row0)
         if(i.ge.n) then
            col = 0._10
            do l = 1,n
               col = MAX(col,buf(l))
            end do
         endif
c
c write out row of error matrix
         if(Line.ge.lintst) then
 
c lintst=lintst-1
            write(Iout,150) Iterat,Heding,Date,Npage
            Npage = Npage + 1
            write(Iout,20) n,words(jfk)
   20       format('0PRODUCT OF THE H MATRIX OF ORDER',i4,
     .          ' AND ITS CALCULATED INVERSE MINUS THE IDENTITY MATRIX'
     .          ,a8)
            jfk  = 2
            Line = 3
   40       format(/i4,1p,16E8.1)
         endif
         if(n.eq.16) then
            write(Iout,40) i,(f(k),k = 1,n)
         else
            write(Iout,60) i,(f(k),k = 1,n)
   60       format(/i4,1p,16E8.1,/(4x,1p,16E8.1))
         endif
         Line = Line + 2 + (n - 1)/16
 
      end do
c
c printout row and column norm
      if(Line.lt.lintst) then
         write(Iout,100) row,col
  100    format(//1x,'ROW NORM=',1pd11.4,2x,'COL NORM=',1pd11.4)
  150    format('1FORM FILTER MATRICES',5x,'(ITERAT=',i2,')',5x,
     .          18A4,1x,2A4,' PAGE',i5)
      else
         write(Iout,150) Iterat,Heding,Date,Npage
         Npage = Npage + 1
         write(Iout,200) n,words(2),row,col
  200    format('0PRODUCT OF THE H MATRIX OF ORDER',i4,
     .        ' AND ITS CALCULATED INVERSE MINUS THE IDENTITY MATRIX',
     .        a8 / '0ROW NORM=',1pd11.4,'  COL NORM=',d11.4)
      endif
 
      return
      end
