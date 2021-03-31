      subroutine FORMU(b)
 
      implicit none
c
c d. white  october 1973  subroutine formu
c
c
c parameter is b matrix into which is read the a priori b
      real*10 b(1)
c
c common
      include 'aprtbf.inc'
      include 'estvec.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
 
c misc. temporary storage
      common/WRKCOM/ Xap(1000),U(1000),Xmx(1000),Sumsq,T,Pointr(1000)
      real*10 U,Xap,Xmx,Sumsq
      integer*2 Pointr
      character*80 T
c local
      integer   i0,i,ij,inc,j,k,mparam
      integer*4 minus1/-1/
c
c read and print title
      read(Ibuf2) T
      write(Iout,100) Iterat,Heding,Date,Npage,T
  100 format('1FORM A PRIORI RHS - ITERAT =',i3,9x,18A4,1x,2A4,
     .       ' PAGE',i5/' TITLE = ',A80/)
      Npage = Npage + 1
c
c skip sumsq residual contribution
      read(Ibuf2)
c
c get mparam and pointer
      read(Ibuf2) mparam,(Pointr(i),i = 1,mparam)
c
c skip right hand side
      read(Ibuf2)
c
c read in a priori b
c fill in lower diagonal half
      i0 = 0
      do i = 1,mparam
         read(Ibuf2) k,(b(i0+j),j=1,i)
         i0 = i0 + i
      end do
c
c read in a priori x
      read(Ibuf2) k,(Xap(i),i=1,k)
      if(k.ne.mparam) call SUICID(
     . 'NO A PRIORI ESTIMATE VECTOR, STOP IN FORMU  ',11)
c
c form difference vector  xap - xbar
      do i = 1,mparam
         Xmx(i) = Xap(i)-Xbar(Pointr(i))
      end do
c
c multiply  u = b * (xap - xbar)
      i0 = 0
      Sumsq = 0._10
      do i = 1,mparam
         U(i)= 0._10
         ij  = i0
         i0  = i0 + i
         inc = 1
         do j = 1,mparam
            ij = ij + inc
            if(j.ge.i) inc = j
            U(i) = U(i) + b(ij)*Xmx(j)
         end do
         Sumsq = Sumsq + Xmx(i)*U(i)
      end do
c
c get ready to write out new rhs contribution, skip title & ptr
      rewind Ibuf2
      write(Ibuf2) T
      write(Ibuf2) minus1,Sumsq
      write(Ibuf2) mparam,(Pointr(i),i=1,mparam)
      write(Ibuf2) mparam,(U(i),i=1,mparam)
      write(Iout,200) Sumsq,(U(i),i=1,mparam)
  200 format(' NEW SUMSQ=',1PD12.5,' NEW RHS CONTRIBUTION'/(7x,10D12.5))
      Line = 5 + (mparam - 1)/10
c
c write b back out
      call FWSIG(Ibuf2,-mparam,b,Xmx)
c
c write xap back
      write(Ibuf2) mparam,(Xap(i),i = 1,mparam)
      rewind Ibuf2
      Itrwnd(Ibuf2) = 0
      return
      end
