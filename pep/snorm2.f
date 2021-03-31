      subroutine SNORM2(f,pvrow,pvcol,pvrwb,b,nparam,itmax,eps,iercod)
 
      implicit none

c m.ash   jan 1972    subroutine snorm2
c iterative cleanup of covariance matrix

c arguments
      real*10 f(1000),pvrow(1000),pvcol(1000),pvrwb(1000)
      real*10 b(1)
      real*4 eps(2)
      integer*4 iercod,nparam
      integer*2 itmax
 
c common
      include 'inodta.inc'

c local
      real*10 finorm, fjnorm, fnormc, fnormr
      integer*4 iter

      call PAGCHK(47,2,0)
      write(Iout,100) nparam
  100 format('0CLEANUP OF INVERSE OF ORDER', i4)
c
c iterative cleanup of covariance matrix
      iercod = 0
      fnormr = 1.0E9_10
      fnormc = 1.0E9_10
      iter   = 0
      do while( .true. )
 
         iter = iter + 1
         if(iter.ne.1) then
            fnormr = finorm
            fnormc = fjnorm
         endif
 
         call MULITR(b,f,pvrow,pvcol,pvrwb,nparam,finorm,fjnorm)
 
         call PAGCHK(60,1,0)
         write(Iout,150) iter,finorm,fjnorm
  150    format(9x,' ITERATION', i4, '  ROW NORM =', 1pd11.4,
     .          '  COLUMN NORM =', 1pd11.4)
c
c check if row or column norm is decreasing
         if(finorm.ge.fnormr .and. fjnorm.ge.fnormc) then
            call PAGCHK(60,1,0)
            write(Iout,160) iter
  160       format(1x,' NORM NON-DECREASING AFTER', i4,
     .             '  ITERATIONS')
c
c check if row or column norm .lt.1
            if(finorm.ge.1.0_10 .and. fjnorm.ge.1.0_10) then
c
c bad error, row and column norm.ge.1
               call PAGCHK(60,2,0)
               write(Iout,170)
  170          format('0BAD ERROR, ROW AND COLUMN NORM.GE.1')
               iercod = 3
c
c check if row or column norm.lt.eps(1)
            else if(finorm.ge.eps(1) .and. fjnorm.ge.eps(1)) then
c
c mild error, row and column norm.ge.eps(1)
               iercod = 2
               call PAGCHK(60,2,0)
               write(Iout,180) eps(1)
  180          format('0MILD ERROR, ROW AND COLUMN NORM.GE.EPS(1)=',
     .                1pe11.3)
            endif
c
c
c check if row or column norm.lt.eps(2)
         else if(finorm.ge.eps(2) .and. fjnorm.ge.eps(2)) then
c
c check if number of iterations.ge.itmax   (from input)
            if(iter.lt.itmax) goto 300
c
c mildest error,number of iterations.ge.itmax
            iercod = 1
            call PAGCHK(60,2,0)
            write(Iout,200) itmax
c
c
c cleanup completed  or  discontinued
  200       format('0MILDEST ERROR, NUMBER OF ITERATIONS.GE.ITMAX=',
     .             i10)
         endif
         call TIMRIT('  ITERATIVE CLEANUP ', 5)
 
         return
  300 end do
      end
