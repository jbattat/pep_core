      subroutine ADDCLN(h,hinv,f,v,buf,n)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 finorm, fjnorm, fnormc, fnormr
      integer   iter
 
c*** end of declarations inserted by spag
 
 
c
c d. white  subroutine addcln  november 1974
c
c         set up for cleanup of inversion
c
c         parameters
      integer*4 n
      real*10 h(n,n),hinv(n,n),f(n),v(n),buf(n)
c
c common
      include 'fcntrl.inc'
      include 'filtim.inc'
      include 'inodta.inc'
c
c local
      integer*2 itmax
      equivalence (Fict(3),itmax)
c
c formats
  100 format('0CLEANUP OF INVERSE OF ORDER', i4)
  200 format(9x,' ITERATION', i4, '  ROW NORM =', 1pd11.4,
     .       '  COLUMN NORM =', 1pd11.4)
  300 format(1x,' NORM NON-DECREASING AFTER', i4, '  ITERATIONS')
  400 format('0BAD ERROR, ROW AND COLUMN NORM.GE.1.0_10')
  500 format('0MILDEST ERROR, NUMBER OF ITERATIONS.GE.ITMAX=', i10)
c
c whats doing?
      if(itmax.gt.0) then
c
c iterative cleanup
         write(Iout,100) n
         Line   = Line + 2
         fnormr = 1.0E9_10
         fnormc = 1.0E9_10
         iter   = 0
         do while( .true. )
 
            iter = iter + 1
            if(iter.ne.1) then
               fnormr = finorm
               fnormc = fjnorm
            endif
 
            call MULIT2(h,hinv,f,v,buf,n,finorm,fjnorm)
 
            write(Iout,200) iter,finorm,fjnorm
            Line = Line + 1
c
c check if row or column norm is decreasing
            if(finorm.ge.fnormr .and. fjnorm.ge.fnormc) then
               write(Iout,300) iter
               Line = Line + 1
c
c check if row or column norm .lt.1
               if(finorm.ge.1.0_10 .and. fjnorm.ge.1.0_10) then
c
c bad error, row and column norm.ge.1
                  write(Iout,400)
                  Line = Line + 2
                  goto 600
               endif
            endif
c
c check if row or column norm.lt.feps(1)
            if((finorm.ge.Feps(1)) .and. (fjnorm.ge.Feps(1)))
     .         then
c
c check if number of iterations.ge.itmax   (from input)
               if(iter.lt.itmax) goto 550
c
c mildest error,number of iterations.ge.itmax
               write(Iout,500) itmax
               Line = Line + 2
            endif
            goto 600
  550    end do
      endif
c
c print out error matrix
  600 if(itmax.eq.-1) return
      call MULER2(h,hinv,f,buf,n)
      return
      end
