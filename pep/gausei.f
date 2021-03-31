      subroutine GAUSEI(a, b, x, n, k2, epst, iersol)
 
      implicit none

c
c        f amuchastegui - november 1968 - subroutine gausei
c
c        subroutine gausei finds the improved solutions of a
c        system of linear equations, using gauss-seidel method
c        it is suposed a priori that we know the first estimation
c        of the solutions - the calling sequence is the following
c
c     call gausei(a,b,x,n,k2,epst,iersol)
c
c        where:
c        a = lower diagonal half of the symmetric coefficient
c            matrix a
c        b = right hand side fo the equations
c        x = first estimations of the solutions and also,
c            final solution
c        n = rank of the matrix a
c        k2 = maximum number of iterations
c        epst= error estimation, program will be stopped
c              when  (x(new)-x(old)) < (epst*x(old))
c        iersol is a warning message so that
c        if the returned value of iersol is different
c        from zero, it means that something was wrong in the
c        clean-up of the solution
c arguments
      real*10 a(1),b(1),x(1)
      real*4 epst
      integer*4 n,iersol
      integer*2 k2

c common
      include 'inodta.inc'

c local
      real*10 psave,sum(2),e,e2,q,rel,qepst
      integer   i,k,l,loop,nk

      e2   = 1.E19_10
      loop = 0
      write(Iout,100) n
  100 format('0GAUSS-SEIDEL CLEAN-UP OF SOLUTION OF ORDER ',i5)
      do while(.true.)
         loop = loop + 1
         if(loop.ne.1) e2 = e
         e   = 0._10
         q   = 0._10
         rel = 0._10
         do i = 1,n
            nk     = (i*(i+1))/2
            sum(1) = 0._10
            sum(2) = 0._10
            do k = 1,n
               if(k.ge.i) then
                  l = i + (k*(k-1))/2
               else
                  l = k + (i*(i-1))/2
               endif
               call XLOAD8(a(l))
               call XDIV8(a(nk))
               call XMUL8(x(k))
               call XADD(sum)
               call XSTORE(sum)
               end do
            call XLOAD8(b(i))
            call XDIV8(a(nk))
            call XSUB(sum)
            call STORND(psave)
            call XADD8(x(i))
            call STORND(x(i))
            q = q + ABS(x(i))
            if(x(i).ne.0._10) rel = rel + ABS(psave/x(i))
            e = e + ABS(psave)
            end do
         write(Iout,150) loop,e,q
  150    format(10x,'ITERATION',i4,'  ERROR',1pd12.2,
     .          '  RELATIVE TO',1pd12.2)
 
c check if error is decreasing
         if(e.lt.e2) then
 
c check if error is.lt.q*epst
            qepst = q*epst
            if(e.ge.qepst) then
 
c check if number of iterations .ge.k2
               if(loop.lt.k2) go to 200
               iersol = 1
               write(Iout,160) loop
  160          format('0MILDEST ERROR, NUMBER OF ITERATIONS .GE. K2 =',
     .                i5)
            endif
         else
            write(Iout,180) loop
  180       format('0 ERROR NON-DECREASING AFTER',i5,'  ITERATIONS')
 
c check if the relative error is .lt.1.0_10
            if(rel.ge.1._10) then
               write(Iout,190)
  190          format('0BAD ERROR, RELATIVE ERROR GREATER THAN 1')
               iersol = 3
            endif
         endif
         return
  200    end do
      end
