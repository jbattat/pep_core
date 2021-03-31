      subroutine XSQRT(a, b)
 
      implicit none
 
c
c ash/forni  february 1967   subroutine xsqrt
c extended precision calcualtion of a= sqrt(b)
c
      real*10 a(2), b(2), c, b8
 
      call XLOAD(b)
      call STORND(b8)
      if(b8.lt.0) then
 
         write(6, 50) b8
   50    format(' NEGATIVE ARGUMENT OF SQUARE ROOT =',
     .          1pd24.16/' STOP IN SUBROUTINE XSQRT')
         stop
      else if(b8.eq.0) then
         a(1) = 0._10
         a(2) = 0._10
      else
 
         c = SQRT(b8)
         call XLOAD(b)
         call XDIV8(c)
         call XADD8(c)
         call XMUL8(0.5_10)
         call XSTORE(a)
      endif
c one iteration away from double precision square root was found
c to give extended precision accuracy
      return
      end
