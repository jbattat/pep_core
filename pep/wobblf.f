      subroutine WOBBLF(jd, fract, x, y)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, jd
      real      x, y
 
c*** end of declarations inserted by spag
 
 
c
c r.w.king    mar 1973    subroutine wobblf
c determine corrections for wobble using analytic models
c
      real*10 fract
      include 'fcntrl.inc'
      include 'dtparm.inc'
      include 'leon.inc'
c
c local
      real*10 a1, a2, h, b1, b2
 
      x  = 0.
      y  = 0.
      I3 = 0
      I4 = 0
      I5 = 0
      I6 = 0
      if( Numdt .ne. 0 ) then
         if( Jddt0 .le. 0 ) then
c
c determine wobble from piecewise linear function
            if( jd .ge. Jddt(1) ) then
               if( jd .lt. Jddt(Numdt) ) then
                  do i = 2, Numdt
                     if( jd .lt. Jddt(i) ) go to 10
                  end do
                  call SUICID(' JD NOT IN DT TABLE, STOP IN WOBBLF 10  '
     .                        , 10)
   10             I1 = i - 1
                  I2 = i
                  I3 = I1 + 200
                  I4 = I2 + 200
                  I5 = I3 + 200
                  I6 = I4 + 200
                  a1 = Dt(I3)
                  a2 = Dt(I4)
                  b1 = Dt(I5)
                  b2 = Dt(I6)
                  h  = Jddt(I2) - Jddt(I1)
                  Ff = ((jd-Jddt(I1)) + fract)/h
                  if( Jct(12) .lt. 0 ) Ff = 0._10
                  x = (a2 - a1)*Ff + a1
                  y = (b2 - b1)*Ff + b1
               endif
            endif
         endif
      endif
c
c
      return
      end
