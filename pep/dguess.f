      subroutine DGUESS(delay,leg,jdarg,utarg)
 
      implicit none

c          r. goldstein     april, 1976
c          this routine returns the first guess for the time delay for
c          a leg of the time delay iter. scheme. it uses a polynomial
c          extrapolation of the delays of up to three previous points
c          to accomplish this.  entries are made to tsav to save the
c          previous points and to dginit to initialize at the beginning
c           of a series.
c
c  inputs:
      real*10 delay,utarg
      integer*4 leg,jdarg
c       delay - the returned value of the extrapolated time delay
c               in seconds
c       leg - the leg of the iteration
c                1 receive to reflect time
c                2 s/c a to s/c b time
c                3 reflect to send time
c           for only 1 spacecraft, use leg=1 and leg=3.
c       jdarg,utarg - time at which to calculate delay
c                (jdarg-days,  utarg-seconds)

c array dimensions
      include 'globdefs.inc'
c
c        commons
      include 'comdat.inc'
      real*10 secday
      equivalence (Comcon(73),secday)
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obscrd.inc'
      real*10 tmdly0
      equivalence (Result,tmdly0)
      include 'spqind.inc'
      include 'yvect.inc'

c local
      real*10 cor,time,ute
      integer   i,j,jctt,jde
      real*10 dxsav(3,3),ssav(2,3),ypsav(2,3),ysav(3)
      integer*4 ns(3),jdsav(3)
      real*10 epsmax/5.E3_10/
      real*10 utsav(3),deltim
c
c
c           ns, the number of stored points is kept track of in
c          the tsav section
c
c           first check if there are any saved points or if the time
c          span is greater than epsmax or if a pseudo-observable,
c          or other unusual conditions.
c
      if(Ncodf.eq.19) then
c
c
         delay = 0._10
         return
      else
         if(ns(leg).gt.0) then
 
            deltim = (jdarg - jdsav(leg))*secday + (utarg - utsav(leg))
            if(ABS(deltim).le.epsmax) then
 
               dxsav(3,leg) = deltim
 
c do extrapolation
               delay = ysav(leg)
               if(ns(leg).gt.1) then
 
c 2nd approx. is linear
                  ssav(2,leg) = dxsav(2,leg) + deltim
 
c be sure times are distinct
                  if(dxsav(2,leg).ne.0._10) then
                     cor   = ypsav(2,leg)*deltim
                     delay = delay + cor
                     if(ns(leg).gt.2) then
 
c 3rd approx. is quadratic
                        if(dxsav(1,leg).ne.0._10 .and. ssav(1,leg)
     .                     .ne.0._10) delay = delay +
     .                      (ssav(2,leg)/ssav(1,leg))
     .                      *(cor - ypsav(1,leg)*deltim)
                     endif
                  endif
               endif
               return
            endif
         endif
c
c
c          the following sections returns a guess for the delay
c          when there are no saved points or if a pseudo-obs,
c          or other unusual conditions.
c
         ns(leg) = 0
         if(leg.eq.2) then
 
            call SUICID('LEG=2 IS NOT IMPLEMENTED IN DGUESS  ', 9)
         else if(leg.ne.3) then
 
            delay = Tguess
c use observed delay/2 unless it is to be demodded
c also insist on positive range and actual range observation
            if(Ncodf.le.2) then
               if(tmdly0.gt.0._10 .and. Nice.le.0 .and. Ict(25).le.0)
     .          delay = tmdly0/2._10
            endif
            if(Dstf(4).gt.0._10) delay = Dstf(4)
            return
         endif
 
         delay = Tmdly2
         return
      endif
c
c
c save the calculated time delays
c
      entry TSAV(time,leg,jde,ute)
      if(ns(leg).lt.3) ns(leg) = ns(leg) + 1
      jdsav(leg)    = jde
      utsav(leg)    = ute
      ssav(1,leg)  = ssav(2,leg)
      dxsav(1,leg) = dxsav(2,leg)
      dxsav(2,leg) = dxsav(3,leg)
      ypsav(1,leg) = ypsav(2,leg)
      if(ns(leg).gt.1 .and. dxsav(2,leg).ne.0._10) ypsav(2,leg)
     .    = (time - ysav(leg))/dxsav(2,leg)
      ysav(leg) = time
      jctt = Jct(6)
      if(mod(jctt,2).ne.0) then
         if(Line.ge.57) call OBSPAG
         write(Iout,50) leg,jde,ute,time,ns(leg)
   50    format(' TSAV: LEG', i2, ', JD.UT=', i8, f10.2, ', DELAY=',
     .          1pd22.15, ', NS=', i2)
         Line = Line + 1
      endif
      return
c
c
c
      entry DGINIT
      do i = 1, 3
         ns(i)    = 0
         jdsav(i) = 0
         utsav(i) = 0._10
         ysav(i)  = 0._10
         do j = 1, 3
            dxsav(i,j) = 0._10
         end do
      end do
 
      return
      end
