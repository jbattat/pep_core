      subroutine HRTRP(jd,fract,nvel,lcntl,icall,yh,x,body,hc,
     .                 jdtb, frtb, p)
 
      implicit none
 
c        subr. hrtrp - j.f.chandler - 1979 nov
c        derived from subr. sctrp - j.f.chandler, 1977 may
c        this is a control routine for interpolating a body from its
c        integration (or n-body or s-body) tape - hermite method
c        input parameters:
c          jd      integer value of the coord. time of the desired pt
c          fract   fractional part of the coordinate time of the point
c          nvel =0 positions determined
c                1 velocities determined
c               -1 velocities determined from position y-vectors
c          lcntl=0 use results of last call if possible
c                1 do not use results (must recalculate y vectors)
c               -1 same as 0
c          icall   indicates which is the calling routine
c                   1 - emtrp,  2 - mntrp,  3 - pltrp,  4 - sbtrp,
c                   5 - sctrp,  6 -(prtrp), 7 -(ertrp)
c          yh      array for interpolation y-vectors
c          x       output coordinate array
c          body    coordinates from integration tape
c          hc      tabular intervals
c          jdtb    tabular julian dates
c          frtb    tabular fractions of a day
c          p       p-q array to be used
c
c        note for future work:
c        the trick of using the results of the last call is good only
c        within the light time iteration since velocity
c        calculations overwrite the y vectors.
c        however, velocities could be interpolated from the same
c        y-vectors, by numerical differentiation
c        (would require new function "dhermt")
c        and settling for (in effect) a lower-order polynomial
c
c        commons
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'trpcom.inc'
 
      real*10 fract, HERMTF, p
      integer*4 icall,j,jd,jh,jo,jy,lcntl,nbt,ndmy,nvel,nvela
c
c names of interpolators
      character*2 name(10)/'EM','MN','PL','SB','SC','PR','ER',3*'  '/
 
c data from integration tape
      real*10 body(6,8,3),frtb(6),hc(6)
      integer*4 jdtb(6)
 
c output
      real*10 yh(8,3),x(6)
c
c
      nvela = iabs(nvel)
      jo    = 3*nvela
      if(nvel.eq.-2) jo = 0
c
c test if routine should attempt to use saved values
      if(lcntl.ne.1) then
c
c
c see if within saved tabular region
         nbt = Nbtrp(icall)
         if(nvel.eq.-2 .and. jd.eq.0) goto 100
         if(nbt.gt.0) then
            p = (jd - jdtb(nbt)) + fract - frtb(nbt)
            if(Idirb(icall).lt.0) then
 
c backward in time
               if(p.le.0._10 .and. p.ge.hc(nbt)) goto 100
 
c forward in time
            else if(p.ge.0._10 .and. p.le.hc(nbt)) then
               goto 100
            endif
         endif
      endif
c
c*  start=1100
c read body tape
      if(icall.eq.1) then
         call EMREED(jd,fract)
      else if(icall.eq.2) then
         call MNREED(jd)
      else if(icall.eq.3) then
         call PLREED(jd)
      else if(icall.eq.4) then
         call SBREED(jd,fract)
      else if(icall.eq.5) then
         call SCREED(jd,fract)
      else if(icall.eq.6) then
         call RTREED(jd,fract,ndmy,10)
      else if(icall.eq.7) then
         call RTREED(jd,fract,ndmy,3)
      else
         call SUICID('BAD ICALL, STOP IN HRTRP', 6)
      endif
 
      if(jd.le.0) return
      Nvels(icall) = 9999
c
c*  start=1500
c find left hand tabular point
      p = (jd - jdtb(3)) + fract - frtb(3)
      if(p*hc(3).ge.0._10) then
         nbt = 3
         if(ABS(p).ge.ABS(hc(3))) then
            p   = p - hc(3)
            nbt = 4
         endif
      else
         p   = p + hc(2)
         nbt = 2
      endif
 
c save tabular region boundaries
      Nbtrp(icall) = nbt
c
c*  start=2500
c set up y vector...warning...mtype=0...no good for partials
  100 jh = nvela + 1
      if(nvela.ne.Nvels(icall)) then
         do j = 1, 3
            call YHERMT(body(1,1,j),yh(1,j),jh,0,hc,nbt)
         end do
         Nvels(icall) = nvela
      endif
c
c*  start=3000
c do interpolation
      do jy = 1, 3
         j    = jy + jo
         x(j) = HERMTF(yh(1,jy),jh,0,p)
      end do
c
c debug printout
      call TRPLST(jd,fract,nvel,lcntl,name(icall)//'TRP(HR)',x)
c
      return
      end
