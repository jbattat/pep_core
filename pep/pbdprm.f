      subroutine PBDPRM(nkp,kp,lplnt,kplnt,ltest,iflag)
      implicit none
 
c originally Ash/Reasonberg 1970 Sep - just for body parameters.
c J.F.Chandler, 1977 Jul
c modular subroutine extracted from routine partl to control logic
c for finding partials on integration tapes w.r.t. solar-system
c parameters or body parameters
c KP is the array of integration control integers (of length NKP),
c specifying which partial derivatives are included.  The first 7 elements
c of KP refer to the motion and initial conditions and are not used here.
c LPLNT is the index of the next element of KP to consider
c KPLNT is the index into the array of integrated partials corresponding
c to LPLNT
c input value of LTEST is abs. value of desired partial
c input IFLAG=1 - desired KP neg. for body parameter(con) or
c                 moon rotation parameter(cond or con)
c input IFLAG=-1 - KP positive for central body parameter(con)
c                  or solar system parameter
c routine searches array KP for desired value (starting at LPLNT)
c and returns updated LPLNT and KPLNT as well as indicator:
c   IFLAG=0  -  value found ok, may call CPARTL
c   IFLAG=1  -  value not found

c parameters
      integer*2 nkp, kp(99),ltest
      integer lplnt,kplnt,iflag

c local variables
      integer i, k, k0

      if(iflag.eq.0) goto 100
      do while( lplnt.le.nkp )
         k = kp(lplnt)
         if(k.eq.0) then

c end of list
            lplnt = nkp + 1
            goto 100
         else if(k*iflag.lt.0) then

c same sign as desired
            i = k*iflag + ltest
            if(i.lt.0) goto 100
         else if(iflag.gt.0) then

c looking for negative, but none left
            lplnt = nkp + 1
            goto 100
         else

c looking for positive kp, skip over neg.
            i = 1
         endif
         lplnt = lplnt + 1
         kplnt = kplnt + 1
         if(i.eq.0) then
c
c returns
            iflag = 0
            return
 
c value not reached, check for harmonic coeff. partials
         else if(k.gt.100) then
            k0 = mod(k,100)
 
c partial is for some planet parameter, see if harmonic
            if(k0.ge.31) then
               if(k0.ne.31) then
 
c tesseral harmonics
                  if(mod(k0,10).gt.4) then
 
c resonant tesserals
                     lplnt = lplnt + 1
                  else
                     kplnt = kplnt + (kp(lplnt+2)*(kp(lplnt+2)-1))/2
     .                 + kp(lplnt+3)
     .                - ((kp(lplnt)*(kp(lplnt)-1))/2 + kp(lplnt+1))
                     lplnt = lplnt + 4
                     goto 50
                  endif
               endif
c note: harmonic coefficient partials violate the usual rule
c of 1 kp integer = 1 partial on tape  --  special handling
c zonal harmonics
               kplnt = kplnt + kp(lplnt+1) - kp(lplnt)
               lplnt = lplnt + 2
            endif
         endif
   50 end do
  100 iflag = 1
      return
      end
