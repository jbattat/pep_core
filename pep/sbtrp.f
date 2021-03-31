      subroutine SBTRP(idum, jd, fract, nvel, lcntl)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 dum, fract
      integer   icall, idum, jd, lcntl, nva, nvel
 
c*** end of declarations inserted by spag
 
 
c
c        sbtrp - dec. 1975 - r. goldstein
c        velocities from numerical differentiation added  (nvel= -1)
c        by j.f.chandler, 1977 may
c        this is a control routine for interpolating a body from its
c        integration tape (or n-body tape)
c        input parameters:
c          jd        integer value of the coord. time of the desired
c                    point
c          fract     fractional part of the coordinate time of the point
c          nvel  = 0 positions determined
c                  1 velocities determined
c                 -1 velocities determined
c              .gt.1 treated as 1
c          lcntl = 0 use results of last call if possible
c                  1 do not use results (must recalculate y vectors)
c                 -1 same as 0, except that nb1 is forced to 2 (this may
c                    require some otherwise unnecessary y-vector setups)
c          lcntl must be set to 1 the first time the routine is called
c        output is put in xsb vector
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdateq.inc'
      include 'emcke.inc'
      include 'emmips.inc'
      include 'sbcom.inc'
      include 'sbdta.inc'
      include 'sbdtavtp.inc'
      include 'spqind.inc'
 
      real*10 xsba(3)
      data icall/4/
 
      if(Ksb(88).gt.0) then
c
c everett interpolation
         call EVTRP(jd, fract, nvel, lcntl, icall, Ysbpos, Xsb, Satprb,
     .              i_mxeqn, Jdsb, Fsb, Psb)
c
c
c hermite interpolation
      else if(Ksb(88).eq.-9) then
c
c quasi-elliptic orbit instead of interpolation
         nva = iabs(nvel)
         call PLIPT((jd-Sbcom(1)-.5_10) + fract, Elptg(1,4), -nva,
     .              Xsb(3*nva+1), dum)
         call TRPLST(jd,fract,nvel,lcntl,'PLIPT(SB)',Xsb)
      else if(Ksb(88).eq.-8) then
c
c elliptic orbit instead of interpolation
         nva = iabs(nvel)
         call JLIPT((jd-Sbcom(1)-.5_10)+fract,Elptg(1,4),-nva,
     .              Xsb(3*nva+1),Rympt(4),Rympt2(4),Rympt3(4),dum)
         call TRPLST(jd,fract,nvel,lcntl,'JLIPT(SB)',Xsb)
      else
         call HRTRP(jd,fract,nvel,lcntl,icall,Ysbpos,Xsb,sprb,
     .              hc2,Jdsb,Fsb,Psb)
      endif
      return
 
      entry SBTRPA(xsba)
c           enter here to get acceleration for same (jd,fract) as last
c           call. for everett method, uses whichever y vectors are
c           available and corresponding differentiating interpolation.
c           for hermite method, will have to calculate new y vector.
c
      if(Ksb(88).gt.0) then
c
c everett acceleration
         call EVTRP(0,0._10,-2,0,icall,Ysbpos,xsba,Satprb,i_mxeqn,
     .              Jdsb,Fsb,Psb)
      else
 
c hermite acceleration
         call HRTRP(0,0._10,-2,0,icall,Ysbpos,xsba,sprb,hc2,Jdsb,
     .              Fsb,Psb)
      endif
 
      return
      end
