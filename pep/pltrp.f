      subroutine PLTRP(idumm,jdt,frt,nvel,lcntl)
 
      implicit none

c        PLTRP - Dec. 1975 - R. Goldstein
c        Velocities from numerical differentiation added  (nvel= -1)
c        by J.F.Chandler, 1977 May
c        This is a control routine for interpolating a planet from its
c        integration tape (or n-body tape)

c input parameters:
      real*10 frt
      integer*4 idumm,jdt,nvel,lcntl
c          IDUMM     ignored
c          JDT       integer value of the coord. time of the desired pt
c          FRT       fractional part of the coordinate time of the point
c          NVEL  = 0 positions determined
c                  1 velocities determined
c                 -1 velocities determined from position y-vectors
c          LCNTL = 0 use results of last call if possible
c                  1 do not use results (must recalculate y vectors)
c                 -1 use old results only if jd.fract is in same tab pt.
c        output is put in the xp(1-6) vector
c
c        at present, only everett interpolation is allowed
c        hermite method would be indicated by kp(88)=0

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'coord.inc'
      include 'emcke.inc'
      include 'emmips.inc'
      include 'fcntrl.inc'
      include 'pemctl.inc'
      include 'pqind.inc'
      include 'tapdta.inc'
      include 'tapdtp.inc'
      include 'yvect.inc'

c local
      real*10 dum,frts
      integer*4 icall/3/,jdts,nva
c
c argument array for acceleration
      real*10 xpla(3)
 
      jdts = jdt
      frts = frt
      if(Kp(88).eq.-8) then
c
c elliptic motion instead of integrated
         nva = iabs(nvel)
         call JLIPT(jdt-Pcom(1)-0.5_10+frt, Elptg(1,3), -nva,
     .              Xp(3*nva+1),Rympt(3),Rympt2(3),Rympt3(3),dum)
         call TRPLST(jdt,frt,nvel,lcntl,'PLIPT(PL)',Xp)
      else
c
c everett interpolation
         call EVTRP(jdt, frt, nvel, lcntl, icall, Yplpos, Xp, Planet,
     .              i_mxplprt+1, Jdp, Fp, Ppl)
      endif
 
      if(Jct(55).ne.0) call PLCMC(jdts, frts, Xp, nvel)
      return
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry PLTRPA(xpla)
c           enter here to get acceleration for same (jdt,frt) as last
c           call. for everett method, uses whichever y vectors are
c           available and corresponding differentiating interpolation.
c
c           everett acceleration
      call EVTRP(0,0._10,-2,0,icall,Yplpos,xpla,Planet,i_mxplprt+1,Jdp,
     .           Fp,Ppl)
 
      if(Jct(55).ne.0) call PLCMC(jdts,frts,xpla,2)
      return
      end
