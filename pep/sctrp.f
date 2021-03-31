      subroutine SCTRP(idum, jd, fract, nvel, lcntl, i12i)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 fract
      integer   i12, i12i, icall, idum, jd, lcntl, nvel
 
c*** end of declarations inserted by spag
 
 
c
c           j.f.chandler, 1977 may
c        this is a control routine for interpolating a body from its
c        integration tape (or n-body tape)
c        input parameters:
c          jd        integer value of the coord. time of the desired pt
c          fract     fractional part of the coordinate time of the point
c          nvel  = 0 positions determined
c                  1 velocities determined
c                 -1 velocities determined from position y-vectors
c          lcntl = 0 use results of last call if possible
c                  1 do not use results (must recalculate y vectors)
c                 -1 same as 0, except that nb1 is forced to 2 (this may
c                    require some otherwise unnecessary y-vector setups)
c          lcntl must be set to 1 the first time the routine is called
c          i12     controls which psc and xsc vectors are used
c        output is put in xsc(*,i12) vector
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'sbcom.inc'
      include 'scdta.inc'
      include 'scdtavtp.inc'
      include 'spqind.inc'
      include 'trpcom.inc'
 
      real*10 xsca(3)
 
      data icall/5/
 
      i12 = i12i
      if(Ksc(88).le.0) then
c
c
c hermite interpolation
         call HRTRP(jd, fract, nvel, lcntl, icall, Yscpos, Xsc(1,i12),
     .              sprc, Gc2, Jdsc, Fsc, Psc(1,i12))
         if(i12.eq.1) Nb1(icall) = Nbtrp(icall)
      else
c
c everett interpolation
         call EVTRP(jd, fract, nvel, lcntl, icall, Yscpos, Xsc(1,i12),
     .              Satprc, i_mxeqn, Jdsc, Fsc, Psc(1,i12))
      endif
      return
 
      entry SCTRPA(xsca)
c           enter here to get acceleration for same (jd,fract) as last
c           call. for everett method, uses whichever y vectors are
c           available and corresponding differentiating interpolation.
c           for hermite method, will have to calculate new y vector.
c           results go into xsca.
c
      if(Ksc(88).gt.0) then
c
c everett acceleration
         call EVTRP(0,0._10,-2,0,icall,Yscpos,xsca,Satprc,i_mxeqn,
     .              Jdsc,Fsc,Psc(1,i12))
      else
 
c hermite acceleration
         call HRTRP(0,0._10,-2,0,icall,Yscpos,xsca,sprc,Gc2,
     .              Jdsc,Fsc,Psc(1,i12))
      endif
 
      return
      end
