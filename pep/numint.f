      subroutine NUMINT(ngo)
 
      implicit none
 
c
c     a.d.rasinski   december 1964   subroutine numint
c           main program for nordsieck integration of n first order
c           ordinary differential equations with modifications
c           such that only the first m equations will control the
c           interval size.  branch for corout,fn added oct 67, nov 69.
c           timer added aug 1967
c           forward-backward added december 1967
c
c
c argument
      integer*4 ngo

c     ngo=1 planet-asteroid integration (obsolete)
c     ngo=2 moon integration (obsolete)
c     ngo=3 satellite-probe integration (numint called by planet)
c     ngo=4 n-body integration (numint called by planet)
c     ngo=5 earth rotation integration
c     ngo=6 moon rotation integration (obsolete)
c     ngo=7 planet rotation integration
c     ngo=8 moon orbit and rotation integration
c
c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'adams.inc'
      include 'fcntrl.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'nmstrt.inc'
      include 'output.inc'
      include 'stint.inc'
 
c local variables
      integer*4 i, icall, ii, istrtt, kstop, l2
      logical*4 ovrlap
c
c timer quantities
      integer*4 ihour(2),iminut(2),itot(2)
      real*4 second(2)
      character*4 loon(2)/'REAL','TASK'/
      character*4 gotdat(2)/'GOT ','DATE'/
      real*10 steps, duratn, stpday, ww(2), www(2)
c
c initialize
      ovrlap = .false.
      icall  = 0
      istrtt = 0
  100 if(Tstop.lt.T0) then
         assign 300 to kstop
      else if(Tstop.eq.T0) then
         goto 900
      else
         assign 400 to kstop
      endif
  200 call NINSTP(l2, ngo)
      if(istrtt.eq.0) then
c extra print of initial accelerations
         if(Kout.gt.0 .and. MOD(Intxpr,2).eq.1)
     .    call KOUTXPR(T0,Dydt0,Neq,6,1)
         call TIMRIT(' STARTING PROCEDURE ', 5)
      endif
      istrtt = 1

c extra print of accelerations
      if(Kout.gt.0 .and. Intxpr.gt.0) call KOUTXPR(T,Dydt,Neq,6,Intxpr)
 
      if(l2.gt.0) goto 900
      if(ngo.eq.3) then
         if(Jct(56).le.0) call SBOUT
         if(Jct(56).gt.0) call SBFOUT(ovrlap, icall)
      else if(ngo.eq.4) then
         call BODOUT
      else if(ngo.eq.5) then
         call ERTOUT
      else if(ngo.eq.7) then
         call PRTOUT
      else if(ngo.eq.8) then
         call MOROUT
      else
c        call COROUT
         call SUICID('FN NOT IMPLEMENTED.  STOP NUMINT', 8)
      endif
      if(L1.lt.0) goto 100
 
      goto kstop
  300 if(T.le.Tstop) goto 500
      goto 200
  400 if(T.lt.Tstop) goto 200
c
c pick up last tabular point
  500 if(ngo.eq.3 .and. Jct(56).gt.0) call SBFOUT(ovrlap, icall)
c
c write completion message
      steps  = K2
      duratn = ABS(T - T0) + Both
      stpday = steps/duratn
      write(Iout, 600) K2, duratn, stpday
  600 format('-NORDSIECK NUMERICAL INTEGRATION HAS BEEN COMPLETED IN',
     .       i8, ' STEPS FOR', 1pd20.13, ' DAYS OR', 1pd12.5,
     .       ' STEPS PER DAY')
      call TIMSET(ihour, iminut, second, itot, gotdat)
      do i = 1, 2
         ww(i)  = itot(i)
         ww(i)  = ww(i)*1E-2_10
         www(i) = ww(i)/duratn
         ww(i)  = ww(i)/steps
      end do
      write(Iout, 700) (loon(i), ihour(i), iminut(i), second(i), ww(i),
     .                 www(i), i = 1, 2)
  700 format(17x, a4, ' TIME ELAPSED=', i3, 'H', i3, 'M', f6.2,
     .       'S  OR', f9.4, 'S  PER STEP  OR', f9.4, 'S PER DAY')
      ww(1) = ww(2)/ww(1)
      write(Iout, 800) ww(1)
  800 format(11x, '(TASK TIME)/(REAL TIME)=', f8.5)
      return
 
  900 write(6, 1000) T0, Tstop
 1000 format(//10x, 'T0=', 1pd25.16, ' TSTOP=', 1pd25.16,
     .       ' PROGRAM STOP')
      return
      end
