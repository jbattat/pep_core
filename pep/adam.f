      subroutine ADAM(ngo)
 
      implicit none
c
c w.b.smith   oct 1968    subroutine adam
c adams-moulton integration subroutine
c
c
c
c     ngo=1 planet-asteroid integration (obsolete)
c     ngo=3 satellite-probe integration (called by subroutine planet)
c     ngo=4 n-body integration (called by subroutine planet)
c     ngo=5 earth rotation integration
c     ngo=7 planet rotation integration
c     ngo=8 moon orbit and rotation integration
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'adams.inc'
c npredt = number of predictor terms
c ncoret = number of corrector terms
c neq    = number of equations
c nh     = positive, step size is nh days
c nh    <= zero, step size is 2**nh days
c inthmx = positive, tabular interval is inthmx days
c inthmx<= zero, tabular interval is 2**inthmx days
      include 'coeff.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'nmstrt.inc'
      include 'output.inc'
      include 'stint.inc'
c
c local
      real*10 ytemp,tstopt,xsign
      real*10 dydtnu(6*i_mxeqn),yyy(6*i_mxeqn)
      integer   i,isteps,istrtt,iter,j,k,
     .          kstop,l,l2,mm,nc,ngo,np,nsigns,nst,nstmid
c timer quantities
      integer*4 ihour(2),iminut(2),itot(2)
      real*4    second(2)
      character*4 loon(2)/'REAL','TASK'/
      character*4 gotdat(2)/'GOT ','DATE'/
      character*36 strtmsg
      real*10 steps,duratn,stpday,ww(2),www(2)
c
c external functions
      real*10 BODFN,ERTFN,PRTFN,SBFN,MORFN
c
c
c call sr ''calcof'' to compute coefficients
c
      call CALCOF
c
c
      istrtt = 0
      np     = Npredt+1
      nc     = Ncoret+1
  100 if(Tstop.lt.T0) then
         assign 400 to kstop
      else if(Tstop.eq.T0) then
         return
      else
         assign 500 to kstop
      endif
c
c set up control constants for starting procedure using nint
c to calculate np derivatives prior to t0
      k = 1
      tstopt = Tstop
      xsign  = Nsign
      if(Nh.le.0) then
         Hmx = 2._10**Nh*xsign
      else
         Hmx = Nh*Nsign
      endif
      Tstop  = T0-np*Hmx
      Hmx    = -Hmx
      nsigns = Nsign
      Nsign  = -Nsign
      do while( .true. )
c
c calculation of the np derivatives
         call NINSTP(l2,ngo)
         if(l2.gt.0) return
         k = k+1
         do l = 1,Neq
            if(ngo.eq.3) then
               Dydt(l,k) = SBFN(l,3,T)
            else if(ngo.eq.4) then
               Dydt(l,k) = BODFN(l,3,T)
            else if(ngo.eq.5) then
               Dydt(l,k) = ERTFN(l,3,T)
            else if(ngo.eq.7) then
               Dydt(l,k) = PRTFN(l,3,T)
            else if(ngo.eq.8) then
               Dydt(l,k) = MORFN(l,3,T)
            else
c              Dydt(l,k) = FN(l,3,T)
               call SUICID('FN NOT IMPLEMENTED. STOP IN ADAM',8)
            endif
            if(.not.(Dydt(l,k).le.0._10.or.Dydt(l,k).gt.0._10)) then
               write(Iout,140) l,k
  140          format(' NUMERIC FAILURE FOR EQUATION',I4,' AT STEP',I3,
     .          ' IN STARTING PROCEDURE.')
               call SUICID('NON-NUMERIC DERIVATIVE IN INTEGRATION STARTI
     .NG PROCEDURE, STOP IN ADAM  ',18)
            endif
         end do
         if(k.ge.np) then
c
c extra print of initial accelerations
            if(Kout.gt.0 .and. MOD(Intxpr,2).eq.1)
     .       call KOUTXPR(T0,Dydt0,Neq,6,1)

c set up of control constants for continuing proceedure using
c adams-moulton method
c
            do l = 1,Neq
               Y(l,4) = V0(l)
            end do
            Hc    = -Hmx
            Nsign = nsigns
            if(Inthmx.le.0) then
               Hmx = 2._10**Inthmx*xsign
            else
               Hmx = Inthmx*Nsign
            endif
            Tstop  = tstopt
            nstmid = Hmx/Hc+1E-15_10
            T      = T0
            if(istrtt.eq.0) then
               write(strtmsg,170) K2
  170          format(' STARTING PROCEDURE OF',I7,' STEPS ')
               call TIMRIT(strtmsg,9)
            endif
            istrtt = 1
            goto 300
         endif
      end do
 
  300 do nst = 1,nstmid
c
c prediction of y(n+1)
         do l = 1,Neq
            ytemp = 0._10
            do k = 1,np
               ytemp = ytemp + C(k)*Dydt(l,k)
            end do
 
c 130 y(l,1)=y(l,4)+hc*ytemp
            call XLOAD8(Hc*ytemp)
            call XADD8(Y(l,4))
            call STORND(Y(l,1))
         end do
         j = 1
         do iter = 1,2
c
c calculation of dy/dt at y(n+1)
            do l = 1,Neq
               if(ngo.eq.3) then
                  dydtnu(l) = SBFN(l,j,T+Hc)
               else if(ngo.eq.4) then
                  dydtnu(l) = BODFN(l,j,T+Hc)
               else if(ngo.eq.5) then
                  dydtnu(l) = ERTFN(l,j,T+Hc)
               else if(ngo.eq.7) then
                  dydtnu(l) = PRTFN(l,j,T+Hc)
               else if(ngo.eq.8) then
                  dydtnu(l) = MORFN(l,j,T+Hc)
               else
c                 dydtnu(l) = FN(l,j,T+Hc)
                  call SUICID('FN NOT IMPLEMENTED. STOP IN ADAM',8)
               endif
               if(.not.(dydtnu(l).le.0._10.or.dydtnu(l).gt.0._10)) then
                  write(Iout,340) l,j,T+Hc
  340             format(' NUMERIC FAILURE FOR EQUATION',I4,' ITERATION'
     .             ,I3,' AT TIME',F20.10)
                  call SUICID(
     .       'NON-NUMERIC DERIVATIVE IN INTEGRATION, STOP IN ADAM ',13)
               endif
            end do
            j = j+1
c
c corrector equation being applied twice
            do l = 1,Neq
               ytemp = D(1)*dydtnu(l)
               if(iter.ne.2) then
                  yyy(l) = 0._10
                  do k = 1,Ncoret
                     yyy(l) = yyy(l) + D(k+1)*Dydt(l,k)
                  end do
               endif
 
c 155 y(l,j)=y(l,4)+hc*(ytemp +yyy(l))
               call XLOAD8(Hc*(ytemp+yyy(l)))
               call XADD8(Y(l,4))
               call STORND(Y(l,j))
            end do
         end do
c
c update things, preparatory to the next step.
c
         do l = 1,Neq
            Y(l,4) = Y(l,3)
            do k = 2,np
               mm = np-k+2
               Dydt(l,mm) = Dydt(l,mm-1)
            end do
            Dydt(l,1) = dydtnu(l)
         end do
         T = T+Hc
      end do
c extra print of accelerations
      if(Kout.gt.0 .and. Intxpr.gt.0) call KOUTXPR(T,Dydt,Neq,6,Intxpr)
c
c call output routine
      if(ngo.eq.3) then
         call SBOUT
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
         call SUICID('FN NOT IMPLEMENTED. STOP IN ADAM',8)
      endif
 
c
c test for end
      goto kstop
 
  400 if(T.le.Tstop) goto 600
      goto 300
  500 if(T.lt.Tstop) goto 300
  600 if(L1.lt.0) goto 100
c
c printout time required for integration
      duratn = ABS(T-T0) + Both
      stpday = 1._10/ABS(Hc)
      steps  = duratn*stpday
      isteps = steps+0.01_10
      write(Iout,700) isteps,duratn,stpday
  700 format(
     .     '-ADAMS-MOULTON NUMERICAL INTEGRATION HAS BEEN COMPLETED IN',
     .     i8,' STEPS FOR',1pd20.13,' DAYS OR',1pd12.5,
     .     ' STEPS PER DAY')
      call TIMSET(ihour,iminut,second,itot,gotdat)
      do i = 1,2
         ww(i)  = itot(i)
         ww(i)  = ww(i)*0.01_10
         www(i) = ww(i)/duratn
         ww(i)  = ww(i)/steps
      end do
      write(Iout,800) (loon(i),ihour(i),iminut(i),second(i),ww(i),
     .                 www(i),i = 1,2)
  800 format(17x,a4,' TIME ELAPSED=',i3,'H',i3,'M',f6.2,
     .       'S  OR',f9.4,'S  PER STEP  OR',f9.4,'S PER DAY')
      ww(1) = ww(2)/ww(1)
      write(Iout,900) ww(1)
  900 format(11x,'(TASK TIME)/(REAL TIME)=',f8.5)
      return
      end
