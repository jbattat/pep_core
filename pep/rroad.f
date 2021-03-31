      subroutine RROAD(ngo)
 
      implicit none
c
c w.b. smith  subroutine rroad    11-15-68
c royal road (second sum) method of numerical integration
c
c argument
      integer*4 ngo

c     ngo=1 planet-asteroid integration (obsolete)
c     ngo=2 moon integration (obsolete)
c     ngo=3 satellite-probe integration (called by subroutine planet)
c     ngo=4 n-body integration (called by subroutine planet)
c     ngo=5 earth rotation integration
c     ngo=6 moon rotation integration (obsolete)
c     ngo=7 planet rotation integration
c     ngo=8 moon orbit and rotation integration
c
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'adams.inc'
      include 'coeff.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'nmstrt.inc'
      include 'output.inc'
      include 'stint.inc'

c local variables
      real*10 dydtnu(3*i_mxeqn),yp(6*i_mxeqn,15),dely(3*i_mxeqn),
     . delynu(3*i_mxeqn)
      real*10 hc2,tstopt,xsign,ydemp,yderiv,ytemp
      integer   i,isteps,istrtt,j,k,kstop,l,
     .          l2,lk,lm,ltest,mm,neqhaf,np,npd
      integer   nsigns,nst,nstmid
c
c timer quantities
      integer*4 ihour(2),iminut(2),itot(2)
      real*4    second(2)
      character*4 loon(2)/'REAL','TASK'/
      character*4 gotdat(2)/'GOT ','DATE'/
      character*36 strtmsg
      real*10 duratn,ww(2),www(2),steps,stpday

c external functions
      real*10 BODFN,ERTFN,PRTFN,SBFN,MORFN

c
c-------note that the changes made to the royal road integrator by --
c-------r. king in may 77 to fix the truncation problem on the 370/470--
c-------machines did not restore completely the lincoln 360/67 precison
c-------and that the adams-moulton integrator should be used----------
c
      istrtt = 0
      np     = Npredt+1
      npd    = Itype
c
c call sr ''calcof'' to compute coefficients
c
      call CALCOF
c
c
  100 if(Tstop.lt.T0) then
         assign 300 to kstop
      else if(Tstop.eq.T0) then
         return
      else
         assign 400 to kstop
      endif
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
c for starting, use nord. to compute d2y/dt2 at np points prior to zero
c
         call NINSTP(l2,ngo)
         if(l2.gt.0) return
         k = k+1
         if(k.eq.2) then
            do l = 1,Neq
               yp(l,1)   = V0(l)
               Dydt(l,1) = Dydt0(l)
            end do
         endif
         do l = 1,Neq
            yp(l,k) = Y(l,3)
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
               call SUICID('FN NOT IMPLEMENTED.  STOP RROAD ',8)
            endif
            if(.not.(Dydt(l,k).le.0._10.or.Dydt(l,k).gt.0._10)) then
               write(Iout,140) l,k
  140          format(' NUMERIC FAILURE FOR EQUATION',I4,' AT STEP',I3,
     .          ' IN STARTING PROCEDURE.')
               call SUICID('NON-NUMERIC DERIVATIVE IN INTEGRATION STARTI
     .NG PROCEDURE, STOP IN RROAD ',18)
            endif
         end do
         if(k.ge.np) then
            ltest  = 3
            neqhaf = Neq/2
            do l = 1,neqhaf
               lk = l+ltest
               lm = lk-3
 
c set up y(n)
               Y(lm,4) = V0(lm)
               Y(lk,4) = V0(lk)
               if(l.eq.ltest) ltest = ltest+3
            end do
c
c redefine dydt and define y(n-1) for royal road coding
c
            Hc    = -Hmx
            hc2   = Hc**2
            ltest = 3
            do l = 1,neqhaf
               lk = l+ltest
               lm = lk-3
               dely(l) = Hc*Dydt0(lm)
               do k = 1,np
                  yp(l,k)   = yp(lm,k)
                  Dydt(l,k) = Dydt(lk,k)
                  dely(l)    = dely(l)-hc2*A(k)*Dydt(l,k)
               end do
               Dydt(l,1) = Dydt0(lk)
               if(l.eq.ltest) ltest = ltest+3
            end do

c extra print of initial accelerations
            if(Kout.gt.0 .and. MOD(Intxpr,2).eq.1)
     .       call KOUTXPR(T0,Dydt0,Neq,3,1)
c
c end of starting procedure. restore variables specialized for same.
            Nsign = nsigns
            if(Inthmx.le.0) then
               Hmx = 2._10**Inthmx*xsign
            else
               Hmx = Inthmx*Nsign
            endif
            Tstop  = tstopt
            nstmid = Hmx/Hc+1E-15_10
 
c redefinition of hc and hc2 moved to precede update of dely
            T = T0
            if(istrtt.eq.0) then
               write(strtmsg,170) K2
  170          format(' STARTING PROCEDURE OF',I7,' STEPS ')
               call TIMRIT(strtmsg,9)
            endif
            istrtt = 1
            goto 200
         endif
      end do
 
 
c
c --------------------begin the integration.------------------
c
  200 do nst = 1,nstmid
c compute y(n+1) and y(n+1)' from predictor eqs.
c
c y(n+1) predict
c
         j     = 1
         ltest = 3
         do l = 1,neqhaf
            ytemp = 0._10
            do k = 1,np
               ytemp = ytemp+Cc(k)*Dydt(l,k)
            end do
c added statements accomplish rounding in the incrementing
c of  y  and  y'
 
c delynu(l)=dely(l)+hc2*ytemp
            call XLOAD8(hc2*ytemp)
            call XADD8(dely(l))
            call STORND(delynu(l))
            lm = l+ltest-3
 
c y(lm,j)=y(lm,4)+delynu(l)
            call XLOAD8(delynu(l))
            call XADD8(Y(lm,4))
            call STORND(Y(lm,j))
            if(l.eq.ltest) ltest = ltest+3
         end do
 
         if(Itype.lt.4) then
c
c y(n+1)' predict
            ltest = 3
            do l = 1,neqhaf
               lk    = l+ltest
               ydemp = 0._10
               do k = 1,np
                  ydemp = ydemp+C(k)*Dydt(l,k)
               end do
 
c y(lk,j)=y(lk,4)+hc*ydemp
               call XLOAD8(Hc*ydemp)
               call XADD8(Y(lk,4))
               call STORND(Y(lk,j))
               if(l.eq.ltest) ltest = ltest+3
            end do
         endif
c
c
c compute second deriv. using y-predict
         ltest = 3
         do l = 1,neqhaf
            lk = l+ltest
            if(ngo.eq.3) then
               dydtnu(l) = SBFN(lk,j,T+Hc)
            else if(ngo.eq.4) then
               dydtnu(l) = BODFN(lk,j,T+Hc)
            else if(ngo.eq.5) then
               dydtnu(l) = ERTFN(lk,j,T+Hc)
            else if(ngo.eq.7) then
               dydtnu(l) = PRTFN(lk,j,T+Hc)
            else if(ngo.eq.8) then
               dydtnu(l) = MORFN(lk,j,T+Hc)
            else
c              dydtnu(l) = FN(lk,j,T+Hc)
               call SUICID('FN NOT IMPLEMENTED.  STOP RROAD ',8)
            endif
            if(l.eq.ltest) ltest = ltest+3
            if(.not.(dydtnu(l).le.0._10.or.dydtnu(l).gt.0._10)) then
               write(Iout,240) lk,j,T+Hc
  240          format(' NUMERIC FAILURE FOR EQUATION',I4,' ITERATION'
     .          ,I3,' AT TIME',F20.10)
               call SUICID(
     .       'NON-NUMERIC DERIVATIVE IN INTEGRATION, STOP IN RROAD',13)
            endif
         end do
 
c update second derivatives
         do l = 1,neqhaf
            do k = 2,np
               mm = np-k+2
               Dydt(l,mm) = Dydt(l,mm-1)
            end do
            Dydt(l,1) = dydtnu(l)
         end do
 
         j = 3
c
c
c compute y(n+1) and y(n+1)' from corrector eqs.
         ltest = 3
         do l = 1,neqhaf
            ytemp = 0._10
            do k = 1,np
               ytemp = ytemp+Dd(k)*Dydt(l,k)
            end do
 
c delynu(l)=dely(l)+hc2*ytemp
            call XLOAD8(hc2*ytemp)
            call XADD8(dely(l))
            call STORND(delynu(l))
            lm = l+ltest-3
c
c y(n+1) correct
 
c y(lm,j)=y(lm,4)+delynu(l)
            call XLOAD8(delynu(l))
            call XADD8(Y(lm,4))
            call STORND(Y(lm,j))
            if(l.eq.ltest) ltest = ltest+3
         end do
 
         if(Itype.lt.4) then
 
c y(n+1)' correct
            ltest = 3
            do l = 1,neqhaf
               lk    = l+ltest
               ydemp = 0._10
               do k = 1,np
                  ydemp = ydemp+D(k)*Dydt(l,k)
               end do
 
c y(lk,j)=y(lk,4)+hc*ydemp
               call XLOAD8(Hc*ydemp)
               call XADD8(Y(lk,4))
               call STORND(Y(lk,j))
               if(l.eq.ltest) ltest = ltest+3
            end do
         endif
c
c compute derivative using numerical differ. formula
         if(Itype.gt.4) then
            ltest = 3
            do l = 1,neqhaf
               lk     = l+ltest
               yderiv = 0._10
               call XLOAD8(yderiv)
               do k = 1,npd
                  call XADD8(Pf(k)*yp(l,k))
               end do
 
c 162 yderiv   =yderiv   +pf(k)*yp(l,k)
               call STORND(yderiv)
               Y(lk,3) = yderiv/Hc
               if(l.eq.ltest) ltest = ltest+3
            end do
         endif
c
c update y(n)   which is equal to y(l,3)
c update y(n)'  which is equal to y(l+2,3)
c
         ltest = 3
         do l = 1,neqhaf
            lm = l+ltest-3
 
c update of y(n-1) deleted
            Y(lm,4) = Y(lm,3)
            dely(l)  = delynu(l)
            if(Itype.ne.4) then
               lk = lm+3
               Y(lk,4) = Y(lk,3)
            endif
 
            if(l.eq.ltest) ltest = ltest+3
         end do
         if(Itype.gt.4) then
c
c update y for numerical differ. eqs.
            ltest = 3
            do l = 1,neqhaf
               lm = l+ltest-3
               do k = 2,npd
                  mm = npd-k+2
                  yp(l,mm) = yp(l,mm-1)
               end do
               yp(l,1) = Y(lm,4)
               if(l.eq.ltest) ltest = ltest+3
            end do
         endif
         T = T+Hc
      end do
c extra print of accelerations
      if(Kout.gt.0 .and. Intxpr.gt.0) call KOUTXPR(T,Dydt,Neq,3,Intxpr)

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
         call SUICID('FN NOT IMPLEMENTED.  STOP RROAD ',8)
      endif
 
      goto kstop
 
  300 if(T.le.Tstop) goto 500
      goto 200
  400 if(T.lt.Tstop) goto 200
  500 if(L1.lt.0) goto 100
c
c printout time required for integration
      duratn = ABS(T-T0)+Both
      stpday = 1._10/ABS(Hc)
      steps  = duratn*stpday
      isteps = steps+0.01_10
      write(Iout,600) isteps,duratn,stpday
  600 format(
     .     '-  ROYAL ROAD  NUMERICAL INTEGRATION HAS BEEN COMPLETED IN',
     .     i8,' STEPS FOR',1pd20.13,' DAYS OR',1pd12.5,
     .     ' STEPS PER DAY')
      call TIMSET(ihour,iminut,second,itot,gotdat)
      do i = 1,2
         ww(i)  = itot(i)
         ww(i)  = ww(i)*1E-2_10
         www(i) = ww(i)/duratn
         ww(i)  = ww(i)/steps
      end do
      write(Iout,700) (loon(i),ihour(i),iminut(i),second(i),ww(i),
     .                 www(i),i = 1,2)
  700 format(17x,a4,' TIME ELAPSED=',i3,'H',i3,'M',f6.2,
     .       'S  OR',f9.4,'S  PER STEP  OR',f9.4,'S PER DAY')
      ww(1) = ww(2)/ww(1)
      write(Iout,800) ww(1)
  800 format(11x,'(TASK TIME)/(REAL TIME)=',f8.5)
      return
      end
