      subroutine NINSTP(l2, ngo)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, jj, kn, l2, ngo
 
c*** end of declarations inserted by spag
 
 
c
c   subroutine ''nint'' (version 4)  a.rasinski, dec.1964 (rev. 1024-68)
c           revised july 1965.  branch for fn added oct 67, nov 69.
c           revised sept 1969   variable tabular interval (itype=0)
c           revised dec  1979   step-size control   rdr/kcl
c
c
c
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

c common
      include 'adams.inc'
      include 'dumdum.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'nmstrt.inc'
      include 'output.inc'
      include 'stint.inc'
 
      real*10 hcx, alph, alphz/2._10/
      real*10 dir, tab, hr
 
c double precision statements for function subroutines
      real*10 SBFN, BODFN, ERTFN, PRTFN, MORFN
c
c
c        l1
c        -1  first step of continuing procedure
c         0  starting procedure
c         1  continuing procedure
c
c        k1      count of steps -- not used?
c        k2      count of steps -- used to find end of
c                starting procedure at k2=19 in eval
c                and used to control accuracy tests
c        l4      if .lt. 0 step size not to be changed
c                used to prevent frequent step size changes
c        l5      used in tabular output mode
c                -1 do not change step size
c                set to -1 in tabular output mode except
c                for first step after each tabular point.
c                this prevents the need for a fractional
c                step to get to the next tabular point.
c        n1, n2  range of equation number for integration.
c                first (1 to m), second(m+1 to n).
c        istp       0  do not use stpctl
c                  -1  make transition to stpctl
c                   1  use stpctl
c
      if(L1.lt.0) then
      else if(L1.eq.0) then
         do while( .true. )
c
c starting procedure
c
            dir = Tstop - T0
            Hc  = SIGN(Hc, dir)
            L4  = -30
            K1  = 0
            K2  = 0
            N1  = 1
            N2  = N
            l2  = 0
            T   = T0
            do i = 1, N
               if(Lthrst.le.0) then
                  A(i)  = 0._10
                  B(i)  = 0._10
                  C(i)  = 0._10
                  D(i)  = 0._10
                  Ee(i) = 0._10
               endif
               Y(i, 3) = V0(i)
               Y(i, 4) = V0(i)
            end do
            do i = 1, N
               if(ngo.eq.3) then
                  F(i, 3) = SBFN(i, 3, T)
               else if(ngo.eq.4) then
                  F(i, 3) = BODFN(i, 3, T)
               else if(ngo.eq.5) then
                  F(i, 3) = ERTFN(i, 3, T)
               else if(ngo.eq.7) then
                  F(i, 3) = PRTFN(i, 3, T)
               else if(ngo.eq.8) then
                  F(i, 3) = MORFN(i, 3, T)
               else
c                 F(i, 3) = FN(i, 3, T)
                  call SUICID('FN NOT IMPLEMENTED.  STOP NINSTP', 8)
               endif
               F0(i) = F(i, 3)
               if(.not.(F(i,3).le.0._10.or.F(i,3).gt.0._10)) then
                  write(Iout,40) i,T
   40             format(' NUMERIC FAILURE FOR EQUATION',I4,
     .             ' ITERATION 3 AT TIME',F20.10)
                  call SUICID(
     . 'NON-NUMERIC DERIVATIVE IN STARTING PROCEDURE, STOP IN NINSTP',
     .             15)
               endif
            end do
c
c integration started with 20 steps, reversing sense after
c every 5 steps
c
            do j = 1, 4
               do jj = 1, 5
                  call EVAL(l2, ngo)
                  if(l2.lt.0) goto 50
                  if(l2.gt.0) goto 900
               end do
               Hc = -Hc
               if(T.eq.T0) then
                  do i = 1, N
                     Y(i, 3) = V0(i)
                     Y(i, 4) = V0(i)
                  end do
               endif
            end do
c
c save starting coefficients
            do i = 1, N
               A0(i)  = A(i)
               B0(i)  = B(i)
               C0(i)  = C(i)
               D0(i)  = D(i)
               Ee0(i) = Ee(i)
            end do
            Hsave = Hc
            goto 100
   50    end do
      else
         goto 200
      endif
 
  100 tab = T0 + Hmx
      L1  = 1
      L5  = 0
      dir = Nsign
      do i = 1, N
         Y(i, 3) = V0(i)
         Y(i, 4) = V0(i)
      end do
      do i = 1, N
         if(ngo.eq.3) then
            F(i, 3) = SBFN(i, 3, T)
         else if(ngo.eq.4) then
            F(i, 3) = BODFN(i, 3, T)
         else if(ngo.eq.5) then
            F(i, 3) = ERTFN(i, 3, T)
         else if(ngo.eq.7) then
            F(i, 3) = PRTFN(i, 3, T)
         else if(ngo.eq.8) then
            F(i, 3) = MORFN(i, 3, T)
         else
c           F(i, 3) = FN(i, 3, T)
                  call SUICID('FN NOT IMPLEMENTED.  STOP NINSTP', 8)
         endif
         if(.not.(F(i,3).le.0._10.or.F(i,3).gt.0._10)) then
            write(Iout,140) i,T
  140       format(' NUMERIC FAILURE FOR EQUATION',I4,
     .       ' ITERATION 3 AT TIME',F20.10)
            call SUICID(
     . 'NON-NUMERIC DERIVATIVE IN INTEGRATION,STOP IN NINSTP',
     .       13)
         endif
         Dydt(i, 1) = F(i, 3)
         Dydt0(i)   = Dydt(i, 1)
      end do
c
c continuing procedure
c itype>0 standard call of nordsieck(const. tab interval)
c itype=0 special probe use of nordsieck with variable tab.int.
  200 if(Itype.le.0) then
         if(dir.le.0) then
            if(T0+Hmx.ge.T) goto 1000
         else if(T-T0.ge.Hmx) then
            goto 1000
         endif
      endif
      hr = 0._10
c
c from here to call to eval is logic to keep
c from passing a tab point
      if(dir.le.0) then
         assign 400 to kn
      else
         assign 500 to kn
      endif
  300 N1 = 1
      l2 = 0
      N2 = M
      if(ngo.eq.8) N2=N
      goto kn
  400 if(tab.lt.T+Hc) goto 800
      if(tab.eq.T+Hc) goto 700
      goto 600
  500 if(tab.gt.T+Hc) goto 800
      if(tab.eq.T+Hc) goto 700
c
c a step of hc will pass the tab point
c if(itype.lt.0)go to 40
c
  600 hr = Hc
      Hc = tab - T
c a step of hc will land on the tab point
c
  700 if(L4.ge.0) L4 = -1
  800 call EVAL(l2, ngo)
      if(l2.lt.0) goto 300
      if(l2.eq.0) then
 
         if(tab.ne.T) goto 300
c
c have just hit a tab point.  undo the special
c conditions
         if(hr.ne.0) then
            Hc = hr
            hr = 0._10
         endif
         tab = tab + Hmx
         if(L5.lt.0) then
            L5 = 1
            return
         else if(L5.eq.0) then
            L5 = -1
         endif
      endif
  900 return
 
 1000 if(Istp.le.0) then
         do while( .true. )
            L5 = 0
            N1 = 1
            l2 = 0
            N2 = M
            call EVAL(l2, ngo)
            if(l2.ge.0) then
               if(Istp.eq.0) return
               call STPCTL(hcx)
               alph = alphz
               if(ABS(Hc).gt.alph*hcx) return
               if(ABS(Hc*alph).lt.hcx) return
               if(L4.lt.0) then
                  alph = L4 + 6
                  alph = 1._10 + (alphz - 1._10)*alph/5._10
                  if(ABS(Hc).gt.alph*hcx) return
                  if(ABS(Hc*alph).lt.hcx) return
               endif
               Istp = 1
               N1   = 1
               N2   = N
               write(Iout, 1010) Hc, hcx, K2
 1010          format(
     .            '  STEP SIZE TRANSITION OCCURRED WITH HC, HCX, K2 = ',
     .            2(d12.5,2x), i10)
               Hc = SIGN(hcx, Hc)
               return
            endif
         end do
      else
         call STPCTL(hcx)
         Hc = SIGN(hcx, Hc)
         K2 = 0
         call EVAL(l2, ngo)
         return
      endif
      end
