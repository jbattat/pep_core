      subroutine EVAL(l2,ngo)
 
      implicit none
 
 
c*** start of declarations inserted by spag
 
c*** end of declarations inserted by spag
 
 
c
c w.b.smith   march 1966     subroutine eval
c branch for fn added oct 1967, nov 1969
c
c arguments
      integer*4 l2,ngo
c
c     ngo=1 planet-asteroid integration (obsolete)
c           planet)
c     ngo=2 moon integration (obsolete)
c     ngo=3 satellite-probe integration (numint called by planet)
c     ngo=4 n-body integration (numint called by subroutine planet)
c     ngo=5 earth rotation integration
c     ngo=6 moon  rotation integration (obsolete)
c     ngo=7 planet rotation integration
c     ngo=8 moon orbit and rotation integration
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'dumdum.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'output.inc'
      include 'stbtst.inc'
 
c double precision statements for function subroutines
      real*10 SBFN, BODFN
      real*10 ERTFN, PRTFN, MORFN
 
c quantities in subroutine
      integer*4 i,iter,j,k,ka,kb
      real*10 amx,bmx,emx,sta,g
      real*10 w(6),h(6),e(6*i_mxeqn),hh(6)
      data hh(2),hh(3),hh(4),hh(5),hh(6)/2._10,3._10,4._10,5._10,6._10/
      data w/3.15591931216919312169E-1_10,
     .     2.28333333333333333333_10,3.75000_10,4.250000_10,
     .     3.00000_10,1.00000_10/
      real*10 fdum
      real*10 whc(6)
      real*10 whc1,whc2,psc(3)
      real*10 ytemp
c
c        intm
c         0  correct 2 times -- old form of integration
c        -1  royal road mod. -- 1 correction
c         1  *.-1 plus pseudo correction with dadx matrix
c             (not tested)
c
c
c           y and f predictions from data computed in last successful
c           step
  100 h(1) = Hc
      whc1 = w(1)*Hc
      whc2 = whc1*whc1
      do i = 2,6
         h(i) = Hc/hh(i)
      end do
      do k = N1,N2
         i = k
         F(i,1) = F(i,3) + h(1)
     .            *(A(i)+h(2)*(B(i)+h(3)*(C(i)+h(4)*(D(i)+h(5)*Ee(i)))))
 
c added statements accomplish rounding in the incrementing of y
         ytemp = h(1)
     .              *(F(i,3)+h(2)*(A(i)+h(3)*(B(i)+h(4)*(C(i)+h(5)*
     .              (D(i)+h(6)*Ee(i))))))
         call XLOAD8(ytemp)
         call XADD8(Y(i,4))
         call STORND(Y(i,1))
      end do
c     y(i,1) = y(i,4) + h(1)*(f(i,3) + h(2)*(a(i) + h(3)*(b(i) + h(4)*
c    1 (c(i)+h(5)*(d(i)+h(6)*ee(i))))))
c
c           y is corrected        iteration 1
c
      do k = N1,N2
         i = k
         if(ngo.eq.3) then
            fdum = SBFN(i,1,T + Hc)
         else if(ngo.eq.4) then
            fdum = BODFN(i,1,T + Hc)
         else if(ngo.eq.5) then
            fdum = ERTFN(i,1,T + Hc)
         else if(ngo.eq.7) then
            fdum = PRTFN(i,1,T + Hc)
         else if(ngo.eq.8) then
            fdum = MORFN(i,1,T + Hc)
         else
c           fdum = FN(i,1,T + Hc)
            call SUICID('FN NOT IMPLEMENTED. STOP IN EVAL',8)
         endif
         if(.not.(fdum.le.0._10.or.fdum.gt.0._10)) then
            iter=1
            goto 130
         endif
         F(i,2)= fdum
         e(i)  = F(i,2) - F(i,1)
         ytemp = whc1*e(i)
         call XLOAD8(ytemp)
         call XADD8(Y(i,1))
         call STORND(Y(i,2))
      end do
      goto 170
c failure
  130 write(Iout,140) i,iter,T+Hc
  140 format(' NUMERIC FAILURE FOR EQUATION',I4,' ITERATION'
     . ,I3,' AT TIME',F20.10)
      call SUICID(
     . 'NON-NUMERIC DERIVATIVE IN INTEGRATION, STOP IN EVAL ',13)

c     y(i,2)=y(i,1)+whc1*e(i)
c
  170 if(Intm.ne.0) then
c make royal road correction
c
         do k = N1,N2,6
            ka = k
            kb = ka + 2
            do i = ka,kb
               ytemp = whc2*e(i+3)
               call XLOAD8(ytemp)
               call XADD8(Y(i,2))
               call STORND(Y(i,2))
 
c     y(i,2)=y(i,2)+whc2*e(i+3)
               e(i) = e(i) + whc1*e(i+3)
            end do
         end do
 
         if(Intm.gt.0) then
c make pseudo correction
c
            call SUICID(
     .       'OPTION KK(88)=1 NOT TESTED OR DEBUGGED, STOP IN EVAL',13)
            do k = N1,N2,6
               ka = k
               kb = ka + 2
               do i = ka,kb
                  psc(i) = 0._10
                  do j = ka,kb
                     psc(i) = psc(i) + e(j)*Dadx(i,j)*whc1
                  end do
               end do
               do i = ka,kb
                  j    = i + k - 1
                  e(j) = e(j) + psc(j)*whc1
                  e(j+3) = e(j+3) + psc(j)
                  ytemp = whc1*e(j)
                  call XLOAD8(ytemp)
                  call XADD8(Y(j,1))
                  call STORND(Y(j,3))
 
c     y(j,3)=y(j,1)+whc1*e(j)
                  ytemp = whc1*e(j+3)
                  call XLOAD8(ytemp)
                  call XADD8(Y(j+3,1))
                  call STORND(Y(j+3,3))

c     y(j+3,3)=y(j+3,1)+whc1*e(j+3)
               end do
            end do
         endif
 
         do i = N1,N2
            F(i,2) = F(i,1) + e(i)
         end do
         goto 200
      endif
c
c y is corrected        iteration 2
c
      do k = N1,N2
         i = k
         if(ngo.eq.3) then
            F(i,2) = SBFN(i,2,T + Hc)
         else if(ngo.eq.4) then
            F(i,2) = BODFN(i,2,T + Hc)
         else if(ngo.eq.5) then
            F(i,2) = ERTFN(i,2,T + Hc)
         else if(ngo.eq.7) then
            F(i,2) = PRTFN(i,2,T + Hc)
         else if(ngo.eq.8) then
            F(i,2) = MORFN(i,2,T + Hc)
         else
c           F(i,2) = FN(i,2,T + Hc)
            call SUICID('FN NOT IMPLEMENTED. STOP IN EVAL',8)
         endif
         if(.not.(F(i,2).le.0._10.or.F(i,2).gt.0._10)) then
            iter=2
            goto 130
         endif
         e(i)  = F(i,2) - F(i,1)
         ytemp = w(1)*e(i)*Hc
         call XLOAD8(ytemp)
         call XADD8(Y(i,1))
         call STORND(Y(i,3))
 
c     y(i,3) = y(i,1) + w(1) * e(i)*hc
      end do
  200 K1 = K1 + 1
 
      if(K2.lt.19) goto 500
      if(K2.gt.19 .and. N1.gt.1) goto 500
c
c interval control logic
c
      g   = Epsi/ABS(Hc)
      amx = ABS(Y(1,3) - Y(1,2))
      bmx = ABS(Y(1,2) - Y(1,1))
      emx = ABS(e(1))
      do i = 2,M
         amx = MAX(amx,ABS(Y(i,3)-Y(i,2)))
         bmx = MAX(bmx,ABS(Y(i,2)-Y(i,1)))
         emx = MAX(emx,ABS(e(i)))
      end do
      amx = amx + 1E-35_10
      bmx = bmx + 1E-30_10
c
c stability test is performed if bmx is g.t. epsi/16
      sta = Epsi/16._10
      if(bmx.gt.sta) then
         if(Iptr3.lt.3) goto 300
         if(amx.gt.(1.25E-1_10*bmx)) goto 300
      endif
c
c accuracy test
c
      if(emx.lt.g) then

c tests satisfied. attempt made to double interval
         if(L5.lt.0) goto 400
         if(L4.lt.0) goto 400
         if(260._10*emx.ge.g) goto 400
         if(Iptr3.ge.3) then
            if(25.6_10*amx.gt.bmx) goto 400
         endif
 
         Hc = Hc*2._10
         L4 = -5
         l2 = -1
         goto 600
      else if(emx.eq.g) then
         goto 400
      endif
c
c test unsatisfied. interval halved.
  300 if(Nnotst.le.0) then
 
         Hc = Hc/2._10
         L4 = -5
         if(ABS(Hc).le.Hmn) then
            call SUICID('INTERVAL TO MINIMUM IN EVAL.  STOP, '//
     .                  'GO TO NEXT INTEGRATION JOB  ',-16)
 
            l2 = 1
         else
            l2 = -1
         endif
         goto 600
      else
 
c skip interval halving for filter span stepover
         write(Iout,350)
  350    format(' INTERVAL HALVING SKIPPED FOR FILTER SPAN STEPOVER')
      endif

c no change in step size, finish equations up to N if not yet done
  400 if(N2.lt.N) then
         N1 = N2 + 1
         N2 = N
         goto 100
      endif

c A,B,C,D,Ee are corrected. Y,F,T are updated. K2 and L4 are
c advanced
  500 T = T + Hc
      L4 = L4 + 1
      K2 = K2 + 1
      if(L5.ne.0) L5 = -1
      whc(1) = 1._10
      do i = 2,6
         whc(1) = whc(1)*Hc
         whc(i) = w(i)/whc(1)
      end do
      do i = 1,N
         F(i,3) = F(i,2)
         Y(i,4) = Y(i,Iptr3)
         A(i)   = A(i)+h(1)*(B(i)+h(2)*(C(i)+h(3)*(D(i)+h(4)*Ee(i))))
     .            + whc(2)*e(i)
         B(i)   = B(i)+h(1)*(C(i)+h(2)*(D(i)+h(3)*Ee(i))) + whc(3)*e(i)
         C(i)   = C(i)+h(1)*(D(i)+h(2)*Ee(i)) + whc(4)*e(i)
         D(i)   = D(i)+h(1)*Ee(i) + whc(5)*e(i)
         Ee(i)  = Ee(i) + whc(6)*e(i)
      end do
 
  600 if(Iptr3.ne.3) then
 
c temp  temp  temp
         do i = 1,N
            Y(i,3) = Y(i,Iptr3)
         end do
      endif
 
c temp  temp  temp
      return
      end
