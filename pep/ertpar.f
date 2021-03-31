      subroutine ERTPAR(kick)
 
      implicit none
c
c m.e.ash    nov 1969    subroutine ertpar
c calculate partials w.r.t. precession,obliquity and et-ut1 paramtrs
c f.a.kreimendahl   nov 1979   cleaned up and reno'd
c partials w.r.t. annual wobble terms added

c array dimensions
      include 'globdefs.inc'
c
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 x(6,3),xo(3)
      equivalence (Xm(1,1),x(1,1)),(Angdum(1),xo(1))
      include 'dtparm.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'leon.inc'
      real*10 smoblq
      equivalence (Sobliq,smoblq)
      include 'ltrapx.inc'
      integer*2 kind
      equivalence (kind,Numpar)
      integer*2 ildt(6)
      equivalence (ildt,Ildt1)
      include 'mtrapx.inc'
      integer*2 imdt(6)
      equivalence (imdt,Imdt1)
      include 'number.inc'
      include 'nutprc.inc'
      real*4 dpsi(2),deps(2)
      equivalence (Nutat(1),dpsi),(Nutat(3),deps)
      include 'obscrd.inc'
      include 'param.inc'
      include 'partcm.inc'
      include 'sitcrd.inc'
c
c local
      real*10 xtemp(3),xtmpep(3,2),ytemp(3),xtempd(3),ytempd(3),
     . dxb(2),dy(6,2)
      real*10 dd(2),ddeps,ddpsi,ddeps2,dpc2,dps2,dtda,DOT
      equivalence (dd(1),ddeps),(dd(2),ddpsi)
      integer i, iflag, ii, iitt, il, ir, ixy, j, j0, kick, l, m, mgo,
     .          mtype, n, ngo, nixon, ntype
c
c*  start=100
      ngo = 1
      m   = 7
      l   = 6
 
      call PCOPS(m,'ER  ', Iabs1)
  100 iflag = 1
      call PCOPY(l,u_nmbod,iflag,1,Lerx,Merx)
      if(iflag.gt.0) goto 1300
c
c           erotat con's 1-21 are shared by two mutually-exclusive sets
c           of parameters.  check to see which applies:
c              jct(29)=0
c                1- 6:  ad hoc small rotations and rates
c                7-10:  open
c               11-13:  polynomial in ct-ut2 or a1-ut1--see jct(34)
c               14-21:  seasonal terms in ct-ut2 or a1-ut1--see ict(34)
c              jct(29)>0
c                1- 2:  free core nutation (nearly-diurnal free wobble)
c                3-10:  in-phase components of adjustable nutation terms
c                       14-day, 183-day, 365-day, & 18.6-yr
c               11-18:  out-of-phase components of adjustable nutation t
c               19-20:  (k/c) for fortnightly & monthly ut1 terms
c                  21:  open
      if(Lerx(l).gt.21) goto 600
      if(Jct(29).le.0) goto 150
c
c Jct(29)>0 ....................
c
      if(Lerx(l).gt.20) call SUICID(
     .    'LER(L)>20 WITH JCT(29)>0, STOP IN ERTPAR', 10)

      if(Lerx(l).gt.18) then
c radar and radio observations for a1-ut1 tidal partials
c ler(l)= 19-20 and jct(29).gt.0
         if(kick.eq.2 .or. kick.eq.3) call SUICID(
     .'CANNOT CALCULATE UT1 FOR OPTICAL OBSERVABLE, STOP IN ERTPAR ',15)
         ntype = Lerx(l) - 18
         dtda  = Ut1par(ntype)*Sidvel/Aultsc
         do n = 1, Mouse
            dy(1,n) = Y(2,n)*dtda
            dy(2,n) = -Y(1,n)*dtda
            dy(3,n) = 0._10
            call CHNCOR(derem(1,n),dy(1,n))
            if(Index.gt.3) then
               dy(4,n) = -dy(2,n)*Sidvel*Secday
               dy(5,n) = dy(1,n)*Sidvel*Secday
               dy(6,n) = 0._10
               do j=1,3
                  derem(j+3,n) =
     .             DOT(Dnutpr(1,j),dy(1,n)) + DOT(Nutpr(1,j),dy(4,n))
               end do
            endif
            if(Ncodf.gt.20) then
               do i = 1, Index
                  Derpr(i,n,2) = Derpr(i,n,1)
               end do
            endif
         end do
         goto 1000
      else if(Lerx(l).le.2) then
c partial is w.r.t. free core nutation parameters
         if(Ercom(2).eq.0._10) goto 1100
         if(Lerx(l).eq.1) then
            ddeps = Cfnut*Nutn/Nutn0
            ddpsi = Sfnut/smoblq*Nutn/Nutn0
            if(Index.gt.3) then
               dpc2  = (-Sfnut*(Dmoblq+Deps(2))*Sobliq/smoblq -
     .          Cfnut*Ps(2))*Nutn/Nutn0
               dps2  = (Sfnut*(Dmoblq+Deps(2))*Cobliq/smoblq +
     .          Cfnut*Pc(2))*Nutn/Nutn0
               ddeps2= 0._10
            endif
         else
            ddeps = -Sfnut*Nutn/Nutn0
            ddpsi = Cfnut/smoblq*Nutn/Nutn0
            if(Index.gt.3) then
               dpc2  = (-Cfnut*(Dmoblq+Deps(2))*Sobliq/smoblq +
     .          Sfnut*Ps(2))*Nutn/Nutn0
               dps2  = (Cfnut*(Dmoblq+Deps(2))*Cobliq/smoblq -
     .          Sfnut*Pc(2))*Nutn/Nutn0
               ddeps2= 0._10
            endif
         endif
         goto 900
      endif

c partial is w.r.t. 14, 183, 365, or 6798 day resonance term
      ir = Lerx(l)
      if(ir.gt.10) then
 
c partial is w.r.t. out-of-phase components
         if(mod(ir,2).eq.0) then
 
c partial is w.r.t. cos (dpsi) terms  ler(l)= 12,14,16,18
            ir    = ir/2 - 5
            ddeps = 0._10
            ddpsi = Carg(ir)*Nutn/Nutn0
            if(Index.gt.3) then
               dpc2  = -Carg(ir)*(Dmoblq+Deps(2))*Sobliq*Nutn/Nutn0
               dps2  = Carg(ir)*(Dmoblq+Deps(2))*Cobliq*Nutn/Nutn0
               ddeps2= 0._10
            endif
         else
 
c partial is w.r.t. sin (deps) terms  ler(l)= 11,13,15,17
            ir    = ir/2 - 4
            ddeps = Sarg(ir)*Nutn/Nutn0
            ddpsi = 0._10
            if(Index.gt.3) then
               dpc2  = -Sarg(ir)*Ps(2)*Nutn/Nutn0
               dps2  = Sarg(ir)*Pc(2)*Nutn/Nutn0
               ddeps2= 0._10
            endif
         endif
      else if(mod(ir,2).eq.0) then
 
c partial is w.r.t. sin terms   (dpsi)   ler(l)= 4,6,8,10
         ir    = ir/2 - 1
         ddeps = 0._10
         ddpsi = Sarg(ir)*Nutn/Nutn0
         if(Index.gt.3) then
            dpc2  = -Sarg(ir)*(Dmoblq+Deps(2))*Sobliq*Nutn/Nutn0
            dps2  = Sarg(ir)*(Dmoblq+Deps(2))*Cobliq*Nutn/Nutn0
            ddeps2= 0._10
         endif
      else
 
c partial is w.r.t. cos terms   (deps)   ler(l)= 3,5,7,9
         ir    = ir/2
         ddeps = Carg(ir)*Nutn/Nutn0
         ddpsi = 0._10
         if(Index.gt.3) then
            dpc2  = -Carg(ir)*Ps(2)*Nutn/Nutn0
            dps2  = Carg(ir)*Pc(2)*Nutn/Nutn0
            ddeps2= 0._10
         endif
      endif
      goto 900
c
c Jct(29)=0 ....................
c
  150 if(Lerx(l).le.6) then
c partials w.r.t. ad hoc small rotation angles and rates,
c lerx(l)=1-6 and jct(29)=0
         nixon = Lerx(l) + 2
         goto 700
      endif
c
c lerx(l)= 7-10 undefined for jct(29)=0
      if(Lerx(l).le.10) call SUICID(
     . 'LER(L)>6 AND <11 WITH JCT(29)=0,STOP IN ERTPAR  ', 12)
c
c partials w.r.t. ut2-ut1 or a1-ut1 parameters
c ler(l)=11-21 and jct(29)=0
      ntype = Lerx(l) - 13
      if(ntype.le.0) then
         dtda = 1._10
         if(ntype.eq.-1) dtda = Ttutf
         if(ntype.eq.0) dtda  = Ttutf**2
      else
         mtype = (ntype + 1)/2
         dtda  = Tcos(mtype)
         if(ntype.eq.2*mtype) dtda = dtda*Ttutf
      endif
c*  start=200
c kick=1-4, explained in partl
      if(kick.eq.2 .or. kick.eq.3) goto 300
 
c radio observables
  200 if(Nice.ge.0) call SUICID(
     .' CANNOT CALCULATE UT1 PARTIALS FOR RATE OBSERVABLE, STOP ERTPAR '
     ., 16)
      dtda = dtda*Sidvel/Aultsc
      do n = 1, Mouse
         dy(1,n) = Y(2,n)*dtda
         dy(2,n) = -Y(1,n)*dtda
         dy(3,n) = 0._10
         call CHNCOR(derem(1,n),dy(1,n))
         if(Ncodf.gt.20) then
            do i = 1, 3
               Derpr(i,n,2) = Derpr(i,n,1)
            end do
         endif
      end do
      Ivze(1) = 1
      call CPARTC(kick)
      if(ngo.eq.1) goto 1200
      if(ngo.eq.2) goto 1400
c
c optical observables
  300 if(Jds.ge.2435490) goto 1100
      if(Jct(34).le.0 .or. kick.ne.3) call SUICID(
     .'CANNOT CALCULATE PARTIAL W.R.T. POLYNOMIAL CT-UT, STOP IN ERTPAR'
     ., 16)
c transit observation partial w.r.t. ut2-ut1 and
c et-ut2 parameters
  400 if(kick.ne.2) then
         Deriv(kind,1) = -dtda
         Deriv(kind,2) = -dtda
         if(ngo.eq.1) goto 1200
         if(ngo.eq.2) goto 1400
      endif
c optical observation partial w.r.t. ut2-ut1 and
c et-ut2 parameters
  500 if(Ncodf.le.4) dtda = dtda/(1._10 +
     . (xo(1)*Xsitep(5,1)-xo(2)*Xsitep(4,1))
     . /(Sidvel*(xo(1)**2+xo(2)**2)))
      do i = 1, 3
         Zp(i) = Xsitep(i+3,1)*dtda
      end do
      call CPARTO
      if(ngo.eq.1) goto 1200
      if(ngo.eq.2) goto 1400
c
c*  start=400
c partial w.r.t. precession, obliquity, and nutation constants
  600 if(Lerx(l).eq.22) then
         nixon = 1
      else if(Lerx(l).eq.23) then
         nixon = 2
      else if(Lerx(l).eq.24) then
c
c partial is w.r.t. nutation constant, lerx(l)=24
         ddeps = deps(1)/Nutn
         ddpsi = dpsi(1)/Nutn
         if(Index.gt.3) then
            dpc2  = Pc(2)/Nutn
            dps2  = Ps(2)/Nutn
            ddeps2= deps(2)/Nutn
         endif
         goto 900
      else
c
c error message for earth rotation partial derivatives
         call SUICID(
     .' CANNOT CALCULATE PARTIAL DERIVATIVES  W.R.T. EARTH ROTATION PARA
     .METERS, STOP IN ERTPAR ', 22)
      endif
  700 if(kick.eq.1 .or. kick.eq.4) then
c
c radar and radio observations for prec,oblq partials
         do n = 1, Mouse
            xtemp(1) = Y(1,n) + (Pc(1)*Y(2,n) + Ps(1)*Y(3,n))
            xtemp(2) = Y(2,n) - (Pc(1)*Y(1,n) - deps(1)*Y(3,n))
            xtemp(3) = Y(3,n) - (Ps(1)*Y(1,n) + deps(1)*Y(2,n))
            j0 = 0
            do while(.true.)
               do j = 1, 3
                  derem(j+j0,n) = DOT(Dhprec(1,j,nixon),xtemp)/Aultsc
                  if(Nplnt2.ne.0) Derpr(j+j0,n,2) = derem(j+j0,n)
               end do
               if(j0.gt.0) goto 750
               if(Index.le.3) goto 750
               j0 = 3
               xtemp(1) = (-Y(2,n)+Pc(1)*Y(1,n))*Sidvel*Secday
               xtemp(2) = (Y(1,n)+Pc(1)*Y(2,n))*Sidvel*Secday
               xtemp(3) = (Ps(1)*Y(2,n)-deps(1)*Y(1,n))*Sidvel*Secday
            end do
  750    end do
         goto 1000
      else if(kick.eq.3) then
         goto 1100
      else
c
c meridian circle observations for prec,oblq partials
         if(Ncodf.gt.4) goto 1100
         call PRODCT(Dhprec(1,1,nixon),x(1,Kst),xtemp,3,3,1)
         Zp(1) = -xtemp(1) + (Pc(1)*xtemp(2) + Ps(1)*xtemp(3))
         Zp(2) = -xtemp(2) - (Pc(1)*xtemp(1) - deps(1)*xtemp(3))
         Zp(3) = -xtemp(3) - (Ps(1)*xtemp(1) + deps(1)*xtemp(2))
 
c partial of nutation w.r.t. obliquity (nixon=2) ignored
         call CPORTO
         goto 1200
      endif

c
c radar and radio observations for nutation partials
  900 do n = 1, Mouse
c
c get dN/dx times Y
c note: Dnutde is logically extended by Dnutdp, as ddeps is by ddpsi
         call PRODCT(Dnutde,Y(1,n),xtmpep,-6,3,1)
         call PRODCT(ddeps,xtmpep,xtemp,1,-2,3)

         ytemp(1) = (-Xb1(1,n)*Sv(n) + Xb1(2,n)*Cv(n))*ddpsi*Cobliq
         ytemp(2) = (Xb1(1,n)*Cv(n) + Xb1(2,n)*Sv(n))*ddpsi*Cobliq
         ytemp(3) = 0._10
         do j = 1, 3
            derem(j,n) = (DOT(Prec(1,j),xtemp) + DOT(Nutpr(1,j),ytemp))
     .       /Aultsc*Convds
            if(Nplnt2.ne.0) Derpr(j,n,2) = derem(j,n)
         end do
         if(Index.gt.3) then
            ytempd(1) = -ytemp(2)*Sidvel
            ytempd(2) = ytemp(1)*Sidvel
            ytempd(3) = 0._10
c velocity vector in mean orientation frame
            xtempd(1) = -Y(2,n)*Sidvel
            xtempd(2) = Y(1,n)*Sidvel
            xtempd(3) = 0._10
c get dN/dx times Ydot
            call PRODCT(Dnutde,xtempd,xtmpep,-6,3,1)
            call PRODCT(ddeps,xtmpep,xtempd,1,-2,3)
c add dNdot/dx times Y
            xtempd(1) = xtempd(1) + dpc2*Y(2,n) + dps2*Y(3,n)
            xtempd(2) = xtempd(2) - dpc2*Y(1,n) + ddeps2*Y(3,n)
            xtempd(3) = xtempd(3) - dps2*Y(1,n) - ddeps2*Y(2,n)
            do j=1,3
               derem(j+3,n) = (DOT(Dprec(1,j),xtemp) + 
     .          DOT(Dnutpr(1,j),ytemp) + DOT(Nutpr(1,j),ytempd) +
     .          DOT(Prec(1,j),xtempd))/Aultsc*Convds*Secday
               if(Nplnt2.ne.0) Derpr(j+3,n,2) = derem(j+3,n)
            end do
         endif
      end do
c
c*  start=1000
 1000 Ivze(1) = 1
      call CPARTC(kick)
      goto 1200
 
 1100 Deriv(kind,1) = 0._10
      Deriv(kind,2) = 0._10
 
 1200 if(l.lt.u_nmbod) goto 100
c
c*  start=1200
c partials w.r.t. et-ut2, a1-ut1, or wobble parameters
 1300 if(Numdtx.le.0) then
c
c count over dt partials on input tape
         if(Mumdtx.gt.0) then
            do i = 1, 6
               if(imdt(i).gt.0) Mind = Mind + 1
            end do
         endif
         return
      endif
      Ildt1 = I1
      Ildt2 = I2
      Ildt3 = I3
      Ildt4 = I4
      Ildt5 = I5
      Ildt6 = I6
      iitt  = Iabs1
      if(Mumdtx.le.0) iitt = 0
      if(iitt.gt.0) then
 
c input obslib has dt parameters - check for match
         do i = 1, 6
            if(imdt(i).ne.ildt(i) .and. imdt(i).ne.0) then
               write(Iout,1310) i,ildt(i),i,imdt(i)
 1310          format('0*** DT TABULAR POINT I', i1, ' =', i4,
     .                ' DOES NOT MATCH IMDT', i1, ' =', i4, ' ***')
               call SUICID('DT MISMATCH, STOP IN ERTPAR ', 7)
            endif
         end do
      endif
c
c no saved partials, set up pointers for computing
      if(Jddt0.lt.0) then
 
c no a1-ut1 or et-ut2
         Ildt1 = 0
         Ildt2 = 0
      endif
      if(Jddt0.gt.0) then
 
c no wobble
         do i = 3, 6
            ildt(i) = 0
         end do
      endif
 
c set flags to zero for non-adjustable tabular points
      do i = 1, 6
         ii = ildt(i)
         if(ii.gt.0 .and. Ldtx(ii).le.0) ildt(i) = 0
      end do
c
c copy or compute partials as needed
      call PCOPS(1,'DT  ', iitt)
      il  = 0
      ngo = 2
 1400 do while( .true. )
 
c*  start=1400
         iflag = 0
         call PCOPY(il,6,iflag,1,ildt,imdt)
         if(iflag.gt.0) return
 
c must compute partial: il=1,2 means a1-ut1. 3-6 means wobble.
         dtda = Ff
         if(mod(il,2).ne.0) dtda = 1._10 - Ff
         mgo = il - 2
         if(mgo.le.0) then
            if(kick.eq.1 .or. kick.eq.4) goto 200
            if(kick.eq.2) goto 500
            if(kick.eq.3) goto 400
         endif
c
c radio, radar, laser, or interferometry observation partial
c w.r.t. wobble parameters
         if(Nice.ge.0) call SUICID(
     .' CANNOT CALCULATE UT1 PARTIALS FOR RATE OBSERVABLE, STOP ERTPAR '
     ., 16)
         ixy    = (mgo + 1)/2
         dxb(1) = 0._10
         dxb(2) = 0._10
         do n = 1, Mouse
            dxb(ixy) = -dtda*Rs(n)
            dy(1,n) = (Cv(n)*dxb(1) + Sv(n)*dxb(2))/Aultsc*Convds
            dy(2,n) = (Sv(n)*dxb(1) - Cv(n)*dxb(2))/Aultsc*Convds
            dy(3,n) = dtda*Xb0(ixy,n)/Aultsc*Convds
            call CHNCOR(derem(1,n),dy(1,n))
            if(Ncodf.gt.20) then
               do i = 1, 3
                  Derpr(i,n,2) = Derpr(i,n,1)
               end do
            endif
         end do
         Ivze(1) = 1
         call CPARTC(kick)
      end do
c
      end
