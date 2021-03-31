      subroutine PARTL1(kick)
 
      implicit none

c
c m.e.ash   feb 1970    subroutine partl1
c calculate radar,optical,transit,interferometer partial derivatives
c continuation of subroutine partl
c

c arguments
      integer*4 kick
c           kick =1 subroutine radar is calling program
c           kick =2 subroutine optic is calling program
c           kick =3 subroutine trnsit is calling program
c           kick =4 subroutine fermtr is calling program
c           partl calculates partial derivatives of observations with
c           respect to parameters to be adjusted

c array dimensions
      include 'globdefs.inc'

c        commons
      include 'atmos.inc'
      include 'comdateq.inc'
      real*10 uhat(3), vhat(3)
      equivalence (Sbcom(4), uhat), (Sbcom(7), vhat)
      include 'coord.inc'
      real*10 tpvbus,spvbus,cpvbus
      equivalence (Angdum(5), tpvbus), (Angdum(6), spvbus),
     .            (Angdum(7), cpvbus)
      real*10 tpvmsb(2), tpvmsc(2)
      equivalence (Angdum, tpvmsb), (Angdum(3), tpvmsc)
      include 'difnct.inc'
      include 'emmips.inc'
      include 'empcnd.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'ltrap.inc'
      integer*2 kind
      equivalence (kind, Numpar)
      include 'mtrap.inc'
      include 'mnsprt.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,xfactr,freqtr(2),ct
      equivalence (ctrecf,Dstf), (fdev,Dstf(10))
      equivalence (Dstf(6),xfactr), (Dstf(7),freqtr), (Dstf(9),ct)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'partcm.inc'
      include 'pemctl.inc'
      include 'radcrd.inc'
      include 'sbdta.inc'
      include 'scdta.inc'
      include 'tapdta.inc'
      include 'trnocc.inc'
      include 'yvect.inc'

c external functions
      real*10 DOT

c local variables 
      real*10 atmdl1,atmdl2,dad11,dad12,dad21,dad22,
     . dlyrad,taufct,tfct,timlpt,tt,xvcor
      integer*4 i,ict66,iflag,iswtch,itfact,j,jdpvm,jj,jptest,
     . l,l1,lembry,lim,ll,lplnt,lprbcd,lrad,
     . lsprb,lsprc,lt,ltg0,ltype,m,m1,mm,mt,npat,nstart
      character*4 amper9/'&&&&'/
      logical*4 nosccn
      real*10 sitfct(6)/3*1._10, 3*8.64E4_10/
      real*10 tfact(2)
      integer*2 lplext,lscl
      real*10 deriv1(296)
      real*10 pvmfct(3)
      integer*2 izr2/0/
 
c
c initialize scale factors
      if(Nk1.le.0) then
         pvmfct(1) = 1._10/(Aultsc*Ltvel)
         pvmfct(2) = Convd/Aultsc
         pvmfct(3) = pvmfct(2)
      endif
c
c-----------------------------------------------------------------------
c*  start=1000
c
c          part.derivatives w.r.t. probe init.conds. and parameters
      if(Klans1.le.0) goto 300
      Ksprc = 1
      lsprc = 2
      m     = 7
      l     = 0
      call PCOPS(m, 'SC  ', Iabs1)
 
c*  start=1100
  100 iswtch = 0
  110 iflag  = 0
c
c search for part.deriv.w.r.t.probe in.conds.from probe tape
      call PCOPY(l, 6, iflag, iswtch, Lsc, Msc)
      if(iflag.le.0) then
         call PBDIC(Kisc, lsprc, Ksprc, l, 'SC', Klans1, iswtch)
         if(iswtch.ge.0) then
            call CPARTL(5, 2, kick)
            goto 200
         else
 
c not found on probe tape, copy old partial after all
            iswtch = -1
            goto 110
         endif
      endif
c
c search for partial w.r.t.probe parameter from probe tape
      if(l.le.6) then
         Ksprc = Lparsc
         lsprc = 8
      endif
      iswtch = 0
  130 iflag  = 1
      call PCOPY(l, 30, iflag, iswtch, Lsc, Msc)
      if(iflag.gt.0) goto 300
      lscl  = Lsc(l, 1)
      iflag = 1
      call PBDPRM(Nkisc, Kisc, lsprc, Ksprc, lscl, iflag)
      if(iflag.le.0) then
c
c*  start=1200
c calculate partial of observation
         call CPARTL(5, 2, kick)
         goto 200
c
c*  start=1300
c calc. partial w.r.t.probe parameter not on probe tape
      else if(lscl.gt.15) then
         if(lscl.le.19) then
            if(kick.ne.4 .or. nintrf.lt.0) goto 180
            tt = Tfrqsc
            if(lscl.eq.16) taufct = 1._10
            if(lscl.eq.17) taufct = tt
            if(lscl.eq.18) taufct = tt*tt/2._10
            if(lscl.eq.19) taufct = tt**3/6._10
            tfct = 0._10
            if(nddiff.lt.0) then
               if(lscl.eq.16) tfct = tt
               if(lscl.eq.17) tfct = tt*tt/2._10
               if(lscl.eq.18) tfct = tt**3/6._10
               if(lscl.eq.19) tfct = tt**4/24._10
            endif
            if(Nk1.ge.0) then
               Deriv(kind, 1) = -(taufct*Difdly(2,2) - tfct -
     .                          deriv1(kind))*xfactr
            else
               deriv1(kind) = taufct*Difdly(1, 2) - tfct
            endif
            goto 200
         else if(lscl.eq.24) then
c
c*  start=1700
c do relativity factor or print error message
            call PRTREL(lscl, Nkisc, Kisc, Lparsc, 5, kick, iswtch)
            if(iswtch.ge.0) goto 200
 
c can be copied from old obslib tape, do so
            iswtch = -1
            goto 130
         else
            if(kick.eq.2 .or. kick.eq.4) goto 180
            if(lscl.ne.23) goto 180
            if(kick.eq.3) then
 
c radius for mutual occultation
               if(Ibtrn.eq.6) then
                  if(Nice.le.0) Deriv(kind, 1)
     .                = (Re(1) - Rmxc(1))/Dfdt3(1)/Aukm
                  if(Nice.ge.0) Deriv(kind, 2)
     .                = (Re(2) - Rmxc(2))/Dfdt3(2)/Aukm
               else if(Ibtrn.lt.8) then
                  goto 180
               else
                  call PHSCRH(2)
                  call CPARTC(kick)
                  if(Nice.ge.0) Deriv(kind, 2) = Deriv(kind, 2)
     .                - Rmrx(2)/Rmrx(1)/Re(2)/Aukm
               endif
            else
               Deriv(kind, 1) = Aultsc
               Deriv(kind, 2) = 0._10
            endif
            goto 200
         endif
      else
         if(kick.eq.2 .or. kick.eq.3) goto 180
 
c see if pioneer-venus probe within atmosphere of planet
         jdpvm = Sccom(1)
         if(jdpvm.ne.2443852 .or. tpvmsc(1).eq.0._10) goto 180
         lprbcd = 1
         if(lscl.gt.5) lprbcd  = 2
         if(lscl.gt.10) lprbcd = 3
         l1     = lscl - 1
         itfact = mod(l1, 5)
         do i = 1, Mouse
            tfact(i) = tpvmsc(i)**itfact
            if(itfact.eq.2) tfact(i) = tfact(i)/2._10
            if(itfact.eq.3) tfact(i) = tfact(i)/6._10
            if(itfact.eq.4) tfact(i) = tfact(i)/24._10
c for differenced n-count observables the transmission times
c may be different for each site.  this difference is accounted
c for in tpvmsc but not in dxspcd, an adequate approximation.
            do j = 1, Index
               Derpr(j, i, 1) = 0._10
               Derpr(j, i, 2) = -Dxspcd(j, lprbcd, 2)
     .            *pvmfct(lprbcd)*tfact(i)
            end do
         end do
 
c might call partvl here
         Ivze(1) = 0
         call CPARTC(kick)
         if(Nice.ge.0 .and. lopler.ne.-1) call SUICID(
     .' CANNOT CALCULATE PIONEER-VENUS PROBE PARTIALS FOR VELOCITY OBSER
     .VABLE, STOP IN PARTL1  ', 22)
         goto 200
      endif
c
c set partials to zero and loop
  180 Deriv(kind, 1) = 0._10
      Deriv(kind, 2) = 0._10
c*  start=1900
  200 if(l.lt.30) goto 100
c
c-----------------------------------------------------------------------
c*  start=2000
c
c           partial derivatives w.r.t. planet initial conditions and
c           parameters
c          (includes central body of the satellite-probe)
  300 if(Nplnt0.le.0 .or. Nplnt0.eq.10) goto 1400
      nosccn = kick.ne.3 .or. Klans1.le.0 .or. Ncs1.ne.Nplnt(Klan)
      do i = 1, Index
         derem(i, 1) = 0._10
      end do
      Ivze(1) = 0
      if(Klan.le.0 .or. Klan.gt.u_mxpl) goto 1100
      jptest = 0
      Ksprb  = Lparsb
      lsprb  = 8
      Kplnt  = 1
      lplnt  = 2
      Kembry = Lparem
      lembry = 8
      ltg0   = 100*Nplnt0
      m = 7
      l = 0
      call PCOPS(m, 'PL  ', Iabs1)
c
c*  start=2100
c search for partial w.r.t. planet init.cond. from planet tape
  400 iflag = 0
      call PCOPY(l, 6, iflag, 1, Lpl, Mpl)
      if(iflag.gt.0) goto 410
      if(Kipl(1).ge.-1) then
 
c initial condition partial interpolated from tape
         call PBDIC(Kipl, lplnt, Kplnt, l, 'PL', Klan, 0)
         if(Klanb.le.0 .and. nosccn) goto 500
c
c planet is central body of observed body
         call CPARTL(3, 1, kick)
 
c look for added partial on satellite tapes
         lplext = 100*Ncp0 + l
         if( .not. nosccn) then
c in case of mutual occultation, must also check for partials on
c occulting body tape
            iflag = -1
            call PBDPRM(Nkisc, Kisc, lsprc, Ksprc, lplext, iflag)
            if(iflag.le.0) then
 
c compute partial and add central body term
               call CPARTL(5, 1, kick)
               do j = 1, 6
                  Dersc(j, 1) = Dersc(j, 1) + Derpl(j)
               end do
            else
 
c partial not found, assume zero
               do j = 1, 6
                  Dersc(j, 1) = Derpl(j)
               end do
               Ivze(5) = 1
            endif
         endif
 
         if(Klanb.gt.0) then
            iflag = -1
            call PBDPRM(Nkisb, Kisb, lsprb, Ksprb, lplext, iflag)
            if(iflag.le.0) then
 
c compute partials of satellite position
               call CPARTL(4, 1, kick)
               do j = 1, Index
                  Derpl(j) = Derpl(j) + Dersb(j)*Cmfct
               end do
            endif
         endif
      else
 
c initial condition partial to be gotten from elliptic orbit formula
         if(jptest.le.0) then
            jptest = 1
            timlpt = (Jd - Jdpl0(Klan))
            timlpt = timlpt + (ctrecf - Dstf(4)/Secday)
            call EMIPT(1, timlpt, 3)
         endif
         do j = 1, Index
            Derpl(j) = Dympt(j, l, 3)
         end do
         call PARTVL(Derpl, 1, kick)
         Ivze(3) = 1
      endif
c correct partial for indirect light-time effect, but only for
c observables transverse to the line of sight
      if(kick.eq.2) then
         xvcor=Aultvl*DOT(Xsitp0,Derpl)/(1._10-DOT(Xsitp0,Xp(4))*Aultvl)
         do j=1,3
            Derpl(j)=Derpl(j)+Xp(j+3)*xvcor
      end do
      endif
      if(kick.ne.3) then
         do i = 1, Mouse
            do j = 1, Index
               derem(j, i) = -Derpl(j)
               if(Nplnt2.eq.Nplnt(Klan)) Derpr(j, i, 2) = derem(j, i)
            end do
         end do
      endif
 
c ivze(1)=1
      call CPARTC(kick)
      goto 900
c
c search for partial w.r.t. planet parameter from planet tape
  410 if(l.le.6) then
         Kplnt = Lparp
         lplnt = 8
      endif
      iswtch = 0
  420 iflag  = 1
      call PCOPY(l, 30, iflag, iswtch, Lpl, Mpl)
      if(iflag.gt.0) goto 1000
      if(iswtch.eq.0) then
 
c partial not found, must compute
         iflag = 1
         call PBDPRM(Nkipl, Kipl, lplnt, Kplnt, Lpl(l), iflag)
         if(iflag.le.0) goto 500
 
c do not copy triaxial shape parameter  partials
      else if(Lpl(l).gt.5) then
         iswtch = -1
         goto 420
      endif
c
c
c*  start=2300
c calculate partial w.r.t. planet parameter not on planet tape
      if(kick.eq.2 .or. kick.eq.3) goto 700
      if(kick.ne.4) then
         if(Nspot.le.0) then
 
            if(Lpl(l).le.1) then
c it is assumed that there are no shape partials without also
c having partial derivative with respect to planet radius
               dlyrad = -2._10
               if(Spotf.eq.amper9) dlyrad =
     .             -(DOT(Xspcd(1,1),Xsitp0(1,1))
     .              + DOT(Xspcd(1,1),Xsitp0(1,2)))/Pradls
               if(npshp.gt.0) then
                  call SHPPR1(Lpl(l), dlyrad)
                  Deriv(kind, 2) = 0._10
               else
c partial derivative w.r.t. planet radius in kilometers
c multiplied by current value of au in kilometers
                  Deriv(kind, 1) = dlyrad*Aultsc
                  Deriv(kind, 2) = 0._10
               endif
               goto 900
c
c partial derivative w.r.t. planet shape parameters
            else if(Lpl(l).le.17) then
               if(Klanb.gt.0) goto 800
               call SHPPR1(Lpl(l), dlyrad)
               Deriv(kind, 2) = 0._10
               goto 900
            endif
         else if(Lpl(l).le.17) then
            ict66 = Ict(66)
            if(mod(ict66,2).eq.1 .and. Nplnt(Klan).eq.4) then
               call SPOTPR(Lpl(l), kick, Lpl, Mpl)
            else
               call PRTPR1(Lpl(l), kick)
            endif
            goto 900
         endif
      endif
c
c partial w.r.t. zenith delay and scale ht. for planetary atm.
      if(Lpl(l).ne.18) then
         if(Lpl(l).ne.19) goto 700
         npat = 1
      else
         npat = 2
      endif
      dad11 = Pdlpat(1, 1, npat)
      dad21 = Pdlpat(2, 1, npat)
      if(nddiff.gt.0) then
         dad12 = Pdlpat(1, 2, npat)
         dad22 = Pdlpat(2, 2, npat)
      endif
      if(kick.eq.2 .or. kick.eq.3) goto 880
      if(kick.eq.4) then
         if(nintrf.ge.0) then
            if(Nk1.lt.0) then
 
c store temporarily partial for beg. of first interval
               deriv1(kind) = dad11 - dad21
               if(nddiff.gt.0) Deriv(kind, 2) = dad12 - dad22
               goto 900
            else if(Nk1.eq.0) then
               deriv1(kind) = deriv1(kind)*freqtr(1)
               if(nddiff.gt.0) deriv1(kind) = deriv1(kind) -
     .                 Deriv(kind, 2)*freqtr(2)
            endif
            Deriv(kind, 1) = (dad11 - dad21)*freqtr(1)
            if(nddiff.gt.0) Deriv(kind, 1) = Deriv(kind, 1) -
     .                (dad12 - dad22)*freqtr(2)
            Deriv(kind, 1) = (Deriv(kind,1)-deriv1(kind))*xfactr
         else
            Deriv(kind, 1) = dad11 - dad21
            Deriv(kind, 2) = 0._10
            if(Nice.ge.0) call SUICID(
     .' PLANET ATM. PARAMETERS NOT CODED FOR DIFF. DELAY RATE, STOP IN P
     .ARTL1  ', 18)
 
c respective partials were zeroed if object not within atm of planet
            if(nddiff.gt.0) Deriv(kind, 1) = Deriv(kind, 1)
     .          - (dad21 - dad22)
         endif
      else
         Deriv(kind, 1) = dad11 + dad21
         Deriv(kind, 2) = 0._10
         if(Nice.ge.0) call SUICID(
     .'PLANET ATM. PARAMETERS NOT CODED FOR INSTANTANEOUS DOPPLER OBSERV
     .ABLE, STOP PARTL1  ', 21)
         if(Ncodf.gt.20) call SUICID(
     .' TWO-OBJECT RADAR-TYPE OBSERVABLE NOT CODED FOR PLANETARY ATM., S
     .TOP IN PARTL1  ', 20)
      endif
      goto 900
c
c*  start=2200
c calculate partial of coordinates
  500 call CPARTL(3, 1, kick)
c
c search for cross partial w.r.t. planet i.c. or parameter
c on embary tape
      if(Jct(14).gt.0) then
         iflag  = -1
         lplext = ltg0 + l
         call PBDPRM(Nkiem, Kiem, lembry, Kembry, lplext, iflag)
         if(iflag.le.0) then
            call CPARTL(1, 1, kick)
            goto 600
         endif
      endif
c
c no cross partial from planet tape
      do i = 1, Mouse
         do j = 1, Index
            derem(j, i) = 0._10
         end do
      end do
  600 if(kick.ne.3) then
         do i = 1, Mouse
            do j = 1, Index
               derem(j, i) = derem(j, i) - Derpl(j)
               if( Ncodf.gt.20 .and. (
     .           Nplnt2.eq.Nplnt0 .or.
     .           (Klanb.gt.0 .and. (Nplnt2.eq.Ncp0.or.Ncs1.eq.Ncp0)) )
     .           ) Derpr(j, i, 2) = derem(j, i)
            end do
         end do
      endif
c
c compute partial of observation
      call CPARTC(kick)
      goto 900
  700 if(Lpl(l).le.17) goto 880
c
c partial w.r.t. radius of planet for transit observation
      if(Lpl(l).ne.23) then
c
c*  start=2700
c do relativity factor if this is it
         iswtch = 1
         call PRTREL(Lpl(l), Nkipl, Kipl, Lparp, 3, kick, iswtch)
         if(iswtch.gt.0) goto 900
      else
         if(kick.ne.3) goto 880
         if(Nplnt(Klan).eq.Nps1) then
 
c radius for mutual occultation
            jj = 2
            if(Ibtrn.eq.6) then
               if(Nice.le.0) Deriv(kind, 1) = (Re(1) - Rmxc(1))
     .             /Dfdt3(1)/Aukm
               if(Nice.ge.0) Deriv(kind, 2) = (Re(2) - Rmxc(2))
     .             /Dfdt3(2)/Aukm
               goto 900
            endif
         else
            jj = 1
            if(Klanb.gt.0) goto 880
            if(Ibtrn.eq.2) then
 
c partial w.r.t. planet radius for stellar occultation
               if(Nice.le.0) Deriv(kind, 1) = -1._10/Dfdt3(1)/Aukm
               if(Nice.ge.0) Deriv(kind, 2) = -1._10/Dfdt3(2)/Aukm
               goto 900
            else if(Ibtrn.eq.4 .or. Ibtrn.eq.6) then
 
c transit or satellite occultation
               if(Nice.le.0) Deriv(kind, 1) = -Re(1)/Dfdt3(1)/Trnrd(1)
               if(Nice.ge.0) Deriv(kind, 2) = -Re(2)/Dfdt3(2)/Trnrd(2)
               goto 900
            endif
         endif
         if(Ibtrn.lt.8) goto 880
         call PHSCRH(jj)
         call CPARTC(kick)
         if(Nice.ge.0) Deriv(kind, 2) = Deriv(kind, 2)
     .       - Rmrx(2)/Rmrx(1)/Re(jj)/Aukm
         goto 900
      endif
c
c see if partial which does not affect central planet motion affects
c orbiter motion
  800 if(Klanb.le.0) goto 870
 
c assume klans1=0
c
c search for partial w.r.t. planet parameter from probe tape
c lplext=+lpl(l) + ncp0 + 6       out
      lplext = +Lpl(l) + 100*Ncp0 + 6
      iflag  = -1
      call PBDPRM(Nkisb, Kisb, lsprb, Ksprb, lplext, iflag)
      if(iflag.gt.0) goto 870
      call CPARTL(4, 2, kick)
      if(Lpl(l).eq.24) call SUICID(
     .'RELATIVITY FACTOR OF PLANET AFFECTS PROBE ONLY, STOP IN PARTL1  '
     ., 16)
      goto 900
 
c print error message
  870 call PRTREL(izr2, izr2, izr2, 0, 3, 0, 0)
 
c set partials to zero
  880 Deriv(kind, 1) = 0._10
      Deriv(kind, 2) = 0._10
c
c*  start=2900
  900 if(l.lt.30) goto 400
c
c-----------------------------------------------------------------------
c
c partials wrt planet shape harmonic coefficients
 1000 call SHPPAR(kick, dlyrad)
c
c           partial derivatives w.r.t. planet gravitational potential
c           harmonic coefficients for observations of probe with
c           planet as central body
c
      if(Klanb.le.0) goto 1400
      Ksprb  = Lparsb
      nstart = 8
 
c zonal harmonics
      call HPARTL(kick, Nczone, Lczhar, Mczhar, Ncp0, 31, nstart, Iabs1)
 
c tesseral cosine harmonics
      call HPARTM(Nctess, Lcchar, Mcchar, 41)
 
c tesseral sine harmonics
      call HPARTM(Nctess, Lcshar, Mcshar, 51)
c
c-----------------------------------------------------------------------
c*  start=3000
c
c          part.derivatives w.r.t. probe init.conds. and parameters
 1100 if(Klanb.le.0) goto 1400
      jptest=0
      Ksprc = Lparsc
      lsprc = 8
      Ksprb = 1
      lsprb = 2
      m     = 7
      l     = 0
      call PCOPS(m, 'SB  ', Iabs1)
 
c*  start=3100
 1200 iswtch = 0
 1210 iflag  = 0
c
c search for part.deriv.w.r.t.probe init. conds from probe tape
      call PCOPY(l, 6, iflag, iswtch, Lsb, Msb)
      if(iflag.gt.0) goto 1220
      if(Kisb(1).lt.-1) then
c get partial from elliptic formula
         if(jptest.eq.0) then
            jptest=1
            timlpt=Jd-Jdpl0(Klanb)
            timlpt=timlpt+(ctrecf-Dstf(4)/Secday)
            call EMIPT(1,timlpt,4)
         end if
         do j=1,index
           Dersb(j)=Dympt(j,l,4)
        end do
      else
         call PBDIC(Kisb, lsprb, Ksprb, l, 'SB', Klanb, iswtch)
         if(iswtch.lt.0) goto 1210
      end if
 
      if(Klans1.gt.0 .and. kick.eq.3) then
 
c determine if occulted body is affected by observed body
         lplext = Nplnt(Klanb)*100 + l
         iflag  = -1
         call PBDPRM(Nkisc, Kisc, lsprc, Ksprc, lplext, iflag)
         if(iflag.le.0) call CPARTL(5, 1, kick)
      endif
      call CPARTL(4, 2, kick)
      goto 1300
c
c search for partial w.r.t.probe parameter from probe tape
 1220 if(l.le.6) then
         Ksprb = Lparsb
         lsprb = 8
      endif
      iswtch = 0
 1230 iflag  = 1
      call PCOPY(l, 30, iflag, iswtch, Lsb, Msb)
      if(iflag.gt.0) goto 1400
      iflag = 1
      call PBDPRM(Nkisb, Kisb, lsprb, Ksprb, Lsb(l), iflag)
      if(iflag.le.0) then
c
c*  start=3200
c calculate partial of observation
         call CPARTL(4, 2, kick)
         goto 1300
      endif
c
c*  start=3300
c calc. partial w.r.t.probe parameter not on probe tape
      if(Lsb(l).le.15) then
         if(kick.eq.2 .or. kick.eq.3) goto 1280
         jdpvm = Sbcom(1)
 
c see if pioneer-venus bus
         if(jdpvm.ne.-2443852) then
 
c see if pioneer-venus probe within atmosphere of planet
            if(jdpvm.ne.2443852) goto 1270
            if(tpvmsb(1).eq.0._10) goto 1280
            lprbcd = 1
            if(Lsb(l).gt.5) lprbcd  = 2
            if(Lsb(l).gt.10) lprbcd = 3
            l1     = Lsb(l) - 1
            itfact = mod(l1, 5)
            do i = 1, Mouse
               tfact(i) = tpvmsb(i)**itfact
               if(itfact.eq.2) tfact(i) = tfact(i)/2._10
               if(itfact.eq.3) tfact(i) = tfact(i)/6._10
               if(itfact.eq.4) tfact(i) = tfact(i)/24._10
c for differenced n-count observables, the transmission times
c may be different for each site.  this difference is accounted
c for in tpvmsb but not in dxspcd, an adequate approximation.
               do j = 1, Index
                  Derpr(j, i, 1) = -Dxspcd(j, lprbcd, 1)
     .               *pvmfct(lprbcd)*tfact(i)
               end do
            end do
         else if(kick.ne.1) then
            goto 1280
         else if(Lsb(l).ne.1) then
            if(Lsb(l).ne.2) then
               if(Lsb(l).ne.3) goto 1270
               do i = 1, Mouse
                  do j = 1, Index
                     Derpr(j, i, 1) = -Pcond(8, Klanb)
     .                  *(cpvbus*uhat(j) - spvbus*vhat(j))/Aultsc
                  end do
               end do
            else
               do i = 1, Mouse
                  do j = 1, Index
                     Derpr(j, i, 1) =
     .                  -(spvbus*uhat(j) + cpvbus*vhat(j))/Aultsc
                  end do
               end do
            endif
         else
            do i = 1, Mouse
               do j = 1, Index
                  Derpr(j, i, 1) = Pcond(8, Klanb)*tpvbus*
     .              (cpvbus*uhat(j) - spvbus*vhat(j))/Aultsc
               end do
            end do
         endif
         do i = 1, Mouse
            do j = 1, Index
               Derpr(j, i, 2) = 0._10
            end do
         end do
         call PARTVL(Derpr, 2, kick)
         Ivze(1) = 1
         if(Nice.ge.0 .and. lopler.ne.-1) call SUICID(
     .' CANNOT CALCULATE PIONEER-VENUS PROBE PARTIALS FOR VELOCITY OBSER
     .VABLE, STOP IN PARTL1  ', 22)
         call CPARTC(kick)
         goto 1300
 
c partial w.r.t. xmtr freq. for counted-cycle vlbi observable
      else if(Lsb(l).le.19) then
         if(kick.ne.4 .or. nintrf.lt.0) goto 1280
         tt = Tfrqsb
         if(Lsb(l).eq.16) taufct = 1._10
         if(Lsb(l).eq.17) taufct = tt
         if(Lsb(l).eq.18) taufct = tt*tt/2._10
         if(Lsb(l).eq.19) taufct = tt**3/6._10
         tfct = 0._10
         if(nddiff.lt.0) then
            if(Lsb(l).eq.16) tfct = tt
            if(Lsb(l).eq.17) tfct = tt*tt/2._10
            if(Lsb(l).eq.18) tfct = tt**3/6._10
            if(Lsb(l).eq.19) tfct = tt**4/24._10
         endif
         if(Nk1.ge.0) then
            Deriv(kind, 1) = (taufct*Difdly(2,1) - tfct -
     .                       deriv1(kind))*xfactr
         else
            deriv1(kind) = taufct*Difdly(1, 1) - tfct
         endif
         goto 1300
      else if(Lsb(l).ne.24) then
         if(Lsb(l).eq.23) then
            if(kick.eq.2 .or. kick.eq.4) goto 1280
            if(kick.ne.3) then
               Deriv(kind, 1) = Aultsc
               Deriv(kind, 2) = 0._10
            else if(Ibtrn.ne.2) then
 
c partial w.r.t. satellite radius for mutual occultation
               if(Ibtrn.eq.6) then
                  if(Nice.le.0) Deriv(kind, 1) = -Re(1)
     .                /Dfdt3(1)/Trnrd(1)
                  if(Nice.ge.0) Deriv(kind, 2) = -Re(2)
     .                /Dfdt3(2)/Trnrd(2)
               else if(Ibtrn.lt.8) then
                  goto 1280
               else
                  call PHSCRH(1)
                  call CPARTC(kick)
                  if(Nice.ge.0) Deriv(kind, 2) =
     .                Deriv(kind, 2) - Rmrx(2)/Rmrx(1)/Re(1)/Aukm
               endif
            else
 
c partial w.r.t. satellite radius for stellar occultation
               if(Nice.le.0) Deriv(kind, 1) = -1._10/Dfdt3(1)/Aukm
               if(Nice.ge.0) Deriv(kind, 2) = -1._10/Dfdt3(2)/Aukm
            endif
            goto 1300
         endif
      endif
c
c*  start=3700
c do relativity factor if this is it
 1270 call PRTREL(Lsb(l), Nkisb, Kisb, Lparsb, 4, kick, iswtch)
      if(iswtch.ge.0) goto 1300
 
c can be copied from old obslib tape, do so
      iswtch = -1
      goto 1230
 
c set partials to zero and loop
 1280 Deriv(kind, 1) = 0._10
      Deriv(kind, 2) = 0._10
c
c*  start=3900
 1300 if(l.lt.30) goto 1200
c
c-----------------------------------------------------------------------
c*  start=4000
c
c           partial derivatives w.r.t. spot coordinates on observed body
 1400 m  = 1
      m1 = 2
      if(Nspot.gt.0) goto 1600
 1500 if(Nspot2.le.0) goto 1700
      m  = 2
      m1 = 1
 1600 l  = 0
      call PCOPS(m, 'SPCR', Iabs1)
      do while( .true. )
         iflag = 0
         call PCOPY(l, 3, iflag, 1, Lspcd(1,m), Mspcd(1,m))
         if(iflag.le.0) then
            if(Ncodf.eq.18) then
 
c partial of pulsar phase observable
               Deriv(kind, 1) = Dptdsp*DOT(Dxspcd(1,l,m), Xessbc)
               Deriv(kind, 2) = 0._10
            else if(Ncodf.eq.6 .and. Nplnt0.lt.0) then
               if(l.eq.2) then
                  Deriv(kind,1)=240._10*Aultsc/Convd
                  Deriv(kind,2)=0._10
               else if(l.eq.3) then
                  Deriv(kind,1)=0._10
                  Deriv(kind,2)=3600._10*Aultsc/Convd
               endif
            else
               do i = 1, Mouse
                  do j = 1, Index
                     Derpr(j, i, m)  = -Dxspcd(j, l, m)*sitfct(j)
                     Derpr(j, i, m1) = 0._10
                  end do
               end do
               Ivze(1) = 2 - m
 
c it is assummed that no spot is involved for transit or occultation
               Lprspt = 1
               call CPARTC(kick)
            endif
         endif
         if(l.ge.3) then
            if(m.le.1) goto 1500
            goto 1700
         endif
      end do
c
c*  start=4500
 1700 if(kick.eq.1) then
      else if(kick.eq.3) then
         goto 1900
      else if(kick.ne.4) then
c     it is assummed that no observation biases are involved for
c     transit or occulation observation
c
c              partial derivatives w.r.t. phase corrections for
c              planetary optical observations
c     we must increment counters even if nphase=0
         if(Ncph.gt.0) then
            i = 0
            call PCOPS(m, 'PHS ', Iabs1)
            lim = Ncph
            do while( .true. )
               iflag = 0
               call PCOPY(i, lim, iflag, 1, Lphs, Mphs)
               if(iflag.gt.0) goto 1750
               if(Nice.le.0) Deriv(kind, 1) = Phal*Cthet(i)
               if(Nice.ge.0) Deriv(kind, 2) = Phde*Cthet(i)
               if(i.ge.Ncph) goto 1750
            end do
         endif
 1750    if(Nrbias.le.0) goto 1900
         call SUICID(
     .' NO RADAR BIASES ALLOWED FOR OPTICAL OBS, STOP IN PARTL1',14)
      endif
c*  start=5000
c
c partial derivatives w.r.t. radar observation biases
c we must increment counters even if nrbias=0
      do i = 1, 2
         if(Lrbs(i).gt.0) then
            kind = kind + 1
            if(Iabs1.gt.0) then
               if(Mrbs(i).gt.0) Mind = Mind + 1
            endif
 
c quasi-vlbi bias is initial sampling error
            Deriv(kind, i) = 1._10
            if(Ncodf.eq.18) Deriv(kind, 1) = Dpphdt
            Deriv(kind, 3-i) = 0._10
         endif
      end do
 
c*  start=5500
      if(Ncph.gt.0) then
c
c partial derivatives w.r.t. atmospheric zenith correction
c for radio observables
         lrad = 3
         do i = 1, Ncph
            if(i.eq.3) lrad = 6
            if(Lphs(i).le.0) goto 1850
            kind = kind + 1
            if(Iabs1.gt.0) then
               if(Mphs(i).gt.0) Mind = Mind + 1
            endif
            atmdl1 = Raddum(lrad + i)
            atmdl2 = Radum2(lrad + i)
            if(nintrf.lt.0) then
               Deriv(kind, 1) = atmdl1
               Deriv(kind, 2) = atmdl2
               if(Ncodf.le.20) goto 1850
               if(Nice.lt.0) goto 1760
               if(Nice.eq.0) then
                  call SUICID(
     .'PARTIAL DERIVATIVES W.R.T. ATM PARAMETER CANNOT BE CALCULATED FOR
     . BOTH DELAY AND RATE OBSERVABLES,STOP IN PARTL1', 28)
               else
                  Deriv(kind, 2) = Deriv(kind, 2) - atmdl2
                  goto 1850
               endif
            endif
            if(Nk1.lt.0) then
 
c store temporarily atmdl1 for beginning of first interval
               deriv1(kind) = atmdl1
               if(Ncodf.le.20) goto 1850
            else if(Nk1.eq.0) then
               deriv1(kind) = deriv1(kind)*freqtr(1)
               goto 1780
            else
               goto 1780
            endif
 1760       if(nintrf.ge.0) goto 1800
            Deriv(kind, 1) = Deriv(kind, 1) - atmdl2
            goto 1850
 1780       Deriv(kind, 1) = atmdl1*freqtr(1)
            if(Ncodf.le.20) then
               Deriv(kind, 1) = Deriv(kind, 1) - deriv1(kind)
               goto 1850
            endif
 1800       if(Nk1.lt.0) then
 
c store temporarily atmdl2 for beginning of first interval
               Deriv(kind, 2) = atmdl2
               goto 1850
            else if(Nk1.eq.0) then
               deriv1(kind) = deriv1(kind) - Deriv(kind, 2)*freqtr(2)
            endif
            Deriv(kind, 1) = Deriv(kind, 1) - atmdl2*freqtr(2)
     .                       - deriv1(kind)
c
c partial derivatives w.r.t. interferometer clock corrections
 1850    end do
      endif
c
c-----------------------------------------------------------------------
c*  start=6000
c
c partial derivatives w.r.t. target body initial conditions,
c parameters, gravitational potential harmonic coefficients
c
 1900 if(Numtar.le.0) goto 2500
c
c for planet observations
      if(Klanb.gt.0) then
c
c for probe observations
         Ksprb = Lparsb
         lt    = 8
         ll    = 1
         mm    = 1
      else
         if(Mnplnt.ne.0) call TRGPIC(kick)
         goto 2500
      endif
 2000 if(ll.gt.Numtar) goto 2500
c
c search for start of target planet partials on probe tape
      ltg0   = Ntrg(ll)*100
      lplext = ltg0
      iflag  = -1
      call PBDPRM(Nkisb, Kisb, lt, Ksprb, lplext, iflag)
c
c see if target body is on observation library tape
      mt = 0
      if(Iabs1.gt.0 .and. mm.le.Mumtar) then
         if(Mtrg(mm).eq.Ntrg(ll)) mt = 1
      endif
c
c start loop for target body init.cond. & parameter partials
      m = 7
      l = 0
      call PCOPS(m, 'TBOD', mt)
 
c*  start=6100
 2100 if(l.ge.6) goto 2300
 
c target body initial condition partial
      iflag = 0
      call PCOPY(l, 6, iflag, 1, Ltbod(1,ll), Mtbod(1,mm))
      if(iflag.gt.0) goto 2300
      ltype = l
c
c search probe tape for partial
 2200 lplext = ltg0 + ltype
      iflag  = -1
      call PBDPRM(Nkisb, Kisb, lt, Ksprb, lplext, iflag)
      if(iflag.gt.0) call SUICID(
     .'TARGET BODY PARTIAL NOT ON PROBE TAPE, STOP IN PARTL1   ', 14)
      call CPARTL(4, 2, kick)
c
c*  start=6900
c end of loop for target body init.cond. & parameter partials
      if(l.ge.30) goto 2400
      goto 2100
c
c target body parameter partial
 2300 iflag = 1
      call PCOPY(l, 30, iflag, 1, Ltbod(1,ll), Mtbod(1,mm))
      if(iflag.le.0) then
         ltype = Ltbod(l, ll) + 6
         goto 2200
      endif
c
c partial derivatives w.r.t. target body gravitational
c potential harmonic coefficients
c zonal harmonics
 2400 call HPARTL(kick, Ntzone(ll), Ltzhar(1,ll), Mtzhar(1,mm), Ntrg(ll)
     .            , 31, lt, mt)
 
c tesseral cosine harmonics
      call HPARTM(Nttess(ll), Ltchar(1,ll), Mtchar(1,mm), 41)
 
c tesseral sine harmonics
      call HPARTM(Nttess(ll), Ltshar(1,ll), Mtshar(1,mm), 51)
c
c end of ll loop on target bodies
      if(mt.gt.0) mm = mm + 1
      ll = ll + 1
      goto 2000
c
c           partial derivatives w.r.t. pulsar parameters
 2500 call PSRPAR(kick)

      end
