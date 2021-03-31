      subroutine RADEXT(nvel,norm)
 
      implicit none

c          r. goldstein     april, 1976
c  calculate differential spot delay-doppler (relative to sub-radar
c  point)
c
c     oct. 1976...the option to do bandwidth (nspot=0, nice=1) is
c                 believed to be inoperative in the old code. it has
c                 therefore been deleted from this cleaned up version.
c
c     this subroutine is called by radctl
c
c arguments
      integer*4 norm,nvel

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 dutrec
      equivalence (dutrec,Angdum(10))
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'ltrapobs.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'mnsprt.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      include 'kobequiv.inc'
      real*10 ctrecf,fdev,reflct,ylong,ylat,zlong,zlat,ztopo
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (Dstf(2),ylong),(Dstf(3),ylat)
      equivalence (Dstf(5),zlong),(Dstf(6),zlat),(Dstf(8),ztopo)
      real*4    acctim
      equivalence (acctim,Estf)
      include 'param.inc'
      include 'plnhar.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'
 
c local variables
      real*10 dopsr,dum,dum1,frcns,sendtm,tdsr,xsr1m,xsr2m
      integer   i,j,lemctl,loop1,loop2,lplctl,mnspt1,nn
      real*10 xsr1(3),xsr2(3),sr1hat(3),sr2hat(3),srh(3,2),
     .          arpl(3)
      real*10 xspth(3,2)
      equivalence (sr1hat,srh(1,1)),(sr2hat,srh(1,2))
c
c     xslcns(1-6,j) = position and velocity (not used) of the center
c                     of mass of the solar system relative to the sun
c                     at receive time (j=1), reflection time (j=2)
c                     and send time (j=3)
c
c       receive time site and earth coods. determined in radctl
c
c     tests on required flags
c
      if(Nspot.le.0) call SUICID(' NSPOT NO GOOD IN RADEXT',6)
      if(Klan.le.0) call SUICID(' KLAN  NO GOOD IN RADEXT',6)
      if(Klanb.ne.0) call SUICID(' KLANB NO GOOD IN RADEXT',6)
c
c------------------setup for interpolation (reflection)-----------------
      mnspt1 = 0
      lplctl = -1
      lemctl = 0
      loop1  = 0
      loop2  = 0
      do while(.true.)
c
c*  start=1100
c reflection time iteration loop
c tmdly1, the first guess, was generated in radctl
         call TIMINC(Jd, ctrecf, Jdx, Fract, -Tmdly1/Secday)
 
         call SPOTCD(Jdx, Fract, mnspt1, nvel, Nplnt0, dum1, dum1, 1)
c a planet is involved so obtain refl. pt. wrt sun
c
         call PLTRP(1, Jdx, Fract, 0, lplctl)
         if(Jdx.le.0) then
c
c*  start=9910
            Jd = 0
            return
         endif
         do j = 1, 3
            Xplsc(j)  = Xp(j)*Aultsc
            Xsbsun(j) = Xspcd(j, 1) + Xplsc(j)
            end do
c
c*  start=1400
         do j = 1, 3
            Xsitep(j, 1) = Xemlsc(j, 1) - Xsbsun(j)
            end do
         if(Nswcns.gt.0) then
            call SOTRP(Jdx, Fract, Xslcns(1,2), 0)
            do j = 1, 3
               Xsitep(j, 1) = Xsitep(j, 1)
     .                        + (Xslcns(j,2) - Xslcns(j,1))*Aultsc
               end do
         endif
c
c------------------calculate delay time for reflection------------------
         call UVECTR(3, Xsitep, Rsitp, Xsitp0, dum)
         Tmdly2 = Rsitp(1)
         loop1  = loop1 + 1
         if(loop1.gt.20)
     .        call SUICID('MORE THAN 20 RECEIVE ITERATIONS ', 8)
         if(ABS(Tmdly2-Tmdly1).le.acctim) goto 100
         Tmdly1 = Tmdly2
         end do
c
c iteration to determine sending site position
  100 Nit(20) = Nit(20) + loop1
      Tguess  = Tmdly1
      Dstf(4) = Tmdly1
      call TSAV(Tmdly2, 1, Jd, Utrec)
      call DGUESS(tmdly, 3, Jd, Utrec)
      tmdly = tmdly + Tmdly2
      do while(.true.)
         Tmdly1 = tmdly
         sendtm = Utrec - Tmdly1*fdev
         Sidtm  = Sidtm0 + Sidvel*sendtm
c in computing sidereal time, tmdly1 converted to wwv seconds as
c approximation to universal time (if tmdly1 left in coordinate
c time,sending site position would be off by 3 centimeters)
         call SITCOR(Sidtm, 2, -1, 0)
         call TIMINC(Jd, ctrecf, Jdy, frcns, -Tmdly1/Secday)
         call ETRP(1, Jdy, frcns, 0, lemctl, 2, 3)
c
c------------------determine  sending  site wrt refl. pt. at  send  time
         do j = 1, 3
            Xsitep(j, 2) = Xemlsc(j, 2) - Xsbsun(j)
            end do
         if(Nswcns.gt.0) then
            call SOTRP(Jdy, frcns, Xslcns(1,3), 0)
            do j = 1, 3
               Xsitep(j, 2) = Xsitep(j, 2) + (Xslcns(j,2) - Xslcns(j,3))
     .                        *Aultsc
               end do
         endif
c
c calculate delay time for send
         call UVECTR(3, Xsitep(1,2), Rsitp(2), Xsitp0(1,2), dum)
         tmdly = Tmdly2 + Rsitp(2)
         loop2 = loop2 + 1
         if(loop2.gt.20)
     .       call SUICID('MORE THAN 20 SEND ITERATIONS', 7)
         if(ABS(tmdly-Tmdly1).le.acctim) goto 2000
         end do
 2000 Nit(19) = Nit(19) + loop2
      call TSAV(tmdly - Tmdly2, 3, Jd, Utrec)
c
c - - - - - - - - - - end of sending iteration  - - - - - - - - - - - -
c
c test for dummy observation below horizon
      if(Idumob.eq.1) then
         nn = -2
         call HORIZN(Xsite, Xsitep, nn, Jd)
 
c test for spot on planet turned away from earth
         call CNTOCC(Xspcd, Xsitep, -2, Jd)
         if(Jd.le.0) return
      endif
c*  start=2500
c
c       general relativity correction to time delay not supported here
c
c     differenced range to spot on planet from sub-radar point
      do j = 1, 3
         xsr1(j) = Xemlsc(j, 1) - Xplsc(j)
         xsr2(j) = Xemlsc(j, 2) - Xplsc(j)
         if(Nswcns.gt.0) then
            xsr1(j) = xsr1(j) + Aultsc*(Xslcns(j,2) - Xslcns(j,1))
            xsr2(j) = xsr2(j) + Aultsc*(Xslcns(j,2) - Xslcns(j,3))
         endif
         end do
      call UVECTR(3, xsr1, xsr1m, sr1hat, dum)
      call UVECTR(3, xsr2, xsr2m, sr2hat, dum)
      tdsr = xsr1m + xsr2m - 2._10*Pradls
 
c calculate sub-radar latitude and longitude
      if(nplng.gt.0) then
         call PRODCT(Rot, sr1hat, arpl, 3, 3, 1)
         ylong = ATAN2(arpl(2), arpl(1))/Convd + Pcom(5)
         if(ylong.lt.0._10) ylong = ylong + 360._10
         ylat = ASIN(arpl(3))/Convd
 
c delay correction due to planet shape (subradar point)
         if(npshp.ge.0) then
 
c branch depending on shape model used (see bodred)
            if(Nshape(Npshp1).eq.0 .or. Nshape(Npshp1).eq.1) then
            else if(Nshape(Npshp1).eq.2) then
 
c altitude grid model (local shape model) nshape=2
               call GRDSHP(ylat, ylong, tdsr, ztopo)
               go to 110
            else
               call SUICID(
     .           'NSHAPE OUTSIDE ALLOWABLE RANGE. STOP IN RADEXT  ', 12)
            endif
 
c spherical harmonic model or fourier series (nshape=0 or 1)
            call HARSHP(ylat, ylong, tdsr, ztopo, 1, 0._10, 1)
         endif
  110    if(nplng.gt.1) then
            if(Line.gt.55) call OBSPAG
            write(Iout, 120) ylong, ylat, ztopo, srh
  120       format(' YLONG=', f17.12, '  YLAT=', f17.12,
     .             '  ZTOPO=', 1pd15.8,
     .             '  XSITP0 (UNIT VECTORS SUBRADAR PT TO SITE):'/
     .             6D22.13)
            Line = Line + 2
         endif
      endif
c
c*  start=2700
c form differential delay
      tmdly = tmdly - tdsr
c
c velocity calculations
c
      call SITVEL(2, nvel, norm)
      call SPOTCV(nvel, Nplnt0, 1)
 
c velocity of site
      if(nvel.gt.0) call ETRP(1,Jdy,frcns,1,0,2,3)
c
c------------------correct time delay for various effects---------------
c*  start=3000
c
c       skip corrections if only doppler
      if(Nice.le.0) then
c
c
c       all media corrections for ncodf=3 are not done
c
c     time delay in atomic time changed to time delay in ut time
         if(ntime.gt.0) tmdly = tmdly*fdev
 
c constant bias in time delay
         if(Nrbias.gt.0) tmdly = tmdly + Rbsx(1)
 
c power series for clock errors
         if(Neqnox.gt.0) tmdly = tmdly + Pnox +
     .       dutrec*(Pequat + dutrec*Plat/2._10)
c
c------------------doppler---------------------------------------------
c*  start=4000
c nice.lt.0  no doppler;            nice.gt.0  no time delay
         if(Nice.lt.0 .and. calcvl.le.0) return
      endif
c
c
c       no phase delay doppler for the ncodf=3 observable.
c
c
c            velocity of refl. pt.
c
      call PLTRP(1, Jdx, Fract, 1, lplctl)
      if(Jdx.le.0) then
         Jd = 0
         return
      endif
      do j = 4, 6
         Xplsc(j)  = Xp(j)*Aultvl
         Xsbsun(j) = Xspcd(j, 1) + Xplsc(j)
         end do
c
c
c relative velocity
c 1st get doppler to sub-radar point
      do i = 1, 6
         xspth(i, 1)  = Xsitp0(i, 1)
         Xsitp0(i, 1) = srh(i, 1)
         end do
      call VLRTRD(Xemlsc(4,1), Xplsc(4), Xemlsc(4,2), 1._10, 7, 0)
      if(Nice.ge.0) then
         call RADREL(2)
         dopsr = dop
      endif
 
c 2nd get doppler to spot
      do i = 1, 6
         Xsitp0(i, 1) = xspth(i, 1)
         end do
      call VLRTRD(Xemlsc(4,1), Xsbsun(4), Xemlsc(4,2), 1._10, 7, 0)
      if(Nice.ge.0) then
         call RADREL(2)
 
c now get differenced doppler
         dop = dop - dopsr
c
c------------------correct doppler shift for various effects------------
c
c       media corrections for ncodf=3 is not supported
c
         dop = Freq*dop
c
c constant bias in doppler shift
         if(Nrbias.gt.0) dop = dop + Rbsx(2)
 
c power series for clock errors
         if(Neqnox.gt.0) dop = dop - (Pequat + dutrec*Plat)*Freq
      endif
 
c*  start=9990
      return
      end
