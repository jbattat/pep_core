      subroutine RADPL(nvel, norm)
 
      implicit none

c          r. goldstein, april,1976
c          interpolations are now done by subroutines mntrp,emtrp,
c          pltrp, and sbtrp
c          polynomial extrapolations are done by dguess
c
c       sept. 1976....derivative of deldop (dldpbg)
c         this routine handles radar from a planet
c           ncodf.eq.1  klanb.eq.0    nspot.ge.0
c arguments
      integer norm,nvel

c array dimensions
      include 'globdefs.inc'

c        common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 dutrec
      equivalence (Angdum(10),dutrec)
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'ltrapobs.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'mnsprt.inc'
      real*10 radoff
      equivalence (Rspot,radoff)
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct,ylong,ylat,zlong,zlat,ztopo,tmdly0,dop0
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (Dstf(2),ylong),(Dstf(3),ylat)
      equivalence (Dstf(5),zlong),(Dstf(6),zlat),(Dstf(8),ztopo)
c        dstf(1) =ctrecf=coordinate time at reception (fraction of day)
c        dstf(2) = longitude of subradar point (deg)
c        dstf(3) = latitude of subradar point (deg)
c        dstf(4) = time delay from receive to reflection (sec)
c        dstf(5) = longitude of observed spot away from subradar point
c                  if spot name='&&&&' and nspot=0 (deg)
c        dstf(6) = latitude of observed spot away from subradar point
c                  if spot name='&&&&' and nspot=0 (deg)
c        dstf(8) = round-trip topography at observed spot (light sec)
c        dstf(10)= fdev =1+fractional frequency offset from atomic
c                        time for unit of time delay measurement
      equivalence (Result,tmdly0),(Result(2),dop0)
      real*4 acctim
      equivalence (acctim,Estf)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'plnhar.inc'
      include 'prpgat.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'

c local 
      character*4 amper9/'&&&&'/,pound9/'####'/
      logical   srpt
      real*10 arpl(3),oldtopo
      real*10 dum,frsend,radofe,sendtm
      integer   i,itcnt,j,lplctl,mnspt1,nvsp

c external functions
      real*10 ANCOR,CTATF
c
c tests on required flags
      if(Klan.le.0) call SUICID('KLAN NO GOOD IN RADPL   ',6)
      srpt = Spotf.ne.pound9 .and. Spotf.ne.amper9 .and. Nspot.le.0
      nvsp = 0
      if(Spotf.eq.pound9) nvsp = -1
 
      mnspt1 = 0
      lplctl = -1
c
c*  start=1100
c
c iteration to determine planet position for reflection
      itcnt = 0
 1100 call TIMINC(Jd,ctrecf,Jdx,Fract,-Tmdly1/Secday)
      call PLTRP(1,Jdx,Fract,0,lplctl)
      if(Jdx.le.0) goto 9910
      do j = 1,3
         Xplsc(j)     = Xp(j)*Aultsc
         Xsitep(j,1) = Xemlsc(j,1) - Xplsc(j)
      end do
      if(.not.srpt) then
         call SPOTCD(Jdx,Fract,mnspt1,nvsp,Nplnt0,0._10,0._10,1)
         if(nvsp.lt.0) then
 
c spot is defined at offset in longitude from sub-radar point
            call UVECTR(3,Xsitep,Rsitp,Xsitp0,dum)
 
c transform line-of-sight vector to planetocentric coordinates
            call PRODCT(Rot,Xsitp0,arpl,3,3,1)
 
c rotate vector by longitude offset
            Yspcd(1,1) = arpl(1)*Csroff - arpl(2)*Ssroff
            Yspcd(2,1) = arpl(1)*Ssroff + arpl(2)*Csroff
            Yspcd(3,1) = arpl(3)*radoff
            call PRODCT(Rot,Yspcd,Xspcd,-3,3,1)
         endif
         do j = 1,3
            Xplsc(j)    = Xplsc(j) + Xspcd(j,1)
            Xsitep(j,1) = Xemlsc(j,1) - Xplsc(j)
         end do
      endif
c
      call UVECTR(3,Xsitep(1,1),Rsitp(1),Xsitp0(1,1),dum)
      Tmdly2 = Rsitp(1)
      if(srpt) Tmdly2 = Tmdly2 - Pradls
      Nit(20) = Nit(20) + 1
      if(ABS(Tmdly2-Tmdly1).gt.acctim) then
         itcnt = itcnt + 1
         if(itcnt.gt.100) call SUICID('TOO MANY REFLECT ITERATIONS ',7)
         Tmdly1 = Tmdly2
         goto 1100
      endif
 
c
c iteration to determine sending site position
c*  start=2000
      call TSAV(Tmdly2,1,Jd,Utrec)
      Tguess = Tmdly2
      call DGUESS(tmdly,3,Jd,Utrec)
      tmdly  = tmdly + Tmdly2
      reflct = Tmdly2
      itcnt = 0
 2000 Tmdly1 = tmdly
      sendtm = Utrec - Tmdly1*fdev
      Sidtm  = Sidtm0 + Sidvel*sendtm
      call SITCOR(Sidtm,2,-1,0)
      call TIMINC(Jd,ctrecf,Jdy,frsend,-Tmdly1/Secday)
      call ETRP(1,Jdy,frsend,0,0,2,3)
      if(Jdy.le.0) goto 9910
      do j = 1,3
         Xsitep(j,2) = Xemlsc(j,2) - Xplsc(j)
      end do
      call UVECTR(3,Xsitep(1,2),Rsitp(2),Xsitp0(1,2),dum)
      tmdly = Tmdly2 + Rsitp(2)
      if(srpt) tmdly = tmdly - Pradls
      Nit(19) = Nit(19) + 1
      itcnt = itcnt + 1
      if(ABS(tmdly-Tmdly1).gt.acctim) then
         itcnt = itcnt + 1
         if(itcnt.gt.100) call SUICID('TOO MANY SEND ITERATIONS',6)
         goto 2000
      endif
 
      call TSAV(tmdly - Tmdly2,3,Jd,Utrec)
c
c test for dummy observation below horizon
      if(Idumob.eq.1) then
         call HORIZN(Xsite,Xsitep,2,Jd)
         if(Jd.le.0) goto 9990
      endif
c
c*  start=2200
c calculation of planetary longitude and latitude
      if(nplng.le.0) goto 2500
      if(srpt) call SPOTCD(Jdx,Fract,mnspt1,-1,Nplnt0,0._10,0._10,1)
      if(nvsp.ge.0) call PRODCT(Rot,Xsitp0,arpl,3,3,1)
      ylong = ATAN2(arpl(2),arpl(1))/Convd + Pcom(5)
      if(ylong.lt.0._10) ylong = ylong + 360._10
      ylat = ASIN(arpl(3))/Convd
      if(Spotf.eq.amper9) goto 10
      zlong = ylong
      zlat  = ylat
      if(Spotf.eq.pound9) then
         zlong = ylong + Result(2)
         if(zlong.gt.360._10) zlong = zlong - 360._10
         if(zlong.lt.0._10) zlong   = zlong + 360._10
 
c correct latitude for inclination of radar equator
         zlat = ylat - 0.5*SIN(Result(2)*Convd)**2*arpl(3)
     .          *SQRT(1._10 - arpl(3)**2)/Convd
      endif
 
c delay correction due to planet shape
      oldtopo = ztopo
      ztopo = 0._10
      if(Nice.le.0 .or. nplng.gt.1) then
 
c -- subradar point or point relative to subradar point
         if(npshp.ge.0 .and. Nspot.le.0) then
 
c branch depending on shape model used (see bodred)
            if(Nshape(Npshp1).eq.0 .or. Nshape(Npshp1).eq.1) then
 
c spherical harmonic model or fourier series (nshape=0 or 1)
               call HARSHP(zlat,zlong,tmdly,ztopo,1,0._10,1)
            else if(Nshape(Npshp1).eq.2) then
 
c altitude grid model (local shape model) nshape=2
               call GRDSHP(zlat,zlong,tmdly,ztopo)
            else if(Nshape(Npshp1).eq.3) then

c external topography model, not recalculated here
               ztopo=oldtopo
               tmdly=tmdly-ztopo
            else
               call SUICID(
     .            ' NSHAPE OUTSIDE ALLOWABLE RANGE. STOP IN RADPL  ',12)
            endif
         endif
      endif
   10 if(nplng.gt.1) then
         if(Line.gt.55) call OBSPAG
         if(.not.srpt) then
            radofe = radoff
            if(Spotf.eq.pound9) radofe = radofe + ztopo*0.5_10
            radofe = radofe*Ltvel
            write(Iout,20) ylong,ylat,zlong,zlat,radofe
   20       format(23x,'SUBRADAR LON,LAT',2F8.3,
     .        37x,'SPOT LON,LAT,R',2F8.3,f10.3)
            Line = Line + 1
         else
            write(Iout,30) ylong,ylat,ztopo,Xsitp0
   30       format(' YLONG=',f17.12,'  YLAT=',
     .        f17.12,'  ZTOPO=',1pd15.8,
     . '  XSITP0 (UNIT VECTORS SUBRADAR PT TO SITE):'/ 6D22.13)
            Line = Line + 2
         endif
      endif
c
c*  start=2500
c determination of velocities
 2500 call SITVEL(2,nvel,norm)
      if(nvel.gt.0) then
         call PLTRP(1,Jdx,Fract,1,0)
         call ETRP(1,Jdy,frsend,1,0,2,3)
 
         do j = 4,6
            Xplsc(j) = Xp(j)*Aultvl
         end do
         if(.not.srpt) then
            call SPOTCV(nvel,Nplnt0,1)
            do j = 4,6
               Xplsc(j) = Xplsc(j) + Xspcd(j,1)
            end do
         endif
c
c fill relative velocity arrays
         call VLRTRD(Xemlsc(4,1),Xplsc(4),Xemlsc(4,2),1._10,7,0)
      endif
c
c*  start=3000
c          corrections
c
c       relativistic correction.  since radrel is used, xsbsun must be
c       set up. however, since vectors to the subradar point are not
c       calculated in the above code, we must use vectors to the c.m.
c       of the planet. this produces negligable errors.
c       for logspt, the calculations are done correctly.
c
      do i = 1,6
         Xsbsun(i) = Xplsc(i)
      end do
 
      call RADREL(1)
c
c skip all other corrections if just doppler
      if(Nice.le.0) then
c
c propagation corrections
         call PROPCO(1,1)
         tmdly = tmdly + Sumcor(1)
 
c antenna corrections
         tmdly = tmdly + ANCOR(1) + ANCOR(2)
c
c coordinate time delay changed to atomic time delay
         tmdly = tmdly + (CTATF(Jd,(Ctrec-tmdly-Ctat)/Secday,
     .           Ntrmct,2) - Ctat)
 
c time delay in atomic time changed to time delay in ut time
         if(ntime.gt.0) tmdly = tmdly*fdev
 
c constant bias in time delay
         if(Nrbias.gt.0) tmdly = tmdly + Rbsx(1)
 
c power series for clock errors uses /eqenox/ quantities
         if(Neqnox.gt.0) tmdly = tmdly + Pnox +
     .       dutrec*(Pequat + dutrec*Plat/2._10)
      endif
c
c*  start=4000
c doppler shift
      if(Nice.ge.0) then
         call RADREL(2)
         call PROPCO(2,1)
         dop = dop + Sumcor(2)
c
c effect of planet shape on doppler shift
         if(npshp.le.0) then
         endif
 
         dop = Freq*dop
c
c constant bias in doppler shift
         if(Nrbias.gt.0) dop = dop + Rbsx(2)
 
c power series for clock errors uses /eqenox/ quantities
         if(Neqnox.gt.0) dop = dop - Freq*(Pequat + dutrec*Plat)
      endif
      goto 9990
c
c*  start=9910
c missing coordinates, signal lost point
 9910 Jd = 0
c
 9990 return
 
      end
