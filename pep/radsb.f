      subroutine RADSB(nvel,norm,npath,idopob)
 
      implicit none

c          r. goldstein     april, 1976
c          interpolations are now done by calls to mntrp,emtrp,
c          pltrp, and sbtrp.
c          polynomial extrapolations are now done by dguess
c
c       oct. 1976..large sections of code replaced by calls to
c         subroutines propco, radrel. and phadop
c       june 1978.. converted to implicit real*8 , interface with plcms
c         and streamlined radrel
c
c
c     this subroutine is called by radctl
c
c arguments
      integer*4 nvel,norm,npath,idopob

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 tpvmsb(2),dutrec
      equivalence (Angdum,tpvmsb)
      equivalence (dutrec,Angdum(10))
      include 'eqnphs.inc'
      include 'inodta.inc'
      include 'ltrapobs.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      include 'kobequiv.inc'
      real*10 ctrecf,fdev,reflct,tmdly0,dop0,tc
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (Result(1),tmdly0),(Result(2),dop0)
      real*4    acctim
      equivalence (acctim,Estf)
      equivalence (Save(28),tc)
      include 'param.inc'
      include 'prpgat.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'

c external functions
      real*10 ANCOR,CTATF 

c local
      real*10 dum,frcns,sendtm,tpvm,zsb(3)
      integer   ipct,j,jdpvm,jpct,lemctl,loop1,loop2,
     .          lplctl,lsbctl,nn,nrclc
c
c     xslcns(1-6,j) = position and velocity (not used) of the center
c                     of mass of the solar system relative to the sun
c                     at receive time (j=1), reflection time (j=2)
c                     and send time (j=3)
c
c           receive time site and earth coordinates determined in radctl
c
c     tests on required flags
      if(Nspot.ne.0) call SUICID(' NSPOT NO GOOD IN RADSB ',6)
c
c       time delay and phase delay doppler not allowed in same
c        series...see phadop for more detail
c
c------------------setup for interpolation (reflection)-----------------
      nrclc  = 0
      lsbctl = -1
      lplctl = -1
      lemctl = 0
      loop1  = 0
      loop2  = 0
      do while(.true.)
c
c          reflection time iteration loop
c       tmdly1, the first guess, was generated in radctl
c*  start=1100
c
         call TIMINC(Jd,ctrecf,Jdx,Fract,-Tmdly1/Secday)
 
         if(Klanb.gt.0) then
 
c for pioneer-venus entry probes, call pvcrd is within atm.
            jdpvm = Sbcom(1)
            if(jdpvm.eq.2443852) then
               tpvm = Jdx + Fract
               if(tpvm.ge.Sbcom(1)) then
                  call PVCRD(Jdx,Fract,1,1,nvel)
                  tpvmsb(2) = tpvmsb(1)
 
c partl1 looks for tpvmsb separately for each site
                  go to 20
               endif
            endif
            tpvmsb(1) = 0._10
            call SBTRP(1,Jdx,Fract,0,lsbctl)
            if(Jdx.le.0) then
               Jd = 0
               go to 300
            endif
c
c test filter stepover for phase delay doppler
            if(lopler.ne.-1 .or. npath.eq.1) Recalc = .false.
            if(Recalc) then
c
c*  start=9900
               write(Iout,10) Jdx,Fract
   10          format(' INSUFFICIENT FILTER STEPOVER NEAR JD.F=',
     .                i8,f9.8)
               Jd = 0
               go to 300
            endif
 
   20       do j = 1,3
               Xsbsun(j) = Xsb(j)*Aultsc
 
c save probe w.r.t. planet for occultation test
               zsb(j) = Xsbsun(j)
            end do
c
c if pioneer-venus multiprobe bus, correct antenna
c coordinates for roll
            if(jdpvm.eq.-2443852) call PVBUS(Fract,Xsbsun,Klanb,1)
         endif
 
         if(Klan.ne.0) then
c
c a planet is involved so obtain probe wrt sun
c
            call PLTRP(1,Jdx,Fract,0,lplctl)
            if(Jdx.le.0) then
               Jd = 0
               go to 300
            endif
            do j = 1,3
               Xplsc(j) = Xp(j)*Aultsc
            end do
            if(Klanb.gt.0) then
c
c probe wrt sun
               do j = 1,3
                  Xsbsun(j) = Xsbsun(j)*Cmfct + Xplsc(j)
               end do
            else
 
c dummy radio tracking to a planet
               do j = 1,3
                  Xsbsun(j) = Xplsc(j)
               end do
            endif
         endif
c
c*  start=1400
c get receive site w.r.t. probe
         do j = 1,3
            Xsitep(j,1) = Xemlsc(j,1) - Xsbsun(j)
         end do
         if(Nswcns.gt.0) then
            call SOTRP(Jdx,Fract,Xslcns(1,2),0)
            do j = 1,3
               Xsitep(j,1) = Xsitep(j,1) + (Xslcns(j,2) - Xslcns(j,1))
     .                        *Aultsc
            end do
         endif
c
c------------------calculate delay time for reflection------------------
         call UVECTR(3,Xsitep(1,1),Rsitp(1),Xsitp0(1,1),dum)
         Tmdly2 = Rsitp(1) - Xpdly
         loop1  = loop1 + 1
         if(loop1.gt.20)
     .        call SUICID('MORE THAN 20 REFLECT ITERATIONS ',8)
         if(ABS(Tmdly2-Tmdly1).le.acctim) then
 
c check filter epoch against current span limit
            if(lopler.eq.-1 .and. npath.eq.1) then
               call SBRED2(Jdx + Fract + tc/1728E2_10)
               if(Recalc) then
                  Recalc = .false.
                  if(nrclc.gt.0) call SUICID(
     .               'REPEATED FILTER BREAK CROSSING, STOP IN RADSB   '
     .               ,12)
                  nrclc = 1
                  go to 1900
               endif
            endif
            goto 2000
         endif
         Tmdly1 = Tmdly2
 1900 end do
 
c reflect time converged successfully
 2000 Nit(20) = Nit(20) + loop1
      Tguess  = Tmdly1
      Dstf(4) = Tmdly1
      call TSAV(Tmdly2,1,Jd,Utrec)
c
c*  start=2000
c iteration to determine sending site position
      call DGUESS(tmdly,3,Jd,Utrec)
      tmdly = tmdly + Tmdly2
      do while( .true. )
         Tmdly1 = tmdly
         sendtm = Utrec - Tmdly1*fdev
         Sidtm  = Sidtm0 + Sidvel*sendtm
 
c in computing sidereal time, tmdly1 converted to wwv seconds
         call SITCOR(Sidtm,2,-1,0)
         call TIMINC(Jd,ctrecf,Jdy,frcns,-Tmdly1/Secday)
         call ETRP(1,Jdy,frcns,0,lemctl,2,3)
         if(Jdy.le.0) then
            Jd = 0
            go to 300
         endif
 
c------------------determine  sending  site wrt probe at  send  time---
         do j = 1,3
            Xsitep(j,2) = Xemlsc(j,2) - Xsbsun(j)
         end do
         if(Nswcns.gt.0) then
            call SOTRP(Jdy,frcns,Xslcns(1,3),0)
            do j = 1,3
               Xsitep(j,2) = Xsitep(j,2) + (Xslcns(j,2)-Xslcns(j,3))
     .                        *Aultsc
            end do
         endif
c
c calculate delay time for send
         call UVECTR(3,Xsitep(1,2),Rsitp(2),Xsitp0(1,2),dum)
         tmdly = Tmdly2 + Rsitp(2) - Xpdly
         loop2 = loop2 + 1
         if(loop2.gt.20)
     .        call SUICID('MORE THAN 20 SEND ITERATIONS',7)
         if(ABS(tmdly-Tmdly1).le.acctim) goto 2100
      end do
 
 2100 Nit(19) = Nit(19) + loop2
      call TSAV(tmdly - Tmdly2,3,Jd,Utrec)
c
c - - - - - - - - - - end of sending iteration  - - - - - - - - - - - -
c
c test for dummy observation below horizon
      if(Idumob.eq.1) then
         nn = 2
         if(Ncp0.gt.0) nn = -2
         call HORIZN(Xsite,Xsitep,nn,Jd)
 
c test for planetary orbiter occulted by central body
         if(Ncp0.gt.0) call CNTOCC(zsb,Xsitep,-2,Jd)
         if(Jd.le.0) go to 300
      endif
c
c site velocities and partials
c*  start=2500
c
      call SITVEL(2,nvel,norm)
c
c calculate pioneer-venus probe velocity and partials
      if(Klanb.gt.0 .and. jdpvm.eq.2443852 .and.
     .   tpvm.ge.Sbcom(1)) call PVCRDV(1,nvel,1)
 
      if(nvel.gt.0) then
c
c velocity of site
         call ETRP(1,Jdy,frcns,1,0,2,3)
      endif
c
c------------------correct time delay for various effects---------------
c*  start=3000
c general relativistic correction to time delay
      call RADREL(1)
c
c       correct for propagation effects by calling propco
c       set up ipct as follows:
c       -2: first leg of phase delay doppler
c       -1: second leg of phase delay doppler
c        1: time delay
c
c       we can skip the corrections if we are only doing instantaneous
c     doppler
c
      ipct = 1
      if(Nice.gt.0) then
c
c doppler
c
         if(lopler.ne.-1) go to 30
         if(npath.eq.1) ipct = -2
         if(npath.eq.2) ipct = -1
      endif
c
c
      call PROPCO(ipct,1)
      if(ipct.eq.1) tmdly = tmdly + Sumcor(1)
c
c correct for antenna offsets
      tmdly = tmdly + ANCOR(1) + ANCOR(2)
c skip the rest of the corrections if only doing phase delay
c doppler. (they cancel out).
c
      if((lopler.ne.-1) .or. (Nice.le.0) ) then
c
c
c coordinate time delay changed to atomic time delay
         tmdly = tmdly + (CTATF(Jd,(Ctrec-tmdly-Ctat)/Secday,Ntrmct,2)
     .            - Ctat)
 
c time delay in atomic time changed to time delay in ut time
         if(ntime.gt.0) tmdly = tmdly*fdev
 
c constant bias in time delay
         if(Nrbias.gt.0) tmdly = tmdly + Rbsx(1)
 
c power series for clock errors
         if(Neqnox.gt.0) tmdly = tmdly + Pnox +
     .       (Pequat + Plat*dutrec/2._10)*dutrec
      endif
c
c
c
c------------------doppler---------------------------------------------
c*  start=4000
c
c
c
c       nice.lt.0  no doppler;            nice.gt.0  no time delay
      if(Nice.lt.0 .and. calcvl.le.0) go to 300
      if(Nice.lt.0) goto 30
c
c    there are three possible doppler observables:
c       lopler=-1  observable is phase
c       lopler= 1  observable is instantaneous frequency
c     loppler=0 has been removed from the code and has been gone
c     since before 1975.
c
      if(lopler.eq.-1) then
 
c phase delay doppler
         call PHADOP(npath,idopob)
c
c
c       all propagation corrections were calculated by propco
c       in the time delay calculations. now apply them
c       by calling propco for the third time
c
         jpct = 0
         if(calcvl.le.0) go to 200
      endif
 
   30 if(jdpvm.ne.2443852 .or. tpvm.le.Sbcom(1)) then
         if(Klanb.le.0) go to 40
         call SBTRP(1,Jdx,Fract,1,0)
      endif
      do j = 4,6
         Xsbsun(j) = Xsb(j)*Aultvl
      end do
   40 if(Klan.ne.0) then
         call PLTRP(1,Jdx,Fract,1,0)
         if(Klanb.gt.0) then
            do j = 4,6
               Xplsc(j)  = Xp(j)*Aultvl
               Xsbsun(j) = Xsbsun(j)*Cmfct + Xplsc(j)
            end do
         else
 
c dummy radio tracking to a planet
            do j = 4,6
               Xplsc(j)  = Xp(j)*Aultvl
               Xsbsun(j) = Xplsc(j)
            end do
         endif
      endif
c relative velocity
c store in xsitep array
      call VLRTRD(Xemlsc(4,1),Xsbsun(4),Xemlsc(4,2),1._10,7,0)
      if(Nice.lt.0) go to 300
      if(lopler.ne.-1) then
c
c general relativistic doppler calculation
         call RADREL(2)
         jpct = 2
      endif
c
c------------------correct doppler shift for various effects------------
c
c
c       correct for propagation effects
  200 if(lopler.ne.-1 .or. npath.eq.1) then
         call PROPCO(jpct,1)
         dop = dop + Sumcor(2)
         dop = Freq*dop
c
c constant bias in doppler shift
         if(Nrbias.gt.0) dop = dop + Rbsx(2)
c           power series for clock errors
c
c       here used to be the truncation for mariner2 and 4 doppler
c       observable..it has been removed and may it rest in peace.
c
         if(Neqnox.gt.0) dop = dop - (Pequat + Plat*dutrec)*Freq
      endif
 
c*  start=9990
  300 Recalc = .false.
      return
      end
