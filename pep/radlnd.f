      subroutine RADLND(nvel,norm,npath,idopob)
 
      implicit none

c          r. goldstein     april, 1976
c          interpolations are now done by calls to mntrp,emtrp,
c          pltrp, and sbtrp.
c          polynomial extrapolations are now done by dguess
c
c       sept. 1976..a derivative of sbdldp (sbbg)
c         this routine handles planetary landers only
c          ncodf.eq.2       nspot.gt.0
c       oct. 1976....large sections of code replaced by calls to
c         subroutines propco, phadop, and radrel
c       june 1978...converted to implicit real*8,
c       fixed interface to trp routines
c
c     this subroutine is called by radctl
c
c arguments
      integer*4 nvel,norm,npath,idopob

c array dimensions
      include 'globdefs.inc'

c        common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 dutrec
      equivalence (Angdum(10),dutrec)
      include 'eqnphs.inc'
      include 'ltrapobs.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'mnsprt.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct,tmdly0,dop0
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (Result,tmdly0),(Result(2),dop0)
      real*4 acctim
      equivalence (acctim,Estf)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'prpgat.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'
c           obscon(2)=transmitter frequency multiplicative factor
c
c     xslcns(1-6,j) = position and velocity (not used) of the center
c                     of mass of the solar system relative to the sun
c                     at receive time (j=1), reflection time (j=2)
c                     and send time (j=3)
c
c           receive time site and earth coordinates determined in radctl
c
c external functions
      real*10 ANCOR,CTATF

c local variables
      real*10 dum,dum1,frcns,sendtm
      integer   ipct,j,jpct,lemctl,loop1,loop2,lplctl,mnspt1,nn

c       tests on required flags
c
      if(Klanb.ne.0) call SUICID(
     .' CANNOT CALCULATE OBSERVATIONS OF A LANDER ON A MOON IN SUBROUTIN
     .E RADLND   ',19)
      if(Klan.le.0) call SUICID(' KLAN NO GOOD IN RADLND ',6)
c
c------------------setup for interpolation (reflection)-----------------
      mnspt1 = 0
      lplctl = -1
      lemctl = 0
      loop1  = 0
      loop2  = 0
      do while( .true. )
c
c*  start=1100
c reflection time iteration loop
c tmdly1, the first guess, was generated in radctl
         call TIMINC(Jd,ctrecf,Jdx,Fract,-Tmdly1/Secday)
c
c spot on planet
         call SPOTCD(Jdx,Fract,mnspt1,nvel,Nplnt0,dum1,dum1,1)
         do j=1,3
            Xsbpl(j)=Xspcd(j,1)
         end do
c
c------------------obtain lander wrt receiv.site at refl. time----------
         call PLTRP(1,Jdx,Fract,0,lplctl)
         if(Jdx.le.0) then
c
c*  start=9910
            Jd = 0
            goto 9990
         endif
         do j = 1,3
            Xplsc(j)  = Xp(j)*Aultsc
            Xsbsun(j) = Xspcd(j,1) + Xplsc(j)
         end do
c
c*  start=1400
         do j = 1,3
            Xsitep(j,1) = Xemlsc(j,1) - Xsbsun(j)
            Xsitec(j,1) = Xemlsc(j,1) - Xplsc(j)
         end do
         if(Nswcns.gt.0) then
            call SOTRP(Jdx,Fract,Xslcns(1,2),0)
            do j = 1,3
               dum = (Xslcns(j,2) - Xslcns(j,1))*Aultsc
               Xsitep(j,1) = Xsitep(j,1) + dum
               Xsitec(j,1) = Xsitec(j,1) + dum
            end do
         endif
c
c------------------calculate delay time for reflection------------------
         call UVECTR(3,Xsitep(1,1),Rsitp(1),Xsitp0(1,1),dum)
         Tmdly2 = Rsitp(1) - Xpdly
         loop1  = loop1 + 1
         if(loop1.gt.20)
     .        call SUICID('MORE THAN 20 REFLECT ITERATIONS ',8)
         if(ABS(Tmdly2-Tmdly1).le.acctim) goto 2000
         Tmdly1 = Tmdly2
      end do
c
c*  start=2000
c iteration to determine sending site position
 2000 Nit(20) = Nit(20) + loop1
      Tguess  = Tmdly1
      Dstf(4) = Tmdly1
      call TSAV(Tmdly2,1,Jd,Utrec)
      call DGUESS(tmdly,3,Jd,Utrec)
      tmdly = tmdly + Tmdly2
      do while( .true. )
 
         Tmdly1 = tmdly
         sendtm = Utrec - Tmdly1*fdev
         Sidtm  = Sidtm0 + Sidvel*sendtm
         call SITCOR(Sidtm,2,-1,0)
         call TIMINC(Jd,ctrecf,Jdy,frcns,-Tmdly1/Secday)
         call ETRP(1,Jdy,frcns,0,lemctl,2,3)
         if(Jdy.le.0) then
            Jd = 0
            goto 9990
         endif
c
c------------determine sending site wrt lander and planet at send time---
         do j = 1,3
            Xsitep(j,2) = Xemlsc(j,2) - Xsbsun(j)
            Xsitec(j,2) = Xemlsc(j,2) - Xplsc(j)
         end do
         if(Nswcns.gt.0) then
            call SOTRP(Jdy,frcns,Xslcns(1,3),0)
            do j = 1,3
               dum = (Xslcns(j,2) - Xslcns(j,3))*Aultsc
               Xsitep(j,2) = Xsitep(j,2) + dum
               Xsitec(j,2) = Xsitec(j,2) + dum
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
         nn = -2
         call HORIZN(Xsite,Xsitep,nn,Jd)
c
c
c test for spot on planet turned away from earth
         call CNTOCC(Xspcd,Xsitep,-2,Jd)
c
c test for Sun alignment in either direction
         call SUNOCC(Xsbsun,Xsitep,2,Jd)
 
         if(Jd.le.0) goto 9990
      endif
c
c*  start=2500
c
c     calculate velocities
c
      call SITVEL(2,nvel,norm)
      call SPOTCV(nvel,Nplnt0,1)
 
c velocity of site
      if(nvel.gt.0) call ETRP(1,Jdy,frcns,1,0,2,3)
 
c
c
c
c------------------correct time delay for various effects---------------
c*  start=3000
c
c     general relativistic correction to time delay
      call RADREL(1)
c correct for propagation effects by calling propco
c    set ipct as follows:
c    -2: first leg of phase delay doppler
c    -1: second leg of phase delay doppler
c     1: time delay
c
c skip corrections if only doing instantaneous doppler
c
      ipct = 1
      if(Nice.gt.0) then
c
c doppler
         if(lopler.ne.-1) goto 4010
         if(npath.eq.1) ipct = -2
         if(npath.eq.2) ipct = -1
      endif
c
c
c
      call PROPCO(ipct,1)
      if(ipct.eq.1) tmdly = tmdly + Sumcor(1)
 
c correct for antenna offsets
      tmdly = tmdly + ANCOR(1) + ANCOR(2)
c
c we can skip the rest of the corrections if we are only
c doing phase delay doppler (they cancel out)
c
      if(lopler.ne.-1 .or. Nice.le.0) then
c
c
c coordinate time delay changed to atomic time delay
         tmdly = tmdly + (CTATF(Jd,(Ctrec-tmdly-Ctat)/Secday,Ntrmct,2)
     .                    - Ctat)
 
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
c------------------doppler---------------------------------------------
c*  start=4000
c
      if(Nice.ge.0) then
c
c    there are three possible doppler observables:
c       lopler=-1  observable is phase
c       lopler= 0  observable is averaged integrated frequency (jpl)
c       lopler= 1  observable is instantaneous frequency
c       lopler=0 has been removed from the code and has been gone
c       since before 1975.
c
         if(lopler.ne.1) then
c
c phase delay doppler
            call PHADOP(npath,idopob)
c
c       all propagation corrections were calculated and stored by
c       propco in the time delay calculations. now call propco for
c       a third time to apply those corrections to the doppler shift.
c
            jpct = 0
            if(calcvl.le.0) goto 4030
         endif
      else if(calcvl.le.0) then
         goto 9990
      endif
c
c------------------determination of velocities--------------------------
c
c     instantaneous doppler
c
c      velocity of lander
 4010 call PLTRP(1,Jdx,Fract,1,lplctl)
      if(Jdx.le.0) then
         Jd = 0
         goto 9990
      endif
      do j = 4,6
         Xplsc(j)  = Xp(j)*Aultvl
         Xsbsun(j) = Xspcd(j,1) + Xplsc(j)
      end do
 
c relative velocity
      call VLRTRD(Xemlsc(4,1),Xsbsun(4),Xemlsc(4,2),1._10,7,0)
      if(Nice.lt.0) goto 9990
      if(lopler.gt.0) then
c
c general relativistic doppler calculation
         call RADREL(2)
         jpct = 2
      endif
c
c correct for propagation effects by calling propco
 4030 if(lopler.ne.-1 .or. npath.eq.1) then
         call PROPCO(jpct,1)
         dop = dop + Sumcor(2)
 
         dop = Freq*dop
c
c constant bias in doppler shift
         if(Nrbias.gt.0) dop = dop + Rbsx(2)
 
c power series for clock errors
         if(Neqnox.gt.0) dop = dop - (Pequat + Plat*dutrec)*Freq
      endif
 
c*  start=9990
 9990 return
      end
