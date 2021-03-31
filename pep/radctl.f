      subroutine RADCTL(idopob)
 
      implicit none
 
c          r. goldstein, april,1976
c          interpolations are now done by subroutines mntrp,emtrp,
c          pltrp, and sbtrp
c          polynomial extrapolations are done by dguess
c
c       sept. 1976...a derivative of deldop(dldpbg).
c          this routine does receive site coordinate calculations
c          and other setup functions then transfers control to
c          one of several routines to continue calculations.
c       june 1978...converted to implicit real*8, fix call to radmn,
c       fixed interface to trp routines
c
c arguments
      integer*4 idopob

c        idopob = 0 when radctl first called for observation
c        idopob = 1 set in radctl for phase delay doppler observable
c                   on first call of radctl
c        idopob =-1 set in radctl on second call (radctl called second
c                   time only if idopob was set 1 on first call)
c

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 dutrec
      equivalence (Angdum(10),dutrec)
      include 'difnct.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'ltrapobs.inc'
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
      include 'kobequiv.inc'
      include 'param.inc'
      include 'plnhar.inc'
      include 'prpgat.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'

c local 
      character*4 amper9/'&&&&'/, pound9/'####'/
      real*10 epsilo/0.01_10/
      real*10 angle,csptln,csptlt,dsec,dum,fract0,pnoff,ssptln,
     . ssptlt
      integer*4 jda,jdpn,lemctl,norm,npath,nvel

c external functions
      real*10 A1UT1,A1WWV,CTATF,UT2UT1

      if(itime.eq.5) then
         call RADNRM
         goto 100
      endif
c
c see if this is second call to radctl for phase delay
      if(idopob.gt.0) then
         idopob = -1
      else
c
c*  start=100
c setup time quantities
         npath = 1
         Jd    = Jds
         Utrec = Ihr*3600 + Imin*60
         dsec  = Sec
         if(Numsav.ge.40 .and. ABS(Save(40)-dsec).lt.1.E-3_10)
     .    dsec = Save(40)
         Utrec = Utrec + dsec
         if(itime.ge.3) Utrec = Utrec + tmdly0
c itime.ge.3 means that input time of observation is send time
c instead of receive time (does not work in dummy mode because
c observed value of time delay does not exist)
         if(Ncodf.eq.2) then
 
c set up default quantities if not from obslib tape
            if(Idumob.ne.-2) then
               Obscon(2) = 1._10
               if(Numsav.lt.29) Numsav = 29
               Save(29) = Freq
               Save(28) = 600._10
               if(Klan.eq.17) Save(28) = 300._10
               if(Klanb.gt.0 .and. Ncp0.gt.0) Save(28) = 60._10
            endif
c
c if only doing phase delay doppler, shift utrec by
c 1/2 the count interval
            if(Nice.gt.0) then
               if(lopler.eq.-1) then
                  Contig = .false.
                  Utrec  = Utrec - Save(28)/2.0_10
                  if(Utrec.lt.0.0_10) then
                     Utrec = Utrec + Secday
                     Jd    = Jd - 1
                  endif
 
c If beginning time of this observation is the same as final time of
c previous observation, simply update quantities and move to final
c time.
                  if(Jd.eq.Jdsav .and. ABS(Utrec-Utrsav).le.epsilo)
     .                call PHDMVE(npath,idopob)
               endif
            endif
         endif
      endif
c
c*  start=200
c
      call DGUESS(Tmdly1,1,Jd,Utrec)
      Jdsav  = Jd
      Utrsav = Utrec
      Fract  = Utrec/Secday
 
c calculate elapsed time since first observations for clock part
      if(Nk1.lt.0) fract0 = Fract
      dutrec = (Fract - fract0)*Secday
      if(dutrec.lt.0._10) dutrec = dutrec + Secday
 
c fract is time signal time (uts) fraction of day
      if(itime.eq.2) then
c atuts,ututs read in

      else if(itime.eq.0) then
c observation time is ut2
         Ututs = -UT2UT1(Jd,Fract)
         Atuts = A1UT1(Jd,Fract) + Ututs

      else
c observation time is wwv
         Atuts = A1WWV(Jd,Fract)
         Ututs = Atuts - A1UT1(Jd,Fract)
      endif
c ututs = ut1 - given observation time
c atuts = a.1 - given observation time
c ctat = 32.15 for determining precession, nutation and site cor
      call TIMINC(Jd,Fract,jda,Frect,Atuts/Secday)
      Ctat  = 32.15_10 +
     . Ctvary*(jda-Prm97-0.5_10+Frect+32.15_10/86400._10)**2*Secday
      Utrec = Utrec + Ututs
      call SIDTIM(Jd,Ctat+Atuts-Ututs,Sidtm0,Sidvel,Dera)
      call PLATMO(Jd)
      Kindnp = ndprec*(Ncode - 1)
      call TIMINC(jda,Frect,Jd,ctrecf,Ctat/Secday)
      Ctrec = ctrecf*Secday
c
c*  start=500
c read moon peripheral data set
      call MNREED(Jd)
      if(Jd.le.0) goto 100
c
c nutation-precession determined for receiving
      pnoff = 0._10
      if(nprec.le.0) pnoff = -Tmdly1/Secday
      call TIMINC(Jd,ctrecf,jdpn,Fract,pnoff)
      call PRCNUT(jdpn,Fract)
 
c calculate diurnal polar motion
      if(mod(Jct(28)/4,2).ne.0) then
         call DIURNP(jdpn,Fract,(Sidtm0+Utrec*Sidvel),Dtheta)
         Sidtm0 = Sidtm0 - Dtheta
      endif
      Sidtm = Sidtm0 + Utrec*Sidvel + Dgst
c
c determination of geocentric receiving site coordinates
      norm = 0
      if(Izctl.ne.0) norm = 1
      if(Jsite(1).eq.2 .or. Jsite(2).eq.2) norm = 1
      nvel = 0
      if(calcvl.gt.0) nvel = 1
      if(Ncodf.eq.18) then
         nvel = 1
         if(Save(29).ne.0._10) Freq = Obscon(1)*Save(29)
      else if(Ncodf.eq.2) then
         nvel = 1
         Freq = Obscon(2)*Save(29)
      else
         if(Nice.ge.0) nvel = 1
      endif
      if(Jct(10).gt.0 .or. Jct(11).gt.0) nvel = 1
      if(Jct(45).ne.0) norm = 1
      if(Jct(10).gt.0 .or. Jct(49).gt.0) norm = 2
c nvel =0 positions determined
c nvel =1 positions and velocities determined
c nvel =2 positions, velocities and accelerations determined
c nvel =3 position,velocity,acceleration and jerk determined
c norm =0 normals to sites not determined
c norm =1 normals to sites determined
c norm =2 normals to site determined, and partials
c         determined for solid body tides even if
c         lsite(1-3)=0
c
c in old code, site normals + accel. were needed for special
c output tape if jct(40).ne.0
c but output is incompatible with radmn
c
      call SITCOR(Sidtm,1,nvel,norm)
 
c determine ctat including diurnal & monthly terms
      Ctat = CTATF(jda,Frect,Ntrmct,1)
      call TIMINC(jda,Frect,Jd,ctrecf,Ctat/Secday)
      Ctrec = ctrecf*Secday
c
c determine if spot away from subradar point is observed
c if so setup planet fixed spot coordinates
      if(Spotf.eq.pound9) then
c spot defined by longitude offset from sub-radar point
c spot latitude = sub-radar latitude to a good approx.
         radoff = Pradls
         angle  = Result(2)*Convd
         Csroff = radoff*COS(angle)
         Ssroff = radoff*SIN(angle)
      else if(Spotf.eq.amper9) then
c+++++++ does this code belong in radpl?
c+++++++ if not, then harshp and grdshp belong in radar link
         radoff = Pradls
         if(npshp.ge.0) then
            if(Npshp1.le.0) call SUICID(
     .'OBSERVATIONS AWAY FROM SUBRADAR POINT NOT ALLOWED FOR ELLIPSOIDAL
     . PLANET, STOP IN RADCTL', 22)
            if(Nshape(Npshp1).eq.2) then
               call GRDSHP(zlat,zlong,dum,ztopo)
            else if(Nshape(Npshp1).eq.3) then
            else
               call HARSHP(zlat,zlong,dum,ztopo,1,0.0_10,1)
            endif
            radoff = 0.5_10*ztopo + radoff
         endif
c radoff =radius to observed spot off the subradar point
c (equivalenced to rspot(1))
         angle  = Convd*Dstf(6)
         csptlt = COS(angle)
         ssptlt = SIN(angle)
         angle  = Convd*Dstf(5)
         csptln = COS(angle)
         ssptln = SIN(angle)
         Yspcd(1,1) = radoff*csptln*csptlt
         Yspcd(2,1) = radoff*ssptln*csptlt
         Yspcd(3,1) = radoff*ssptlt
      endif
c
c*  start=1000
c
c determine if planet is moon or probe in cislunar space
      if(Ncp0.ne.3 .and. Ncp0.ne.10) then
c
c
c determination of pos. of rec. site rel. to sun
         lemctl = -1
         call ETRP(1,Jd,ctrecf,0,lemctl,1,3)
         if(Jd.le.0) goto 100
         if(nvel.gt.0) call ETRP(1,Jd,ctrecf,1,0,1,3)
 
c get offset to solar system barycenter
         if(Nswcns.gt.0) then
            call SOTRP(Jd,ctrecf,Xslcns,0)
            if(Jd.le.0) goto 100
            if(Nvlcns.gt.0) call SOTRP(Jd,ctrecf,Xslcns,1)
         endif
c
c determination of nutation-precession for sending
         if(nprec.gt.0 .and. nddiff.ge.0) then
            call TIMINC(Jd,ctrecf,jdpn,Fract,-Tmdly1/Secdy2)
            call PRCNUT(jdpn,Fract)
         endif
         Sidtm0 = Sidtm0 + Dgst
c
c transfer control to one of 6 routines
c
         if(Ncodf.eq.18) then
 
c pulsar phase observations
            call RADPSR(nvel,norm)
            return
 
c determine if planet is sun
         else if(Nplnt0.le.0) then
            call SDELDP(nvel,norm)
 
         else if(Ncodf.eq.2) then
c
c ncodf=2, active target
c now decide if orbiter or lander
c
            if(Nspot.gt.0) call RADLND(nvel,norm,npath,idopob)
            if(Nspot.eq.0) call RADSB(nvel,norm,npath,idopob)
         else if(Ncodf.eq.3) then
c
c ncodf=3, differential radar and/or doppler
c
            call RADEXT(nvel,norm)
         else
c
c ncodf=1, passive target
c now decide if radar from a planet or a moon of a planet
c
            if(Klanb.gt.0) call RADNS(nvel,norm)
            if(Klanb.eq.0) call RADPL(nvel,norm)
         endif
      else if(Ncp0.eq.10) then
         call RADMN(nvel,norm,npath,idopob)
 
c decide if observable is one-way
      else if(Nsite2.ne.0 .or. Ncodf.eq.19) then
         call RADMN(nvel,norm,npath,idopob)
      else
         call RADGPS(nvel,norm,npath,idopob)
      endif
c
c*  start=9900
  100 Tmdly1 = Deriv(2,1)
      return
      end
