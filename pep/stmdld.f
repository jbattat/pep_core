      subroutine STMDLD

      implicit none

c
c     m.e.ash   aug 1970    subroutine stmdld
c     compute delay and doppler for signal
c     sent from a ground station to a satellite over to another
c     satellite and then down to another ground station
c        calculate theoretical value of radar or radio transponder
c        observation of one satellite by another (both with same
c        central body)
c        interpolation through sbtrp,sctrp added 1977 oct - j.f.chandler
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 der3(6),rsc1
      equivalence (Angdum,der3),(Angdum(7),rsc1)
      real*10 rsc2,xsd(3,2),der4(3)
      equivalence (Angdum(8),rsc2),(xsd,Raddum(5)),(der4,der3(4))
      include 'eqnphs.inc'
      include 'ltrapobs.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'namtim.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct,tmdly0,dop0
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (Result,tmdly0),(Result(2),dop0)
      real*4    acctim
      equivalence (acctim,Estf)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'

c external functions
      real*10 A1UT1,A1WWV,CTATF,DOT,UT2UT1

c local variables
      real*10 fr1,fr2,fr3,tmdly4,tmdly5,xpdly0,xpdly1
      integer*4 i,iter1,iter2,iter3,iter4,j,jd3,lsw,norm,nvel,lmnctl
c
c determine time quantities
      lsw    = 1
      lmnctl = 1
      Xpdly  = 0._10
      xpdly0 = 0._10
      Utrec  = Ihr*3600 + Imin*60
      Utrec  = Utrec + Sec
      Fract  = Utrec/Secday

c fract is time signal time (uts) fraction of day
      if(itime.eq.2) then
c atuts,ututs read in

      else if(itime.eq.0) then
c observation time is UT2
         Ututs = -UT2UT1(Jd,Fract)
         Atuts = A1UT1(Jd,Fract) + Ututs

      else
c observation time is UTC
         Atuts = A1WWV(Jd,Fract)
         Ututs = Atuts - A1UT1(Jd,Fract)
      endif
c ututs = ut1 - given observation time
c atuts = a.1 - given observation time
      Ctat  = CTATF(Jd,Fract + Atuts/Secday,1,1)
      Ctrec = Utrec + (Ctat + Atuts)
      Utrec = Utrec + Ututs
      call SIDTIM(Jd,Ctat+Atuts-Ututs,Sidtm0,Sidvel,Dera)
      call PLATMO(Jd)
      Kindnp = ndprec*(Ncode - 1)
      do while(Ctrec.lt.0)
         Jd    = Jd - 1
         Ctrec = Ctrec + Secday
      end do
      do while( Ctrec.ge.Secday )
         Jd    = Jd + 1
         Ctrec = Ctrec - Secday
      end do
      ctrecf = Ctrec/Secday
      Tmdly1 = tmdly0*0.5_10
c
c
c nutation-precession determined for receiving
      Fract = ctrecf
      if(nprec.le.0) Fract = ctrecf - Tmdly1/Secday

c read moon tape for precession-nutation data
      call MNREED(Jd)
      if(Jd.le.0) goto 999
      call PRCNUT(Jd,Fract)
      Sidtm = Sidtm0 + Utrec*Sidvel + Dgst
c
c determination of geocentric receiving site coordinates
      norm = 0
      nvel = 0
      if(Nice.ge.0) nvel  = 1
      if(neatm.ge.0) norm = 1
      if(neatm.eq.0 .and. Nice.gt.0) norm = 0
c           nvel =0 positions determined
c           nvel =1 positions and velocities determined
c           nvel =2 positions, velocities and accelerations determined
c           nvel = 3 position,velocity,acceleration and jerk determined
c           norm =0 normals to sites not determined
c           norm =1 normals to sites determined
      call SITCOR(Sidtm,1,nvel,norm)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           iteration to calculate position of first satellite
c
c           determine reflection time from first satellite
      iter3  = 0
      tmdly4 = 0._10

  140 call TIMINC(Jd,ctrecf,Jdx,fr1,-tmdly4/Secday)
      call SCTRP(1,Jdx,fr1,0,lsw,1)
      if(Jdx.le.0) goto 999
      do i = 1, 3
         Xscsun(i,1) = Xsc(i,1)*Aultsc
         end do
      lsw = -1
c
c decide if lunar satellite instead of earth satellite
      if(Npcent(Klans1).eq.10) then
         call MNTRP(1,Jdx,fr1,0,lmnctl,1)
         if(Jdx.le.0) goto 999
         do i=1,3
            Xscsun(i,1)=Xscsun(i,1) + Xm(i,1)*Mnltsc
         end do
         lmnctl = -1
      endif
c
c determine time delay from first site to first satellite
      do i = 1, 3
         Xsitep(i,1) = Xsite(i,1) - Xscsun(i,1)
      end do
      Rsitp(1) = SQRT(DOT(Xsitep(1,1),Xsitep(1,1)))
      tmdly5   = Rsitp(1) - xpdly0
      iter3    = iter3 + 1
      if(ABS(tmdly5-tmdly4).gt.acctim) then
         if(iter3.gt.1000) call SUICID(
     .    ' MORE THAN 1000 ITERATIONS, STOP IN STMDLD  ', 11)
         tmdly4 = tmdly5
         goto 140
      endif
      Nit(20) = Nit(20) + iter3
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           iteration to calculate position of second satellite
c
c           determine reflection time
      iter1   = 0
      Tmdly1  = tmdly5*2._10
      lsw     = 1

  150 call TIMINC(Jd,ctrecf,Jdy,fr2,-Tmdly1/Secday)
c
c see if central body is observed
      if(Ncs1.ne.Nplnt0) then
c
c see if planet or moon is observed
         if(Klanb.le.0) call SUICID(
     .' RADAR OBSERVATIONS OF PLANETS FROM BODY NOT THE EARTH ARE NOT CO
     .DED, STOP IN STDLDP', 21)

         call SBTRP(1,Jdy,fr2,0,lsw)
         if(Jdy.le.0) goto 999
         do i = 1, 3
            Xsbsun(i) = Xsb(i)*Aultsc
            end do
         lsw = -1
c
c decide if lunar satellite instead of earth satellite
         if(Npcent(Klanb).eq.10) then
            call MNTRP(1,Jdy,fr2,0,lmnctl,2)
            if (Jdy.le.0) goto 999
            do i=1,3
               Xsbsun(i)= Xsbsun(i) + Xm(i,2)*Mnltsc
               end do
            lmnctl = -1
         endif
      else
         do i = 1, 3
            Xsbsun(i) = 0._10
c logic to go here for observing subradar point
c or spot on central body
            end do
      endif
c
c determine reflection to receive time delay, decide
c whether to re-iterate for reflection time
      do i = 1, 3
         xsd(i,1) = Xscsun(i,1) - Xsbsun(i)
         end do
      rsc1   = SQRT(DOT(xsd(1,1),xsd(1,1)))
      Tmdly2 = rsc1 + tmdly5 - (Xpdly + xpdly0)
      iter1  = iter1 + 1
      if(ABS(Tmdly2-Tmdly1).gt.acctim) then
         if(iter1.gt.1000) call SUICID(
     .' MORE THAN 1000 REFLECTION TIME ITERATIONS, STOP IN STDLDP  ',15)
         Tmdly1 = Tmdly2
         goto 150
      endif
      Nit(19) = Nit(19) + iter1
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c which 2-spacecraft time delay observable is this
      xpdly1 = Xpdly
      if(Ncodf.eq.14) then
         xpdly1 = xpdly0
c
c start iteration back to first spacecraft
         iter4  = 0
         Tmdly1 = Tmdly2 + (Tmdly2 - tmdly5)
         lsw    = 0
  160    call TIMINC(Jd,ctrecf,jd3,fr3,-Tmdly1/Secday)
         call SCTRP(1,jd3,fr3,0,0,2)
         if(jd3.le.0) goto 999
         do i = 1, 3
            Xscsun(i,2) = Xsc(i,2)*Aultsc
            end do
c
c decide if lunar satellite instead of earth satellite
         if(Npcent(Klans1).eq.10) then
            call MNTRP(1,jd3,fr3,0,lmnctl,1)
c
c write over moon from previous pass through this satellite
            if(jd3.le.0) goto 999
            do i = 1,3
               Xscsun(i,2) = Xscsun(i,2) + Xm(i,1)*Mnltsc
               end do
         endif
c
c determine send to reflection time delay, decide
c whether to re-iterate for send time
         do i = 1, 3
            xsd(i,2) = Xscsun(i,2) - Xsbsun(i)
            end do
         rsc2  = SQRT(DOT(xsd(1,2),xsd(1,2)))
         tmdly = Tmdly2 + rsc2 - (Xpdly + xpdly0)
         iter4 = iter4 + 1
         if(ABS(tmdly-Tmdly1).gt.acctim) then
            if(iter4.gt.1000) call SUICID(
     .        ' MORE THAN 1000 SEND TIME ITERATIONS, STOP IN STMDLD',13)
            Tmdly1 = tmdly
            goto 160
         endif
         Tmdly2  = tmdly
         Nit(17) = Nit(17) + iter4
      endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           iteration to calculate geocentric sending site coordinates
c
      iter2  = 0
      Tmdly1 = Tmdly2 + tmdly5

  250 Sidtm = Sidtm0 + Sidvel*(Utrec-Tmdly1) + Dgst
      call SITCOR(Sidtm,2,-1,0)
c
c determine time delay, decide whether to reiterate for
c send time
      do i = 1, 3
         if(Ncodf.ne.14) then
            Xsitep(i,2) = Xsite(i,2) - Xsbsun(i)
         else
            Xsitep(i,2) = Xsite(i,2) - Xscsun(i,2)
         endif
         end do
      Rsitp(2) = SQRT(Xsitep(1,2)**2 + Xsitep(2,2)**2 + Xsitep(3,2)**2)
      tmdly    = Tmdly2 + Rsitp(2) - xpdly1
      iter2    = iter2 + 1
      if(ABS(tmdly-Tmdly1).gt.acctim) then
         if(iter2.gt.1000) call SUICID(
     .        ' MORE THAN 1000 SEND TIME ITERATIONS, STOP IN STDLDP',13)
         Tmdly1 = tmdly
         goto 250
      endif
      Nit(18) = Nit(18) + iter2
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c test for dummy observation below horizon
      if(Idumob.eq.1) then
         call HORIZN(Xsite,Xsitep,2,Jd)
c
c test for dummy observation occulted by central body
         call CNTOCC(Xsbsun,xsd(1,1),-1,Jd)
         if(Jd.le.0) return
      endif
c
c calculate unit vectors
      do j = 1, 3
         der3(j)     = 0._10
         der3(j + 3) = 0._10
         do i = 1, 2
            Xsitp0(j,i) = Xsitep(j,i)/Rsitp(i)
            end do
         end do
c
c correct time delay for various effects
      if(Nice.le.0) then
c
c time delay in atomic time changed to time delay in ut time
         if(ntime.gt.0) tmdly = tmdly*fdev
c
c constant bias in time delay
         if(Nrbias.gt.0) tmdly = tmdly + Rbsx(1)
      endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c calculate velocites and doppler shift
      if(Nice.ge.0) then
c
c calculate doppler observable
         call SUICID(
     .'CANNOT CALCULATE 2-SPACECRAFT DOPPLER OBSERVABLE, STOP IN STMDLD'
     ., 16)
      endif

      return
c
c missing ephemeris data
  999 Jd = 0
      return
      end
