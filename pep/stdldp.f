      subroutine STDLDP
 
      implicit none
c
c m.e. ash - may 1970 - subroutine stdldp
c calculate theoretical value of radar or radio transponder
c observation of one satellite by another (both with same
c central body)
c
c interpolation through sbtrp,sctrp added 1977 oct - j.f.chandler

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'coord.inc'
      include 'eqnphs.inc'
      include 'ltrapobs.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'mnsprt.inc'
      include 'number.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct,tmdly0,dop0
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (Result(1),tmdly0),(Result(2),dop0)
      real*4    acctim
      equivalence (Estf,acctim)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'
c
c external functions
      real*10 A1WWV,CTATF,DOT

c local
      real*10 ctuts,fr2,fr3
      integer   i,iter1,iter2,j,k,lsw,mnspt1,nvel
c
c determine time quantities
      nvel = 0
      if(Nice.ge.0) nvel = 1
      lsw   = 1
      Xpdly = 0._10
      Utrec = Ihr*3600 + Imin*60
      Utrec = Utrec + Sec
      Fract = Utrec/Secday
      Atuts = A1WWV(Jd,Fract)
      Ctat  = CTATF(Jd,Fract + Atuts/Secday,1,1)
      Ututs = 0.
      ctuts = Ctat + Atuts
      Ctrec = Utrec + ctuts
      do while(.true.)
         if(Ctrec.lt.0) then
            Jd    = Jd - 1
            Ctrec = Ctrec + Secday
         else if(Ctrec.lt.Secday) then
            go to 100
         else
            Jd    = Jd + 1
            Ctrec = Ctrec - Secday
         endif
      end do
  100 ctrecf = Ctrec/Secday
      Tmdly1 = tmdly0/2._10
      if(Dstf(4).gt.0._10) Tmdly1 = Dstf(4)
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           calculate position of observing satellite at
c           receive time
c
      call SCTRP(1,Jd,ctrecf,0,lsw,1)
      if(Jd.le.0) then
 
c missing ephemeris data
         Jd = 0
         return
      endif
      do i = 1,3
         Xscsun(i,1) = Xsc(i,1)*Aultsc
      end do
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           iteration to calculate position of observed satellite
c           at reflection time
c
c           determine reflection time
      iter1  = 0
      mnspt1 = 0
  200 call TIMINC(Jd,ctrecf,Jdx,fr2,-Tmdly1/Secday)
c
c see if central body is observed
      if(Ncs1.ne.Nplnt0) then
c
c see if planet or moon is observed
         if(Klanb.le.0) call SUICID(
     .' RADAR OBSERVATIONS OF PLANETS FROM BODY NOT THE EARTH ARE NOT CO
     .DED, STOP IN STDLDP',21)
c
c read probe peripheral data set
         call SBTRP(1,Jdx,fr2,0,lsw)
         if(Jdx.le.0) then
            Jd = 0
            return
         endif
         do i = 1,3
            Xsbsun(i) = Xsb(i)*Aultsc
         end do
         lsw = -1
      else if(Nspot.gt.0) then
         call SPOTCD(Jdx,fr2,mnspt1,0,Nplnt0,ctuts,Atuts,1)
         if(Jdx.le.0) then
            Jd = 0
            return
         endif
         do i = 1,3
            Xsbsun(i) = Xspcd(i,1)
         end do
      else
         do i = 1,3
            Xsbsun(i) = 0._10
 
c change this for observing subradar point instead of mass center
         end do
      endif
c
c determine reflection to receive time delay, decide
c whether to re-iterate for reflection time
      do i = 1,3
         Xsitep(i,1) = Xscsun(i,1) - Xsbsun(i)
      end do
      Rsitp(1) = SQRT(Xsitep(1,1)**2 + Xsitep(2,1)**2 + Xsitep(3,1)**2)
      Tmdly2   = Rsitp(1) - Xpdly
      iter1    = iter1 + 1
      if(ABS(Tmdly2-Tmdly1).gt.acctim) then
         if(iter1.gt.100) call SUICID(
     .'MORE THAN 100 REFLECTION TIME ITERATIONS, STOP IN STDLDP',14)
         Tmdly1 = Tmdly2
         goto 200
      endif
      Nit(20) = Nit(20) + iter1
      Dstf(4) = Tmdly2
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           iteration to calculate position of observing satellite
c           at send time
c
c           determine send time
      iter2  = 0
      Tmdly1 = 2._10*Tmdly2
      do while(.true.)
         call TIMINC(Jd,ctrecf,Jdy,fr3,-Tmdly1/Secday)
         call SCTRP(1,Jdy,fr3,0,0,2)
         if(Jdy.le.0) then
            Jd = 0
            return
         endif
         do i = 1,3
            Xscsun(i,2) = Xsc(i,2)*Aultsc
         end do
c
c determine send to reflection time delay, decide
c whether to re-iterate for send time
         do i = 1,3
            Xsitep(i,2) = Xscsun(i,2) - Xsbsun(i)
         end do
         Rsitp(2) = SQRT(Xsitep(1,2)**2 + Xsitep(2,2)**2 + Xsitep(3,
     .              2)**2)
         tmdly    = Tmdly2 + Rsitp(2) - Xpdly
         iter2    = iter2 + 1
         if(ABS(tmdly-Tmdly1).le.acctim) goto 410
         if(iter2.gt.100) call SUICID(
     .     ' MORE THAN 100 SEND TIME ITERATIONS, STOP IN STDLDP ',13)
         Tmdly1 = tmdly
      end do
 
  410 Nit(19) = Nit(19) + iter2
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c test for dummy observation occulted by central body
      if(Idumob.eq.1) then
         if(Ncs1.eq.Nplnt0) then
            if(Nspot.le.0) go to 420
         endif
         call CNTOCC(Xsbsun,Xsitep,2,Jd)
         if(Jd.le.0) return
      endif
c
c calculate unit vectors
  420 do j = 1,3
         do i = 1,2
            Xsitp0(j,i) = Xsitep(j,i)/Rsitp(i)
         end do
      end do
c
c calculate spot partials and velocity
      if(Nspot.gt.0) call SPOTCV(nvel,Nplnt0,1)
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
c calculate velocity of observing satellite at receive
c and send times
         call SCTRP(1,Jd,ctrecf,1,0,1)
         call SCTRP(1,Jdy,fr3,1,0,2)
         do k = 1,2
            do i = 4,6
               Xscsun(i,k) = Xsc(i,k)*Aultvl
            end do
         end do
c
c see if central body observed
         if(Ncs1.ne.Nplnt0) then
c
c see if planet or moon is observed
            if(Klanb.le.0) call SUICID(
     .' RADAR OBSERVATIONS OF PLANETS FROM BODY NOT THE EARTH ARE NOT CO
     .DED, STOP IN STDLDP',21)
c
c calculate velocity of observed satellite at reflection time
            call SBTRP(1,Jdx,fr2,1,0)
            do i = 4,6
               Xsbsun(i) = Xsb(i)*Aultvl
            end do
         else if(Nspot.gt.0) then
            do i = 4,6
               Xsbsun(i) = Xspcd(i,1)
            end do
         else
            do i = 4,6
               Xsbsun(i) = 0._10
            end do
         endif
c
c this routine could use the same code as the radar link for
c computing dopler:  vlrtrd + radrel
c calculate relative velocities
         do j = 1,2
            do i = 4,6
               Xsitep(i,j) = Xscsun(i,j) - Xsbsun(i)
            end do
         end do
c
c classical calculation of doppler shift
         do i = 1,2
            Quan1(i) = DOT(Xsitep(4,i),Xsitp0(1,i))
         end do
         Frect = DOT(Xsitp0(1,1),Xsbsun(4))
         Fract = DOT(Xsitp0(1,2),Xscsun(4,2))
         dop   = ((Quan1(2)*Fract-Quan1(1)*Frect) + Quan1(1)
     .           *Quan1(2)) - Quan1(1) - Quan1(2)
c
c multiply by frequency
         dop = Freq*dop
c
c constant bias in doppler shift
         if(Nrbias.gt.0) dop = dop + Rbsx(2)
      endif
 
      return
      end
