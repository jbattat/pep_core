      subroutine RADGPS(nvel,norm,npath,idopob)
 
      implicit none

c     r.w.king  oct 1981  subroutine radgps
c     calculation of one-way time delay and doppler shift for
c     observations of earth satellites with speed of light assumed
c     infinite (no falsi iteration required)

c arguments
      integer*4 nvel,norm,npath,idopob

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 dutrec
      equivalence (Angdum(10),dutrec)
      include 'eqnphs.inc'
      include 'fcntrl.inc'
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
      real*4    acctim
      equivalence (acctim,Estf)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'plnhar.inc'
      include 'prpgat.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'
c
c quantities internal to this routine
      real*10 xemrf(6),dum
      integer*4 etide1
      integer*4 i,ipct,j,jpct,kspt,lemctl,lmnctl,lsbctl,nn,nvlem

c external functions
      real*10 ANCOR,CTATF
c
c
c-------------initialization--------------------------------------------
c
      lsbctl = -1
      lmnctl = -1
      lemctl = 0
      kspt   = 1
      etide1 = 0
c earth vel. not needed if no g.r., no tides, and
c no correction to solar system barycenter
      nvlem = nvel
      if(ndop.le.0 .and. Nswcns.le.0 .and. Jct(10).le.0 .and.
     .    Jct(11).le.0) nvlem = 0
c
c
c-----------calculate earth, moon, sun vectors at observation time------
c
c           earth-moon barycenter w.r.t. sun  (stored in xem(i,j) in au)
c*  start=1000
      if(Nswcns.gt.0 .or. ntmdly.gt.0 .or. ndop.gt.0 .or.
     .    Jct(10).gt.0 .or. Jct(11).gt.0) then
c if calculations are geocentric with no relativity and no solid
c body tides, then the e-m barycenter and solar system barycenter
c postions are not needed.
         call ETRP(1,Jd,ctrecf,0,lemctl,1,2)
         if(Jd.le.0) return
         if(nvlem.gt.0) call ETRP(1,Jd,ctrecf,1,0,1,2)
c
c sun w.r.t. ssbc
         if(Nswcns.gt.0) call SOTRP(Jd,ctrecf,Xslcns,0)
c
c correct site coordinates for solid body tides
         if(Jct(10).gt.0) call ETIDE(etide1,nvel,1,kspt)
      endif
c
c
c-----------calculate satellite position at observation time------------
c
c
c     transmit time same as receive time
      Jdx   = Jd
      Fract = ctrecf
 
      call SBTRP(1,Jdx,Fract,0,lsbctl)
      if(Jdx.le.0) then
c
c*  start=9910
         Jd = 0
         return
      else
         do j = 1, 3
            Xsbsun(j) = Xsb(j)*Aultsc
            Xplsc(j)  = Xsbsun(j)
         end do
c
c determination of range
c*  start=1400
         do j = 1, 3
            Xsitep(j,1) = Xsite(j,1) - Xplsc(j)
         end do
         if(Nswcns.gt.0 .or. ntmdly.gt.0 .or. ndop.gt.0) then
c moon coordinates not previously calculated if nplnt0 is
c earth satellite
            if(Klanb.ne.0 .and. Ncp0.eq.3)
     .          call MNTRP(1,Jdx,Fract,0,lmnctl,1)
            if(Jdx.le.0) then
               Jd = 0
               return
            else
               call EMTRP(1,Jdx,Fract,0,lemctl,1)
               if(Jdx.le.0) then
                  Jd = 0
                  return
               else
                  if(Nswcns.gt.0)
     .                call SOTRP(Jdx,Fract,Xslcns(1,2),0)
                  do i = 1, 3
 
c save earth w.r.t. sun at reflect for relativity
                     xemrf(i) = (Xem(i,1) - Xm(i,1)*Mnfct)*Aultsc
                     if(Nswcns.gt.0) Xsitep(i,1) = Xsitep(i,1)
     .                   + (Xemlsc(i,1) - xemrf(i))
     .                   - (Xslcns(i,1) - Xslcns(i,2))*Aultsc
                  end do
               endif
            endif
         endif
 
         call UVECTR(3,Xsitep(1,1),Rsitp(1),Xsitp0(1,1),dum)
         Tmdly2  = Rsitp(1) - Xpdly
         Nit(20) = Nit(20) + 1
         Tmdly1  = Tmdly2
         tmdly   = Tmdly2
         call TSAV(tmdly - Tmdly2,3,Jd,Utrec)
c
c-----------------------------------------------------------------------
c
c           calculations after iterations are completed
c
c           test for dummy observation below horizon
c*  start=2200
         if(Idumob.eq.1) then
            nn = 1
            call HORIZN(Xsite,Xsitep,nn,Jd)
c
c
            if(Jd.le.0) return
         endif
c
c
c determine if velocities are needed
         if(nvel.gt.0) then
c
c obtain probe velocity relative to central body
            call SBTRP(1,Jdx,Fract,1,lsbctl)
            do j = 4, 6
               Xsbsun(j) = Xsb(j)*Aultvl
               Xplsc(j)  = Xsbsun(j)
            end do
c
c determine velocity of observing site relative to probe
            if(Nswcns.gt.0 .or. ntmdly.gt.0 .or. ndop.gt.0)
     .          then
 
c obtain velocities in sun-centered frame
               do j = 4, 6
                  Xsbsun(j)    = xemrf(j) + Xplsc(j)
                  Xemlsc(j,1) = Xemlsc(j,1) + Xsite(j,1)
               end do
               if(Nswcns.gt.0) then
                  call VLRTRD(Xemlsc(4,1),Xsbsun(4),Xemlsc(4,2),1._10,
     .                        1,0)
                  goto 50
               endif
            endif
 
c use earth-centered frame
            call VLRTRD(Xsite(4,1),Xplsc(4),Xsite(4,2),1._10,7,0)
c
c
c------------------correct time delay for various effects---------------
c*  start=3000
         endif
c
c set up vectors for radrel
   50    do j = 1, 3
            Xsbsun(j)    = xemrf(j) + Xplsc(j)
            Xemlsc(j,1) = Xemlsc(j,1) + Xsite(j,1)
         end do
c     call radrel(1)
c     skip radrel call for now since doesn't work for nsite2=0
c
c        correct for propagation effects
c        set up ipct as follows:
c          -2: 1st leg of phase delay doppler
c          -1: 2nd leg of phase delay doppler
c           1: time delay
c
c        skip the corrections if doing only instantaneous doppler
c
         ipct = 1
         if(Nice.gt.0) then
            if(lopler.ne.-1) goto 100
            if(npath.eq.1) ipct = -2
            if(npath.eq.2) ipct = -1
         endif
         call PROPCO(ipct,1)
         if(ipct.eq.1) tmdly = tmdly + Sumcor(1)
 
c apply antenna correction
         tmdly = tmdly + ANCOR(1) + ANCOR(2)
 
c convert coordinate time delay to atomic
         tmdly = tmdly + (CTATF(Jd,(Ctrec-tmdly-Ctat)/Secday,Ntrmct,2)
     .           - Ctat)
 
c time delay in atomic time changed to time delay in utc time
         if(ntime.gt.0) tmdly = tmdly*fdev
c
c constant bias in time delay
         if(Nrbias.gt.0) tmdly = tmdly + Rbsx(1)
 
c power series for clock errors uses /eqenox/ quantities
         if(Neqnox.gt.0) tmdly = tmdly + Pnox +
     .                               dutrec*(Pequat + dutrec*Plat/2._10)
c
c doppler shift
c
c nice.gt.0 in time delay calculations means phase delay
         if(Nice.gt.0) then
c
c------------------doppler phase observable----------------------------
c*  start=4000
c
c    there are three possible doppler observables:
c          lopler=-1  observable is phase
c          lopler= 0  observable is averaged integrated frequency (jpl)
c                     (no longer implemented)
c          lopler= 1  observable is instantaneous frequency
c
            if(lopler.ne.1) then
               call PHADOP(npath,idopob)
               if(npath.ne.1) return
 
c now get corrections for phase delay doppler
               jpct = 0
               goto 200
            endif
         endif
      endif
c
c-----------compute doppler shift for instantaneous frequency-----------
c
  100 if(Nice.lt.0) return
      call RADREL(2)
      jpct = 2
c
c propagation corrections
  200 call PROPCO(jpct,1)
      dop = dop + Sumcor(2)
 
      dop = Freq*dop
c
c constant bias in doppler shift
      if(Nrbias.gt.0) dop = dop + Rbsx(2)
 
c power series for clock errors uses /eqenox/ quantities
      if(Neqnox.gt.0) dop = dop - Freq*(Pequat + dutrec*Plat)
 
c*  start=9990
      return
      end
