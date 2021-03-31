      subroutine RADMN(nvel, norm, npath, idopob)
 
      implicit none

c     king  june 1977  subroutine radmn
c     cleaned up by j.f.chandler - 1978 june
c     calculation of time delay and doppler shift for observations
c     of the moon or of artificial space probes in cislunar space with
c     central body either the earth or moon or of spots on the moon.
c     this is a revised version of old pep subroutine mdeldp written
c     written and modified by ash/slade/becker/king  oct 67 - jun 75
c        ncodf=1 radar or laser ranging observable
c        ncodf=2 radio tracking observable
c        ncodf=19 doppler count observable for les-8/9
c
c arguments
      integer*4 nvel,norm,npath,idopob

c array dimensions
      include 'globdefs.inc'

c        common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 dutrec
      equivalence (Angdum(10), dutrec)
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'ltrapobs.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1), tmdly), (Deriv(2,2), dop)
      include 'mnsprt.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct,tmdly0,dop0
      equivalence (ctrecf, Dstf), (fdev, Dstf(10)), (reflct, Dstf(4))
      equivalence (Result, tmdly0), (Result(2), dop0)
      real*4 acctim
      equivalence (acctim, Estf)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'plnhar.inc'
      include 'prpgat.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'

c local variables 
      real*10 bdotx1,bdotx2,bdotxm,uem,xemrf(6),mnrf(6)
      real*10 dum,fresn,sendtm,sidtm1
      integer*4 etide1,i,ipct,iprsp,j,jpct,kspt,lemctl,lmnctl,
     .          lsbctl,mnspt1,mtide1,nn,nvlem
      character*4 amper9/'&&&&'/

c external functions
      real*10 ANCOR, CTATF, DOT

c-------------initialization--------------------------------------------
c
c receive time site coordinates were determined in radctl
      sidtm1 = Sidtm
 
      lsbctl = -1
      lmnctl = -1
      lemctl = 0
      mnspt1 = 0
      kspt   = 1
      etide1 = 0
      mtide1 = 0
      call ZFILL(Xemlsc, 16*12)
      call ZFILL(xemrf, 16*6)
c Earth vel. not needed if no g.r., no tides, and
c no correction to solar system barycenter.
c However, Moon vel. is needed if there is Doppler, so do it anyhow.
      nvlem = nvel
c
c set up probe/spot logic:  iprsp= -1,0,1 if spot, neither, probe
      iprsp = 0
      if(Spotf.eq.amper9 .or. Nspot.gt.0) iprsp = -1
      if(Klanb.gt.0) iprsp = 1
c
c-----------calculate earth, moon, sun vectors at receive time(s)-------
c
c           earth-moon barycenter w.r.t. sun  (stored in xem(i,j) in au)
c*  start=1000
      if(Nswcns.gt.0 .or. ntmdly.gt.0 .or. ndop.gt.0 .or.
     .    Jct(10).gt.0 .or. Jct(11).gt.0 .or. Jct(61).ge.0) then
c if calculations are geocentric with no relativity and no solid
c body tides, then the e-m barycenter and solar system barycenter
c postions are not needed.
         call ETRP(1, Jd, ctrecf, 0, lemctl, 1, 2)
         if(Jd.le.0) return
         if(nvlem.gt.0) call ETRP(1, Jd, ctrecf, 1, 0, 1, 2)
c
c sun w.r.t. ssbc
         if(Nswcns.gt.0) call SOTRP(Jd, ctrecf, Xslcns, 0)
c
c correct site coordinates for solid body tides
         if(Jct(10).gt.0) call ETIDE(etide1, nvel, 1, kspt)
c
c correct site coordinates for fluid displacements
         if(Jct(49).gt.0) call FLURED(Jd,ctrecf,1)
      endif
c
c-----------iteration to determine send time from observed object-------
c-----------given receive time at first site----------------------------
c*  start=1100
c
c           determine reflection time
  100 call TIMINC(Jd, ctrecf, Jdx, Fract, -Tmdly1/Secday)
c
c obtain probe position at reflection time
      if(Klanb.gt.0) then
         call SBTRP(1, Jdx, Fract, 0, lsbctl)
         if(Jdx.le.0) then
c
c*  start=9910
            Jd = 0
            return
         endif
         do j = 1, 3
            Xsbsun(j) = Xsb(j)*Aultsc
            Xplsc(j)  = Xsbsun(j)
         end do
         if(Ncp0.eq.3) goto 200
      endif
c
c obtain moon position at reflection time
      call MNTRP(1, Jdx, Fract, 0, lmnctl, 1)
      if(Jdx.le.0) then
         Jd = 0
         return
      endif
      do j = 1, 3
         Xplsc(j) = Xm(j, 1)*Mnltsc
      end do
c
c determine position of spot on surface of moon
c * * * warning * * *
c velocities not recalculated for moon at reflect time
      if(iprsp.lt.0) then
         call MNSPT(Jdx, Fract, mnspt1, nvel, kspt, nlibpr)
         if(Jdx.le.0) then
            Jd = 0
            return
         endif
         if(Jct(11).gt.0) then
            call MNSPTV(nvel, kspt, nlibpr)
            call MTIDE(mtide1, nvel, 1)
            mtide1 = 1
         endif
         do j = 1, 3
            Xsbsun(j) = Xspcd(j, 1)
         end do
      endif
c
c increment xsitep for probe or spot
      if(iprsp.ne.0) then
         do i = 1, 3
            Xplsc(i) = Xplsc(i) + Xsbsun(i)
         end do
      endif
c
c determination of reflection to receive time delay, decide
c whether to re-iterate for reflection time
c*  start=1400
  200 do j = 1, 3
         Xsitep(j,1) = Xsite(j,1) - Xplsc(j)
      end do
      if(Nswcns.gt.0 .or. ntmdly.gt.0 .or. ndop.gt.0) then
c moon coordinates not previously calculated if nplnt0 is
c earth satellite
         if(Klanb.ne.0 .and. Ncp0.eq.3)
     .       call MNTRP(1, Jdx, Fract, 0, lmnctl, 1)
         if(Jdx.le.0) then
            Jd = 0
            return
         endif
         call EMTRP(1, Jdx, Fract, 0, lemctl, 2)
         if(Jdx.le.0) then
            Jd = 0
            return
         endif
c save earth w.r.t. sun at reflect for relativity
c and moon velocity for lorentz contraction
         do i=1,3
            xemrf(i) = (Xem(i,2) - Xm(i,1)*Mnfct)*Aultsc
            mnrf(i)  = xemrf(i) + Xm(i,1)*Mnltsc
         end do
         if(nvlem.gt.0) then
            call MNTRP(1, Jdx, Fract, 1, lmnctl, 1)
            call EMTRP(1, Jdx, Fract, 1, lemctl, 2)
            do i=4,6
               xemrf(i) = (Xem(i,2) - Xm(i,1)*Mnfct)*Aultvl
               mnrf(i)  = xemrf(i) + Xm(i,1)*Mnau*Aultvl
            end do
         endif

c convert to solar-system barycentric frame, including lorentz contraction
         if(Nswcns.gt.0) then
            call SOTRP(Jdx, Fract, Xslcns(1,2), 0)
            bdotx1=0._10
            uem=0._10
            if(nvlem.gt.0 .and. Jct(17).eq.0) then
               bdotx1=DOT(Xsite(1,1),Xemlsc(4,1))
               if(iprsp.lt.0) then
                  uem = Gmc2/SQRT(DOT(mnrf,mnrf))
                  bdotxm = DOT(Xspcd(1,1),mnrf(4))
                  do i=1,3
                     Xplsc(i)    = Xplsc(i) - 0.5_10*bdotxm*mnrf(i+3)
     .                             - uem*Xspcd(i,1)
                     Xsitep(i,1) = Xsite(i,1) - Xplsc(i)
                  end do
               endif
               uem = Gmc2/SQRT(DOT(Xemlsc,Xemlsc))
            endif
            do i = 1,3
               Xsitep(i,1) = Xsitep(i,1)
     .          + (Xemlsc(i,1) - xemrf(i))
     .          - (Xslcns(i,1) - Xslcns(i,2))*Aultsc
     .          - 0.5_10*bdotx1*Xemlsc(i+3,1) - uem*Xsite(i,1)
            end do
c            write(6,1111) (Xemlsc(i,1),i=1,6),(Xsite(i,1),i=1,3),
c     .       (-0.5_10*bdotx1*Xemlsc(i+3,1),i=1,3),
c     .       (-uem*Xsite(i,1),i=1,3),

         endif
      endif
 
      call UVECTR(3, Xsitep(1,1), Rsitp(1), Xsitp0(1,1), dum)
      Tmdly2  = Rsitp(1) - Xpdly
      Nit(20) = Nit(20) + 1
      if(ABS(Tmdly2-Tmdly1).le.acctim) goto 220
      Tmdly1 = Tmdly2
      goto 100
c
c-----------------------------------------------------------------------
c
c           iteration to determine send time site coordinates
c*  start=2000
  220 reflct = Tmdly2
      Tguess = Tmdly2
      call TSAV(Tmdly2, 1, Jd, Utrec)
 
      call DGUESS(tmdly, 3, Jd, Utrec)
      tmdly = tmdly + Tmdly2
c
c is this one-way les-8/9 doppler count observable
      if(Ncodf.ne.19) then
         do while(.true.)
c
c precession-nutation not evaluated at send time for cislunar observ
c
c start iteration loop
            Tmdly1 = tmdly
            sendtm = Utrec - Tmdly1*fdev
            Sidtm  = Sidtm0 + Sidvel*sendtm + Dgst
            call SITCOR(Sidtm, 2, nvel, norm)
            if(Nswcns.gt.0 .or. ntmdly.gt.0 .or. ndop.gt.0 .or.
     .         Jct(10).gt.0 .or. Jct(11).gt.0 .or. Jct(61).ge.0) then
               call TIMINC(Jd, ctrecf, Jdy, fresn, -Tmdly1/Secday)
               call ETRP(1, Jdy, fresn, 0, lemctl, 2, 2)
               if(Jdy.le.0) then
                  Jd = 0
                  return
               endif
 
c need earth velocity for tide model
               if(nvlem.gt.0) call ETRP(1, Jdy, fresn, 1, 0, 2, 2)
               if(Nswcns.gt.0) call SOTRP(Jdy, fresn, Xslcns(1,3), 0)
c
c correct second site coordinates for solid body tides
               if(Jct(10).gt.0) call ETIDE(etide1, nvel, 2, kspt)
               etide1 = 1
c
c correct second site coordinates for fluid displacements
               if(Jct(49).gt.0) call FLURED(Jdy,fresn,2)
            endif
 
            bdotx2=0._10
            if(nvlem.gt.0 .and. Jct(17).eq.0) then
               bdotx2=DOT(Xsite(1,2),Xemlsc(4,2))
               uem = Gmc2/SQRT(DOT(Xemlsc(1,2),Xemlsc(1,2)))
            endif
            do i = 1, 3
               Xsitep(i,2) = Xsite(i,2) - Xplsc(i)
               if(Nswcns.gt.0) Xsitep(i,2) = Xsitep(i,2)
     .             + (Xemlsc(i,2) - xemrf(i))
     .             - (Xslcns(i,3) - Xslcns(i,2))*Aultsc
     .             - 0.5_10*bdotx2*Xemlsc(i+3,2) - uem*Xsite(i,2)
            end do
            call UVECTR(3, Xsitep(1,2), Rsitp(2), Xsitp0(1,2), dum)
            tmdly   = Rsitp(2) - Xpdly + Tmdly2
            Nit(19) = Nit(19) + 1
            if(ABS(tmdly-Tmdly1).le.acctim) goto 250
         end do
      endif
c
c-----------------------------------------------------------------------
c
c        calculations after iterations are completed
c
c        test for dummy observation below horizon
c*  start=2200
  250 call TSAV(tmdly - Tmdly2, 3, Jd, Utrec)
      if(Idumob.eq.1) then
         nn = 2
         if(Ncodf.eq.19) nn = 1
         if(Ncp0.eq.10) nn  = -nn
         call HORIZN(Xsite, Xsitep, nn, Jd)
c
c test for spot on moon occulted by moon as seen from sites
c test for lunar orbiter occulted by moon
         if(Ncp0.eq.10 .or. iprsp.lt.0)
     .    call CNTOCC(Xsbsun, Xsitep, -iabs(nn), Jd)
c
c test for sun above the horizon at observing site(s)
         call SUNHOR(Xsite, Xemlsc, nn, Jd)
c
c test for sun above the horizon at observed spot
         if(iprsp.lt.0) call SUNHRSP(Xspcd, mnrf, Jd)
         if(Jd.le.0) return
      endif
c
c determine lunar spot velocity and partials
c*  start=2500
      if(iprsp.lt.0) then
 
c already determined if solid body tides were calculated
         if(Jct(11).le.0) call MNSPTV(nvel, kspt, nlibpr)
         do j = 4, 6
            Xsbsun(j) = Xspcd(j, 1)
         end do
      endif
c
c determine if velocities are needed
      if(nvel.le.0) goto 400
 
      if(Jct(10).le.0 .and. nvlem.gt.0)
     .    call ETRP(1, Jdy, fresn, 1, 0, 2, 2)
 
c obtain probe velocity relative to central body
      if(Klanb.gt.0) then
         call SBTRP(1, Jdx, Fract, 1, lsbctl)
         do j = 4, 6
            Xsbsun(j) = Xsb(j)*Aultvl
            Xplsc(j)  = Xsbsun(j)
         end do
         if(Ncp0.eq.3) goto 300
      endif
c
c obtain moon velocity relative to earth
      do j = 4, 6
         Xplsc(j) = Xm(j, 1)*Mnltsc/Secday
      end do
c
c increment xplsc for probe or spot
      if(iprsp.ne.0) then
         do j = 4, 6
            Xplsc(j) = Xplsc(j) + Xsbsun(j)
         end do
      endif
c
c        determine velocity of observing site relative to
c        observed thing (moon,spot,probe)
c
c        obtain velocities in sun-centered frame
c        (xemlsc and xemrf are zero if earth-centered frame used)
  300 do j = 4, 6
         Xsbsun(j) = xemrf(j) + Xplsc(j)
         do i = 1, 2
            Xemlsc(j, i) = Xemlsc(j, i) + Xsite(j, i)
         end do
      end do
      if(Nswcns.le.0) then
 
c use earth-centered frame
         call VLRTRD(Xsite(4,1), Xplsc(4), Xsite(4,2), 1._10, 7, 0)
      else
         call VLRTRD(Xemlsc(4,1), Xsbsun(4), Xemlsc(4,2), 1._10, 7, 0)
      endif
c
c write out efg tape
      if(Jct(30).ne.0) call EFGOUT(sidtm1)
c
c------------------correct time delay for various effects---------------
c*  start=3000
c
c           set up vectors for radrel
  400 do j = 1, 3
         Xsbsun(j) = xemrf(j) + Xplsc(j)
         do i = 1, 2
            Xemlsc(j, i) = Xemlsc(j, i) + Xsite(j, i)
         end do
      end do
      call RADREL(1)
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
         if(lopler.ne.-1) goto 500
         if(npath.eq.1) ipct = -2
         if(npath.eq.2) ipct = -1
      endif
      call PROPCO(ipct, 1)
      if(ipct.eq.1) tmdly = tmdly + Sumcor(1)
 
c apply antenna correction
      tmdly = tmdly + ANCOR(1) + ANCOR(2)
 
c convert coordinate time delay to atomic
      tmdly = tmdly + (CTATF(Jd,(Ctrec-tmdly-Ctat)/Secday,Ntrmct,2)
     .        - Ctat)
 
c time delay in atomic time changed to time delay in utc time
      if(ntime.gt.0) tmdly = tmdly*fdev
c
c constant bias in time delay
      if(Nrbias.gt.0) tmdly = tmdly + Rbsx(1)
 
c power series for clock errors uses /eqenox/ quantities
      if(Neqnox.gt.0) tmdly = tmdly + Pnox +
     .                            dutrec*(Pequat + dutrec*Plat/2._10)
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
            call PHADOP(npath, idopob)
            if(npath.ne.1) return
 
c now get corrections for phase delay doppler
            jpct = 0
            goto 600
         endif
      endif
c
c-----------compute doppler shift for instantaneous frequency-----------
c
  500 if(Nice.lt.0) return
      call RADREL(2)
      jpct = 2
c
c propagation corrections
  600 call PROPCO(jpct, 1)
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
