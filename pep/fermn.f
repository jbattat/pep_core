      subroutine FERMN(nvel,norm,npath,kspt)
 
      implicit none

c arguments
      integer*4 nvel,norm,npath,kspt
c
c           r. king   june 1977  subroutine fermn
c           calculate interferometry observables for spots on moon or
c           artificial space probes in cislunar space with central body
c           either the earth or moon.
c           this is a revised version of old pep subroutine mnterf,
c           written originally by m. ash and r. preston, feb 1970, and
c           modified dec 1971 by r.king to calculate differential
c           n-count (quasi-vlbi) observable - either accumulated
c           differential cycle count (difnct) or differential delay
c           increment (ddi)
c           fermn is analogous to subroutine radmn
c
c
c array dimensions
      include 'globdefs.inc'
 
c common 
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 ctat2,dfract
      equivalence (ctat2,Angdum(9)),(dfract,Angdum(10))
      include 'difnct.inc'
      real*10 tc
      equivalence (tc,Cnttim)
      include 'fcntrl.inc'
      include 'mnsprt.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'nutprc.inc'
      real*4 dpsi(2),deps(2)
      equivalence (dpsi,Nutat),(deps,Nutat(3))
      include 'obscrd.inc'
      real*10 ctrecf,reflct,fdev
      equivalence (ctrecf,Dstf),(reflct,Dstf(4)),(fdev,Dstf(10))
      real*4 delit
      equivalence (delit,Estf)
      integer*2 nintrf, nddiff, nlibpr, ntime
      equivalence (nintrf,Kob(12)),(nddiff,Kob(13)),
     .            (nlibpr,Kob(14)),(ntime,Kob(18))
      include 'param.inc'
      include 'radcrd.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'

c external functions
      real*10 CTATF,DOT
c
c quantities internal to this routine
      real*10 xemrf(6,2)
      integer*4 etide1
      real*10 ctrecv,dfdly1,dfdly2,dum,fractr,fresn
      integer*4 i,icount,idum,j,j1,jdesn,klb,lemctl,lmnctl,lsbctl,
     . lscctl,mnspt1,mtide1,n,n2,nn,nspt,nswesn,numsit,nvlesn
c
c           if two objects are observed, coordinates for the second are
c           determined first and saved in /radcrd/ and /difnct/
c
c-----------initialization----------------------------------------------
c
      lsbctl = 1
      lscctl = 1
      lemctl = 1
      lmnctl = 1
      call DGUESS(Tmdly1,1,Jds,Utrec)
      mnspt1 = 0
      etide1 = 0
      mtide1 = 0
      numsit = 2
      if(nddiff.lt.0) numsit = 1
      icount = 0
      nswesn = 0
      if(Ict(27).ne.0) nswesn = 1
      nvlesn = nvel
      if(prmter(81).gt.0._10) nvlesn = 1
c if two objects are observed, coordinates of the second
c are determined first and stored in /radcrd/
c
      nspt = Nspot
      klb  = Klanb
      if(kspt.ne.1) then
         nspt = Nspot2
         klb  = Klans1
      endif
c
c
c-----------calculate earth, moon, sun vectors at receive time(s)-------
c
c           earth-moon barycenter w.r.t. sun  (stored in xem(i,j) in au)
      if(nswesn.gt.0 .or. Reldel.gt.0._10 .or. Jct(10).gt.0 .or.
     .    Jct(11).gt.0) then
c if calculations are geocentric with no relativity and no solid
c body tides, then the e-m barycenter and solar system barycenter
c postions are not needed.
         call EMTRP(1,Jd,ctrecf,0,lemctl,1)
         if(Jd.le.0) return
         lemctl = 0
         if(nvlesn.gt.0) call EMTRP(1,Jd,ctrecf,1,lemctl,1)
c
c moon w.r.t. earth   (stored in xm(1-6,1), in a.u.)
         call MNTRP(1,Jd,ctrecf,0,lmnctl,1)
         if(Jd.le.0) return
         lmnctl = -1
         if(nvlesn.gt.0) call MNTRP(1,Jd,ctrecf,1,lmnctl,1)
c
c sun w.r.t. ssbc
         if(nswesn.gt.0) then
            call SOLCNT(Jd,ctrecf,Xslcns,nvlesn)
c
c earth w.r.t. s.s.b.c
            do i = 1, 3
               Xemlsc(i,1) = (Xem(i,1) - Xm(i,1)*Mnfct)*Aultsc
               if(nvlesn.gt.0) Xemlsc(i + 3,1)
     .             = (Xem(i+3,1) - Xm(i+3,1)*Mnfct)*Aultvl
            end do
         endif
c
c for counted-cycle vlbi observable, calculate earth w.r.t.
c sun and sun w.r.t. s.s.b.c. for the receive
c time at the second site
         if(nintrf.ge.0 .and. numsit.ne.1) then
            call MNTRP(1,Jd,Ctrec2/Secday,0,lmnctl,2)
            if(nvlesn.gt.0) call MNTRP(1,Jd,Ctrec2/Secday,1,
     .          lmnctl, 2)
            call EMTRP(1,Jd,Ctrec2/Secday,0,lemctl,2)
            if(nvlesn.gt.0) call EMTRP(1,Jd,Ctrec2/Secday,1,
     .          lemctl, 2)
            do i = 1, 3
               Xemlsc(i,2) = (Xem(i,2) - Xm(i,2)*Mnfct)*Aultsc
               if(nvlesn.gt.0) Xemlsc(i + 3,2)
     .             = (Xem(i+3,2) - Xm(i+3,2)*Mnfct)*Aultvl
            end do
            if(nswesn.gt.0) then
               if(Ctrec2.ne.Ctrec)
     .             call SOLCNT(Jd,Ctrec2/Secday,Xslcns(1,3),nvlesn)
            endif
         endif
c
c correct site coordinates for solid body tides
         if(Jct(10).gt.0) then
            call ETIDE(etide1,nvel,1,kspt)
            if(nintrf.ge.0 .and. numsit.eq.2)
     .          call ETIDE(etide1,nvel,2,kspt)
         endif
      endif
c
c
c-----------iteration to determine send time from observed noise source-
c-----------given receive time at first site ( n = 1 )------------------
c
      icount = 0
      n = 1
      ctrecv = Ctrec
  100 call TIMINC(Jd,ctrecv/Secday,Jdx,fractr,-Tmdly1/Secday)
c
c determine if space probe is observed body (either with earth
c or moon as central body)
      if(klb.gt.0) then
c
c obtain probe position at send time
         if(kspt.eq.2) then
            call SCTRP(idum,Jdx,fractr,0,lscctl,n)
            if(Jdx.le.0) then
c     since npath always = 1 for conventional vlbi observable,
c     store diff. delay rate in second position of difdly array.
c
c
c
               Jd = 0
               return
            else
               lscctl = -1
               do i = 1, 3
                  Xsbsun(i) = Xsc(i,n)*Aultsc
               end do
            endif
         else
            call SBTRP(1,Jdx,fractr,0,lsbctl)
            if(Jdx.le.0) then
               Jd = 0
               return
            else
               lsbctl = -1
               do i = 1, 3
                  Xsbsun(i) = Xsb(i)*Aultsc
               end do
            endif
         endif
c
c test to see if moon is to be interpolated for a moon
c centered probe (lunar orbiter) or an earth centered probe
c which hits the moon (ranger)
         if(Npcent(klb).eq.3) then
            do i = 1, 3
               Xsitep(i,n) = Xsbsun(i)
            end do
            goto 200
         endif
      endif
c
c obtain moon position at transmit time
      call MNTRP(1,Jdx,fractr,0,lmnctl,1)
      if(Jdx.le.0) then
         Jd = 0
         return
      else
         do i = 1, 3
            Xsitep(i,n) = Xm(i,1)*Mnltsc
         end do
         lmnctl = -1
c
c determine position of spot on surface of moon
         if(nspt.gt.0) then
            call MNSPT(Jdx,fractr,mnspt1,nvel,kspt,nlibpr)
            if(Jdx.le.0) then
               Jd = 0
               return
            else
c
c correct spot coordinates for solid body tides on moon
               if(Jct(11).gt.0) then
                  call MNSPTV(nvel,kspt,nlibpr)
                  call MTIDE(mtide1,nvel,kspt)
                  mtide1 = 1
               endif
c
c increment xsitep by spot coordinates
               do i = 1, 3
                  Xsitep(i,n) = Xsitep(i,n) + Xspcd(i,kspt)
               end do
            endif
c
c increment xsitep if observed body is moon centered probe
         else if(klb.gt.0) then
            do i = 1, 3
               Xsitep(i,n) = Xsitep(i,n) + Xsbsun(i)
            end do
         endif
      endif
c
c form vector pointing from source to site
  200 do i = 1, 3
         Xsitep(i,n) = Xsite(i,n) - Xsitep(i,n)
      end do
      if(nswesn.gt.0) then
         n2 = 2*n
c moon coordinates not previously calculated if nplnt0 is
c earth satellite
         if(klb.ne.0 .and. Npcent(klb).eq.3) then
            call MNTRP(1,Jdx,fractr,0,lmnctl,1)
            if(nvlesn.gt.0) call MNTRP(1,Jdx,fractr,1,lmnctl,
     .          1)
         endif
         call EMTRP(1,Jdx,fractr,0,lemctl,2)
         if(nvlesn.gt.0) call EMTRP(1,Jdx,fractr,1,lemctl,2)
         call SOLCNT(Jdx,fractr,Xslcns(1,n2),nvlesn)
 
c save earth w.r.t. sun at reflect for relativity
         do i = 1, 3
            xemrf(i,n) = (Xem(i,2) - Xm(i,1)*Mnfct)*Aultsc
            if(nvlesn.gt.0) xemrf(i + 3,n)
     .          = (Xem(i+3,2) - Xm(i+3,1)*Mnfct)*Aultvl
         end do
         do i = 1, 3
            Xsitep(i,n) = Xsitep(i,n) + (Xemlsc(i,n) - xemrf(i,n))
     .                     - (Xslcns(i,n) - Xslcns(i,n2))*Aultsc
         end do
      endif
c
c determination of send to receive time delay, decide
c whether to re-iterate for source send time
      call UVECTR(3,Xsitep(1,n),Rsitp(n),Xsitp0(1,n),dum)
      Tmdly2 = Rsitp(n)
      icount = icount + 1
      if(ABS(Tmdly2-Tmdly1).gt.delit) then
         Tmdly1 = Tmdly2
         if(icount.le.15) goto 100
         call SUICID('MORE THAN 15 DELAY ITERATIONS,STOP IN FERMN ',11)
      endif
      Nit(n + 18) = Nit(n + 18) + icount
      icount = 0
 
c save tmdly1 for initial guess in next observation
      Dstf(4) = Tmdly1
c
c            end of iteration loop for first site position
c           (also second site position if counted-cycle vlbi)
c           for counted-cycle observable determine send time
c           given receive time at second site  ( n = 2 )
      if(nintrf.lt.0) then
c
c-----------------------------------------------------------------------
c
c           for conventional vlbi observable  determine receive time
c           at second site given  source send time
         reflct = Tmdly2
         dfdly1 = 0.0_10
      else
         if(n.eq.2 .or. numsit.eq.1) goto 500
         n = 2
         ctrecv = Ctrec2
         goto 100
      endif
  300 Sidtm2 = Sidtm0 + Sidvel*(Utrec-dfdly1) + Dgst
      call SITCOR(Sidtm2,2,nvlesn,norm)
      if(nswesn.gt.0 .or. prmter(81).gt.0 .or. Jct(10).gt.0 .or.
     .    Jct(11).gt.0) then
         call TIMINC(Jd,ctrecf,jdesn,fresn,-dfdly1/Secday)
         call EMTRP(1,jdesn,fresn,0,lemctl,2)
         if(nvlesn.gt.0) call EMTRP(1,jdesn,fresn,1,lemctl,2)
         call MNTRP(1,jdesn,fresn,0,lmnctl,2)
         if(nvlesn.gt.0) call MNTRP(1,jdesn,fresn,1,lmnctl,2)
 
         if(nswesn.gt.0) then
            call SOLCNT(jdesn,fresn,Xslcns(1,2),nvlesn)
            do i = 1, 3
               Xemlsc(i,2) = (Xem(i,2) - Xm(i,2)*Mnfct)*Aultsc
            end do
         endif
c
c correct second site coordinates for solid body tides
         if(Jct(10).gt.0) call ETIDE(etide1,nvel,2,kspt)
         etide1 = 1
      endif
 
      if(klb.gt.0) then
         if(Npcent(klb).eq.3) then
            do i = 1, 3
               Xsitep(i,2) = Xsbsun(i)
            end do
            goto 400
         endif
      endif
      do i = 1, 3
         Xsitep(i,2) = Xm(i,1)*Mnltsc
         if(nspt.gt.0) Xsitep(i,2) = Xsitep(i,2) - Xspcd(i,kspt)
         if(klb.gt.0) Xsitep(i,2)  = Xsitep(i,2) - Xsb(i)
      end do
  400 do i = 1, 3
         Xsitep(i,2) = Xsite(i,2) - Xsitep(i,2)
         if(nswesn.gt.0) Xsitep(i,2) = Xsitep(i,2)
     .       + (Xemlsc(1,2) - xemrf(i,1)) - (Xslcns(i,3) - Xslcns(i,2))
     .       *Aultsc
      end do
 
      call UVECTR(3,Xsitep(1,2),Rsitp(2),Xsitp0(1,2),dum)
      icount = icount + 1
      dfdly2 = Tmdly2 - Rsitp(2)
      if(ABS(dfdly2-dfdly1).gt.delit) then
         dfdly1 = dfdly2
         if(icount.le.15) goto 300
         call SUICID('MORE THAN 15 DELAY ITERATIONS, STOP IN FERMN',11)
      endif
      Difdly(1,kspt) = dfdly2
      Nit(20) = Nit(20) + icount
c
c
c           end of iteration loop for conventional vlbi
c           second site postion
c
c---------test for dummy observation below horizon or s/c---------------
c           occulted by moon
c
  500 if(Idumob.eq.1) then
         nn = numsit
         if(klb.gt.0 .and. Npcent(klb).eq.10) nn = -numsit
         call HORIZN(Xsite,Xsitep,nn,Jd)
         if(klb.ne.0 .and. Npcent(klb).eq.10)
     .       call CNTOCC(Xsbsun,Xsitep,-numsit,Jd)
         if(Jd.le.0) return
      endif
c
c
c-----------calculate velocities, and partials----------
c
c
c           determine lunar spot velocity and partials
      if(nspt.gt.0) then
 
c already determined if solid body tides were calculated
         if(Jct(11).le.0) call MNSPTV(nvel,kspt,nlibpr)
      endif
c
c determine if velocities are needed
      if(nvel.gt.0) then
c note:  send times assumed same for both sites (cf. fersb)
c
c determine probe velocity relative to central body
         if(klb.gt.0) then
            if(kspt.eq.2) then
               call SCTRP(idum,Jdx,fractr,1,lscctl,1)
 
c note: puts coords. in xsb(i,1) for either site
               do i = 4, 6
                  Xsbsun(i) = Xsc(i,1)*Aultvl
               end do
            else
               call SBTRP(1,Jdx,fractr,1,lsbctl)
               do i = 4, 6
                  Xsbsun(i) = Xsb(i)*Aultvl
               end do
            endif
c
c see if this is earth satellite rather than lunar orbiter
            if(Npcent(klb).eq.3) then
               do j = 1, numsit
                  do i = 4, 6
                     Xsitep(i,j) = Xsbsun(i)
                  end do
               end do
               goto 550
            endif
         endif
c
c determine velocity of spot on moon relative to earth
         if(nspt.gt.0) then
            do j = 1, numsit
               do i = 4, 6
 
c error:  xm(i,1) is moon pos. at transmit time for site 2 only
                  Xsitep(i,j) = Xm(i,1)*Mnltsc/Secday + Xspcd(i,kspt)
               end do
            end do
c
c determine lunar orbiter velocity relative to earth
         else if(klb.gt.0) then
            do j = 1, numsit
               do i = 4, 6
                  Xsitep(i,j) = Xm(i,1)*Mnltsc/Secday + Xsbsun(i)
               end do
            end do
         endif
c
c determine velocity of observing site relative to
c observed thing (moon,spot,probe)
  550    do j = 1, numsit
            do i = 4, 6
               Xsitep(i,j) = Xsite(i,j) - Xsitep(i,j)
            end do
         end do
c
c correct for velocity of site and source w.r.t. s.s.b.c.
         if(nswesn.gt.0) then
            do j = 1, numsit
               j1 = 1
               if(j.eq.2) j1 = 3
               do i = 4, 6
                  Xsitep(i,j) = Xsitep(i,j)
     .                           + (Xemlsc(i,j) - xemrf(i,j))
     .                           - (Xslcns(i,j1) - Xslcns(i,2*j))*Aultvl
               end do
            end do
         endif
c
c calculate beta (=v/c)
         do i = 1, numsit
            Beta(i,kspt) = DOT(Xsitep(4,i),Xsitp0(1,i))
         end do
         if(Nice.gt.0) then
c
c---------------compute differential delay rate-------------------------
c
            Difdly(2,kspt) = Beta(1,kspt) - Beta(2,kspt)
            return
         endif
      endif
c
c
c-----------corrections to differential delay---------------------------
c
c
c
c          convert atomic time delay to utc
      if(ntime.gt.0) then
         do i = 1, numsit
            Rsitp(i) = Rsitp(i)*fdev
         end do
      endif
c
c calculate differential delay in atomic time
c ctat2 for counted-cycle observable calculated in ferctl
      if(nintrf.lt.0) ctat2 = CTATF(jdesn,fresn - Ctat/Secday,3,
     .                            2)
      if(numsit.eq.2) Difdly(npath,kspt) = Rsitp(1) - Rsitp(2)
     .    - Ctat + ctat2
c     npath=1 for conventional interferometry observable
c
c      clock terms added in ferctl
c
c           save quantities for partials and differential
c           n-count calculations
      do i = 1, numsit
         Rsave(npath,i,kspt) = Rsitp(i)
      end do
      if(kspt.eq.2) then
         do i = 1, numsit
            do j = 1, 3
               Ysitep(j,i) = Xsitep(j,i)
               if(nvel.gt.0) Ysitep(j + 3,i) = Xsitep(j + 3,i)
               Ysitp0(j,i) = Xsitp0(j,i)
            end do
            Rsitp2(i) = Rsitp(i)
         end do
      endif
c calculate beta from phase deleay differences if velocities
c not used (applies to uncondensed alsep ddr observations)
      if(nintrf.lt.0) then
         Difdly(2,kspt) = Beta(1,kspt) - Beta(2,kspt)
      else if(nvel.le.0) then
         if(npath.eq.2) then
            do i = 1, numsit
               Beta(i,kspt) = (Rsave(2,i,kspt) - Rsave(1,i,kspt))/tc
            end do
         endif
      endif
      return
      end
