      subroutine FERSB(nvel,norm,npath,kspt)
 
      implicit none
 
 
c*** start of declarations inserted by spag
 
c*** end of declarations inserted by spag
 
c
c     r.king/s.gourevitch  june 1977   subroutine fersb
c     calculate vlbi observables for probes and spots not in cis-lunar
c     space. this is a revised version of old pep subroutine sbterf,
c     written by g.slater/r.cappallo/r.king  jul 72  - aug 74
c     this routine computes delay and delay rate vlbi observables
c     for probes and spots not in cis-lunar space.
c
c arguments
      integer*4 nvel,norm,npath,kspt
c
c           nvel= 0  positions only calculated
c               = 1  positon and velocity calculated
c           norm= 0  no site normals calculated
c               = 1  site normals calulated
c           npath=1  for a counted cycle observable, fersb is
c                    called at the beginning of the interval
c                =2  fersb is called at the end of the interval
c           kspt= 1  only one observed probe for fersb to calculate
c               = 2  two observed probes for fersb to calculate;
c                    do both on single call and return kspt=1
c

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 tpvmsb(2),tpvmsc(2),ctat2,dutrec
      equivalence (Angdum,tpvmsb),(Angdum(3),tpvmsc)
      equivalence (Angdum(9),ctat2),(Angdum(10),dutrec)
      include 'difnct.inc'
      include 'fcntrl.inc'
      include 'ltrapobs.inc'
      real*10 dfdly,ddr
      equivalence (dfdly,Deriv(2,1)),(ddr,Deriv(2,2))
      include 'mnsprt.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct
      real*4    acctim
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (acctim,Estf)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'radcrd.inc'
      include 'sbcom.inc'
      include 'scdta.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'
c
c quantities internal to this routine
      real*10 ctrecv,dfdlya,dum,dum1,fr,fract2,fruct,frv,tpvm,tpvm0
      integer   i,idum,j,jdpvm,jduck,loop1,loop2,lplctl,
     .          lsbctl,lscctl,mnspt1,n,n2,nn,numsit,nvlesn
      integer*2 npl, nspt, klb
      real*10 bdotx(2)
c external functions
      real*10 CTATF,DOT

c           xslcns(1-6,j)= position and velocity (not used) of center
c           of mass of the solar system relative to the sun
c              1= receive time for site 1
c              2= send time for site 1
c              3= receive time for site 2
c              4= send time for site 2(same as for site 1 if conven-
c                 tional vlbi observable, but may be different for
c                 n-count observable)
c           this convention makes the times analogous to those used
c           in the radar routines
c
c     receive time at site 1 and site coordinates determined
c         in ferctl
c
c
c------------------------initialization---------------------------------
      mnspt1 = 0
      lsbctl = -1
      lscctl = -1
      lplctl = -1
      loop1  = 0
      loop2  = 0
      numsit = 2
      if(nddiff.lt.0) numsit = 1
      call DGUESS(Tmdly1,1,Jd,Utrec)
 
c need velocites for ct-at correction even if delay only observable
      nvlesn = 1
      Nswcns = 0
      Nvlcns = 0
      if(Ict(27).gt.0) Nswcns = 1
      if(Ict(27).gt.1) Nvlcns = 1
c
c if two objects are observed, coordinates of the second
c are determined first and stored in /radcrd/ and /difcnt/
c
      npl   = Nplnt0
      nspt  = Nspot
      klb   = Klanb
      tpvm0 = Sbcom(1)
      if(kspt.ne.1) then
         npl   = Nplnt2
         nspt  = Nspot2
         klb   = Klans1
         tpvm0 = Sccom(1)
      endif
c
c
c
c----------iteration to determine send time from observed source--------
c----------given receive time at first site (n=1)-----------------------
c
      n  = 1
      n2 = 1
      ctrecv = Ctrec
 
  100 frv = ctrecv/Secday
c
c obtain receiving site wrt solar system barycenter
      call ETRP(1,Jd,frv,0,0,n,2)
      if(Jd.le.0) return
      if(nvlesn.gt.0) call ETRP(1,Jd,frv,1,0,n,2)
      if(Nswcns.gt.0) then
c
c apply solar system barycenter offset
         call SOTRP(Jd,frv,Xslcns(1,n2),0)
         if(Jd.le.0) return
         if(Nvlcns.gt.0) call SOTRP(Jd,frv,Xslcns(1,n2),1)
         do i = 1, 3
            Xemlsc(i,n) = Xemlsc(i,n) - Xslcns(i,n2)*Aultsc
            if(Nvlcns.gt.0) Xemlsc(i + 3,n) = Xemlsc(i + 3,n)
     .          - Xslcns(i + 3,n2)*Aultvl
         end do
      endif
      bdotx(n) = DOT(Xsite(1,n),Xemlsc(4,n))
      do i = 1, 3
         Xemlsc(i,n) = Xemlsc(i,n) + .5*bdotx(n)*Xemlsc(i + 3,n)
     .                  + Xsite(i,n)
      end do
c
c get transmit time
  200 call TIMINC(Jd,frv,Jdx,fr,-Tmdly1/Secday)
c determine if orbiting probe is observed body (either with
c sun or planet as central body)
      if(klb.le.0) goto 600
c
c obtain probe position at send time
c for pioneer-venus entry probes, call pvcrd if within atmosphere
      jdpvm = tpvm0
      if(jdpvm.eq.2443852) then
         tpvm = Jdx + fr
         if(tpvm.ge.tpvm0) then
            call PVCRD(Jdx,fr,n,kspt,nvel)
            if(kspt.eq.1) goto 300
            if(kspt.eq.2) goto 400
         endif
      endif
      if(kspt.eq.2) then
         call SCTRP(idum,Jdx,fr,0,lscctl,n)
         tpvmsc(n) = 0._10
         if(Jdx.gt.0) goto 400
c     since npath always = 1 for conventional vlbi observable,
c     store diff. delay rate in second position of difdly array
c
c
c
         Jd = 0
         return
      else
         tpvmsb(n) = 0._10
         call SBTRP(idum,Jdx,fr,0,lsbctl)
         if(Jdx.le.0) then
            Jd = 0
            return
         endif
      endif
  300 do j = 1, 3
         Xsbsun(j) = Xsb(j)*Aultsc
      end do
      goto 500
  400 do j = 1, 3
         Xsbsun(j) = Xsc(j,n)*Aultsc
      end do
 
c if pioneer-venus bus, correct antenna coordinates for roll
  500 if(jdpvm.eq.-2443852) call PVBUS(fr,Xsbsun,klb,4)
c
c test to see if planet is to be interpolated for a planet-
c centered orbiter or a planetary lander (spot on planet)
  600 if(Klan.ne.0) then
c
c obtain planet position at send time
         call PLTRP(1,Jdx,fr,0,lplctl)
         if(Jdx.le.0) then
            Jd = 0
            return
         else
            do j = 1, 3
               Xplsc(j) = Xp(j)*Aultsc
            end do
c
c increment by spot position
            if(nspt.gt.0) then
               call SPOTCD(Jdx,fr,mnspt1,nvel,Nplnt0,dum1,dum1,
     .                     kspt)
               do i = 1, 3
                  Xsbsun(i) = 0._10
                  Xplsc(i)  = Xplsc(i) + Xspcd(i,kspt)
               end do
            endif
c
c obtain probe wrt sun
            do j = 1, 3
               Xsbsun(j) = Xsbsun(j) + Xplsc(j)
            end do
c
c obtain probe wrt solar system barycenter
            if(Nswcns.gt.0) then
               call SOTRP(Jdx,fr,Xslcns(1,n2+1),0)
               do j = 1, 3
                  Xsbsun(j) = Xsbsun(j) - Xslcns(j,n2 + 1)*Aultsc
               end do
            endif
         endif
      endif
c
c------------------obtain receiving site wrt probe----------------------
      do j = 1, 3
         Xsitep(j,n) = Xemlsc(j,n) - Xsbsun(j)
      end do
c
c------------------calculate delay time---------------------------------
      Rsitp(n) = SQRT(Xsitep(1,n)**2 + Xsitep(2,n)**2 + Xsitep(3,n)**2)
      Tmdly2   = Rsitp(n)
      loop1    = loop1 + 1
      if(loop1.gt.15) call SUICID(
     .                   'MORE THAN 15 SEND ITERATIONS, STOP IN FERSB '
     .                   , 11)
      if(ABS(Tmdly2-Tmdly1).le.acctim) then
c
c------------------iteration converged -go to second site--------------
c
         Nit(20) = Nit(20) + loop1
         call TSAV(Tmdly2,1,Jd,Utrec)
         Tmdly1  = Tmdly2
         Dstf(4) = Tmdly1
         if(n.eq.1) Fract  = fr
         if(n.eq.2) fract2 = fr
c
c for differential n-count observable, determine send time
c given receive time at second site
         if(nintrf.lt.0) then
c
c-----------------------------------------------------------------------
c
c           for conventional interferometry observable, determine
c           receive time at second site given source send time
            dfdly = 0._10
            do while( .true. )
               dfdlya = dfdly
               Sidtm2 = Sidtm0 + Sidvel*(Utrec-dfdlya) + Dgst
               call SITCOR(Sidtm2,2,nvlesn,norm)
               call TIMINC(Jd,Ctrec/Secday,jduck,fruct,
     .                     -dfdlya/Secday)
c
c obtain earth w.r.t. ssbc
               call ETRP(1,jduck,fruct,0,0,2,2)
               if(jduck.le.0) then
                  Jd = 0
                  return
               else
                  if(nvlesn.gt.0)
     .                call ETRP(1,jduck,fruct,1,0,2,2)
                  if(Nswcns.gt.0) then
c
c apply solar system barycenter offset
                     call SOTRP(jduck,fruct,Xslcns(1,3),0)
                     if(Nvlcns.gt.0)
     .                   call SOTRP(jduck,fruct,Xslcns(1,3),1)
                     do i = 1, 3
                        Xemlsc(i,2) = Xemlsc(i,2) - Xslcns(i,3)
     .                                 *Aultsc
                        if(Nvlcns.gt.0) Xemlsc(i + 3,2)
     .                      = Xemlsc(i + 3,2) - Xslcns(i + 3,3)
     .                      *Aultvl
                     end do
                  endif
                  bdotx(2) = DOT(Xsite(1,2),Xemlsc(4,2))
                  do i = 1, 3
                     Xemlsc(i,2) = Xemlsc(i,2) + .5*bdotx(2)
     .                              *Xemlsc(i + 3,2) + Xsite(i,2)
                  end do
c
c obtain second site wrt probe
                  do j = 1, 3
                     Xsitep(j,2) = Xemlsc(j,2) - Xsbsun(j)
                  end do
c
c calculate differential delay
                  Rsitp(2) = SQRT(Xsitep(1,2)**2 + Xsitep(2,2)**2 +
     .                       Xsitep(3,2)**2)
                  dfdly    = Rsitp(1) - Rsitp(2)
                  loop2    = loop2 + 1
                  if(loop2.gt.15) call SUICID(
     .                'MORE THAN 15 RECV ITERATIONS, STOP IN FERSB ',
     .                11)
                  if(ABS(dfdly-dfdlya).le.acctim) then
                     Nit(19) = Nit(19) + loop2
                     goto 650
                  endif
               endif
            end do
         else if(n.ne.2 .and. numsit.ne.1) then
            n = 2
            ctrecv = Ctrec2
            loop1  = 0
            n2     = 3
            goto 100
         endif
c
c-----------calculations after iterations are complete------------------
c
c calculate unit vectors
  650    do n = 1, numsit
            call UVECTR(2,Xsitep(1,n),Rsitp(n),Xsitp0(1,n),dum)
         end do
c
c test for dummy observation below horizon
         if(Idumob.eq.1) then
            nn = numsit
            if(Klan.gt.0) nn = -numsit
            call HORIZN(Xsite,Xsitep,nn,Jd)
c
c test for planetary orbiter occulted by central body
c test for spot on planet turned away from earth
            if(Klan.gt.0) call CNTOCC(Xsb,Xsitep,-numsit,Jd)
            if(Jd.le.0) return
         endif
c
c determine geocentric velocitiy and partials for second site
c determined in ferctl for n-count observables
         if(nintrf.lt.0) call SITVEL(2,nvlesn,norm)
c
c determine spot velocity and partials
         if(nspt.gt.0) call SPOTCV(nvel,npl,kspt)
 
         do n = 1, numsit
            fr = Fract
            n2 = 2
            if(n.ne.1 .and. nintrf.ge.0) then
               fr = fract2
               n2 = 4
            endif
c
c determine pioneer-venus probe velocity and partials
            if(klb.le.0) goto 720
            if(jdpvm.eq.2443852) then
               if(tpvm.ge.tpvm0) then
                  call PVCRDV(n,nvel,kspt)
                  if(nvel.le.0) goto 750
c observed body partials (spot or pvm probe) assumed the same
c for both sites, so exit loop after n=1
                  if(kspt.eq.1) goto 660
                  if(kspt.eq.2) goto 680
               endif
            endif
c
c obtain probe velocity
            if(nvel.le.0) goto 750
            if(kspt.eq.2) then
               call SCTRP(idum,Jdx,fr,1,lscctl,n)
               goto 680
            else
               call SBTRP(idum,Jdx,fr,1,lsbctl)
            endif
  660       do i = 4, 6
               Xsbsun(i) = Xsb(i)*Aultvl
            end do
            goto 700
  680       do i = 4, 6
               Xsbsun(i) = Xsc(i,n)*Aultvl
            end do
  700       if(jdpvm.eq.-2443852) call PVBUSV(Xsbsun)
c
c obtain velocity of site
  720       if(nvel.le.0) goto 750
            do i = 4, 6
               Xemlsc(i,n) = Xemlsc(i,n) + Xsite(i,n)
            end do
c
c obtain velocity of planet
            if(Klan.ne.0) then
               call PLTRP(idum,Jdx,fr,1,lplctl)
               do i = 4, 6
                  Xplsc(i) = Xp(i)*Aultvl
c
c increment by spot velocity
                  if(nspt.gt.0) then
                     Xsbsun(i) = 0._10
                     Xplsc(i)  = Xplsc(i) + Xspcd(i,kspt)
                  endif
c
c obtain probe velocity wrt sun
                  Xsbsun(i) = Xsbsun(i) + Xplsc(i)
               end do
            endif
c
c obtain probe velocity wrt solar system barycenter
            if(Nvlcns.gt.0) then
               call SOTRP(Jdx,fr,Xslcns(1,n2),1)
               do j = 4, 6
                  Xsbsun(i) = Xsbsun(i) - Xslcns(i,n2)*Aultvl
               end do
            endif
c
c relative velocity
            do i = 4, 6
               Xsitep(i,n) = Xemlsc(i,n) - Xsbsun(i)
            end do
c
c calculate beta (=v/c)
            Beta(n,kspt) = DOT(Xsitep(4,n),Xsitp0(1,n))
 
         end do
         if(Nice.gt.0) then
c
c---------------compute differential delay rate-----------------------
c
            Difdly(2,kspt) = Beta(1,kspt) - Beta(2,kspt)
            return
         endif
c
c
c-----------corrections to differential delay---------------------------
c
c
c           convert atomic time delay to utc
  750    if(ntime.gt.0) then
            do i = 1, numsit
               Rsitp(i) = Rsitp(i)*fdev
            end do
         endif
c
c calculate differential delay in atomic time
         do i = 1, numsit
            Rsitp(i) = Rsitp(i) - bdotx(i)
         end do
         if(nintrf.lt.0) ctat2 = CTATF(jduck,fruct - Ctat/Secday,
     .                               3, 2)
         if(numsit.eq.2) Difdly(npath,kspt) = Rsitp(1) - Rsitp(2)
     .       - Ctat + ctat2
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
         if(nvel.le.0) then
            if(nintrf.lt.0) then
               Difdly(2,kspt) = Beta(1,kspt) - Beta(2,kspt)
            else if(npath.eq.2) then
               do i = 1, numsit
                  Beta(i,kspt) = (Rsave(2,i,kspt) - Rsave(1,i,kspt))
     .                            /Cnttim
               end do
            endif
         endif
      else
         Tmdly1 = Tmdly2
         goto 200
      endif
      return
      end
