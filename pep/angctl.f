      subroutine ANGCTL(kathy)
 
      implicit none
 
c subroutine angctl - j.f.chandler - 1979 jan
c set up for optical observations and call appropriate subroutine

c arguments
      integer*4 kathy
c          kathy=-1  azimuth,elevation observation (topocentric)
c          kathy= 0  meridian circle observation (geocentric referred
c                    to true equinox and equator of date)
c          kathy= 1  photographic observation (topocentric without
c                    elliptic aberration removed referred to the mean
c                    equinox and equator of ref. epoch) (astrometric)
c          kathy= 2  photographic observation (topocentric with
c                    all aberration removed referred to the mean
c                    equinox and equator of ref. epoch) (astrographic)

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdateq.inc'
      include 'coord.inc'
      include 'coordxeq.inc'
      include 'coordoeq.inc'
      include 'fcntrl.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*4    acctim,accdst,accprc
      equivalence (acctim,Estf),(accdst,Estf(2)),
     .            (accprc,Estf(3))
      real*10 dist
      equivalence (dist,Dstf(4))
      include 'kobequiv.inc'
      include 'jdfrequv.inc'
      include 'param.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'
 
c local
      real*10 vcentr(3)/3*0._10/
      integer*4 j,k01,lsw,norm

c external functions
      real*10 A1UT1,A1WWV,CTATF,CTUTF
 
      if(dist.le.0._10) dist = Tguess
      if(Dstf(5).le.0._10) Dstf(5) = dist
 
c initialize aberration velocity
      do j = 1,3
         Ves(j) = 0._10
      end do
c
c determine time quantities
      Utrec = 6E1_10*(60*Ihr+Imin)+Sec
      Fract = Utrec/Secday
      Ututs = 0.
 
      if(Jds.lt.2435490) then

c before atomic time - use coarse interpolation
         Ctat  = CTUTF(Jds,Fract)
         Atuts = 0._10

      else if(kathy.eq.0 .or. (Jds.lt.2441318 .and. kathy.gt.0)) then
 
c not azimuth-elevation
         Ctat = 32.15_10
         if(itime.le.1) Atuts = A1UT1(Jds,Fract)

c photographic obs after 1.0 jan 1972 assumed to have utc epoch
      else
         Atuts = A1WWV(Jds,Fract)
         Ututs = Atuts-A1UT1(Jds,Fract)
         Ctat  = CTATF(Jds,Fract+Atuts/Secday,1,1)
      endif
      Ctut  = (Ctat+Atuts)/Secday
c
c*  start=1000
c
c determine sidereal time quantities
      call SIDTIM(Jds,Ctat+Atuts-Ututs,Sidtm0,Sidvel,Dera)
      call PLATMO(Jds)
 
c read moon tape for precession-nutation data
      Jd = Jds
      call MNREED(Jd)
      if(Jd.gt.0) then
 
         lsw = -1
         call TIMINC(Jds,Fract,Jd,fr(1),Ctut)
         Jdx   = Jd
         fr(2) = fr(1)
 
c determine nutation-precession
         call PRCNUT(Jdx,fr(2))
         Sidtm = Sidtm0+Dgst-Longr(1)
         if(Kst1.ne.2) then
c
c get moon position and velocity at receive time
            call MNTRP(1,Jdx,fr(2),0,lsw,1)
            if(Jdx.le.0) then
 
c*  start=9000
               Jd = 0
               return
            else
               call MNTRP(1,Jdx,fr(2),-1,0,1)
c
c read earth-moon barycenter tape, perform interpolation
c to determine sun relative to earth (vel. vice-versa)
               call EMTRP(1,Jd,fr(1),0,lsw,1)
               if(Jd.le.0) return
               call EMTRP(1,Jd,fr(1),-1,0,1)
               do j = 1,3
                  x(j,2)   = Mnfct*Xm(j,1)-Xem(j,1)
                  x(j+3,2) = -Mnfct*Xm(j+3,1)+Xem(j+3,1)
 
c decide if annual or only elliptic aberration is needed
                  if(kathy.gt.0) then
                     if(Ncodf.ne.25) then
                        Ves(j) = Velipt(j)
                        goto 210
                     endif
                  endif
                  Ves(j) = (x(j+3,2)-vcentr(j))*Aultvl
  210          end do
            endif
         endif
c
c topocentric correction
         if(kathy.ne.0) then
            k01  = 1
            norm = 0
            if(Ncodf.gt.20 .or. kathy.eq.-1) norm = 1
            call SITCOR(Sidtm0+Dgst+Sidvel*(Utrec+Ututs),1,k01,norm)
 
c get site coordinates in au
            do j = 1,3
               Xsitau(j) = Xsite(j,1)/Aultsc
            end do
 
c include site velocity for az.-el. aberration
            if(kathy.lt.0) then
               do j = 1,3
                  Ves(j) = Ves(j)+Xsite(j+3,1)
               end do
            endif
         endif
c
c*  start=2000
         if(kathy.gt.2) then
            Jd = 0
         else
            if(Ncodf.gt.20) then
               call ANGDIF(kathy)
            else if(kathy.eq.-1) then
               call ANGAZ(kathy)
            else if(kathy.eq.0) then
               call ANGMER(kathy)
            else
               call ANGEQ(kathy)
            endif
 
c any final stuff common to all types goes here
            Tguess = dist
         endif
 
      endif
 
c*  start=9900
      return
      end
