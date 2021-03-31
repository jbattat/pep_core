      subroutine HORIZN(xsite,xsitep,nn,jd)
 
      implicit none
c
c m.e.ash    oct 1970    subroutine horizn
c test for observed body below horizon in dummy observation mode
c
c arguments to HORIZN and SUNHOR
      real*10 xsite(9,2),xsitep(6,2)
      integer*4 nn,jd

c arguments to CNTOCC and SUNOCC (in addition to xsitep,nn,jd)
      real*10 xsb(6)

c  horizn called from rad...,fer...,ang...
c     nn         = negative, entry point cntocc will be called
c                  directly after this call to horizn
c     nn         = positive, no such call to cntocc
c     n=iabs(nn) = number of sites (1 or 2)
c     xsite(.,i) = position and velocity of site i relative to
c                  center of earth (i=1,2)
c     xsitep(.,i)= position and velocity of site i relative to
c                  observed body (i=1,2)
c     (note that only position coordinates 1,2,3 are used in this
c     subroutine. coordinates referred to mean equinox and equator of
c     ref. epoch and units in light seconds, but these do not matter)
c     if body below horizon at one of the n sites, set jd=-10 if
c     ict(37) so indicates
c ict(37)=-1 do not skip dummy observation which is below horizon
c            of observing site
c ict(37)= 0 skip such a dummy observation
c ict(37).gt.0 same as ict(37)=0 plus increment time of observation by
c            ict(37) minutes for first occurance
c     iswtch =0 do not increment time by ict(37) minutes
c               (set zero by entry point hrzswt called from subroutine
c               obsred after new dummy card read or set zero after time
c               incremented by ict(37) minutes)
c     iswtch =1 all right to increment time by ict(37) minutes
c               (set one if found to be above horizon
c               as signal to next point)

c array dimensions
      include 'globdefs.inc'

c common
      include 'fcntrl.inc'
      include 'param.inc'
      include 'redobs.inc'
      include 'statsrad.inc'

c external functions
      real*10 DOT

c local
      real*10 coszen,dstcnt,fdysc1,fdysc2,fdysc3,fdysc4,q,rsb2,rsitp
      integer*4 i,iswtch,jswtch,kswtch,n

c
c eps(13)= sine of elevation angle below which observed body is below
c          horizon at transmitting site or 2nd observing site, if any,
c          for dummy observations in subroutine horizn
c eps(17)= sine of elevation angle below which observed body is below
c          horizon at observing site for dummy observations in
c          subroutine horizn
c eps(18)= impact distance in kilometers below which observed body is
c          occulted by its central body for dummy observations in
c          subroutine horizn (entry point cntocc)
c
c     fdyscp = time in seconds of current dummy observation
c
c jct(35) = hour of start of dummy observations each day
c jct(36) = minute of start
c jct(37) = hour of  end  of dummy observations each day
c jct(38) = minute of end
c if jct(35 to 38) are zero, every time within day is allowed
c dummy observations done just between jct(35&36) and jct(37&38)
c if jct(37&38).lt.jct(35&36), then period of dummy obs include 0 hr utc
c these time limits could be the start and end of an optical observatory
c observing evening. this feature not needed for radar predictions
c if jct(35 to 38) nonzero, best to have ict(37) & ict(38) .le.0
c
c           test for correct time
      if(fdysc3.lt.0) then
         if(Fdyscp.gt.fdysc2 .or. Fdyscp.lt.fdysc1) then
            kswtch = 1
            goto 100
         endif
      else if(fdysc3.eq.0) then
         goto 100
      else if(Fdyscp.lt.fdysc2 .and. Fdyscp.gt.fdysc1) then
         kswtch = 1
         goto 100
      endif
      jd = -10
      if(kswtch.gt.0) then
         if(Intdyc.eq.0) Fdyscp = Fdyscp + fdysc4
c forward in time dummy observation series assumed
c ict(37) and ict(38) had better be non-positive
         kswtch = 0
      endif
      return
c
c should we test for below horizon
  100 if(Ict(37).ge.0) then
         n = iabs(nn)
c
c calculate zenith angle
         do i = 1,n
            coszen = -DOT(xsite(1,i),xsitep(1,i))
     .               /SQRT(DOT(xsite(1,i),xsite(1,i))
     .               *DOT(xsitep(1,i),xsitep(1,i)))
c     coszen = cosine of zenith angle at site i = sine of elevation
c     note that horizontal plane is taken as normal to vector
c     from center of earth to site and not normal to geoid
c     Eps(17) is limit for receive site, Eps(13) for send site
c
c           are we below horizon at ith site
            if(coszen.le.Eps(21-4*i)) then
c
c we are below horizon
               jd = -10
c
c should we increment time for first below horizon
               if(Ict(37).le.0) return
               if(iswtch.le.0) return
               if(nn.gt.0) goto 200
               if(Ict(38).le.0) goto 200
               if(jswtch.le.0) return
               goto 200
             endif
         end do
c
c we are above the horizon
         iswtch = 1
      endif
      return
c we cannot increment time if already incremented for occultation
c******** interaction of incrementing logic in two cases to be fixed
c
c increment time by ict(37) minutes
  200 Fdyscp = Fdyscp + 60*Ict(37)
      iswtch = 0
      return
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           entry point to set switch for starting new dummy obs.card
c           called from subroutine obsred after reading dummy card
      entry HRZSWT
      iswtch = 0
      jswtch = 0
      kswtch = 0
      fdysc1 = 60*(Jct(35)*60 + Jct(36))
      fdysc2 = 60*(Jct(37)*60 + Jct(38))
      fdysc3 = fdysc2 - fdysc1
      fdysc4 = 8.64E4_10 - fdysc3
      if(fdysc3.lt.0.0_10) fdysc4 = ABS(fdysc3)
      return
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c test for central body occultation in dummy observation mode
      entry CNTOCC(xsb,xsitep,nn,jd)
c  cntocc called from stdldp,stangl,stmdld,rad...,fer...,ang...
c     nn         = negative, entry point horizn was called
c                  just before this call to cntocc
c     nn         = positive, no such call to horizn
c     n=iabs(nn) = number of sites (1 or 2)
c     xsb(.)     = position and velocity of observed body relative
c                  to central body
c     xsitep(.,i)= position and velocity of site i relative to
c                  observed body (i=1,2)
c     (note that only position coordinates 1,2,3 are used in this
c     subroutine. coordinates referred to mean equinox and equator of
c     reference epoch and units in light seconds)
c     if observed body occulted by central body as seen at one of the
c     n sites, set jd=-10 if ict(38) so indicates
c ict(38)=-1 do not skip dummy observation which is is occulted by
c            central body of observed body
c ict(38)= 0 skip such a dummy observation
c ict(38).gt.0 same as ict(38)=0 plus increment time of observation by
c            ict(38) minutes for first occurance
c     jswtch =0 do not increment time by ict(38) minutes
c               (set zero by entry point hrzswt called from subroutine
c               obsred after new dummy card read or set zero after time
c               incremented by ict(38) minutes)
c     iswtch =1 all right to increment time by ict(38) minutes
c               (set one if found to be not occulted
c               as signal to next point)
c
c           should we test for occultation
      if(Ict(38).lt.0) return
      n = iabs(nn)
c
c calculate impact distances of vectors from observed body
c to observing sites relative to center of central body
c of observed body
      rsb2 = DOT(xsb,xsb)
      do i = 1,n
         rsitp = SQRT(DOT(xsitep(1,i),xsitep(1,i)))
         q     = DOT(xsb,xsitep(1,i))/rsitp
c
c does the central body occult the observed body
         if(q.lt.0._10) then
            dstcnt = SQRT(MAX(0._10,rsb2-q**2))*Ltvel
            if(dstcnt.le.Eps(18)) then
c
c we are occulted
               jd = -10
c
c should we increment time for first occulted
               if(Ict(38).le.0) return
               if(jswtch.le.0) return
               if(nn.gt.0) goto 300
               if(Ict(37).le.0) goto 300
               if(iswtch.le.0) return
               goto 300
            endif
         endif
      end do
c
c we are not occulted
      jswtch = 1
      return
c we cannot increment time if already incremented for below horizon
c******** interaction of incrementing logic in two cases to be fixed
c
c increment time by ict(38) minutes
  300 Fdyscp = Fdyscp + 60*Ict(38)
      jswtch = 0
 
      return
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c test for Sun occultation in dummy observation mode
      entry SUNOCC(xsb,xsitep,nn,jd)
c sunocc called from stdldp,stangl,stmdld,rad...,fer...,ang...
c   nn         = number of sites (1 or 2)
c   xsb(.)     = position and velocity of observed body relative to Sun
c   xsitep(.,i)= position and velocity of site i relative to
c                observed body (i=1,2)
c   (Note that only position coordinates 1,2,3 are used in this
c   subroutine. Coordinates referred to mean equinox and equator of
c   reference epoch and units in light seconds)
c   If observed body occulted by Sun or transiting Sun as seen at
c   one of the 'n' sites, set jd=-10
c   Note: size of Sun avoidance specified by Prmter(95) is also used for
c   calculating (real) transit observations
c
c is observation already skipped
      if(jd.le.0) return
c
c calculate impact distances of vectors from observed body
c to observing sites relative to center of Sun
      rsb2 = DOT(xsb,xsb)
      do i = 1,nn
         rsitp = SQRT(DOT(xsitep(1,i),xsitep(1,i)))
         q     = DOT(xsb,xsitep(1,i))/rsitp
c
c does the Sun occult the observed body or vice versa
         dstcnt = SQRT(MAX(0._10,rsb2-q**2))*Ltvel
         if(dstcnt.le.Prmter(95)) then
c
c we are occulted or transiting
            jd = -10
            return
         endif
      end do

      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c test for Sun above horizon in dummy observation mode
      entry SUNHOR(xsite,xsitep,nn,jd)
c sunhor called from rad...,fer...,ang...
c   nn         = number of sites (1 or 2)
c   xsite(.,i) = position and velocity of site i relative to
c                center of earth (i=1,2)
c   xsitep(.,i)= position and velocity of site i relative to
c                Sun (i=1,2)
c   (Note that only position coordinates 1,2,3 are used in this
c   subroutine. Coordinates referred to mean equinox and equator of
c   reference epoch and units in light seconds)
c   If Sun is above the horizon as seen at
c   one of the 'n' sites, set jd=-10
c
c is observation already skipped
      if(jd.le.0 .or. Jct(61).lt.0) return
      n = iabs(nn)
c
c calculate zenith angle
      do i = 1,n
         coszen = -DOT(xsite(1,i),xsitep(1,i))
     .    /SQRT(DOT(xsite(1,i),xsite(1,i))
     .    *DOT(xsitep(1,i),xsitep(1,i)))
c     coszen = cosine of zenith angle at site i = sine of elevation
c     note that horizontal plane is taken as normal to vector
c     from center of earth to site and not normal to geoid
c
c           are we below horizon at ith site
         if(coszen.ge.Eps(18+i)) then
c
c sun is above horizon
            jd = -10
            return
         endif
      end do
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c test for Sun above horizon at observed spot in dummy observation mode
      entry SUNHRSP(xsite,xsitep,jd)
c sunhrsp called from rad...,fer...,ang...
c   xsite = position and velocity of spot relative to
c           center of observed body
c   xsitep= position and velocity of spot relative to sun
c   (Note that only position coordinates 1,2,3 are used in this
c   subroutine. Coordinates referred to mean equinox and equator of
c   reference epoch and units in light seconds)
c   If Sun is above the horizon as seen at the observed spot,
c   set jd=-10
c
c is observation already skipped
      if(jd.le.0 .or. Jct(62).lt.0) return
c
c calculate zenith angle
      i=1
      coszen = -DOT(xsite(1,i),xsitep(1,i))
     . /SQRT(DOT(xsite(1,i),xsite(1,i))
     . *DOT(xsitep(1,i),xsitep(1,i)))
c     coszen = cosine of zenith angle at spot = sine of elevation
c     note that horizontal plane is taken as normal to vector
c     from center of body to spot and not normal to surface
c
c           are we below horizon
      if(coszen.ge.Eps(21)) then
c
c sun is above horizon
         jd = -10
         return
      endif
      return
      end
