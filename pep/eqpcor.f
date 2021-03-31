      subroutine EQPCOR(kathy)
 
      implicit none

c k.m.becker    october 1967    subroutine eqpcor
c moon phase correction inserted feb 1970
c new formulation added 1979 sep - j.f.chandler
c calculates equator-equinox and phase corrections to right
c ascension and declination

c arguments
      integer*4 kathy

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdat.inc'
      real*10 mnfcte,mnfct
      equivalence (Comcon(130),mnfcte),(Comcon(84),mnfct)
      include 'coord.inc'
      real*10 x(6,3)
      equivalence (Xm(1,1),x(1,1))
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'ltrapobs.inc'
      include 'number.inc'
      include 'obscrd.inc'
      include 'sitcrd.inc'

c external functions
      real*10 DOT

c local
      real*10 cc,cnsq,cof1,cof2,cog1,cog2,g,pn,rerp,
     .          resq,rn,rnsq,rpsq,s,ss,sthet,xdum
      integer   i
      real*10 xes(3),xsun(3),xps(3)
c xsun(1-3)=mean sun relative to embary in equatorial coordinates
c xes(1-3) =vector from sun to earth
c xps(1-3) =vector from sun to observed body
      real*10 u(3),umag
c
c           nice=-1  observation of right ascension only
c           nice= 0  observation of right ascension and declination
c           nice= 1  observation of declination only
c
c           equinox-equator corrections
      if(Neqnox.gt.0) then
         if(Sitf(1)(1:4).eq.'@REF') then

c pure rotation biases
            if(Nice.ne.0 .or. Jct(39).eq.-1 .or. Jct(39).eq.2)
     .       call SUICID(
     .       'MUST HAVE NORMAL RA+DEC FOR PURE ROTATION BIASES',12)
            u(1) = -Xsitp0(1,1) + Plat*Convds*Xsitp0(2,1)
     .             - Pequat*Convds*Xsitp0(3,1)
            u(2) = -Plat*Convds*Xsitp0(1,1) - Xsitp0(2,1)
     .             + Pnox*Convds*Xsitp0(3,1)
            u(3) = Pequat*Convds*Xsitp0(1,1) - Pnox*Convds*Xsitp0(2,1)
     .             - Xsitp0(3,1)
            do i=1,3
               Xsitp0(i,2) = -u(i)
            end do
            umag = SQRT(u(1)**2+u(2)**2)
            Deriv(2,1) = ATAN2(u(2),u(1))
            if(Deriv(2,1).lt.0._10) Deriv(2,1) = Deriv(2,1)+Twopi
            Deriv(2,1) = Deriv(2,1)/Convhs
            Deriv(2,2) = ATAN2(u(3),umag)/Convds
         else

c equinox/equator/latitude biases
            if(Nice.le.0) Deriv(2,1) = Deriv(2,1)
     .       + Pnox - Pequat*Calph*Tdelt/15._10
            if(Nice.ge.0) Deriv(2,2) = Deriv(2,2) + Pequat*Salph + Plat
         endif
      endif
 
      entry DPHCOR(kathy)
      if(Nphase.gt.0 .and. Ncph.gt.0) then
c planetary phase corrections
c observed body is moon
         if(Nplnt0.eq.10) then
            call MENSUN(Jd,Dstf(1),xsun,0)
            if(kathy.gt.0) then
               do i = 1,3
                  xps(i) = -xsun(i) + Xm(i,1)*mnfcte
                  xes(i) = -xsun(i) - Xm(i,1)*mnfct
               end do
            else
               call CORCHN(xps,Xm)
               call CORCHN(xes,xsun)
               do i = 1,3
                  xdum   = -xes(i) + xps(i)*mnfcte
                  xes(i) = -xes(i) - xps(i)*mnfct
                  xps(i) = xdum
               end do
            endif
 
c observed body is planet
         else if(kathy.gt.0) then
            do i = 1,3
               xes(i) = -x(i,2)
               xps(i) = Xp(i)
            end do
         else
            call CORCHN(xps,Xp)
            call CORCHN(xes,x(1,2))
            do i = 1,3
               xes(i) = -xes(i)
            end do
         endif
 
c calculate quantities for phase correction
         resq = DOT(xes,xes)
         rpsq = DOT(xps,xps)
         rerp = DOT(xes,xps)
         cof1 = rerp - rpsq
         cof2 = rerp - resq
         rnsq = -(cof1 + cof2)
         rn   = SQRT(rnsq)
         s    = SQRT(resq*rpsq - rerp**2)
         g    = rn*s
         cog1 = cof1/g
         cog2 = cof2/g
         do i = 1,3
            Xppxe(i) = xps(i) - xes(i)
            Xrho(i)  = cog1*xes(i) + cog2*xps(i)
         end do
         pn = SQRT(rnsq*rpsq)
         cc = cof1/pn
         ss = s/pn
         Cthet(1) = D
         Drho     = D*Aphs(1)
         if(Ncph.gt.1) then
 
c compute cosine terms
            sthet = 0._10
            do i = 2,Ncph
               Cthet(i) = Cthet(i - 1)*cc - sthet*ss
               sthet    = Cthet(i - 1)*ss + sthet*cc
               Drho     = Drho + Aphs(i)*Cthet(i)
            end do
         endif
 
c compute phase correction to right ascension
         cnsq = rnsq - Xppxe(3)**2
         if(Nice.le.0) then
            Phal = (Xrho(2)*Xppxe(1) - Xrho(1)*Xppxe(2))/cnsq/Convhs
            if(Ncodf.gt.20) Phal = Phal*15._10
            Deriv(2,1) = Deriv(2,1) + Drho*Phal
         endif
 
c correction to declination
         if(Nice.ge.0) then
            Phde = Xrho(3)/SQRT(cnsq)/Convds
            Deriv(2,2) = Deriv(2,2) + Drho*Phde
         endif
      endif
      return
      end
