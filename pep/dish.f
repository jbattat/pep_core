      subroutine DISH(eps,ad,icode,ndish)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      az, bmu, bnu, cosmga, deps, factor, fmu, fnu, omeg, 
     .          phi0, phiz, pi2, pimz, pvar, pvmu, pvnu, pvomeg, s2phi0,
     .          secmga
      real      sphi0, tanmga
      integer   nd
 
c*** end of declarations inserted by spag
 
 
      real*4    eps
      real*10 ad(3)
      integer*4 icode, ndish
c
c        r.b. goldstein dec. 1978
c        based on  1. code written by rdr in feb 76
c                  2. tm 391-240 by c. christensen
c                  3. notes prepared by david band  - uses shield factor
c
c     ndish - set up from kkp(78) by plnorb
c             0 no dish
c             1 mvm dish
c             2 pvo dish
c
c
c        eps  earth-probe-sun angle in radians
c        ad   effective area of dish in 3 directions
c        icode:
c        .le.0     error
c        .gt.0     number of coordinate system used
c        1         earth, s/c, sun
c        2         s/c (internal)
c        note - presently only icode=0 and 1 are used
c
      include 'dish1.inc'
 
      logical   old/.false./
      integer*2 rgn
      real*4    pi/3.14159265/
      real*4    mvmu, mvnu, mvar, mvomeg, eta1, eta2, alph, fzf, fzb, 
     .          fyf, fyb, aa(3)
      data mvmu/.075/, pvmu/.1875/, mvnu/.035/, pvnu/.1291/,
     .     mvar/1.17E4/, pvar/9.37E3/, mvomeg/36.87/, pvomeg/39.94/
 
      data factor/22160./
c
c
c setup on first call
c
      if(.not. (old)) then
         old = .true.
         nd  = ndish + 1
         if(nd.eq.2) then
 
c mvm
            fmu  = mvmu
            bmu  = mvmu
            fnu  = mvnu
            bnu  = mvnu
            az   = mvar
            omeg = mvomeg
         else if(nd.eq.3) then
 
c pvo
            fmu  = pvmu
            bmu  = pvmu
            fnu  = pvnu
            bnu  = pvnu
            az   = pvar
            omeg = pvomeg
         else
            icode = 0
            return
         endif
 
         tanmga = tan(pi*omeg/180.)
         phiz   = (90. - omeg)*pi/180.
         pimz   = pi - phiz
c
c
         secmga = tanmga*tanmga + 1
         secmga = sqrt(secmga)
         cosmga = 1.0/secmga
         aa(1)  = 2*cosmga/(1. + cosmga)
         aa(2)  = (1. - cosmga)*(2. + secmga)/(3.*(1+cosmga))
         aa(3)  = alog(cosmga)/(tanmga*tanmga)
         A(1)   = aa(1) + aa(2)
         A(2)   = -0.5 - 3.*aa(3)
         A(3)   = aa(1) - aa(2)
         A(4)   = -1.5 - 5.0*aa(3)
         A(5)   = 0.5 - aa(3)
         pi2    = pi/2.
      endif
c
c
c
      ad(1) = 0.
c
c        find the region (rgn)
c             1    0<eps<phiz
c             2    phiz< eps <pi/2
c             3    pi/2 < eps < pi-phiz
c             4    pi-phiz < eps < pi
      rgn = 0
      if((eps.ge.0) .and. (eps.le.phiz)) rgn   = 1
      if((eps.gt.phiz) .and. (eps.lt.pi2)) rgn = 2
      if((eps.ge.pi2) .and. (eps.lt.pimz)) rgn = 3
      if((eps.ge.pimz) .and. (eps.le.pi)) rgn  = 4
      if(rgn.eq.0) call SUICID('EPS OUT OF RANGE IN DISH', 6)
      fzf  = 0.
      fzb  = 0.
      fyf  = 0.
      fyb  = 0.
      eta1 = 0.
      eta2 = 0.
c
c branch to correct s/c
c
      if(ndish.ne.2) then
c
c mvm - use pvo for 'simple' regions
c use fits to numerical quadrature for other regions
c
         deps = abs(eps - pi2)
         if(rgn.ne.1) then
            if(rgn.eq.3) then
               ad(2) = -(0.14 + 0.55*deps)
            else if(rgn.eq.4) then
               goto 100
            else
               ad(2) = -0.14 + (1.53 - 0.86*deps)*deps
            endif
            ad(3) = 0.08 + 0.26*deps + 0.2*deps*deps
            ad(2) = ad(2)*factor/10.
            ad(3) = ad(3)*factor
            goto 300
         endif
c
c
c pvo
c
      else if(rgn.eq.2 .or. rgn.eq.3) then
         sphi0 = 0.
         if(abs(eps-pi2).ge.5.E-5) sphi0 = 1./(tan(eps)*tanmga)
         phi0   = asin(sphi0)
         s2phi0 = sin(2.*phi0)
         eta1   = (2.*phi0 + s2phi0)/pi
         eta2   = 0.5 - eta1/2.
         alph   = eps - pi
         call DISHFC(alph,bmu,bnu,fzb,fyb)
         if(rgn.eq.2) call DISHFC(eps,fmu,fnu,fzf,fyf)
         goto 200
      else if(rgn.eq.4) then
         goto 100
      endif
      eta1 = 1.
      eta2 = 0.
      call DISHFC(eps,fmu,fnu,fzf,fyf)
      goto 200
  100 alph = eps - pi
      eta1 = 0.
      eta2 = 1.
      call DISHFC(alph,bmu,bnu,fzb,fyb)
c
c
  200 ad(2) = az*(eta1*fyf + eta2*fyb)
      ad(3) = az*(eta1*fzf + eta2*fzb)
 
  300 icode = 1
      return
      end
