      real*10 function BODFN(k,j,s)
 
      implicit none

c
c subroutine ''bodfn''     w.b.smith   10-24-68
c evaluate right sides of n-body equations of motion
c m.e.ash   may 1974   asteroid perturbations added
c
c arguments
      integer*4 k,j
      real*10 s
c k = equation number
c j = iteration number
c s = time (julian date + 0.5)

c array dimensions
      include 'globdefs.inc'

c common 
      include 'aprtbf.inc'
      include 'bdctrl.inc'
      include 'bddtaint.inc'
      include 'bodstf.inc'
      include 'ellips.inc'
      include 'funcon.inc'
      include 'intstf.inc'
      include 'orblun.inc'
      include 'output.inc'
      include 'param.inc'
      include 'prtcod.inc'

c local variables
      integer*4 jx,jy,jz,jzt,kkk,km,kx,ky,kz,l,ndex,npl,nppl
      real*10 rjp(15,15), r(15), r3(15), xsp3(3,15)
      real*10 xjp(3,15,15), xes(3), xms(3),
     .       rjp3(15,15),x(6)
      real*10 va(3,15),vhp(3,15),sumr(3),sumrm(3)
      real*10 plj2(3,15),sunj2(3)
      real*10 alf,astsm,g,harter,omega,prefxs,q,r2
      real*10 rjp2,sum,suma,sumb,savemnfn

c external functions
      real*10 DOT,MORFN

      kkk = (k-1)/6
      km  = k - 6*kkk
      if(km.lt.4) then
         BODFN = Y(k+3, j)
         return
      endif
 
c km = equation number for a given body
c kz = index of the current body
c npl= planet number of the current body
      l    = km - 3
      kz   = kkk + 1
      npl  = Nplbdy(kz)
      Mp   = Mass(npl)
      if(k.eq.4) then
         if(j.ne.2) then
 
c calculations done once per step
            Jd    = s
            Fract = s - Jd

            if(Orbint) then
               if(Kbdy(19).lt.0) then
               else if(Kbdy(19).eq.0) then
c determine elliptic orbit coordinates for encke's method
                  Tlpt = s - Tlpt0
                  call ELIPT(0,Tlpt)
               else
c determine brown mean lunar orbit coordinates
                  call LUNORB(Jd,Fract,0)
               endif
            endif

c setup once per step for asteroids
            if(Astflg) call BDAST0(s)
 
c read moon or nutation from perturbing planet tape or moon tape
            if(Lpert.gt.0) then
               if(Ipert2.eq.0) then
                  call PRTCRD(Jd,Fract)
               else
                  call MNPCRD(Jd,Fract)
               endif
            endif
            if(Kbdy(22).ge.0) then
               Tvary = s-Tvary0
               Gamat = Gama*(1._10 + Gmvary*Tvary)
            endif
         endif
 
c calculations done once per iteration
 
c get heliocentric positions and velocities of planets
c if moon is being integrated simultaneously, the integrated quantity
c is the geocentric moon position, and we still have the heliocentric
c position of the embary as the "planet" position being integrated
         do kx=1,3
            ky=kx
            do jy=1,Nbodyp
               nppl=Nplbdy(jy)
               if(nppl.le.9) then
                  Xpert(kx,nppl)=Y(ky,j)
                  Xpert(kx+3,nppl)=Y(ky+3,j)
               endif
               xjp(kx,jy,Nsun)=Y(ky,j)
               xjp(kx,Nsun,jy)=-Y(ky,j)
               vhp(kx,jy)=Y(ky+3,j)
               ky=ky+6
            end do
c if moon is an integrated body, copy that as well
            do jy=Nbodyp+1,Nbody
               if(jy.eq.Jmoon) then
                  if(Kbdy(19).lt.0) then
                     Xpert(kx,10)=Y(ky,j)
                     Xpert(kx+3,10)=Y(ky+3,j)
                  else if(Kbdy(19).eq.0) then
                     Xpert(kx,10)=Y(ky,j) + Ylpt(kx)
                     Xpert(kx+3,10)=Y(ky+3,j) + Ylpt(kx+3)
                  else
                     Xpert(kx,10)=Y(ky,j) + Ylun(kx)
                     Xpert(kx+3,10)=Y(ky+3,j) + Ylun(kx+3)
                  endif
               else
c don't need lunar libration in the planetary part, just orbit
               endif
               ky=ky+6
            end do
         end do
c
         if(Ksepem.gt.0) then
c this block of code used when earth and moon
c effects on other bodies are considered separately
c here, we set the moon position and velocity vectors and replace the
c embary position and velocity with the corresponding earth quantities
            do kx=1,3
               xes(kx)=xjp(kx,Jembry,Nsun)-Mass(10)*Xpert(kx,10)
               xms(kx)=xjp(kx,Jembry,Nsun)+Masse*Xpert(kx,10)
               xjp(kx,Jembry,Nsun)=xes(kx)
               xjp(kx,Nsun,Jembry)=-xes(kx)
               xjp(kx,Jmoonx,Nsun)=xms(kx)
               xjp(kx,Nsun,Jmoonx)=-xms(kx)
               vhp(kx,Jmoonx)=vhp(kx,Jembry)+Masse*Xpert(kx+3,10)
               vhp(kx,Jembry)=vhp(kx,Jembry)-Mass(10)*Xpert(kx+3,10)
            end do
         endif

         if(Kbdy(36).ge.0) then
c set up velocity arrays for new g.r. formalism
c
c get sun velocity from reaction
            do kx=1,3
               sum=0._10
               do jy=1,Nbodyx
                  sum=sum-vhp(kx,jy)*Mass1(jy)
               end do
               va(kx,Nsun)=sum/Mascnt
            end do
c get barycentric velocities of planets
            do kx=1,3
               do jy=1,Nbodyx
                  va(kx,jy)=vhp(kx,jy)+va(kx,Nsun)
               end do
            end do
         endif
 
c find distance and distance**3 and (x,y,z)/r**3 for planet w.r.t. sun
         do jy = 1, Nbodyx
 
c wrt sun
            r2     = DOT(xjp(1,jy,Nsun),xjp(1,jy,Nsun))
            r(jy)  = SQRT(r2)
            r3(jy) = r(jy)*r2
            rjp(jy,Nsun)=r(jy)
            rjp(Nsun,jy)=r(jy)
            do kx=1,3
               xsp3(kx,jy) = Epfacs*xjp(kx,jy,Nsun)/r3(jy)
            end do
            jzt = jy + 1
 
c between bodies
            do jz = jzt, Nbodyx
               do kx = 1, 3
                  xjp(kx,jz,jy) = xjp(kx,jz,Nsun)-xjp(kx,jy,Nsun)
                  xjp(kx,jy,jz) = -xjp(kx,jz,jy)
               end do
               rjp2 = DOT(xjp(1,jz,jy),xjp(1,jz,jy))
               rjp(jz,jy)  = SQRT(rjp2)
               rjp(jy,jz)  = rjp(jz,jy)
               rjp3(jz,jy) = rjp2*rjp(jz,jy)
               rjp3(jy,jz) = rjp3(jz,jy)
            end do
         end do
c
         if(Kbdy(36).ge.0) call BODGR(xjp,rjp,va,Nsun,sumr)
c
c set up for solar j2 effects
         if(Kbdy(23).ge.0) then
            do kx=1,3
               sunj2(kx)=0._10
            end do
            do jy = 1, Nbodyx            
               g = DOT(xjp(1,jy,Nsun),C3)/r(jy)
               prefxs = Shar2*Gamat/r(jy)/r3(jy)
               do kx=1,3
                  harter = prefxs*(xjp(kx,jy,Nsun)*(7.5_10*g**2-1.5_10)
     .             /r(jy) - 3._10*C3(kx)*g)
                  plj2(kx,jy) = harter*Epfacp(jy)
                  sunj2(kx)   = sunj2(kx) + harter*Epfacs*Mass1(jy)
               end do
            end do
         endif
c
c moon-related quantities
c set acceleration of sun due to non-spherical gravity fields of earth,moon
         do kx=1,3
            Asunmor(kx)=0._10
            Aembmor(kx)=0._10
         end do
         if(Jmoon.gt.0) then
            savemnfn = MORFN(Jmoon*6-2,j,s)
         endif
      endif

      if(kz.gt.Nbodyp) then
         if(k.eq.Jmoon*6-2) then
            BODFN=savemnfn
         else
            BODFN=MORFN(k,j,s)
         endif
         return
      endif

c calculations done once per iteration for each body (km = 4)
      if(l.eq.1) then
 
c set up position-velocity vector  x(6)
         ndex=6*(kz-1)
         do kx = 1, 6
            x(kx)= Y(ndex+kx,j)
         end do
      endif
 
c point-mass effects of sun and planets
      omega = 0._10
 
      if(npl.eq.3 .and. Ksepem.gt.0) then
c
c earth-moon
         q = -Gamat*(Mface*xes(l)/r3(Jembry)+Mfacm*xms(l)/r3(Jmoonx))
         do jx = 1, Nbodyp
            if(jx.ne.Jembry) omega = omega + Mass1(jx)*
     .       ((Epfacp(Jembry)*Masse*xjp(l,jx,Jembry)/rjp3(jx,Jembry) +
     .       Epfacp(Jmoonx)*Mass(10)*xjp(l,jx,Jmoonx)/rjp3(jx,Jmoonx))
     .                  - xsp3(l,jx))
         end do
      else
c non earth moon case or embary as a point mass
c acceleration due to sun
c@@      q = -Gamat*mfac*xsp3(l,kz)
         q = -Gamat*Mfac(kz)*x(l)/r3(kz)
 
c acceleration wrt sun due to planets (and moon if separate)
         do jx = 1, Nbodyx
            if(jx.ne.kz) omega = omega + Mass1(jx)
     .       *(Epfacp(kz)*xjp(l,jx,kz)/rjp3(jx,kz) - xsp3(l,jx))
         end do
      endif
      BODFN = q + omega*Gamat
 
c extra forces
c           general relativity
      if(Kbdy(21).ge.0) then
         if(Kbdy(36).ge.0) then
            if(l.eq.1) then
               call BODGR(xjp,rjp,va,kz,sumr)
               if(npl.eq.3 .and. Ksepem.gt.0) then
                  call BODGR(xjp,rjp,va,Jmoonx,sumrm)
                  do kx=1,3
                     sumr(kx)=Masse*sumr(kx)+Mass(10)*sumrm(kx)
                  end do
               endif
            endif
         else
            if(l.eq.1) then
               alf = Gamat/Cvel2
               suma = DOT(x(4),x(4))
               sumb = DOT(x(1),x(4))
            endif
            sumr(l)= ( x(l)*(4._10*alf/r(kz)-suma/Cvel2)
     .                 +4._10*x(l+3)*sumb/Cvel2 ) / r3(kz)
         endif
         BODFN=BODFN+Relftb(kz)*Gamat*sumr(l)
      endif
c           sun j2
c treating embary as a point mass
      if(Kbdy(23).ge.0) then
         harter = plj2(l,kz)+sunj2(l)
         BODFN  = BODFN + harter
      endif
 
c add accelerations due to asteroids (individual and ring)
      if(Astflg) then
         call BDASTF(l, x, r(kz), r3(kz), astsm)
         BODFN = BODFN + astsm*Gamat
      endif
c sun acceleration due to harmonic expansion of earth,moon gravity fields
c as calculated in the moon integration code (otherwise 0)
      BODFN=BODFN-Asunmor(l)
      if(npl.eq.3) BODFN=BODFN+Aembmor(l)
      return
      end
