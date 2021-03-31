      subroutine MORHAR(ncall,j)
 
      implicit none

c
c     r.king and r.cappallo   july 1977   subroutine morhar
c     evaluation of the effect of the gravity fields of the earth and
c     moon in the right side of the equations for the lunar orbit and
c     rotation.  called by morfn, this routine incorporates (at present)
c     m.slade's coding from monfn for the orbital effect and king and
c     cappallo's coding from mrtfn for the rotation effect. the two sets
c     of code should be standardized soon and the effects of 4th degree
c     forces and torques added.
c
c
c
      integer ncall,j
c         ncall= -2 setup once per iteration of a given step for partial
c              = -1 setup once per iteration of a given step for motion
c              =  0 setup once per step
c              =  k evaluate acceleration for equation k.gt.0
c         j    = pointer into integration output array
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'harmor.inc'
      include 'intstf.inc'
      include 'metuna.inc'
      include 'mnrtlb.inc'
      include 'morcrd.inc'
c
c     ecor  = earth relative to sun (re is distance)
c     mcor  = moon relative to sun (rm is distance)
c     mecor = moon relative to earth (rem is distance)
c     mcor etc. are equatorial coordinates
c     smcor etc. are selenodetic coordinates
c
c
      include 'morstf.inc'
c use summh for both He and Hs, sump for both Dhedx and Dhsdx
      real*10 summh(3,2),sump(3,3,2)
      equivalence (He,summh),(Dhedx,sump)
      include 'nutces.inc'
      include 'output.inc'
      include 'param.inc'

c
c quantities internal to this routine
      real*10 slat,clat,slat1(3),clatr,rh1,rh(19),slng(19,3),
     .          clng(19,3),lng1(3),leg(19),leg1(19),leg2(19),
     .          gleg(54),gleg1(54),gleg2(54),zone(3,19,0:2),
     .          tess1(3,54,0:2),tess2(3,54,0:2),snh1(54),snh2(54),di2,
     .          djh,sumt,clap2(3,3),slap2(3,3),
     .          sumpt,dcnm
      real*10 ezone(3,19),szone(3,19),eslng(19),eclng(19),
     .          sslng(19),sclng(19),tesse1(3,54),tesse2(3,54),
     .          tesss1(3,54),tesss2(3,54)
      equivalence (zone(1,1,1),ezone),(zone(1,1,2),szone),
     . (slng(1,1),eslng),(slng(1,2),sslng),
     . (clng(1,1),eclng),(clng(1,2),sclng),
     . (tess1(1,1,1),tesse1),(tess1(1,1,2),tesss1),
     . (tess2(1,1,1),tesse2),(tess2(1,1,2),tesss2)
      real*10 cor(3),r,r2,ddw1dq,ddw2dq,ddw3dq
      real*10 a(3),dadx(3,3),x(3),asum(3),dsumh,asun(3),
     .          dadxsm(3,3),dum
      real*10 drrdpsi(3),drrdth(3),drrdphi(3),drr(3,3),drrx(3,3),
     . drrf(3,3),dady(3,3)
      equivalence (drr,drrdpsi),(drr(1,2),drrdth),(drr(1,3),drrdphi)
      integer i,i1,icmkkk,ih,ii,inh,ip,is,jh,jj,l
c
c kronecker delta function
      real*10 delta(3,3)/1._10,3*0._10,1._10,3*0._10,1._10/
 
      if(ncall.gt.0) then
 
         if(.not.Orbint .or. ncall.gt.Korb) then
c
c-----------------------------------------------------------------------
c
c           effect of harmonics on rotation partials
c
            icmkkk = Icrtrl(Kkk)
c
c check if parameter is zonal harmonic coefficient
            if(icmkkk.ne.Jzone) then
c
c check if parameter is tesseral cosine harmonic coefficient
               if(icmkkk.eq.Jcos) then
                  ii = 0
                  dcnm = 1._10
                  do i = 1, Ntesp
                     i1 = i + 1
                     do jh = 1, i1
                        ii = ii + 1
                        if(Icos(ii).eq.Kkk) go to 100
                     end do
                  end do
c
c check if parameter is tesseral sine harmonic coefficient
               else if(icmkkk.eq.Jsin) then
                  ii = 0
                  do i = 1, Ntesp
                     i1 = i + 1
                     do jh = 1, i1
                        ii = ii + 1
                        if(Isin(ii).eq.Kkk) then
                           djh    = jh
                           ddw1dp = -Gamat3*Mma*(Smecor(2)*(eslng(jh)*
     .                              tesse1(3,ii)+eclng(jh)*tesse2(3,ii)
     .                              *djh) - Smecor(3)
     .                              *(eslng(jh)*tesse1(2,ii)+eclng(jh)
     .                              *tesse2(2,ii)*djh))
                           ddw2dp = -Gamat3*Mmb*(Smecor(3)*(eslng(jh)*
     .                              tesse1(1,ii)+eclng(jh)*tesse2(1,ii)
     .                              *djh) - Smecor(1)
     .                              *(eslng(jh)*tesse1(3,ii)+eclng(jh)
     .                              *tesse2(3,ii)*djh))
                           ddw3dp = -Gamat3*Mmc*(Smecor(1)*(eslng(jh)*
     .                              tesse1(2,ii)+eclng(jh)*tesse2(2,ii)
     .                              *djh) - Smecor(2)
     .                              *(eslng(jh)*tesse1(1,ii)+eclng(jh)
     .                              *tesse2(1,ii)*djh))
                           if(Kmr(81).ge.1) then
                              ddw1dp = ddw1dp -
     .                                 Gamat*Mma*(Smcor(2)*(sslng(jh)
     .                                 *tesss1(3,ii)+sclng(jh)
     .                                 *tesss2(3,ii)*djh) - Smcor(3)
     .                                 *(sslng(jh)*tesss1(2,ii)
     .                                 +sclng(jh)*tesss2(2,ii)*djh))
                              ddw2dp = ddw2dp -
     .                                 Gamat*Mmb*(Smcor(3)*(sslng(jh)
     .                                 *tesss1(1,ii)+sclng(jh)
     .                                 *tesss2(1,ii)*djh) - Smcor(1)
     .                                 *(sslng(jh)*tesss1(3,ii)
     .                                 +sclng(jh)*tesss2(3,ii)*djh))
                              ddw3dp = ddw3dp -
     .                                 Gamat*Mmc*(Smcor(1)*(sslng(jh)
     .                                 *tesss1(2,ii)+sclng(jh)
     .                                 *tesss2(2,ii)*djh) - Smcor(2)
     .                                 *(sslng(jh)*tesss1(1,ii)
     .                                 +sclng(jh)*tesss2(1,ii)*djh))
                           endif
                           return
                        endif
                     end do
                  end do
               endif
            else
               do i = 1, Nzonp
                  if(Izone(i).eq.Kkk) then
                     ddw1dp = Gamat3*Mma*(Smecor(2)*ezone(3,i)
     .                        - Smecor(3)*ezone(2,i))
                     ddw2dp = Gamat3*Mmb*(Smecor(3)*ezone(1,i)
     .                        - Smecor(1)*ezone(3,i))
                     ddw3dp = Gamat3*Mmc*(Smecor(1)*ezone(2,i)
     .                        - Smecor(2)*ezone(1,i))
                     if(Kmr(81).ge.1) then
                        ddw1dp = ddw1dp +
     .                           Gamat*Mma*(Smcor(2)*szone(3,i) -
     .                           Smcor(3)*szone(2,i))
                        ddw2dp = ddw2dp +
     .                           Gamat*Mmb*(Smcor(3)*szone(1,i) -
     .                           Smcor(1)*szone(3,i))
                        ddw3dp = ddw3dp +
     .                           Gamat*Mmc*(Smcor(1)*szone(2,i) -
     .                           Smcor(2)*szone(1,i))
                     endif
                     return
                  endif
               end do
            endif
 
c rhs for the kth orbital equation
         else if(Kkk.gt.0) then
 
c variational equations rhs
            icmkkk = Icmtrl(Kkk)
            dsumh  = 0._10
            do jj = 1, 3
               dsumh = dsumh + dadxsm(Kk, jj)*Dmcor(jj, Kkk)
            end do
            if(icmref(kkk).gt.0) then
               do jj=1,3
                  dsumh=dsumh+dady(Kk,jj)*Y(Korb+6*Icmref(Kkk)+jj,j)
               end do
            endif
c also pick up contribution due to g*m changes
            if(icmkkk.eq.3) then
               dsumh  = dsumh + Gamat/Gamtem*asum(Kk)
            else if(icmkkk.eq.32) then
               dsumh = dsumh + Gamem/Gamtem*Tvary*asum(Kk)
     .          + Gama/Gamat*Tvary*asun(Kk)
 
c check if parameter is harmonic coefficient
            else if(icmkkk.eq.Jzone) then
               do i = 1, Nzonp
                  if(Imzone(i).eq.Kkk) dsumh=dsumh-zone(Kk,i,0)
               end do
            else if(icmkkk.eq.Jcos) then
               ii = 0
               do i = 1, Ntesp
                  i1 = i + 1
                  do jh = 1, i1
                     ii = ii + 1
                     if(Imcos(ii).eq.Kkk) dsumh=dsumh-tess1(kk,ii,0)
                  end do
               end do
            else if(icmkkk.eq.Jsin) then
               ii = 0
               do i = 1, Ntesp
                  i1 = i + 1
                  do jh = 1, i1
                     ii = ii + 1
                     if(Imsin(ii).eq.Kkk) dsumh=dsumh-tess2(kk,ii,0)
                  end do
               end do
            endif

c include indirect contribution of J2, beta, and gamma to C22
            if(Icmref(Kkk).gt.0) then
               dsumh = dsumh - Dc22(Icmref(Kkk))*tess1(Kk,2,0)
            endif
            Fn(1) = Fn(1) + dsumh
         else
 
c moon acc.rel. to earth same (after scaling mass) as sun's
            Fn(1) = Fn(1) + asum(Kk) + asun(Kk)
         endif
         return
      endif

c-----------------------------------------------------------
c setup operations

      do i=1,3
         Asunmor(i)=0._10
         Aembmor(i)=0._10
      end do

c harmonic forces calculated once/iteration and saved
      if(Orbint) then
         ip = 0
         if(Iparm.gt.1) ip = 1000
         do i = 1, 3
            asum(i) = 0._10
            asun(i) = 0._10
            do jj = 1, 3
               dadxsm(i, jj) = 0._10
            end do
         end do
         if(Keh.gt.0) then
c
c earth zonal field interacting with the moon
c
            call GFIELD(Mecor, Nutprc, Erzhar, dum, dum, Gamtem,
     .       Erad, Keh, 0, ip+10, asum, dadxsm, dum,dum,dum)
            if(Km(33).ge.0) then
 
c now get solar-earth j2 interaction
               do i = 1, 3
                  x(i) = -Ecor(i)
               end do
               call GFIELD(x, Nutprc, Erzhar, dum, dum, Gamat, Erad,
     .          2, 0, 10, asun, dum,dum,dum,dum)
               do i=1,3
                  Asunmor(i)=Asunmor(i)+asun(i)*Merth
               end do
            endif
         endif
c ignoring effect of sol-ej2 in varl. eqns.
c
c lunar zonal and tesseral fields acting on earth
c zone,tess1,tess2 as returned by GFIELD include the GM factor
c
         if(Kmh.gt.0) then
            do i = 1, 3
               x(i) = -Mecor(i)
            end do
            call GFIELD(x, Mrotlb, Zhar, Char, Shar, Gamtem, Mrad,
     .       Kmh, Kmh, ip+110, a, dadx, zone, tess1, tess2)
            if(Kkm(60).gt.0)
     .       write(6,99111) (zone(i,1,0),i=1,3),(tess1(i,2,0),i=1,3),x
99111       format('zone j2',1p3d20.12/'tess1 c22 ',1p3d20.12/
     .       'x',3d20.12)
            do i = 1, 3
               asum(i) = asum(i) - a(i)
               if(Iparm.gt.1) then
                  do jj = 1, 3
                     dadxsm(i, jj) = dadxsm(i, jj) + dadx(i, jj)
                  end do
               endif
            end do
            if(Iparmr.gt.1 .and. Iparm.gt.1) then
               drrdpsi(1)=0._10
               drrdpsi(2)=0._10
               drrdpsi(3)=1._10
               drrdth(1)=Cpsi
               drrdth(2)=Spsi
               drrdth(3)=0._10
               drrdphi(1)=Stheta*Spsi
               drrdphi(2)=-Stheta*Cpsi
               drrdphi(3)=Ctheta
               do jj=1,3
                  call CROSS(drr(1,jj),x,drrx(1,jj))
                  call CROSS(a,drr(1,jj),drrf(1,jj))
               end do
               call PRODCT(dadx,drrx,dady,3,3,3)
               do i=1,3
                  do jj=1,3
                     dady(i,jj)=dady(i,jj)+drrf(i,jj)
                  end do
               end do
            endif
c 
c lunar zonal and tesseral fields acting on sun
            do i = 1, 3
               x(i) = -Mcor(i)
            end do
            call GFIELD(x, Mrotlb, Zhar, Char, Shar, Gamat, Mrad,
     .       Kmh, Kmh, ip+10, a, dadx, dum, dum, dum)
            do i=1,3
               asum(i)=asum(i)-a(i)
               Asunmor(i)=Asunmor(i)+a(i)*Mmoon
               if(Iparm.gt.1) then
                  do jj=1,3
                     dadxsm(i,jj)=dadxsm(i,jj)+dadx(i,jj)*Masse
                  end do
               endif
            end do
         endif
         if(Keh.gt.0 .or. Kmh.gt.0) then
            do i=1,3
               Aembmor(i)=-Asunmor(i)/Mass(3)
            end do
         endif
      endif
c
c           determine quantities for harmonic terms in rotation
c
c           is=1     forcing terms due to earth
c           is=2     forcing terms due to sun
c           loop begins at statement 113, ends at 180
c
      if(.not.Rotint) return
      is = 1
      do i = 1, 3
         cor(i) = Smecor(i)
      end do
      r  = Srem
      r2 = Srem2

      do while(.true.)
 
c latitude terms
         slat = cor(3)/r
         clat = SQRT(1._10 - slat**2)
         do i = 1, 3
            slat1(i) = -cor(i)*slat/r
         end do
         slat1(3) = 1._10 + slat1(3)
         clatr    = clat*r
         rh1   = r/Mrad
         rh(1) = rh1**2
         if(Ntop.ge.2) then
            do i = 2, Ntop
               rh(i) = rh(i - 1)*rh1
            end do
         endif
 
c longitude terms
         if(Ntess.gt.1) then
            slng(1,is) = cor(2)/clatr
            clng(1,is) = cor(1)/clatr
            do i = 2, Ntess
               slng(i,is) = slng(i-1,is)*clng(1,is)
     .          + clng(i-1,is)*slng(1,is)
               clng(i,is) = clng(i-1,is)*clng(1,is)
     .          - slng(i-1,is)*slng(1,is)
            end do
            lng1(1) = -slng(1,is)/clat
            lng1(2) = clng(1,is)/clat
            lng1(3) = 0._10
         endif
c
c determine p(n), p'(n,h), p(n,h)
         call LEGNDR(slat, clat, Nzone, Ntess, leg, leg1, gleg, gleg1)
c
c zonal harmonic terms for equations of motion
         do jj = 1, 3
            summh(jj,is) = 0._10
            do i = 1, Nzon1
               zone(jj,i,is) = leg1(i)*slat1(jj)/rh(i)/r2
               summh(jj,is) = summh(jj,is) + Zhar(i)*zone(jj,i,is)
            end do
c
c tesseral harmonic terms for equations of motion
            ih = 0
            do i = 1, Ntes1
               sumt = 0._10
               i1   = i + 1
               do jh = 1, i1
                  djh = jh
                  ih  = ih + 1
                  tess1(jj,ih,is) = gleg1(ih)*slat1(jj)/rh(i)/r2
                  tess2(jj,ih,is) = gleg(ih)*lng1(jj)/rh(i)/r2
                  snh1(ih)= Char(ih)*clng(jh,is) + Shar(ih)*slng(jh,is)
                  snh2(ih)= -Char(ih)*slng(jh,is) + Shar(ih)*clng(jh,is)
                  sumt = sumt + snh1(ih)*tess1(jj,ih,is) +
     .             djh*snh2(ih)*tess2(jj,ih,is)
               end do
               summh(jj,is) = summh(jj,is) - sumt
            end do
         end do
c
c
c determine harmonic quantitites for partials
c
         if(Iparmr.gt.1) then
            if(is.eq.2 .and. Kmr(81).lt.1) return
c
c determine p(n), p'(n), p(n,h), p'(n,h)
            call LEGND2(slat, clat, Nzone, Ntess, leg, leg1, leg2, gleg,
     .                  gleg1, gleg2)
c
c partials of zonal harmonic terms wrt coordinates
            do jj = 1, 3
               do l = 1, 3
                  sump(jj,l,is) = 0._10
                  clap2(jj,l)=((delta(2,jj)*slng(1,is)
     .             +delta(1,jj)*clng(1,is))
     .             *(delta(1,l)*slng(1,is)-delta(2,l)*clng(1,is))
     .             + (delta(2,l)*slng(1,is)+delta(1,l)*clng(1,is))
     .             *(delta(1,jj)*slng(1,is)-delta(2,jj)*clng(1,is)))
     .             /r2/clat**2
                  slap2(jj, l) = (3._10*cor(3)*cor(jj)*cor(l)/r2 -
     .                           delta(3,jj)*cor(l) - delta(3,l)*cor(jj)
     .                           - delta(jj,l)*cor(3))/r2/r
                  do i = 1, Nzon1
                     if(i.gt.Nzonp) go to 10
                     di2 = i + 2
                     sump(jj,l,is) = sump(jj,l,is)
     .                             + (r*leg1(i)*slap2(jj,l) + leg2(i)
     .                             *slat1(jj)*slat1(l)/r - di2*leg1(i)
     .                             *slat1(jj)*cor(l)/r2)*Zhar(i)/rh(i)
     .                             /r2
                  end do
c
c partials of tesseral harmonic terms wrt coordinates
   10             if(Ntess.gt.1) then
                     inh = 0
                     do i = 1, Ntes1
                        if(i.gt.Ntesp) go to 20
                        sumpt = 0._10
                        di2   = i + 2
                        i1    = i + 1
                        do ih = 1, i1
                           djh   = ih
                           inh   = inh + 1
                           sumpt = sumpt +
     .                             (di2*cor(l)/r2*(snh1(inh)*gleg1(inh)
     .                             *slat1(jj)+djh*snh2(inh)*gleg(inh)
     .                             *lng1(jj))
     .                             - (snh1(inh)*(gleg1(inh)*slap2(jj,l)
     .                             +(gleg2(inh)*slat1(jj)*slat1(l)
     .                             -djh**2*gleg(inh)*lng1(jj)*lng1(l))
     .                             /r2)+snh2(inh)
     .                             *djh*(gleg1(inh)*slat1(jj)*lng1(l)
     .                             /r2+gleg(inh)*clap2(jj,l)+gleg1(inh)
     .                             *lng1(jj)*slat1(l)/r2))*r)/rh(i)/r2
                        end do
                        sump(jj,l,is) = sump(jj,l,is) + sumpt
                     end do
                  endif
 
   20          end do
            end do
         endif
         if(is.eq.2 .or. Kmr(81).lt.0) return
         is = 2
         do i = 1, 3
            cor(i) = Smcor(i)
         end do
         r  = Srm
         r2 = Srm2
      end do
c
c pick up indirect c22 contribution to beta, gamma, or j2 partial
      entry MORHR2
      jh     = 2
      ii     = 2
      dcnm   = Dc22(Kkk)
c at this point, JH is the order of the cosine tesseral term,
c II is the index counting from C21, and
c DCNM is the partial of the specified coefficient w.r.t. the
c  parameter of interest, i.e., either 1 (for the coefficent itself)
c  or the appropriate partial (for beta, gamma, or J2 indirect terms)
  100 djh    = jh
      ddw1dq = -Gamat3*Mma*(Smecor(2)*(eclng(jh)*tesse1(3,ii)-eslng(jh)*
     .         tesse2(3,ii)*djh) - Smecor(3)
     .         *(eclng(jh)*tesse1(2,ii)-eslng(jh)*tesse2(2,ii)*djh))
      ddw2dq = -Gamat3*Mmb*(Smecor(3)*(eclng(jh)*tesse1(1,ii)-eslng(jh)*
     .         tesse2(1,ii)*djh) - Smecor(1)
     .         *(eclng(jh)*tesse1(3,ii)-eslng(jh)*tesse2(3,ii)*djh))
      ddw3dq = -Gamat3*Mmc*(Smecor(1)*(eclng(jh)*tesse1(2,ii)-eslng(jh)*
     .         tesse2(2,ii)*djh) - Smecor(2)
     .         *(eclng(jh)*tesse1(1,ii)-eslng(jh)*tesse2(1,ii)*djh))
      if(Kmr(81).ge.1) then
         ddw1dq = ddw1dq -
     .            Gamat*Mma*(Smcor(2)*(sclng(jh)*tesss1(3,ii)-sslng(jh)
     .            *tesss2(3,ii)*djh) - Smcor(3)
     .            *(sclng(jh)*tesss1(2,ii)-sslng(jh)*tesss2(2,ii)*djh))
         ddw2dq = ddw2dq -
     .            Gamat*Mmb*(Smcor(3)*(sclng(jh)*tesss1(1,ii)-sslng(jh)
     .            *tesss2(1,ii)*djh) - Smcor(1)
     .            *(sclng(jh)*tesss1(3,ii)-sslng(jh)*tesss2(3,ii)*djh))
         ddw3dq = ddw3dq -
     .            Gamat*Mmc*(Smcor(1)*(sclng(jh)*tesss1(2,ii)-sslng(jh)
     .            *tesss2(2,ii)*djh) - Smcor(2)
     .            *(sclng(jh)*tesss1(1,ii)-sslng(jh)*tesss2(1,ii)*djh))
      endif
 
c correct for indirect effect of beta, gamma, and J2 on C22
c or simply include effect of tesseral cosine harmonics other than C22
      ddw1dp = ddw1dq*dcnm + ddw1dp
      ddw2dp = ddw2dq*dcnm + ddw2dp
      ddw3dp = ddw3dq*dcnm + ddw3dp
c
c
      return
      end
