      subroutine SOLPRB(ncall)
 
      implicit none

c
c        ash / amuchastegui - october 1969 - subroutine solprb
c        computation of forces acting on sun-centered space probe
c     l.friedman's original coding in sbfn put into solprb.
c
c     low thrust forces (gas leakage and radiation pressure) and
c     asteroid ring (inside) from l.friedman's mit thesis.
c
c arguments
      integer*4 ncall
c        ncall=-2 setup once per iteration of a given step for partials
c        ncall=-1 setup once per iteration of a given step for motion
c        ncall= 0 setup once per step
c        ncall= k evaluate acceleration for equation k.gt.0
c
c
c array dimensions
      include 'globdefs.inc'
c        common
      include 'astroi.inc'
      include 'intstf.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'sbrot.inc'
      include 'sbstuf.inc'
      include 'sbthng.inc'
      include 'smlbdy.inc'
      include 'smlstf.inc'
      include 'trghar.inc'
      include 'xprcom.inc'

c local
      real*10 dum,yscor(3),ys3(3),ysc3,dadxsm(3,3),td,tg,tg1,f1,f2
      real*10 rysml,rysml2,rysml3
 
      real*10 sumh(3),sumhp(3,10),sumhs(3),sump(3),sumpp(3),sumps(3)
      real*10 ltftm,ltftm2,termc,termr,dadxj2s(3,3),gf,gf2

      real*10 dist,pfact,pfact2,psumh,psumlt,psumr,
     . q,rysb2,rysb3,sumpkk,sunlam,fn1sav
      real*10 ngcut,ngfct,dngfct,ngfct1,ngfct2,thhat(3),
     . dadxng(3,3),dadxngv(3,3)
      integer*4 i,iden,ikkk,j,k,l,ntarg

c external functions
      real*10 DOT

      if(ncall.eq.0) then
c
c      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c set-up once per step
c   general realtivity setup
         if(Kp(61).eq.0 .or. Kp(61).eq.1) call SOLGR(ncall)
 
c compute limited asteroid coordinates
         if(Kp(30).ge.0 .and. Numsml.gt.0) then
            do i = 1,Numsml
               call JLIPT((Jd-Jdsml0(i))+Fract,Elptsm(1,i),0,
     .          Ysml(1,i),rysml,rysml2,rysml3,dum)
               do j = 1, 3
                  Ysml3(j,i) = Ysml(j,i)/rysml3
               end do
            end do
         endif

c effect of second harmonic of sun due to other planets
         if(Kp(63).ge.0) then
            do j=1,3
               sumhs(j)=0._10
            end do
            do l=1,10
               if(Kpb(l)) then
                  gf  = DOT(C3,Pccor(1,l))/Rpc(l)
                  gf2 = (7.5_10*gf**2 - 1.5_10)/Rpc(l)
                  do j=1,3
                     sumhp(j,l)=(Pccor(j,l)*gf2-3._10*gf*C3(j))/Rpc(l)
     .                /Rpc3(l)
                     sumhs(j)=sumhs(j)+Mass1(l)*sumhp(j,l)
                  end do
               endif
            end do
         endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      else if(ncall.eq.-1) then
c
c set-up once per iteration of a given step for motion
c
c determine if relativity factor is included
         if(Kp(61).ge.0) then
            if(Kp(61).lt.2) call SOLGR(ncall)
 
c central body sun
            Rvsb2 = DOT(Sbcor(4), Sbcor(4))
            Sbvsb = DOT(Sbcor, Sbcor(4))
c gamapm is the generalized metric parameter gamma
c betapm is the generalized metric parameter beta
c einstein general relativity values are 1
            Rfact1 = 0.5_10*(Gamapm+Betapm)*(Alph4/Rsb)
     .       - Gamapm*(Rvsb2/Cvel2)
            Rfact2 = 0.5_10*(Gamapm+Betapm)*(Alph16/Rsb3)
     .       - Gamapm*3._10*Rvsb2/(Cvel2*Rsb2)
            Rfact3 = 0.5_10*(1._10 + Gamapm)*3._10*Sbvsb/Rsb2
            Rfact4 = 0.5_10*(1._10 + Gamapm)*4._10*Sbvsb/Cvel2
         endif
c
c determine if second harmonic of sun is included
         if(Kp(63).ge.0) then
            Gfact  = DOT(C3, Sbcor)/Rsb
            Gfact2 = (7.5_10*Gfact**2 - 1.5_10)/Rsb
         endif
c
c determine if second harmonic of integrated planet included
         if(Nplnt.le.30) then
            if(Kp(82).ge.0) then
               pfact  = (Bodrot(3,1)*Sbcor(1) + Bodrot(3,2)*Sbcor(2)
     .          + Bodrot(3,3)*Sbcor(3))/Rsb
               pfact2 = (7.5_10*pfact**2 - 1.5_10)/Rsb
            endif
         endif
c
c determine if low-thrust forces are included
         if(Kp(81).ge.1) then
            call CROSS(Pccor(1,3), Sbcor, Nvect)
            q = SQRT(DOT(Nvect,Nvect))
            Nvect(1) = Nvect(1)/q
            Nvect(2) = Nvect(2)/q
            Nvect(3) = Nvect(3)/q
            Tvect(1) = (Nvect(2)*Sbcor(3) - Nvect(3)*Sbcor(2))/Rsb
            Tvect(2) = (Nvect(3)*Sbcor(1) - Nvect(1)*Sbcor(3))/Rsb
            Tvect(3) = (Nvect(1)*Sbcor(2) - Nvect(2)*Sbcor(1))/Rsb
c lambda=1,    radiation pressure force  no shadow
c lambda=0, no radiation pressure force     shadow
            Lambda = 1._10
c the shadow coding has not been tested --since it is programmed
c as a discontinuous effect (lamda=0 or 1) it has been found to
c lead to numerical integration difficulties - hence it is
c bypassed for now until more work is done with it ----l.f.
            if( .false. ) then
               do i = 1, Numtar
                  ntarg = Ntrg(i)
                  if(ntarg.le.10 .and. Rpb(ntarg).le.0.005_10) then
                     sunlam = DOT(Sbcor, Pccor(1,ntarg))/Rpc(ntarg)
                     dist   = sunlam - Rpc(ntarg)
                     if(dist.gt.0._10 .and.
     .                Rsb2-sunlam**2.lt.Trad(i)**2) then
                        Lambda = 0._10
                        go to 10
                     endif
                  endif
               end do
            endif
   10       ltftm  = (Jd - Jdp0) + Fract
            ltftm2 = ltftm**2
            Ltf1   = (Con(3) - Ltf11*ltftm - Ltf12*ltftm2)
     .       /Rsb + Lambda*Con(6)/Rsb3
            Ltf2   = Con(4) - Ltf21*ltftm - Ltf22*ltftm2 +
     .       Lambda*Con(7)/Rsb2
            Ltf3   = Con(5) - Ltf31*ltftm - Ltf32*ltftm2 +
     .       Lambda*Con(8)/Rsb2
         endif
c
c determine if distributed asteroidal pertubation is included
         if(Nbelt.gt.0) call ASBFNC(ncall)
c
c determine if limited asteroids are included
         if(Kp(30).ge.0 .and. Numsml.gt.0) then
            do j = 1, 3
               Sumsml(j) = 0._10
               do iden = 1, 5
                  Dsmsml(j,iden) = 0._10
               end do
               do k=1,3
                  dadxsm(j,k) = 0._10
               end do
            end do
            do i = 1, Numsml
               do j = 1, 3
                  yscor(j) = Ysml(j,i) - Sbcor(j)
               end do
               rysb2 = DOT(yscor, yscor)
               rysb3 = rysb2*SQRT(rysb2)
               td    = -Smlmas(i)/rysb3
               tg1   = -3._10*td/rysb2
               do j = 1, 3
                  ysc3     = yscor(j)/rysb3
                  ys3(j)   = ysc3 - Ysml3(j,i)
                  Sumsml(j)= Sumsml(j) + Smlmas(i)*ys3(j)
                  tg       = tg1*yscor(j)
                  dadxsm(j,j)=dadxsm(j,j)+td
                  do k=1,3
                     dadxsm(j,k) =dadxsm(j,k)+tg*yscor(k)
                  end do
               end do
               iden = Denpts(i) - 33
               if(iden.gt.0 .and. iden.le.5) then
                  do j = 1, 3
                     Dsmsml(j,iden) = Dsmsml(j,iden) + Smlvol(i)*ys3(j)
                  end do
               endif
            end do
            do j=1,3
               do k=1,3
                  dadxsm(j,k)= Gamat*dadxsm(j,k)
               end do
            end do
         endif
c
c determine if cometary nongravitational force is included
         if(Kp(83).ge.0) then
            if(Kp(61).lt.0) then
c these quantities are needed now, even if no relativity
               Rvsb2 = DOT(Sbcor(4), Sbcor(4))
               Sbvsb = DOT(Sbcor, Sbcor(4))
            endif
            ngcut=  1._10 + (Rsb/Ngr0)**Ngn
            ngfct=  Ngalph*(Rsb/Ngr0)**(-Ngm) * ngcut**(-Ngk)
            ngfct2= Rsb2*Rvsb2 - Sbvsb**2
            ngfct1= SQRT(ngfct2)*Rsb
            do i=1,3
               thhat(i)= (Rsb2*Sbcor(3+i) - Sbvsb*Sbcor(i))/ngfct1
            end do
         endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else if(ncall.lt.-1) then
c
c set-up once per iteration of a given step for partials
c
c   determine if general relativity effect is included in partials
         if(Kp(61).ge.1) then
            do i = 1, Kount
               Vsbdsb(i) = DOT(Sbcor(4), Dsbcor(4,i))
               Sdsb(i)   = DOT(Sbcor(1), Dsbcor(4,i)) 
     .          + DOT(Dsbcor(1,i), Sbcor(4))
            end do
         endif
c
c determine if second harmonic of sun is included in partials
         if(Kp(63).ge.1) then
            Gfact5 = 15._10*Gfact/Rsb2
            do i=1,3
               termr = -5._10*Sbcor(i)*Gfact2/Rsb2 + Gfact5*
     .          (C3(i)-Gfact*Sbcor(i)/Rsb)
               termc = Gfact5*Sbcor(i) - 3._10*C3(i)/Rsb
               do j=1,3
                  dadxj2s(j,i)=termr*Sbcor(j)+termc*C3(j)
                  if(i.eq.j) dadxj2s(j,i)=dadxj2s(j,i)+Gfact2
               end do
            end do
         endif
c
c determine if low-thrust force parameters are included in partials
         if(Kp(81).ge.1) then
            Aquad  = 1._10 - Con(1)*ltftm - Con(2)*ltftm2
            Ffact1 = Con(3)*ltftm
            Ffact2 = Con(4)*ltftm
            Ffact3 = Con(5)*ltftm
            Ffact4 = Con(3)*ltftm2
            Ffact5 = Con(4)*ltftm2
            Ffact6 = Con(5)*ltftm2
         endif
c
c determine if distributed asteroids is included in partials
c (calculations already done above along with motion expressions)
c
c determine if cometary nongravitational force is included for partials
         if(Kp(83).gt.0) then
            dngfct= -(Ngm + Ngn*Ngk*(1._10-1._10/ngcut))*ngfct/Rsb
c   symmetric part of grav.grad. matrices, both position and velocity
            do i=1,3
               do j=1,i
                  dadxng(i,j)=(Nga1*(dngfct-ngfct/Rsb)
     .             + Nga2*ngfct*Sbvsb/ngfct1)
     .             *Sbcor(i)*Sbcor(j)/Rsb2 
     .             + Nga2*ngfct*Sbvsb/ngfct1 *thhat(i)*thhat(j)
                  dadxng(j,i)=dadxng(i,j)
                  dadxngv(i,j)=-Nga2*ngfct/ngfct1*(Sbcor(i)*Sbcor(j)
     .             + thhat(i)*thhat(j)*Rsb2)
                  dadxngv(j,i)=dadxngv(i,j)
               end do
               dadxng(i,i)=dadxng(i,i) + Nga1*ngfct/Rsb
     .          - Nga2*ngfct*Sbvsb/ngfct1
               dadxngv(i,i)=dadxngv(i,i) + Nga2*ngfct*Rsb2/ngfct1
            end do
c   assymmetric part of grav.grad. matrix (only position)
            do i=1,3
               do j=1,3
                  dadxng(i,j)=dadxng(i,j) + Nga2*
     .             (dngfct*Sbcor(i)*thhat(j)/Rsb
     .             -ngfct*thhat(i)*Sbcor(j)/Rsb2)
               end do
            end do
         endif
c
c      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      else if(Kkk.eq.0) then
c
c accelerations for motion
c
c effect of violation of principle of equivalence on motion
         if(Jpoev.gt.0) then
c for perturbing planets
            sumpp(Kk)=0._10
            sumps(Kk)=0._10
            do l=1,10
               if(Kpb(l)) then
                  sumpp(Kk)=sumpp(Kk) + Mass1(l)*Pbcor3(Kk,l)
                  sumps(Kk)=sumps(Kk) - Mass1(l)*Pccor3(Kk,l)
               endif
            end do
c for central force
            sumpp(Kk) = sumpp(Kk)-Sbforc(Kk)
            sumps(Kk) = sumps(Kk)-Massp*Sbforc(Kk)
            sump(Kk)  = Dltbod*sumpp(Kk)+Dltsun*sumps(Kk)
            Sbfpoe(Kk)= Gamat*(Con1(10)*sumpp(Kk)+Sunpoe*sumps(Kk))
            Fn(1)=Fn(1) + Gamat*sump(Kk)
         endif
c
c effect of general relativity
         if(Kp(61).ge.0) then
            if(Kp(61).eq.2) Sumr(Kk)
     .          = (Rfact1*Sbcor(Kk)+Rfact4*Sbcor(ncall))/Rsb3
            Relacc(Kk)=Rlfact*Gamat*Sumr(Kk)
            Fn(1) = Fn(1)+Relacc(Kk)
         endif
c
c effect of second harmonic of sun
         if(Kp(63).ge.0) then
            sumh(Kk) = (Sbcor(Kk)*Gfact2-3._10*Gfact*C3(Kk))/Rsb4
            Fn(1)    = Fn(1) - Shar2*(Gamat3*sumh(Kk)-Gamat*sumhs(Kk))
         endif
c
c effect of second harmonic of integrated planet
         if(Nplnt.le.30) then
            if(Kp(82).ge.0) then
               sumpkk = (Sbcor(Kk)*pfact2-3._10*pfact*Bodrot(3,Kk))/Rsb4
               Fn(1)  = Fn(1)-Gamat3*Phar2*sumpkk
            endif
         endif
c
c effect of low-thrust forces
         if(Kp(81).ge.1) then
            Sumltf = Ltf1*Sbcor(Kk) + Ltf2*Tvect(Kk) + Ltf3*Nvect(Kk)
            Fn(1)  = Fn(1) + Sumltf
         endif
c
c effect of distributed asteroidal perturbation
         if(Nbelt.gt.0) call ASBFNC(ncall)
c
c effect of limited asteroids
         if(Kp(30).ge.0 .and. Numsml.gt.0) then
            Fn(1) = Fn(1) + Gamat*Sumsml(Kk)
         endif
c
c effect of cometary nongravitational forces
         if(Kp(83).ge.0) then
            Fn(1) = Fn(1) + ngfct*(Nga1*Sbcor(Kk)/Rsb + Nga2*thhat(Kk))
         endif
c
c      - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      else
c
c compute accelerations for partials
         ikkk  = Icntrl(Kkk)
c
c   effect of violation of the principle of equivalence
         if(Jpoev.gt.0) then

c
c check if alpha is time variable gravitational constant
            if(ikkk.eq.32) then
               Fn(1)=Fn(1) + Gama*Tvary*sump(Kk)
c check if alpha is perturbing planet mass or integrated planet mass
            else if(ikkk.ge.1 .and. ikkk.le.10 .and. Kp(ikkk+30).ge.0)
     .          then
               if(ikkk.eq.Nplnt) then
                  Fn(1)=Fn(1) - Gamat*Dltsun*Sbforc(Kk)
               else if(Nplnt.ne.3 .or. ikkk.ne.10) then
                  Fn(1)=Fn(1) + Gamat*(Dltbod*Pbcor3(Kk,ikkk)
     .                              - Dltsun*Pccor3(Kk,ikkk))
               endif
            else if(ikkk.eq.40) then
c deltas
               Fn(1)=Fn(1) + Gamat*sumps(Kk)
            else if(ikkk.eq.-20) then
c deltap
               Fn(1)=Fn(1) + Gamat*sumpp(Kk)
            endif

         endif
c
c   effect of general relativity
         if(Kp(61).ge.1) then
            fn1sav = Fn(1)
            if(Kp(61).eq.2) then
               psumr = (Dsbcor(Kk,Kkk)*Rfact1 + Rfact4*Dsbcor(Kk+3,Kkk)
     .          - Sbcor(Kk)*(Rsbdsb(Kkk)*Rfact2
     .          + 2._10*Gamapm*Vsbdsb(Kkk)/Cvel2)
     .          - 4._10*Sbcor(Kk+3)*(Rfact3*Rsbdsb(Kkk)
     .          - 0.5_10*(1._10+Gamapm)*Sdsb(Kkk))/Cvel2)
     .          /Rsb3
               if(Kkk.eq.1) then
                  f1 = (Dsbcor(Kk,Kkk)*Rfact1
     .             - Sbcor(Kk)*Rsbdsb(Kkk)*Rfact2
     .             - 4._10*Sbcor(Kk+3)*Rfact3*Rsbdsb(Kkk)/Cvel2)/Rsb3
                  f2 = (Rfact4*Dsbcor(Kk+3,Kkk)
     .             - Sbcor(Kk)*2._10*Gamapm*Vsbdsb(Kkk)/Cvel2
     .             + 2._10*Sbcor(Kk+3)*(1._10+Gamapm)*Sdsb(Kkk)/Cvel2)
     .             /Rsb3
               endif
            else
               f1 = DOT(Dsbcor(1,Kkk),Dadxr(1,Kk))
               f2 = DOT(Dsbcor(4,Kkk),Dadvr(1,Kk))
               psumr = f1 + f2
            endif
            psumr = Gamat*Rlfact*psumr
            Fn(1) = Fn(1) + psumr
            if(Kkk.eq.1) then
               Indpsav(Kk,1,2)=Gamat*Rlfact*f1
               Indpsav(Kk+3,1,2)=Gamat*Rlfact*f2
            endif
c
c check if alpha is a planet mass
            if(ikkk.ge.1 .and. ikkk.le.9 .and. Kp(61).lt.2) then
               Fn(1) = Fn(1) + Gamat*Rlfact*Sumrp(Kk,ikkk)
c
c check if alpha is relativity factor
            else if(ikkk.eq.31) then
               Fn(1) = Fn(1) + Gamat*Sumr(Kk)
c
c check if alpha is time variable gravitational constant
            else if(ikkk.eq.32) then
               Fn(1) = Fn(1)+ Gama*Tvary*Rlfact*Sumr(Kk)
c
c check if alpha is metric param beta
            else if(ikkk.eq.41 .or. ikkk.eq.43) then
               if(Kp(61).eq.2) Sumrb(Kk) = 0.5_10*Sbcor(Kk)*Alph4/Rsb4
               Fn(1) = Fn(1) + Gamat*Rlfact*Sumrb(Kk)
 
c beta also gives principle of equivalence violation
               if(ikkk.eq.43) Fn(1) = Fn(1)+ 4._10*Sbfpoe(Kk)
c
c check if alpha is metric param gamma
            else if(ikkk.eq.42 .or. ikkk.eq.44) then
               if(Kp(61).eq.2) Sumrg(Kk)= (Sbcor(Kk)*0.5_10*Alph4/Rsb +
     .            (2._10*Sbcor(Kk+3)*Sbvsb-Sbcor(Kk)*Rvsb2)/Cvel2)/Rsb3
               Fn(1) = Fn(1) + Rlfact*Gamat*Sumrg(Kk)
 
c gamma also gives principle of equivalence violation
               if(ikkk.eq.44) Fn(1) = Fn(1) - Sbfpoe(Kk)
            endif
            if(Kkk.eq.1) Relpar(Kk)=Fn(1)-fn1sav
         endif
c
c effect of second harmonic of sun
         if(Kp(63).ge.1) then
            psumh = DOT(dadxj2s(1,Kk),Dsbcor(1,Kkk))/Rsb4
            Fn(1) = Fn(1) - Gamat3*Shar2*psumh
c
c check if alpha is second harmonic of sun
            if(ikkk.eq.33) then
               Fn(1) = Fn(1) - Sunrd2*(Gamat3*sumh(Kk)-Gamat*sumhs(Kk))
c
c check if alpha is time variable gravitational constant
            else if(ikkk.eq.32) then
               Fn(1) = Fn(1) - Tvary*Shar2*(Gama3*sumh(Kk)
     .          -Gama*sumhs(Kk))
c
c check if alpha is mass of integrated planet
            else if(ikkk.eq.Nplnt .and. Nplnt.le.30) then
               Fn(1) = Fn(1) + Gamat*Shar2*sumh(Kk)
c
c check if alpha is mass of another planet
            else if(ikkk.gt.0 .and. ikkk.le.10 .and. Kpb(ikkk)) then
               termc = Gamat*Shar2*sumhp(Kk,ikkk)
               if(ikkk.eq.10) then
                  termc = termc*Mass(3)
               else if(ikkk.eq.3 .and. Kpb(10)) then
                  termc = termc*Masse
               endif
               Fn(1) = Fn(1) + termc
            endif
         endif
c
c effect of low-thrust forces on partial derivatives
c the effect on position and velocity is not included
         if(Kp(81).ge.1) then
            if(ikkk.lt.0) then
               psumlt = 0._10
               if(ikkk.eq.-1) psumlt = psumlt -
     .             Ffact1*Sbcor(Kk)/Rsb - Ffact2*Tvect(Kk)
     .             - Ffact3*Nvect(Kk)
               if(ikkk.eq.-2) psumlt = psumlt -
     .             Ffact4*Sbcor(Kk)/Rsb - Ffact5*Tvect(Kk)
     .             - Ffact6*Nvect(Kk)
               if(ikkk.eq.-3) psumlt = psumlt +
     .             Aquad*Sbcor(Kk)/Rsb
               if(ikkk.eq.-4) psumlt = psumlt + Aquad*Tvect(Kk)
               if(ikkk.eq.-5) psumlt = psumlt + Aquad*Nvect(Kk)
               if(Lambda.ne.0._10) then
                  if(ikkk.eq.-6) psumlt = psumlt + Sbcor(Kk)/Rsb3
                  if(ikkk.eq.-7) psumlt = psumlt + Tvect(Kk)/Rsb2
                  if(ikkk.eq.-8) psumlt = psumlt + Nvect(Kk)/Rsb2
               endif
               Fn(1) = Fn(1) + psumlt
            endif
         endif
c
c effect of distributed asteroidal perturbation on partials
c the effect on position and velocity is not included
         if(Nbelt.gt.0) call ASBFNC(ncall)
c
c effect of limited asteroids
         if(Kp(30).gt.0 .and. Numsml.gt.0) then
            Fn(1)= Fn(1)+ DOT(Dsbcor(1,Kkk),dadxsm(1,Kk))
c
c check if alpha is asteroid class density
            iden = ikkk - 33
            if(iden.gt.0 .and. iden.le.5) then
               Fn(1) = Fn(1) + Gamat*Dsmsml(Kk,iden)
c
c check if alpha is time variable gravitational constant
            else if(ikkk.eq.32) then
               Fn(1) = Fn(1)+ Gama*Tvary*Sumsml(Kk)
            endif
         endif
c
c effect of cometary nongravitational forces
         if(Kp(83).gt.0) then
            Fn(1)=Fn(1) + DOT(Dsbcor(1,Kkk),dadxng(1,Kk)) 
     .       + DOT(Dsbcor(4,Kkk),dadxngv(1,Kk))
c check if alpha is radial coefficient
            if(ikkk.eq.-14) Fn(1)=Fn(1)+ngfct*Sbcor(Kk)/Rsb
c check if alpha is tangential coefficient
            if(ikkk.eq.-15) Fn(1)=Fn(1)+ngfct*thhat(Kk)
         endif
      endif
 
      return
      end
