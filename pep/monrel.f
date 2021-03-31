      subroutine MONREL(ncall)
 
      implicit none
c
c     j.f.chandler - sep 2011
c
c     compute general relativistic correction to force on moon.
c     based on eq. (6.34) of will, theory and experiment in
c     gravitational physics, cambridge univ. press, 1981.
c     (actually this equation is in isotropic rather than harmonic
c     coordinates, but these are the same to 1st post-newtonian order).
c
c     this incorporates parts of the 1977 version of King and Cappallo,
c     stemming from Slade's monfn routine
c
c arguments
      integer ncall
c
c        ncall=-2 setup once per iteration of a given step for partials
c                 (not used - already done with motion)
c        ncall=-1 setup once per iteration of a given step for motion
c        ncall= 0 setup once per step
c        ncall= k evaluate acceleration for equation k.gt.0
c
c array dimensions
      include 'globdefs.inc'

c commons
c
c
      include 'bddtaint.inc'
      include 'intstf.inc'
      include 'metuna.inc'
      include 'morcrd.inc'
c rvec(4-6) not filled in, except for earth-sun-moon combinations
c ecor  = earth relative to sun (re is distance) = rvec(.,3,10)
c mcor  = moon relative to sun (rm is distance) = rvec(.,11,10)
c mecor = moon relative to earth (rem is distance) = rvec(.,11,3)
c mcor etc. are equatorial coordinates
c smcor etc. are selenodetic coordinates
      include 'morstf.inc'
      include 'param.inc'
      include 'prtcod.inc'
      include 'tapdtplp.inc'
      include 'xprcom.inc'
      include 'yvectplp.inc'

c external functions
      real*10 DOT,DOTN,VECMTPRD
c
c local
      real*10 bp1,f1,f2,f3,f4,f5,f6,f7,gp1,sum,tsum,tsumx,f4r2,df1,df2
      real*10 fx,tgv,tgva,tgvb,tgvx,
     . vbrab,sumr(3),sumrg(3),sumrb(3),tsump,sumrp(3,11),dfdp,
     . rbc3,vcrab,vavc,dva(3),dvb(3),mf,tsuma,df1d,df2d,tgab,tgbc,tgac,
     . tguu,tduu,df1dc,df2dc,df3dc,xij(11,11),mfx,vab(3),
     . dfdxb1(3),dfdxb2(3),dfdxb3(3),tgpb,tgpc,tgp,tdp,
     . tgp1,tgp3,sumt
      real*10 sumrx(42)
      equivalence (sumrx(1),sumr),(sumrx(4),sumrg),(sumrx(7),sumrb),
     . (sumrx(10),sumrp)
      integer a,b,c,i,j,jj,k,kt,kkk1,avals(2)/3,11/,loop,icmkkk
c
c     note on subscripts in this routine:
c         a,b,c (integers) match formula (6.34) or (6.78) in will
c         a is equivalenced to nplnt
c         1-9 refer to planets (3 = earth)
c         10  is sun, not moon
c         11  is integrated body (moon)
c
c     However, after computing sumrp (partials w.r.t. masses of these
c     bodies) in terms of the above subscript convention, the results
c     for 3 and 11 are converted into partials w.r.t. MASS(3) and
c     MASS(10) and stored in the 3rd and 10th slots.  Similarly,
c     the Dadxrp and Dadvrp arrays (partials w.r.t. coordinates) are
c     converted for 3 and 11 (Earth and Moon vs Sun) into 3 and 10
c     (Embary vs Sun and Moon vs Earth).
c
      if(ncall.lt.-1) then
         return
      else if(ncall.lt.0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     set up once per iteration of a given step for motion
c
         do i = 1,3
            v(i,11)= Mcor(i+3) + v(i,10)
            v(i,3) = Ecor(i+3) + v(i,10)
         end do

         do i = 1,10
            if(Kgr(i).ge.0) then
               vab2(i,11) = DOT(v(1,i),v(1,11))
               vab2(11,i) = vab2(i,11)
               vab2(i,3) = DOT(v(1,i),v(1,3))
               vab2(3,i) = vab2(i,3)
            endif
         end do

         do i=1,11
            if(Kgr(i).ge.0) then
               do j=1,i
                  if(Kgr(j).ge.0) then
                     xij(i,j)= DOT(Rvec(1,i,j),Mecor)
                     if(j.lt.i) xij(j,i)= -xij(i,j)
                  endif
               end do
            endif
         end do
  
         vab2(11,11) = DOT(v(1,11),v(1,11))
 
         call ZFILL(sumrx,16*42)
         call ZFILL(Dadxrp,16*198)

         do loop=1,2
            a=avals(loop)
            do b = 1,11
               if(b.ne.a .and. Kgr(b).ge.0) then
                  f4   = Xmass(b)/Rab(a,b)
                  f5   = Xmass(a)/Rab(a,b)
                  tsum = B2g2*f4 + B2g21*f5
                  f4r2 = f4/Rab(a,b)**2
c     also set up for partials
                  bp1  = 2._10*(f4 + f5)*A44
                  gp1  = bp1 -(vab2(a,a)-2._10*vab2(a,b)+vab2(b,b))
     .             /Cvel2
                  do i=1,3
                     dfdxb1(i)=0._10
                     dfdxb2(i)=0._10
                     dfdxb3(i)=0._10
                  end do
                  f1   = 0._10
                  f2   = 0._10
                  f3   = 0._10
                  tgab = 0._10
                  tguu = 0._10
                  do c = 1,11
                     if(Kgr(c).ge.0 .and. c.ne.a .and. c.ne.b) then
                        df1= Xmass(c)/Rab(b,c)
                        df2= Xmass(c)/Rab(a,c)
                        f7 = DOT(Rvec(1,a,b),Rvec(1,b,c))/Rab(b,c)**2
                        f1 = f1 + df1
                        f2 = f2 + df2
                        f3 = f3 + df1*f7
                        df1d = df1/Rab(b,c)**2
                        df2d = df2/Rab(a,c)**2
                        df1dc= (B2m1-1.5_10*f7)*df1d
                        df2dc= B2g2*df2d
                        df3dc= 0.5_10*df1d
                        tgac=-df2dc
                        tgbc=-df3dc
                        tgpc=df3dc-df1dc
                        if(b.eq.3 .or. b.eq.11) then
                           tgbc=tgbc-df1dc+df3dc
                           tgab=tgab-df3dc
                        endif
                        if(c.ne.10) then
                           dfdp=((B2m1-f7/2._10)/Rab(b,c) +
     .                      B2g2/Rab(a,c))*A44*f4r2
                           do i=1,3
                              sumrp(i,c)= sumrp(i,c) + dfdp*Rvec(i,a,b)
                           end do
                           if(c.eq.3 .or. c.eq.11) then
                              tgac=tgac+df2dc
                              tgbc=tgbc+df1dc
                              tgab=tgab+df3dc
                           endif
                        endif
                        tguu= tguu + tgac*xij(a,c) + tgbc*xij(b,c)
                        tdp=G7*A44*f4*df1d
                        do i=1,3
                           dfdxb1(i)=dfdxb1(i)+tgpc*Rvec(i,b,c)-
     .                      df3dc*Rvec(i,a,b)
                           dfdxb2(i)=dfdxb2(i)+df1d*Rvec(i,b,c)
                           dfdxb3(i)=dfdxb3(i)+df3dc*Rvec(i,b,c)+
     .                      df2dc*Rvec(i,a,c)
                           tgp=(df1dc*Rvec(i,b,c)+df3dc*Rvec(i,a,b)+
     .                      df2dc*Rvec(i,a,c))*f4r2*A44
                           tgpb=3._10*tdp*Rvec(i,b,c)/Rab(b,c)**2
                           Dadxrp(i,i,c)=Dadxrp(i,i,c)+tdp
                           Dadxrp(i,i,b)=Dadxrp(i,i,b)-tdp
                           do j=1,3
                              Dadxrp(i,j,c)=Dadxrp(i,j,c)+
     .                         tgp*Rvec(j,a,b)-tgpb*Rvec(j,b,c)
                              Dadxrp(i,j,b)=Dadxrp(i,j,b)+
     .                         tgpb*Rvec(j,b,c)
                           end do
                        end do
                     endif
                  end do

                  vbrab= DOT(v(1,b),Rvec(1,a,b))/Rab(a,b)
                  f6   = B2m1*f1 + B2g2*f2 - f3/2._10
                  f7   = (G22*vab2(a,b)-G11*vab2(b,b)-Gamapm*vab2(a,a)
     .             + 1.5_10*vbrab**2) / Cvel2
                  tsumx= (tsum + f6)*A44 + f7
                  tsump= f4r2*(B2g2*A44/Rab(a,b)) + tsumx/Rab(a,b)**3
                  tsuma= f4r2*(B2g21*A44/Rab(a,b))
                  tdp=f4r2*tsumx
                  do i=1,3
                     Dadxrp(i,i,b)=Dadxrp(i,i,b)-tdp
                     Dadxrp(i,i,a)=Dadxrp(i,i,a)+tdp
                     tgp=((A44*tsum+3._10*tsumx)*Rvec(i,a,b)/Rab(a,b)**2
     .                 - 3._10*vbrab*(v(i,b)-vbrab*Rvec(i,a,b)/Rab(a,b))
     .                /Cvel2/Rab(a,b))*f4r2
                     tgp1=tgp+A44*f4r2*dfdxb1(i)
                     tgp3=tgp+A44*f4r2*dfdxb3(i)
                     do j=1,3
                        Dadxrp(i,j,b)=Dadxrp(i,j,b)+tgp1*Rvec(j,a,b)
     .                   - G7*A44*f4r2*dfdxb2(j)*Rvec(i,a,b)
                        Dadxrp(i,j,a)=Dadxrp(i,j,a)-tgp3*Rvec(j,a,b)
     .                   + G7*A44*f4r2*dfdxb2(j)*Rvec(i,a,b)
                     end do
                  end do
                  tguu = (tguu + tgab*xij(a,b))*A44
                  if(b.ne.3 .and. b.ne.11) tguu= tguu -
     .             (3._10*tsumx+tsum*A44)*xij(a,b)/Rab(a,b)**2
                  tguu = tguu*f4r2/Mass(3)
                  tduu = tsumx*f4r2/Mass(3)
                  do i=1,3
                     sumrp(i,3) = sumrp(i,3)+ tguu*Mass(10)*Rvec(i,a,b)
                     sumrp(i,11)= sumrp(i,11)- tguu*Masse*Rvec(i,a,b)
                     if(b.ne.3 .and. b.ne.11) then
                        sumrp(i,3) = sumrp(i,3)+ tduu*Mass(10)*Mecor(i)
                        sumrp(i,11)= sumrp(i,11)- tduu*Masse*Mecor(i)
                     endif 
                  end do
                  fx   = -2._10*f4r2/Mascnt/Cvel2
                  do c=1,11
                     if(c.ne.10 .and. Kgr(c).ge.0) then
                        vavc =vab2(a,c)
                        if(c.eq.3) then
                           mf=-Mass(10)/Mass(3)
                        else if(c.eq.11) then
                           mf=Masse/Mass(3)
                        endif
                        do i=1,3
                           dvb(i)= v(i,c)
                           if(c.eq.3 .or. c.eq.11) then
                              dvb(i)=Xpert(i+3,3)+v(i,10)
                              if(b.eq.3 .or. b.eq.11) dvb(i)= dvb(i) +
     .                         Mecor(i+3)*mf*Mascnt
                           endif
                        end do
                        vcrab=DOT(dvb,Rvec(1,a,b))/Rab(a,b)
                        if(c.eq.3 .or. c.eq.11) then
                           vavc = Mass(10)*vab2(a,11)+Masse*vab2(a,3) +
     .                      mf*Mascnt*(vab2(a,11)-vab2(a,3))
                           if(b.ne.3 .and. b.ne.11) then
                              vcrab= vcrab + Mascnt*mf/Rab(a,b) *
     .                         (DOT(v(1,b),Mecor)
     .                         -vbrab*xij(a,b)/Rab(a,b))
                              vavc = vavc - G11*Mascnt*mf*
     .                         (vab2(a,11)-vab2(a,3)
     .                         +vab2(b,3)-vab2(b,11))
                           endif
                        endif
                        dfdp = fx*(vavc+1.5_10*vbrab*vcrab)
                        do i=1,3
                           sumrp(i,c)= sumrp(i,c)+dfdp*Rvec(i,a,b)
                        end do
                     endif
                  end do

                  do i = 1,3
                     sumr(i)  = sumr(i) + tsumx*f4r2*Rvec(i,a,b)
                     if(b.ne.10) sumrp(i,b) = sumrp(i,b)
     .                + tsump*Rvec(i,a,b)
                     sumrp(i,a)= sumrp(i,a) + tsuma*Rvec(i,a,b)
                     sumrb(i) = sumrb(i) + (bp1 + 2._10*(f1+f2)*A44)
     .                *f4r2*rvec(i,a,b)
                     sumrg(i) = sumrg(i) + (gp1 + 2._10*f2*A44)
     .                *f4r2*rvec(i,a,b)
                  end do

                  f6 = f4*A44
                  tduu= 0._10
                  do c = 1,11
                     if(Kgr(c).ge.0 .and. c.ne.a .and. c.ne.b) then
                        rbc3 = Rab(b,c)**3
                        f7 = f6*Xmass(c)/rbc3
                        do i = 1,3
                           sumr(i)  = sumr(i) - f7*rvec(i,b,c)*G7
                           sumrg(i) = sumrg(i) - f7*rvec(i,b,c)*2._10
                        end do
                        if(b.eq.3 .or. b.eq.11) then
                           tduu= tduu - f7
                        else if(c.eq.3 .or. c.eq.11) then
                           tduu= tduu + f7
                        endif
                        tguu= 0._10
                        if(b.ne.3 .and. b.ne.11) tguu= tguu+
     .                   xij(a,b)/Rab(a,b)**2
                        if(c.ne.10) then
                           dfdp=-f6*G7/rbc3
                           do i=1,3
                              sumrp(i,c)= sumrp(i,c)+dfdp*Rvec(i,b,c)
                           end do
                           if(c.eq.3 .or. c.eq.11) tguu=tguu-3._10*
     .                      xij(b,c)/Rab(b,c)**2
                        endif
                        if(b.ne.10) then
                           dfdp=-A44/Rab(a,b)*G7*Xmass(c)/rbc3
                           do i=1,3
                              sumrp(i,b)= sumrp(i,b)+dfdp*Rvec(i,b,c)
                           end do
                           if(b.eq.3 .or. b.eq.11) then
                              tguu=tguu+3._10*xij(b,c)/Rab(b,c)**2
                           endif
                        endif
                        tguu= tguu*f7*G7/Mass(3)
                        do i=1,3
                           sumrp(i,3) = sumrp(i,3)+ tguu*Mass(10)
     .                      *Rvec(i,b,c)
                           sumrp(i,11)= sumrp(i,11)- tguu*Masse
     .                      *Rvec(i,b,c)
                        end do
                     endif
                  end do
                  tduu= tduu*G7/Mass(3)
                  do i=1,3
                     sumrp(i,3) = sumrp(i,3)+tduu*Mass(10)*Mecor(i)
                     sumrp(i,11)= sumrp(i,11)-tduu*Masse*Mecor(i)
                  end do
               endif
 
            end do
 
            do b = 1,11
               if(b.ne.a .and. Kgr(b).ge.0) then
                  vbrab= DOT(v(1,b),Rvec(1,a,b))/Rab(a,b)
                  f1= DOT(Rvec(1,a,b),v(1,a))
                  f2= DOT(Rvec(1,a,b),v(1,b))
                  f6= Rab(a,b)**3*Cvel2
                  fx= Xmass(b)/f6
                  f3= (G22*f1-G21*f2)*fx
                  f4= 2._10*(f1 - f2)*fx
                  do i = 1,3
                     vab(i)= v(i,a)-v(i,b)
                     sumr(i)= sumr(i) + f3*vab(i)
                     sumrg(i)= sumrg(i) + f4*vab(i)
                  end do
                  do i=1,3
                     tgp= 3._10*f3*Rvec(i,a,b)/Rab(a,b)**2 -
     .                fx*(G22*vab(i)+v(i,b))
                     vcrab = 3._10*vbrab*Rvec(i,a,b)/Rab(a,b)
                     tgvb  = (G22*vab(i) + vcrab)*fx
                     tgva= (-G22*vab(i) + 2._10*v(i,a))*fx
                     tgv= -G21*Rvec(i,a,b)*fx
                     tgvx= G22*Rvec(i,a,b)*fx
                     Dadvrp(i,i,b)=Dadvrp(i,i,b)-f3
                     Dadvrp(i,i,a)=Dadvrp(i,i,a)+f3
                     do j=1,3
                        Dadxrp(i,j,b)=Dadxrp(i,j,b)+vab(j)*tgp
                        Dadxrp(i,j,a)=Dadxrp(i,j,a)-vab(j)*tgp
                        Dadvrp(i,j,b)=Dadvrp(i,j,b)+Rvec(j,a,b)*tgvb+
     .                   vab(j)*tgv
                        Dadvrp(i,j,a)=Dadvrp(i,j,a)+Rvec(j,a,b)*tgva+
     .                   vab(j)*tgvx
                        do c=1,11
                           if(Kgr(c).ge.0 .and. c.ne.10)
     .                      Dadvrp(i,j,c)=Dadvrp(i,j,c) -
     .                      (Rvec(j,a,b)*(2._10*v(i,a) + vcrab) +
     .                      vab(j)*Rvec(i,a,b))*Xmass(c)/Mascnt*fx
                        end do
                     end do
                  end do
                  if(b.ne.10) then
                     dfdp= (G22*f1-G21*f2)/f6
                     do i=1,3
                        sumrp(i,b)= sumrp(i,b)+dfdp*vab(i)
                     end do
                  endif
                  do c=1,11
                     if(c.ne.10 .and. Kgr(c).ge.0) then
                        do i=1,3
                           dva(i)= v(i,c)
                        end do
                        if(c.eq.3 .or. c.eq.11) then
                           if(c.eq.3) then
                              mf=-Mass(10)/Mass(3)
                           else
                              mf=Masse/Mass(3)
                           endif
                           if(b.ne.3 .and. b.ne.11) then
                              mfx= G22*mf*Mascnt
                           else
                              mfx= mf*Mascnt
                           endif
                           dfdp= -f3*mf
                           do i=1,3
                              dva(i)= Xpert(i+3,3)+v(i,10)
                              dva(i)= dva(i) + Mecor(i+3)*mfx
                              if(b.ne.3 .and. b.ne.11) then
                                 sumrp(i,c)= sumrp(i,c)+dfdp*Mecor(i+3)
                              endif
                           end do
                        endif
                        dfdp=-fx*DOT(Rvec(1,a,b),dva)/Mascnt
                        if((c.eq.3 .or. c.eq.11) .and.
     .                   b.ne.3 .and. b.ne.11)
     .                   dfdp= dfdp + fx*mf*(3._10*(G22*f1-G21*f2)*
     .                   xij(a,b)/Rab(a,b)**2 -
     .                   (G22*DOT(Mecor,v(1,a))-G21*DOT(Mecor,v(1,b))))
                        do i=1,3
                           sumrp(i,c)= sumrp(i,c)+dfdp*vab(i)
                        end do
                     endif
                  end do
               endif
            end do
c
c subtract acceleration of earth
c
            if(loop.eq.1) then
               do i=1,3
                  sumr(i)= -sumr(i)
                  sumrg(i)= -sumrg(i)
                  sumrb(i)= -sumrb(i)
                  do b=1,11
                     sumrp(i,b)= -sumrp(i,b)
                  end do
                  do j=1,3
                     do b=1,11
                        if(Kgr(b).ge.0) then
                           Dadxrp(i,j,b) = -Dadxrp(i,j,b)
                           Dadvrp(i,j,b) = -Dadvrp(i,j,b)
                        endif
                     end do
                  end do
               end do
            else
               do i=1,3
                  sumrp(i,10)=Mass(3)*(sumrp(i,11)-sumrp(i,3))
                  sumrp(i,3) =Masse*sumrp(i,3)+Mass(10)*sumrp(i,11)
                  do j=1,3
                     do b=1,11
                        if(Kgr(b).ge.0) then
                           Dadxrp(i,j,b) = Gamat*Relmon*Dadxrp(i,j,b)
                           Dadvrp(i,j,b) = Gamat*Relmon*Dadvrp(i,j,b)
                        endif
                     end do
                     Dadxr(i,j) = -Mass(10)*Dadxrp(i,j,3)+
     .                Masse*Dadxrp(i,j,11)
                     Dadvr(i,j) = -Mass(10)*Dadvrp(i,j,3)+
     .                Masse*Dadvrp(i,j,11)
                     Dadxrp(i,j,3) = Dadxrp(i,j,3)+Dadxrp(i,j,11)
                     Dadvrp(i,j,3) = Dadvrp(i,j,3)+Dadvrp(i,j,11)
                  end do
               end do

               if(Kkm(7).gt.0) then
c save quantities for extra printout, starting with rel. accel.
                  do i=1,3
                     do jj=1,3
                        Dadxsav(i,jj,1,1,2)=Dadxr(i,jj)
                        Dadxsav(i,jj,2,1,2)=Dadvr(i,jj)
                        Dadxsav(i,jj,1,3,2)=Dadxrp(i,jj,3)
                        Dadxsav(i,jj,2,3,2)=Dadvrp(i,jj,3)
                        if(Numtar.gt.0) then
                           b=Nplpt(1)
                           Dadxsav(i,jj,1,2,2)=Dadxrp(i,jj,b)
                           Dadxsav(i,jj,2,2,2)=Dadvrp(i,jj,b)
                        endif
                     end do
                  end do
               endif
            endif
         end do
         return
      else if(ncall.eq.0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c     setup once per step
c         note: xpert(i,nplnt) is filled only once per step
c
c     compute velocity of sun relative to cm of solar system
c     shortcut: combine earth and moon by using xpert(.,3) and mass(3)
c
         do i = 1,3
            sum = 0._10
            do j = 1,9
               if(Kgr(j).ge.0) sum = sum - Xpert(i+3,j)*Mass(j)
            end do
            v(i,10) = sum/Mascnt
         end do
 
         do k = 1,9
            if(k.ne.3 .and. Kgr(k).ge.0) then
               do j = k,9
                  if(j.ne.3 .and. Kgr(j).ge.0) then
                     do i = 1,3
                        Rvec(i,j,k) = Xpert(i,j) - Xpert(i,k)
                        Rvec(i,k,j) = -Rvec(i,j,k)
                     end do
                     Rab(j,k) = SQRT(DOT(Rvec(1,j,k),Rvec(1,j,k)))
                     Rab(k,j) = Rab(j,k)
                  endif
               end do
            endif
         end do
 
         do j = 1,9
            if(j.ne.3 .and. Kgr(j).ge.0) then
               do i = 1,3
                  v(i,j) = Xpert(i+3,j) + v(i,10)
                  Rvec(i,j,10) = Xpert(i,j)
                  Rvec(i,10,j) = -Xpert(i,j)
               end do
               Rab(j,10) = SQRT(DOT(Rvec(1,j,10),Rvec(1,j,10)))
               Rab(10,j) = Rab(j,10)
            endif
         end do

         do i = 1,10
            if(i.ne.3 .and. Kgr(i).ge.0) then
               do j = i,10
                  if(j.ne.3 .and. Kgr(j).ge.0) then
                     vab2(i,j) = DOT(v(1,i),v(1,j))
                     if(i.ne.j) then
                        vab2(j,i) = vab2(i,j)
                     endif
                  endif
               end do
            endif
         end do
         return
      else
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c effect of general relativity on the motion of the moon
         if(Kkk.le.0) then
            Relacc(Kk)=Gamat*Relmon*sumr(Kk)
            Fn(1) = Fn(1) + Relacc(Kk)
            return
         endif
c
c
c effect of general relativity on partial derivatives
         f1=DOT(Dmcor(1,Kkk),Dadxr(1,Kk))
         f2=DOT(Dmcor(4,Kkk),Dadvr(1,Kk))
         sumt=f1+f2
         if(Kkk.eq.1) then
            Indpsav(Kk,1,2)=f1
            Indpsav(Kk+3,1,2)=f2
         endif
         do kt=0,Numtar
            b=Nplpt(kt)
            kkk1=Kpt(Kkk,kt)-1
            if(b.le.9 .and. Kgr(b).ge.0 .and. kkk1.gt.0) then
               f1=DOT(Dtcor(1,kkk1,kt),Dadxrp(1,Kk,b))
               f2=DOT(Dtcor(4,kkk1,kt),Dadvrp(1,Kk,b))
               sumt=sumt+f1+f2
               if(Kkk.eq.1) then
                  if(kt.eq.0) then
                     Indpsav(Kk,3,2)=f1
                     Indpsav(Kk+3,3,2)=f2
                  else if(kt.eq.1) then
                     Indpsav(Kk,2,2)=f1
                     Indpsav(Kk+3,2,2)=f2
                  endif
               endif
            endif
         end do

         icmkkk=Icmtrl(Kkk)
c partial derivative w.r.t. RELFCT
         if(icmkkk.eq.31) then
            sumt = sumt + Gamat*sumr(Kk)
c partial derivative w.r.t. GMVARY
         else if(icmkkk.eq.32) then
            sumt = sumt + Gama*Tvary*sumr(Kk)
c
c partial derivative w.r.t. metric parameter beta
         else if(icmkkk.eq.41 .or. icmkkk.eq.43) then
            sumt = sumt + Gamat*Relmon*sumrb(Kk)
c
c partial derivative w.r.t. metric parameter gamma
         else if(icmkkk.eq.42 .or. icmkkk.eq.44) then
            sumt = sumt + Gamat*Relmon*sumrg(Kk)
c
c partial derivative w.r.t. planet mass
         else if(icmkkk.ge.1 .and. icmkkk.le.10) then
            sumt = sumt + Gamat*Relmon*sumrp(Kk,icmkkk)
         endif
         Fn(1)=Fn(1)+sumt
         if(Kkk.eq.1) Relpar(Kk)=sumt
c
         return
      endif
      end
