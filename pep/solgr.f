      subroutine SOLGR(ncall)
 
      implicit none
c
c     r.w.babcock - april 1984
c
c     compute general relativistic correction to force on sun-centered
c     body.  based on eq. (6.34) of will, theory and experiment in
c     gravitational physics, cambridge univ. press, 1981.
c     (actually this equation is in isotropic rather than harmonic
c     coordinates, but these are the same to 1st post-newtonian order).
c
c        ncall=-2 setup once per iteration of a given step for partials
c        ncall=-1 setup once per iteration of a given step for motion
c        ncall= 0 setup once per step
c        ncall= k evaluate acceleration for equation k.gt.0
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'bddtaint.inc'
      include 'intstf.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'prtcod.inc'
      include 'sbembr.inc'
      include 'sbstuf.inc'
      include 'sbthng.inc'
      include 'xprcom.inc'
c
c local
      real*10 bp1,DOT,f1,f2,f3,f4,f5,f6,f7,gp1,sum,tsum
      real*10 rab(0:11,0:11),rvec(3,0:11,0:11),v(3,0:11),
     . vab2(0:11,0:11),xpcop(6,10),df1,df1d,df1dc,df2,df2d,df2dc,df3dc,
     . df6dx(3),df7ab,df7bc,dfx,dfv,dfdp,dva(3),dvb(3),f4r2,fx,fxa,rbc3,
     . sumrt(3,3),sumrbt(3,3),sumrgt(3,3),sumrpt(3,11,3),
     . dadxrt(3,3,3),dadvrt(3,3,3),
     . tdv,tdx,tg,tga,tgb,tgc,tgv,tgva,tgvb,tgvx,tgx,tgxv,tsuma,
     . tsump,tsumx,vab(3),vavc,vbrab,vcrab
      integer  a,b,c,i,j,k,lbd,ncall,nbods, avals(3)/0,11,10/
      logical doterm
      real*10 arrays(180)
      equivalence (arrays(1),sumrt),(arrays(10),sumrbt),
     . (arrays(19),sumrgt),(arrays(28),sumrpt),(arrays(127),dadxrt),
     . (arrays(154),dadvrt)
c
c     note on subscripts in this routine:
c         a,b,c (integers) match formula (6.34) or (6.78) in will
c          0  is sun
c         1-9 refer to planets (3 = embary, or earth if kgr(10).ge.0)
c         10  is moon if kgr(10).ge.0, else ignored
c         11  is integrated body (will duplicate if a planet)
c
      if(ncall.lt.0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c* start=500
c
c     set up once per iteration of a given step for motion
c
c if embary is integrated, both earth and moon need to be calculated
c each iteration of a given step,
c but the earth-moon vector is just once per step
         if(Nplnt.eq.3 .and. Kgr(10).ge.0) then
            nbods=3
            do i = 1,3
               v(i,11) = Sbcor(i+3)-Mass(10)*Xpert(i+3,10) + v(i,0)
               v(i,10) = v(i,11)+Xpert(i+3,10)
               rvec(i,11,0)= Xes(i)
               rvec(i,10,0)= Xms(i)
               rvec(i,0,11)= -rvec(i,11,0)
               rvec(i,0,10)= -rvec(i,10,0)
               do j = 1,9
                  if(Kgr(j).ge.0) then
                     rvec(i,j,11)= Xpes(i,j)
                     rvec(i,j,10)= Xpms(i,j)
                     rvec(i,11,j)= -rvec(i,j,11)
                     rvec(i,10,j)= -rvec(i,j,10)
                  endif
               end do
            end do
            rab(11,0)= Res
            rab(0,11)= Res
            rab(10,0)= Rms
            rab(0,10)= Rms
            do j = 1,9
               if(Kgr(j).ge.0) then
                  rab(j,11)= Rpes(j)
                  rab(11,j)= Rpes(j)
                  rab(j,10)= Rpms(j)
                  rab(10,j)= Rpms(j)
               endif
            end do
         else
            nbods=2
            do i = 1,3
               v(i,11) = Sbcor(i+3) + v(i,0)
               rvec(i,11,0) = Sbcor(i)
               rvec(i,0,11) = -Sbcor(i)
               do j = 1,10
                  if(Kgr(j).ge.0) then
                     rvec(i,j,11) = Pbcor(i,j)
                     rvec(i,11,j) = -rvec(i,j,11)
                  endif
               end do
            end do
            rab(11,0)= Rsb
            rab(0,11)= Rsb
            do j = 1,10
               if(Kgr(j).ge.0) then
                  rab(j,11)= Rpb(j)
                  rab(11,j)= Rpb(j)
               endif
            end do
         endif
 
         do i = 0,10
            if(Kgr(i).ge.0) then
               vab2(i,11) = DOT(v(1,i),v(1,11))
               vab2(11,i) = vab2(i,11)
            endif
         end do
         vab2(11,11) = DOT(v(1,11),v(1,11))
 
         call ZFILL(arrays,16*180)
c calculate acceleration for 2 or 3 bodies: sun+planet or sun+earth+moon
         do lbd=1,nbods
            a=avals(lbd)
            do b = 0,11
               if(b.ne.a .and. Kgr(b).ge.0) then
                  f4   = Xmass(b)/rab(a,b)
                  f5   = Xmass(a)/rab(a,b)
                  tsum = B2g2*f4 + B2g21*f5
                  f4r2 = f4/rab(a,b)**2
c     also set up for partials
                  bp1  = 2._10*(f4 + f5)*A44
                  gp1  = bp1 -(vab2(a,a)-2._10*vab2(a,b)+vab2(b,b))
     .             /Cvel2
                  do i=1,3
                     df6dx(i)=0._10
                  end do
                  f1   = 0._10
                  f2   = 0._10
                  f3   = 0._10
                  do c = 0,11
                     if(Kgr(c).ge.0 .and. c.ne.a .and. c.ne.b) then
                        df1= Xmass(c)/rab(b,c)
                        df2= Xmass(c)/rab(a,c)
                        f7 = DOT(rvec(1,a,b),rvec(1,b,c))/rab(b,c)**2
                        f1 = f1 + df1
                        f2 = f2 + df2
                        f3 = f3 + df1*f7
                        df1d = df1/rab(b,c)**2
                        df2d = df2/rab(a,c)**2
                        df1dc= (B2m1-1.5_10*f7)*df1d
                        df2dc= B2g2*df2d
                        df3dc= 0.5_10*df1d
                        if(c.ne.0) then
                           dfdp=((B2m1-f7/2._10)/rab(b,c) +
     .                      B2g2/rab(a,c))*A44*f4r2
                           do i=1,3
                              sumrpt(i,c,lbd)= sumrpt(i,c,lbd) +
     .                         dfdp*rvec(i,a,b)
                           end do
                        endif
                        doterm=.false.
                        tgc=0._10
                        if(a.eq.11 .or. a.eq.10) then
                           tg =-df3dc
                           tgx=-df2dc
                           doterm=.true.
                        else
                           tg =0._10
                           tgx=0._10
                        endif
                        if(b.eq.11 .or. (Nplnt.eq.3 .and. b.eq.10)) then
                           tg =tg-df1dc+df3dc
                           tgc=tgc-df3dc
                           doterm=.true.
                        endif
                        if(c.eq.11 .or. (Nplnt.eq.3 .and. c.eq.10)) then
                           tg =tg+df1dc
                           tgx=tgx+df2dc
                           tgc=tgc+df3dc
                           doterm=.true.
                        endif
                        if(doterm) then
                           do i=1,3
                              df6dx(i)=df6dx(i)+tg*rvec(i,b,c)+
     .                         tgx*rvec(i,a,c)+tgc*rvec(i,a,b)
                           end do
                        endif
                     endif
                  end do

                  vbrab= DOT(v(1,b),rvec(1,a,b))/rab(a,b)
                  f6   = B2m1*f1 + B2g2*f2 - f3/2._10
                  f7   = (G22*vab2(a,b)-G11*vab2(b,b)-Gamapm*vab2(a,a)
     .             + 1.5_10*vbrab**2) / Cvel2
                  tsumx= (tsum + f6)*A44 + f7
                  tsump= f4r2*(B2g2*A44/rab(a,b)) + tsumx/rab(a,b)**3
                  tsuma= f4r2*(B2g21*A44/rab(a,b))
                  fx   = -2._10*f4r2/Mascnt/Cvel2
                  do c=1,11
                     if(Kgr(c).ge.0) then
                        if(c.eq.3 .or. c.eq.10 .or.
     .                   (c.eq.11.and.Nplnt.eq.3)) then
                           do i=1,3
                              dvb(i)=Xpert(i+3,3)+v(i,0)
                           end do
                           vavc=DOT(v(1,a),dvb)
                        else
                           do i=1,3
                              dvb(i)= v(i,c)
                           end do
                           vavc=vab2(a,c)
                        endif
                        vcrab=DOT(dvb,rvec(1,a,b))/rab(a,b)
                        dfdp = fx*(vavc+1.5_10*vbrab*vcrab)
                        do i=1,3
                           sumrpt(i,c,lbd)= sumrpt(i,c,lbd)
     .                      +dfdp*rvec(i,a,b)
                        end do
                     endif
                  end do

                  dfx = tsum*A44 + 3._10*(tsumx+vbrab**2/Cvel2)
                  dfv = 3._10*vbrab/rab(a,b)
                  tgx = 0._10
                  tgxv= 0._10
                  tdx = 0._10
                  tgva= -2._10*Xmsb
                  tgvb= 0._10
                  tgvx= -dfv*Xmsb
                  if(a.eq.11 .or. a.eq.10) then
                     tgx = -dfx
                     tgxv= dfv
                     tdx = tsumx*f4r2
                     tgva= tgva-2._10*Gamapm
                     tgvb= G22
                  endif
                  if(b.eq.11 .or. (Nplnt.eq.3.and.b.eq.10)) then
                     tgx = tgx+dfx
                     tgxv= tgxv-dfv
                     tdx = tdx-tsumx*f4r2
                     tgva= tgva+G22
                     tgvb= tgvb-G22
                     tgvx= tgvx+dfv
                  endif
                  tgxv= tgxv*f4r2/Cvel2
                  tgva= tgva*f4r2/Cvel2
                  tgvb= tgvb*f4r2/Cvel2
                  tgvx= tgvx*f4r2/Cvel2
                  tgx=tgx*f4r2/rab(a,b)**2
                  tgc=A44*f4r2
                  do i=1,3
                     dadxrt(i,i,lbd)=dadxrt(i,i,lbd)+tdx
                     tg  = tgx*rvec(i,a,b)+tgxv*v(i,b)+tgc*df6dx(i)
                     tgv = tgva*v(i,a)+tgvb*v(i,b)+tgvx*rvec(i,a,b)
                     do j=1,3
                        dadxrt(i,j,lbd)=dadxrt(i,j,lbd)+tg*rvec(j,a,b)
                        dadvrt(i,j,lbd)=dadvrt(i,j,lbd)+tgv*rvec(j,a,b)
                     end do
                  end do

                  do i = 1,3
                     sumrt(i,lbd)=sumrt(i,lbd) + tsumx*f4r2*rvec(i,a,b)
                     if(b.ne.0) sumrpt(i,b,lbd) = sumrpt(i,b,lbd)
     .                + tsump*rvec(i,a,b)
                     if(a.ne.0) sumrpt(i,a,lbd) = sumrpt(i,a,lbd)
     .                + tsuma*rvec(i,a,b)
                     sumrbt(i,lbd)=sumrbt(i,lbd)+(bp1+2._10*(f1+f2)*A44)
     .                *f4r2*rvec(i,a,b)
                     sumrgt(i,lbd)=sumrgt(i,lbd)+(gp1+2._10*f2*A44)
     .                *f4r2*rvec(i,a,b)
                  end do

                  f6 = f4*A44
                  do c = 0,11
                     if(Kgr(c).ge.0 .and. c.ne.a .and. c.ne.b) then
                        rbc3 = rab(b,c)**3
                        f7 = f6*Xmass(c)/rbc3
                        do i = 1,3
                           sumrt(i,lbd)=sumrt(i,lbd) - f7*rvec(i,b,c)*G7
                           sumrgt(i,lbd)=sumrgt(i,lbd) -
     .                      f7*rvec(i,b,c)*2._10
                        end do
                        if(c.ne.0) then
                           dfdp=-f6*G7/rbc3
                           do i=1,3
                              sumrpt(i,c,lbd)= sumrpt(i,c,lbd)+
     .                         dfdp*rvec(i,b,c)
                           end do
                        endif
                        if(b.ne.0) then
                           dfdp=-A44/rab(a,b)*G7*Xmass(c)/rbc3
                           do i=1,3
                              sumrpt(i,b,lbd)= sumrpt(i,b,lbd)+
     .                         dfdp*rvec(i,b,c)
                           end do
                        endif
                        df7ab= f7*G7/rab(a,b)**2
                        df7bc= 3._10*f7*G7/rab(b,c)**2
                        tgc=0._10
                        tdx=0._10
                        if(a.eq.11 .or. a.eq.10) then
                           tgx= df7ab
                        else
                           tgx= 0._10
                        endif
                        if(b.eq.11 .or. (Nplnt.eq.3.and.b.eq.10)) then
                           tgx= tgx-df7ab
                           tgc= df7bc
                           tdx= -f7*G7
                        endif
                        if(c.eq.11 .or. (Nplnt.eq.3.and.c.eq.10)) then
                           tgc= tgc-df7bc
                           tdx= tdx+f7*G7
                        endif
                        do i=1,3
                           dadxrt(i,i,lbd)=dadxrt(i,i,lbd)+tdx
                           tg=tgx*rvec(i,a,b)+tgc*rvec(i,b,c)
                           do j=1,3
                              dadxrt(i,j,lbd)=dadxrt(i,j,lbd)+
     .                         tg*rvec(j,b,c)
                           end do
                        end do
                     endif
                  end do
               endif

            end do
 
            do b = 0,11
               if(b.ne.a .and. Kgr(b).ge.0) then
                  f1= DOT(rvec(1,a,b),v(1,a))
                  f2= DOT(rvec(1,a,b),v(1,b))
                  f6= rab(a,b)**3*Cvel2
                  fx= Xmass(b)/f6
                  f3= (G22*f1-G21*f2)*fx
                  f4= 2._10*(f1 - f2)*fx
                  do i = 1,3
                     vab(i)= v(i,a)-v(i,b)
                     sumrt(i,lbd)=sumrt(i,lbd) + f3*vab(i)
                     sumrgt(i,lbd)=sumrgt(i,lbd) + f4*vab(i)
                  end do
                  if(b.ne.0) then
                     dfdp= (G22*f1-G21*f2)/f6
                     do i=1,3
                        sumrpt(i,b,lbd)= sumrpt(i,b,lbd)+dfdp*vab(i)
                     end do
                  endif
                  do c=1,11
                     if(Kgr(c).ge.0) then
                        if(c.eq.3 .or. c.eq.10 .or.
     .                   (c.eq.11.and.Nplnt.eq.3)) then
                           do i=1,3
                              dva(i)=Xpert(i+3,3)+v(i,0)
                           end do
                        else
                           do i=1,3
                              dva(i)= v(i,c)
                           end do
                        endif
                        dfdp=-fx*DOT(rvec(1,a,b),dva)/Mascnt
                        do i=1,3
                           sumrpt(i,c,lbd)= sumrpt(i,c,lbd)+dfdp*vab(i)
                        end do
                     endif
                  end do
                  tga = 0._10
                  tgb = 0._10
                  tgx = 0._10
                  tgvx= -Xmsb
                  tdv = 0._10
                  fxa = 3._10*f3/rab(a,b)**2
                  if(a.eq.11 .or. a.eq.10) then
                     tga = G22
                     tgb = -G21
                     tgx = -fxa
                     tgvx= tgvx+G22
                     tdv = f3
                  endif
                  if(b.eq.11 .or. (Nplnt.eq.3.and.b.eq.10)) then
                     tga = tga-G22
                     tgb = tgb+G21
                     tgx = tgx+fxa
                     tgvx= tgvx-G21
                     tdv = tdv-f3
                  endif
                  tga = tga*fx
                  tgb = tgb*fx
                  tgvx= tgvx*fx
                  do i=1,3
                     tg =tga*v(i,a)+tgb*v(i,b)+tgx*rvec(i,a,b)
                     tgv=tgvx*rvec(i,a,b)
                     dadvrt(i,i,lbd) = dadvrt(i,i,lbd)+tdv
                     do j=1,3
                        dadxrt(i,j,lbd)= dadxrt(i,j,lbd)+tg*vab(j)
                        dadvrt(i,j,lbd)= dadvrt(i,j,lbd)+tgv*vab(j)
                     end do
                  end do
               endif
            end do

c copy partial w.r.t. integrated planet mass into the range 1-9
            if(Nplnt.le.9) then
               do i=1,3
                  sumrpt(i,Nplnt,lbd)=sumrpt(i,11,lbd)
               end do
            endif

c combine earth and moon into embary if necessary, both for the mass
c used as a parameter and for the accelerations whose partial
c derivatives are being taken

            if(Kgr(10).ge.0) then
               do i=1,3
                  sumrpt(i,3,lbd) = sumrpt(i,3,lbd) +
     .             Mass(10)*(sumrpt(i,10,lbd)-sumrpt(i,3,lbd))
               end do
            endif
            if(lbd.eq.3) then
               do i=1,3
                  sumrt(i,2) = sumrt(i,2) +
     .             Mass(10)*(sumrt(i,3)-sumrt(i,2))
                  sumrbt(i,2)= sumrbt(i,2) +
     .             Mass(10)*(sumrbt(i,3)-sumrbt(i,2))
                  sumrgt(i,2)= sumrgt(i,2) +
     .             Mass(10)*(sumrgt(i,3)-sumrgt(i,2))
                  do b=1,9
                     if(Kgr(b).ge.0 .or. b.eq.3) sumrpt(i,b,2)=
     .                sumrpt(i,b,2) + Mass(10)*
     .                (sumrpt(i,b,3)-sumrpt(i,b,2))
                  end do
                  do j=1,3
                     dadxrt(i,j,2)= dadxrt(i,j,2) +
     .                Mass(10)*(dadxrt(i,j,3)-dadxrt(i,j,2))
                     dadvrt(i,j,2)= dadvrt(i,j,2) +
     .                Mass(10)*(dadvrt(i,j,3)-dadvrt(i,j,2))
                  end do
               end do
            endif
         end do
c
c subtract off acceleration of sun
c
         do i=1,3
            Sumr(i) = sumrt(i,2)-sumrt(i,1)
            Sumrb(i)= sumrbt(i,2)-sumrbt(i,1)
            Sumrg(i)= sumrgt(i,2)-sumrgt(i,1)
            do b=1,9
               if(Kgr(b).ge.0 .or. b.eq.Nplnt)
     .          Sumrp(i,b)=sumrpt(i,b,2)-sumrpt(i,b,1)
            end do
            do j=1,3
               Dadxr(i,j)=dadxrt(i,j,2)-dadxrt(i,j,1)
               Dadvr(i,j)=dadvrt(i,j,2)-dadvrt(i,j,1)
               Dadxsav(i,j,1,1,2)=Dadxr(i,j)
               Dadxsav(i,j,2,1,2)=Dadvr(i,j)
            end do
         end do
         return
      else if(ncall.eq.0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c* start=100
c
c     setup once per step
c         note: xpert(i,nplnt) is filled only once per step
c
c     compute velocity of sun relative to cm of solar system
c
         do i = 1,3
            sum = 0._10
            do j = 1,9
               if(Kgr(j).ge.0 .or. j.eq.Nplnt) sum = sum
     .          - Xpert(i+3,j)*Mass(j)
            end do
            v(i,0) = sum/Mascnt
         end do

         do k=1,9
            do i=1,6
               xpcop(i,k)=Xpert(i,k)
            end do
         end do
c get separate earth and moon coordinates if needed
         if(Kgr(10).ge.0) then
            do i=1,6
               xpcop(i,3)=xpcop(i,3)-Mass(10)*Xpert(i,10)
               xpcop(i,10)=xpcop(i,3)+Xpert(i,10)
            end do
            if(Nplnt.eq.3) then
               do i=1,3
                  rvec(i,10,11) = Xpert(i,10)
                  rvec(i,11,10) = -Xpert(i,10)
               end do
               rab(10,11) = SQRT(DOT(rvec(1,10,11),rvec(1,10,11)))
               rab(11,10) = rab(10,11)
            endif
         endif

         do k = 1,10
            if(Kgr(k).ge.0) then
               do j = k,10
                  if(Kgr(j).ge.0) then
                     do i = 1,3
                        rvec(i,j,k) = xpcop(i,j) - xpcop(i,k)
                        rvec(i,k,j) = -rvec(i,j,k)
                     end do
                  endif
               end do
            endif
         end do
 
         do j = 1,10
            if(Kgr(j).ge.0) then
               do i = 1,3
                  v(i,j) = xpcop(i+3,j) + v(i,0)
                  rvec(i,j,0) = xpcop(i,j)
                  rvec(i,0,j) = -xpcop(i,j)
               end do
            endif
         end do
 
         do i = 0,10
            if(Kgr(i).ge.0) then
               do j = i,10
                  if(Kgr(j).ge.0) then
                     vab2(i,j) = DOT(v(1,i),v(1,j))
                     if(i.ne.j) then
                        vab2(j,i) = vab2(i,j)
                        rab(i,j)  = SQRT(DOT(rvec(1,i,j),rvec(1,i,j)))
                        rab(j,i)  = rab(i,j)
                     endif
                  endif
               end do
            endif
         end do
         return
      else
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c* start=1000
c
c     gr acceleration for motion and partials is calculated
c     in setup call and added to fn in solprb
c
         return
      endif
      end
