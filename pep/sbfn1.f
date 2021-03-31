      subroutine SBFN1(k)
 
      implicit none
c
c     ash/amuchastegui - june 1969 - subroutine sbfn1
c     revised 1977 april - j.f.chandler
c        evaluate right side of motion and partials
c           principle of equivalence violation is programmed only for
c          sun-centered bodies.  extra terms would appear otherwise.

c arguments
      integer*4 k

c array dimensions
      include 'globdefs.inc'
c commons
      include 'cnthar.inc'
      include 'drgprt.inc'
      include 'ellips.inc'
      include 'fcntrl.inc'
      include 'fmmips.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'mascn1.inc'
      include 'namtimq.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'prtcod.inc'
      include 'rtpars.inc'
      include 'sbembr.inc'
      include 'sbrot.inc'
      include 'sbstuf.inc'
      real*10 dadxx(3,3,6)
      equivalence (dadxx,Dadxc)
      logical*1 lggfg(8),lggfgs(16)
      equivalence (lggfg,Ggfg),(Ggfgs,lggfgs)
      include 'sbthng.inc'
      include 'tapdtplp.inc'
      include 'trghar.inc'
      include 'xprcom.inc'
      include 'yvectplp.inc'

c external functions
      real*10 DOT

c local
      real*10 di2,di3,di4,djh,fnfb,psum,psump,psumt1,psumt2,
     . psunp,sum,sumast,sumda,sumrot,td,termd,termg,
     . tg,tgaa,tgar,tgax,tgrr,tgrrp,tgrx,tgxx,tm1,tm2,
     . tq,tqp,zd,zgaa,zgax,zgxx
      integer   i,i1,icn,ib,ih,ikkk,iq,itg,ix,jb,jh,jq,kkk1,
     .          kq,kt,l,ll,nt,ntt,ntz
      character*2 mtxid(19)/'  ','C ','H ','S ','P ','A ','T1',
     .          'T2','T3','T4','T5','T6','T7','T8','T9','TT',
     .          'TE','L ','M '/
      real*10 dsmast(3,12),sums(3),sump(3),dsump(3,10),
     .       tmprot(3),talng1(3,i_mxtrg),calng1(3),dydm(3)
      real*10 tesst1(3,i_mxtrg,5),tesst2(3,i_mxtrg,5),sumtc,sumpar,
     . tzone(3,i_mxtrg,4),sumht(i_mxtrg,3),tsnh1(i_mxtrg,5),
     . tsnh2(i_mxtrg,5),sumhc(3),tessc1(3,(u_mxtes*(u_mxtes+1))/2-1),
     . tessc2(3,(u_mxtes*(u_mxtes+1))/2-1),czone(3,u_mxzon-1),
     . csnh1((u_mxtes*(u_mxtes+1))/2-1),
     . csnh2((u_mxtes*(u_mxtes+1))/2-1)
      real*10 mscfrc
c
c
      if(Kkk.gt.0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=2000
c
c           determine right side of equations for partial derivatives
c           with respect to parameter alpha(kkk)
c           (kkk goes from 1 to number of partials)
c
c           effect of perturbing planets on partial derivatives
         psump = 0._10
         ikkk  = Icntrl(Kkk)
         if(dflgp) then

c calculate gravity gradient w.r.t. coordinates of integrated body
            do iq = 1,9
               Dadxp(iq,1) = 0._10
            end do
            do l = 1,10
               if(Kpb(l) .and. Kp(l+30).gt.0) then
                  if(Nplnt.eq.3 .and. Kp(40).ge.0) then
 
c for earth-moon barycenter
                     tg    = Mass1(l)*Gamat
                     termd =-tg*(Masse/Rpes3(l)+Mass(10)/Rpms3(l))
                     tg    = 3._10*tg
                     do iq = 1,3
                        do jq = 1,3
                           termg = tg*(
     .                      Xpes(iq,l)*Xpes(jq,l)/Rpes5(l)*Masse +
     .                      Xpms(iq,l)*Xpms(jq,l)/Rpms5(l)*Mass(10))
                           Dadxp(iq,jq) = Dadxp(iq,jq) + termg
                           if(iq.eq.jq) goto 10
                           Dadxp(jq,iq) = Dadxp(jq,iq) + termg
                        end do
   10                   Dadxp(iq,iq) = Dadxp(iq,iq) + termd
                     end do
                  else
 
c for everything but earth-moon barycenter
                     termd = -Mass1(l)*Gamat/Rpb3(l)
                     tg    = -termd*3._10/Rpb2(l)
                     do iq = 1,3
                        do jq = 1,3
                           termg = tg*Pbcor(iq,l)*Pbcor(jq,l)
                           Dadxp(iq,jq) = Dadxp(iq,jq) + termg
                           if(iq.eq.jq) goto 20
                           Dadxp(jq,iq) = Dadxp(jq,iq) + termg
                        end do
   20                   Dadxp(iq,iq) = Dadxp(iq,iq) + termd
                     end do
                  endif
               endif
   40       end do
         endif
c*  start=2100
c
c check if alpha is mass of perturbing planet l
         if(ikkk.ge.1 .and. ikkk.le.10) then
            l = ikkk
            if(Kpb(l)) then
 
c psunp = pbcor3(kk,l)-pccor3(kk,l)
               psunp = dsump(Kk,l)
               if(l.ne.3) then
                  if(l.eq.10) then
                     if(Ncentr.eq.3) then
                        psunp = Mass(3)*(psunp - Sbcor(Kk)/Rsb3)
                     else
                        psunp = Mass(3)*(psunp - dsump(Kk,3))
                     endif
                  endif
               else if(Ncentr.eq.10) then
                  psunp = Masse*psunp - Mass(10)*Sbcor(Kk)/Rsb3
               else if(Kp(40).gt.0) then
                  psunp = Masse*psunp + Mass(10)*dsump(Kk,10)
               endif
               psump = psump + psunp
            else if(Nplnt.eq.3 .and. l.eq.10 .and. Kp(40).gt.0) then
c special case for partial derivative of EMbary orbit w.r.t. Moon mass
               psunp=0._10
               do i=1,9
                  if(Kp(i+30).gt.0 .and. i.ne.3) psunp=psunp+Mass1(i)*
     .             (Xpms(Kk,i)/Rpms3(i)-Xpes(Kk,i)/Rpes3(i))
               end do
               psunp=psunp+Massp1*(Xes(Kk)/Res3-Xms(Kk)/Rms3)
               psump=psump+psunp
            endif
         endif
c*  start=2200
c
c effect of perturbing asteroids or satellites on partial derivatives
         if(Numast.gt.0) then
            if(dflga) then
c calculate gravity gradient w.r.t. coordinates of integrated body
               do iq = 1,9
                  Dadxa(iq,1) = 0._10
               end do
               do l = 1,Numast
                  if(Kpast(l).gt.0) then
                     termd = -Astmas(l)*Gamat/Rastb3(l)
                     tg    = -termd*3._10/Rastb2(l)
                     do iq = 1,3
                        do jq = 1,3
                           termg = tg*Yastb(iq,l)*Yastb(jq,l)
                           Dadxa(iq,jq) = Dadxa(iq,jq) + termg
                           if(iq.eq.jq) goto 50
                           Dadxa(jq,iq) = Dadxa(jq,iq) + termg
                        end do
   50                   Dadxa(iq,iq) = Dadxa(iq,iq) + termd
                     end do
                  endif
               end do
            endif
c check if alpha is mass of asteroid/satellite or central body thereof
            do l = 1,Numast
               if(ikkk.eq.Nplast(l)) then
                  psump = psump + dsmast(Kk,l)*Astmac(l)
               else if(ikkk.eq.Ncnast(l)) then
                  psump = psump + dsmast(Kk,l)*Astmab(l)
               endif
            end do
         endif

c*  start=2300
         Fn(1) = Fn(1) + Gamat*psump
c
c check if alpha is time variable gravitational constant
         if(ikkk.eq.32) Fn(1) = Fn(1) + Gama*Tvary*sump(Kk)

c*  start=2500

c add indirect terms in partial w.r.t. parameter affecting target bodies
         do kt=1,Numtar
            iq=Kpt(Kkk,kt)
            if(iq.gt.0) then
               tm1=DOT(Dtcor(1,iq-1,kt),Dadp(1,Kk,kt))
               Fn(1)=Fn(1)+tm1
               if(Kkk.eq.1) then
                  nt=Nplpt(kt)
                  if(kt.eq.1) then
                     Indpsav(Kk,2,1)=tm1
                  else if(Nplnt.eq.3 .and. Kp(40).ge.0 .and. nt.eq.10)
     .                then
                     Indpsav(Kk,3,1)=tm1
                  endif
               endif
            endif
         end do
c
c effect of target  body harmonics on partial derivatives
         do ll = 1,Numtar
            if(Ntopt(ll).gt.0) then
               sumpar = 0._10
c
c effect of target body  zonal harmonics on partials
               if(Ntzone(ll).gt.1) then
                  ntz = Ntzon1(ll)
                  if(Rtb(ll).gt.Hsitb(ll)) ntz = 1
                  if(dflgt .or. ikkk.eq.Jtzone(ll))then
                     if(dflgt) then
                        tg   = Masst(ll)*Gamat/Rtb3(ll)
                        zgaa = 0._10
                        zgax = 0._10
                        zgxx = 0._10
                        zd   = 0._10
                        di3  = 3._10
                        di4  = 4._10
                     endif
                     do i = 1,ntz
                        if(dflgt .and. i.le.Ntzonp(ll)) then
                           di2  = di3
                           di3  = di4
                           di4  = i + 4
                           tm2  = Tzhar(ll,i)/Rpbh(ll,i)
                           zgaa = zgaa + Tleg2(i,ll)*tm2
                           tm1  = tm2*(di2*Tleg(i,ll) + Tslat(ll)
     .                      *Tleg1(i,ll))
                           tm2  = tm2*(di3*Tleg1(i,ll) + Tslat(ll)
     .                      *Tleg2(i,ll))
                           zgax = zgax + tm2
                           zgxx = zgxx + di4*tm1 + Tslat(ll)*tm2
                           zd   = zd + tm1
                        endif
 
c check if alpha is zonal harmonic coefficient
                        if(ikkk.eq.Jtzone(ll)) then
                           if(Itzone(ll,i).eq.Kkk)
     .                      sumpar = sumpar + tzone(Kk,ll,i)
                        endif
                     end do
                     if(dflgt) then
                        zgaa = zgaa*tg
                        zgax = zgax*tg/Rtb(ll)
                        zgxx = zgxx*tg/Rtb2(ll)
                        zd   = -zd*tg
                        do iq = 1,3
                           do jq = 1,3
                              termg = zgxx*Tbcor(iq,ll)*Tbcor(jq,ll)
     .                         + zgax*(Tbcor(iq,ll)*Trgrot(3,jq,ll)
     .                         + Tbcor(jq,ll)*Trgrot(3,iq,ll))
     .                         + zgaa*Trgrot(3,iq,ll)
     .                         *Trgrot(3,jq,ll)
                              Dadxt(iq,jq,ll) = termg
                              if(iq.eq.jq) goto 110
                              Dadxt(jq,iq,ll) = termg
                           end do
  110                      Dadxt(iq,iq,ll) = termg + zd
                        end do
                     endif
                  endif
               endif
c*  start=2600
c
c effect of target body tesseral harmonics on partials
               if(Nttess(ll).gt.1 .and. Rtb(ll).le.Hsitb(ll)) then
                  if(dflgt .or. ikkk.eq.Jtcos(ll)
     .             .or. ikkk.eq.Jtsin(ll)) then
                     if(dflgt) then
                        if(Ntzone(ll).le.1) then
                           do iq = 1,9
                              Dadxt(iq,1,ll) = 0._10
                           end do
                        endif
                        tgxx  = 0._10
                        tgax  = 0._10
                        tgrx  = 0._10
                        tgrr  = 0._10
                        tgar  = 0._10
                        tgaa  = 0._10
                        tgrrp = 0._10
                        td    = 0._10
                        do iq = 1,3
                           talng1(iq,ll)=-Trgrot(2,iq,ll)*Tslng(ll,1)
     .                      - Trgrot(1,iq,ll)*Tclng(ll,1)
                        end do
                        di3 = 3._10
                        di4 = 4._10
                     endif
                     ih  = 0
                     ntt = Nttes1(ll)
                     do i = 1,ntt
                        di2 = di3
                        di3 = di4
                        di4 = i + 4
                        i1  = i + 1
                        do jh = 1,i1
                           djh = jh
                           ih  = ih + 1
                           if(i.le.Nttesp(ll) .and. dflgt) then
                              tq    = tsnh1(ll,ih)/Rpbh(ll,i)
                              tqp   = djh*tsnh2(ll,ih)/Rpbh(ll,i)
                              tm1   = di2*Tgleg(ih,ll)
     .                         + Tslat(ll)*Tgleg1(ih,ll)
                              tm2   = di3*Tgleg1(ih,ll)
     .                         + Tslat(ll)*Tgleg2(ih,ll)
                              tgar  = tgar + Tgleg1(ih,ll)*tqp
                              tgrrp = tgrrp + Tgleg(ih,ll)*tqp
                              tgaa  = tgaa + Tgleg2(ih,ll)*tq
                              tgrx  = tgrx + tqp*tm1
                              tgax  = tgax + tq*tm2
                              tgxx  = tgxx +
     .                         tq*(di4*tm1 + Tslat(ll)*tm2)
                              tgrr  = tgrr +
     .                         djh*djh*tq*Tgleg(ih,ll)
                              td    = td + tq*tm1
                           endif
 
c check if alpha is tess. harmonic cosine coefficient
                           if(ikkk.eq.Jtcos(ll) .and.
     .                      Itcos(ll,ih).eq.Kkk) then
                              sumpar = sumpar +
     .                         Tclng(ll,jh)*tesst1(Kk,ll,ih)
     .                         - Tslng(ll,jh)*tesst2(Kk,ll,ih)*djh
 
c check if alpha is tess. harmonic sine coefficient
                           else if(ikkk.eq.Jtsin(ll) .and.
     .                         Itsin(ll,ih).eq.Kkk) then
                              sumpar = sumpar +
     .                         Tslng(ll,jh)*tesst1(Kk,ll,ih)
     .                         + Tclng(ll,jh)*tesst2(Kk,ll,ih)*djh
                           endif
                        end do
                     end do
                     if(dflgt) then
                        tgar  = tgar*tg
                        tgrrp = tgrrp*tg/Tclat(ll)
                        tgaa  = tgaa*tg
                        tgrx  = tgrx*tg/Rtb(ll)
                        tgax  = tgax*tg/Rtb(ll)
                        tgxx  = tgxx*tg/Rtb2(ll)
                        tgrr  = -tgrr*tg
                        td    = -td*tg
                        do iq = 1,3
                           do jq = 1,3
                              termg = tgxx*Tbcor(iq,ll)*Tbcor(jq,ll) +
     .                         tgax*(Tbcor(iq,ll)*Trgrot(3,jq,ll) +
     .                         Tbcor(jq,ll)*Trgrot(3,iq,ll))
     .                         + tgrx*(Tbcor(iq,ll)*Tlng1(jq,ll)
     .                         + Tbcor(jq,ll)*Tlng1(iq,ll))
     .                         + tgar*(Trgrot(3,iq,ll)*Tlng1(jq,ll)
     .                         + Trgrot(3,jq,ll)*Tlng1(iq,ll))
     .                         + tgaa*Trgrot(3,iq,ll)*Trgrot(3,jq,ll)
     .                         + tgrrp*(Tlng1(jq,ll)*talng1(iq,ll)
     .                         + Tlng1(iq,ll)*talng1(jq,ll))
     .                         + tgrr*Tlng1(iq,ll)*Tlng1(jq,ll)
                              Dadxt(iq,jq,ll)= Dadxt(iq,jq,ll) + termg
                              if(iq.eq.jq) goto 130
                              Dadxt(jq,iq,ll)= Dadxt(jq,iq,ll) + termg
                           end do
  130                      Dadxt(iq,iq,ll) = Dadxt(iq,iq,ll) + td
                        end do
                     endif
                  endif
               endif
c*  start=2800
c
               Fn(1) = Fn(1) + Gamat*Masst(ll)*sumpar

               kt = Itgast(ll)
c check if alpha is mass of target planet
               if(ikkk.eq.Ntrg(ll)) then
                  sumpar = Gamat*Masstc(ll)
c check if alpha is mass of target planet's center
               else if(kt.gt.0 .and. ikkk.eq.Ncnast(kt)) then
                  sumpar = Gamat*Astmab(kt)
c check if alpha is time variable grav. constant
               else if(ikkk.eq.32) then
                  sumpar = Gama*Tvary*Masst(ll)
               else
                  goto 160
               endif
               Fn(1) = Fn(1) + sumpar*sumht(ll,Kk)
            endif
  160    end do
c*  start=2900
c
         if(Ncentr.le.0) goto 400
c
c effect of sun on partial derivatives
         if(Kp(Ncentr+30).gt.0) then
            if(dflgs) then
 
               termd = -Gamat/Rb3
               tg    = Gamat*3._10/Rb5
               do iq = 1,3
                  do jq = 1,3
                     termg = tg*Bcor(iq)*Bcor(jq)
                     Dadxs(iq,jq) = termg
                     if(iq.eq.jq) goto 170
                     Dadxs(jq,iq) = termg
                  end do
  170             Dadxs(iq,iq) = Dadxs(iq,iq) + termd
               end do
            endif
c
c check if alpha is time variable grav. constant
            if(ikkk.eq.32) Fn(1) = Fn(1) + Gama*Tvary*sums(Kk)
         endif
c*  start=3000
c
c effect of central body zonal harmonics on part. deriv.
         sumpar = 0._10
         tg     = Gamat3/Rsb3
         if(Nczone.gt.1) then
            if(dflgh .or. ikkk.eq.Jczone) then
               if(dflgh) then
                  zgaa = 0._10
                  zgax = 0._10
                  zgxx = 0._10
                  zd   = 0._10
                  di3  = 3._10
                  di4  = 4._10
               endif
               do i = 1,Nczon1
                  if(i.le.Nczonp .and. dflgh) then
                     di2  = di3
                     di3  = di4
                     di4  = i + 4
                     tm2  = Czhar(i)/Rsbh(i)
                     zgaa = zgaa + Cleg2(i)*tm2
                     tm1  = tm2*(di2*Cleg(i) + Cslat*Cleg1(i))
                     tm2  = tm2*(di3*Cleg1(i) + Cslat*Cleg2(i))
                     zgax = zgax + tm2
                     zgxx = zgxx + di4*tm1 + Cslat*tm2
                     zd   = zd + tm1
                  endif
 
c check if alpha is zonal harmonic coefficient
                  if(ikkk.eq.Jczone) then
                     if(Iczone(i).eq.Kkk) sumpar = czone(Kk,i)
                  endif
               end do
 
               if(dflgh) then
                  zgaa = zgaa*tg
                  zgax = zgax*tg/Rsb
                  zgxx = zgxx*tg/Rsb2
                  zd   = -zd*tg
 
                  do iq = 1,3
                     do jq = 1,3
                        termg = zgxx*Sbcor(iq)*Sbcor(jq)
     .                          - zgax*(Sbcor(iq)*Cntrot(3,jq)
     .                          + Sbcor(jq)*Cntrot(3,iq))
     .                          + zgaa*Cntrot(3,iq)*Cntrot(3,jq)
                        Dadxh(iq,jq) = termg
                        if(iq.eq.jq) goto 180
                        Dadxh(jq,iq) = termg
                     end do
  180                Dadxh(iq,iq) = Dadxh(iq,iq) + zd
                  end do
               endif
            endif
         endif
c*  start=3100
c effect of central body tesseral harmonics on partials
         if(Nctess.gt.1) then
            if( .not.(.not.dflgh .and. ikkk.ne.Jcsin .and.
     .          ikkk.ne.Jccos) ) then
               if(dflgh) then
                  if(Nczone.le.1) then
                     do iq = 1,9
                        Dadxh(iq,1) = 0._10
                     end do
                  endif
                  do iq = 1,3
                     calng1(iq) = -Cntrot(2,iq)*Cslng(1)
     .                            - Cntrot(1,iq)*Cclng(1)
                  end do
                  tgar  = 0._10
                  tgrrp = 0._10
                  tgaa  = 0._10
                  tgrx  = 0._10
                  tgax  = 0._10
                  tgxx  = 0._10
                  tgrr  = 0._10
                  td    = 0._10
                  di3   = 3._10
                  di4   = 4._10
               endif
               ih = 0
               do i = 1,Nctes1
                  di2 = di3
                  di3 = di4
                  di4 = i + 4
                  i1  = i + 1
                  do jh = 1,i1
                     djh = jh
                     ih  = ih + 1
                     if(dflgh) then
                        if(i.le.Nctesp) then
                           tq    = csnh1(ih)/Rsbh(i)
                           tqp   = djh*csnh2(ih)/Rsbh(i)
                           tm1   = di2*Cgleg(ih) + Cslat*Cgleg1(ih)
                           tm2   = di3*Cgleg1(ih) + Cslat*Cgleg2(ih)
                           tgar  = tgar + Cgleg1(ih)*tqp
                           tgrrp = tgrrp + Cgleg(ih)*tqp
                           tgaa  = tgaa + tq*Cgleg2(ih)
                           tgrx  = tgrx + tqp*tm1
                           tgax  = tgax + tq*tm2
                           tgxx  = tgxx + tq*(di4*tm1 + Cslat*tm2)
                           tgrr  = tgrr + djh*djh*tq*Cgleg(ih)
                           td    = td + tq*tm1
                        endif
                     endif
 
c check if alpha is tess. harmonic cosine coefficient
                     if(ikkk.eq.Jccos) then
                        if(Iccos(ih).eq.Kkk) then
                           sumpar = Cclng(jh)*tessc1(Kk,ih) - Cslng(jh)
     .                              *tessc2(Kk,ih)*djh
                           goto 190
                        endif
                     endif
 
c check if alpha is tess. harmonic sine coefficient
                     if(ikkk.eq.Jcsin) then
                        if(Icsin(ih).eq.Kkk) sumpar = Cslng(jh)
     .                      *tessc1(Kk,ih) + Cclng(jh)*tessc2(Kk,ih)
     .                      *djh
                     endif
  190             end do
               end do
               if(dflgh) then
                  tgar  = tgar*tg
                  tgrrp = tgrrp*tg/Cclat
                  tgaa  = tgaa*tg
                  tgrx  = -tgrx*tg/Rsb
                  tgax  = -tgax*tg/Rsb
                  tgxx  = tgxx*tg/Rsb2
                  tgrr  = -tgrr*tg
                  td    = -td*tg
                  do iq = 1,3
                     do kq = 1,3
                        termg = tgxx*Sbcor(iq)*Sbcor(kq)
     .                          + tgax*(Sbcor(iq)*Cntrot(3,kq)
     .                          + Sbcor(kq)*Cntrot(3,iq))
     .                          + tgrx*(Sbcor(iq)*Clng1(kq) + Sbcor(kq)
     .                          *Clng1(iq))
     .                          + tgar*(Cntrot(3,iq)*Clng1(kq)
     .                          + Cntrot(3,kq)*Clng1(iq))
     .                          + tgaa*Cntrot(3,iq)*Cntrot(3,kq)
     .                          + tgrrp*(Clng1(kq)*calng1(iq)
     .                          + Clng1(iq)*calng1(kq)) + tgrr*Clng1(iq)
     .                          *Clng1(kq)
                        Dadxh(iq,kq) = Dadxh(iq,kq) - termg
                        if(iq.eq.kq) goto 200
                        Dadxh(kq,iq) = Dadxh(kq,iq) - termg
                     end do
  200                Dadxh(iq,iq) = Dadxh(iq,iq) - td
                  end do
               endif
            endif
         endif
 
         Fn(1) = Fn(1) - Gamat3*sumpar
c*  start=3300
c
c check for partials derivatives w.r.t. central body rotation
c matrix parameters  con 6,7,8,9
         ix = 0
         do while(.true.)
            ix = ix + 1
            if(ix.gt.Icrot) goto 400
            if(Kkk.eq.Icrotp(ix,1)) then
               sumrot = Gamat3*(Pcnrot(Kk,1,ix)*sumhc(1)
     .                  + Pcnrot(Kk,2,ix)*sumhc(2) + Pcnrot(Kk,3,ix)
     .                  *sumhc(3))
               if(Kk.eq.1) call PRODCT(Pcnrot(1,1,ix),Sbcor,
     .             tmprot,3,3,1)
               sumrot = sumrot + DOT(Dadxh(1,Kk),tmprot(1))
               Fn(1)  = Fn(1) + sumrot
               goto 400
            endif
         end do
      endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=1000
c
c           determine right side of equations of motion
c
c           effect of perturbing planets on motion of satellite or probe
      sump(Kk) = 0._10
      do l=1,10
         if(Kpb(l)) then
            dsump(Kk,l) = Pbcor3(Kk,l) - Pccor3(Kk,l)
            sump(Kk) = sump(Kk) + Mass1(l)*dsump(Kk,l)
         endif
      end do
c
c*  start=1100
c
c effect of perturbing asteroids or satellites on motion
      sumast = 0._10
      do l = 1,Numast
         dsmast(Kk,l) = Yastb(Kk,l)/Rastb3(l) - Yastc3(Kk,l)
         sumast = sumast + Astmas(l)*dsmast(Kk,l)
      end do
      sump(Kk) = sump(Kk) + sumast
 
      Fn(1) = Fn(1) + Gamat*sump(Kk)
c
c*  start=1200
c effect of target  body  harmonics on equations of motion
      do ll = 1,Numtar
c
c effect of target body  zonal   harmonics on
c equations of motion
         sumht(ll,Kk) = 0._10
         if(Ntopt(ll).gt.0) then
            if(Ntzone(ll).gt.1) then
               ntz = Ntzon1(ll)
               if(Rtb(ll).gt.Hsitb(ll)) ntz = 1
               do i = 1,ntz
                  di2 = i + 2
                  tzone(Kk,ll,i) = (-di2*Tbcor(Kk,ll)*Tleg(i,ll)/
     .                               Rtb(ll) - Tleg1(i,ll)
     .                               *Tslat1(Kk,ll))/Rpbh(ll,i)
     .                               /Rtb2(ll)
                  sumht(ll,Kk)    = sumht(ll,Kk) + Tzhar(ll,i)
     .                               *tzone(Kk,ll,i)
c substract effect of target body harmonic on central body
c in appropiate orbiter routine (we ignore effect on sun
c of target body harmonic for sun centered probe)
               end do
            endif
c
c effect of target body tesseral harmonic on eq. of motion
            if(Nttess(ll).gt.1) then
               ih  = 0
               ntt = Nttes1(ll)
               if(Rtb(ll).le.Hsitb(ll)) then
                  do i = 1,ntt
                     di2   = i + 2
                     sumtc = 0._10
                     i1    = i + 1
                     do jh = 1,i1
                        djh = jh
                        ih  = ih + 1
                        tesst1(Kk,ll,ih)
     .                     = (di2*Tbcor(Kk,ll)*Tgleg(ih,ll)/Rtb(ll)
     .                     + Tgleg1(ih,ll)*Tslat1(Kk,ll))/Rpbh(ll,i)
     .                     /Rtb2(ll)
                        tesst2(Kk,ll,ih) = Tgleg(ih,ll)
     .                     *Tlng1(Kk,ll)/Rpbh(ll,i)/Rtb2(ll)
                        if(Kk.le.1) then
 
c only once per iteration
                           tsnh1(ll,ih) = Tchar(ll,ih)
     .                        *Tclng(ll,jh) + Tshar(ll,ih)
     .                        *Tslng(ll,jh)
     .                        *Tslng(ll,jh) + Tshar(ll,ih)
     .                        *Tclng(ll,jh)
                        endif
                        sumtc = sumtc + tsnh1(ll,ih)
     .                          *tesst1(Kk,ll,ih)
     .                          + djh*tsnh2(ll,ih)
     .                          *tesst2(Kk,ll,ih)
c substract effect of target body harmonic on central body
c in appropiate orbiter routine (we ignore effect on sun
c of target body harmonic for sun centered probe)
                     end do
                     sumht(ll,Kk) = sumht(ll,Kk) + sumtc
                  end do
               endif
            endif
            Fn(1) = Fn(1) + Gamat*Masst(ll)*sumht(ll,Kk)
         endif
      end do
c
c*  start=1300
      if(Ncentr.gt.0) then
c
c effect of the sun on the motion of the satellite or probe
         if(Kp(Ncentr+30).ge.0) then
            sums(Kk) = Ccor3(Kk) - Bcor(Kk)/Rb3
            Fn(1)    = Fn(1) + Gamat*sums(Kk)
         endif
c
c*  start=1400
c effect of central body zonal harmonics on equations of motion
         sumhc(Kk) = 0._10
         if(Nczone.gt.1) then
            do i = 1,Nczon1
               di2 = i + 2
               czone(Kk,i) = (di2*Sbcor(Kk)*Cleg(i)/Rsb - Cleg1(i)*
     .                        Cslat1(Kk))/Rsbh(i)/Rsb2
               sumhc(Kk)    = sumhc(Kk) + Czhar(i)*czone(Kk,i)
            end do
         endif
c
c effect of central body tess. harmonics on equations of motion
         if(Nctess.gt.1) then
            ih = 0
            do i = 1,Nctes1
               di2   = i + 2
               sumtc = 0._10
               i1    = i + 1
               do jh = 1,i1
                  djh = jh
                  ih  = ih + 1
                  tessc1(Kk,ih) = (-di2*Sbcor(Kk)*Cgleg(ih)/Rsb +
     .                             Cgleg1(ih)*Cslat1(Kk))/Rsbh(i)/Rsb2
                  tessc2(Kk,ih) = Cgleg(ih)*Clng1(Kk)/Rsbh(i)/Rsb2
                  if(Kk.le.1) then
 
c only once per iteration
                     csnh1(ih) = Cchar(ih)*Cclng(jh) + Cshar(ih)
     .                           *Cslng(jh)
                     csnh2(ih) = -Cchar(ih)*Cslng(jh) + Cshar(ih)
     .                           *Cclng(jh)
                  endif
                  sumtc = sumtc + csnh1(ih)*tessc1(Kk,ih)
     .                    + djh*csnh2(ih)*tessc2(Kk,ih)
               end do
               sumhc(Kk) = sumhc(Kk) + sumtc
            end do
         endif
         Fn(1) = Fn(1) - Gamat3*sumhc(Kk)
      endif
c
c*  start=1600
c effect of central force on motion of the satellite or probe
      if(Nplnt.ne.3) then
 
c for everything but earth-moon barycenter
         Sbforc(Kk) = Sbcor(Kk)/Rsb3
      else if(Kp(40).lt.0) then
         Sbforc(Kk) = Sbcor(Kk)/Rsb3
      else
 
c for earth-moon barycenter
         Sbforc(Kk) = Masse*Xes(Kk)/Res3 + Mass(10)*Xms(Kk)/Rms3
      endif
c
c*  start=1500
c additional accelerations for sun centered probe
      if(lcentr.le.0) then
         call SOLPRB(k)
 
c additional accelerations for earth satellite
      else if(lcentr.eq.3) then
         call ERTORB(k)
c
c additional accelerations for lunar orbiter
      else if(lcentr.eq.10) then
         call MONORB(k)
c
c additional accelerations for planetary orbiter
c (artificial space probe)
      else if(Nplnt.le.30) then
c
c additional accelerations for natural planetary satellites
         call PLNSAT(k)
      else
         call PLNORB(k)
      endif

      if (Nmmsc1.le.0) then
         sum = Gamat3*Sbforc(Kk)
      else
c for possibility of mascons altering central body force
         sum = MSCFRC(Kk,Gamat3,Sbforc)
      endif
      if(Kp(100).lt.0) then
c no change for true orbit integration
c     else if(Kp(100).gt.0) then
c subtract mean orbit acceleration (use osculating orbit for now)
      else
c subtract elliptic orbit acceleration
         sum = sum - Gamt30*Ylpt(Kk)/Rylpt3
      endif
 
c add central force to other perturbing forces
      Fn(1) = Fn(1) + sum
      return
c*  start=3400
c
c effect of additional accelerations for solar probe on partials
  400 if(lcentr.le.0) then
         call SOLPRB(k)
c
c effect of additional accelerations for earth satel on partials
      else if(lcentr.eq.3) then
         call ERTORB(k)
c
c effect of additional accelerations for lunar orbit on partials
      else if(lcentr.eq.10) then
         call MONORB(k)
c
c effect of additional accelerations for planet orb. on partials
c (artificial space probe)
      else if(Nplnt.le.30) then
c
c effect of additional accel. for planet natural satellites
         call PLNSAT(k)
      else
         call PLNORB(k)
      endif
c*  start=3500
c
c check if alpha is initial condition of target planet
      itg   = (ikkk - 1)/100
      psump = 0._10
      if(itg.gt.0) then
         kkk1 = ikkk - itg*100
         if(kkk1.le.6) then
            i=Npltg(itg)
            if(i.gt.0) then
               if(Ntopt(i).gt.0) Fn(1) = Fn(1)
     .          - DOT(Dtcor(1,kkk1,i),Dadxt(1,Kk,i))
               goto 500
            endif
c*  start=3600
c
c check if alpha is init. cond. of central body
            if(itg.eq.Ncentr) then
 
c effect of sun
               if(Kp(Ncentr+30).gt.0) then
 
c psumt1= dot( bcor(1),dccor(1,kkk1) )
                  psumt2 = DOT(Ccor(1),Dccor(1,kkk1))
                  Fn(1)  = Fn(1) + DOT(Dccor(1,kkk1),Dadxs(1,Kk))
                  psump  = psump + (Dccor(Kk,kkk1) - 3._10*Ccor(Kk)
     .                     *psumt2/Rc2)/Rc3
               endif
 
c effect of perturbing planets
               do l = 1,10
                  if(Kp(l+30).gt.0) then
                     if(l.ne.Ncentr) then
                        psumt2 = DOT(Pccor(1,l),Dccor(1,kkk1))
                        psump  = psump + Mass1(l)
     .                           *((Dccor(Kk,kkk1)-3._10*Pccor(Kk,l)
     .                           *psumt2/Rpc2(l))/Rpc3(l))
                     endif
                  endif
               end do
               Fn(1) = Fn(1) + DOT(Dccor(1,kkk1),Dadxp(1,Kk))
c           effect of perturbing asteroids
c     if(numast.le.0) goto 2328
c     do 2327 l=1,numast
c     if(kpast(l).le.0) goto 2327
c           eliminate concentric satellites
c     if(ncnast(l).eq.ncentr) goto 2327
c2327 continue
c           effect of target body harmonics
               if(Numtar.gt.0) then
                  do iq = 1,3
                     sumpar = 0._10
                     do l = 1,Numtar
                        if(Ntopt(l).gt.0) then
 
c eliminate concentric targets
                           kt = Ktrg(l)
                           if(kt.le.u_mxpl) then
                              if(Npcent(kt).ne.Ncentr)
     .                            sumpar = sumpar + Dadxt(iq,Kk,l)
                           endif
                        endif
                     end do
                     Fn(1) = Fn(1) + sumpar*Dccor(iq,kkk1)
                  end do
               endif
               Fn(1) = Fn(1) + Gamat*psump
            endif
         endif
      endif
c*  start=4000
c effect of central force on partial derivatives
  500 if(dflgc) then
         if(Nplnt.eq.3) then
            if(Kp(40).ge.0) then
 
c for earth-moon barycenter
               termd = Gamat3*(Masse/Res3 + Mass(10)/Rms3)
               tg    = -Gamat3*3._10
               do iq = 1,3
                  do jq = 1,iq
                     termg = tg*(Masse*Xes(iq)*Xes(jq)/Res5 + Mass(10)
     .                       *Xms(iq)*Xms(jq)/Rms5)
                     Dadxc(iq,jq) = termg
                     if(iq.eq.jq) goto 510
                     Dadxc(jq,iq) = termg
                  end do
  510             Dadxc(iq,iq) = Dadxc(iq,iq) + termd
               end do
               goto 600
            endif
         endif
 
c for everything but earth-moon barycenter
         termd = Gamat3/Rsb3
         tg    = -termd*3._10/Rsb2
         do iq = 1,3
            do jq = 1,iq
               termg = tg*Sbcor(iq)*Sbcor(jq)
               Dadxc(iq,jq) = termg
               if(iq.eq.jq) goto 520
               Dadxc(jq,iq) = termg
            end do
  520       Dadxc(iq,iq) = Dadxc(iq,iq) + termd
         end do
      endif
c*  start=4500
c form dadx each iteration
  600 if(k.eq.10) then
         ll = 5 + Numtar
         do iq = 1,9
            sumda = 0._10
            do l = 1,ll
               sumda = sumda + dadxx(iq,1,l)
            end do
            Dadx(iq,1) = sumda
            Dadxsav(iq,1,1,1,1) = sumda
         end do
         if(Pntflg) then
            do jq = 1,ll
               iq = min0(6,jq)
               if(lggfg(iq)) then
                  if(Line.ge.52) then
                     call NEWPGT(Iout,Npage,24000)
                     Line = 2
                  endif
                  Line = Line + 4
                  write(Iout,610) mtxid(jq + 1),
     .                             (dadxx(iq,1,jq),iq = 1,9)
               endif
            end do
            write(Iout,610) mtxid(1),(Dadx(iq,1),iq = 1,9)
            if(Posflg) then
               do jb = 1,2
                  if(Line.ge.52) then
                     call NEWPGT(Iout,Npage,24000)
                     Line = 2
                  endif
                  write(Iout,610) mtxid(jb+10),
     .                             (Dadxl(ib,3*jb-2),ib = 1,9)
                  Line = Line + 4
               end do
  610          format('0   IN SBFN1  DADX',a2,(t26,1p,3D25.15))
            endif
            Line = Line + 6
            if(Line.lt.59) write(Iout,620)
  620       format('0')
            Pntflg = .false.
         endif
      endif
 
c note - dadx is symmetric except for dadxl (drag contribution)
      fnfb = DOT(Dadx(1,Kk),Dsbcor(1,Kkk))
      if(Kkk.eq.1) then
         Parsav(Kk,1)=Dsbcor(Kk,1)
         Parsav(Kk+3,1)=Dsbcor(Kk+3,1)
         Indpsav(Kk,1,1)=fnfb
      endif
      if(Posflg) fnfb = fnfb + DOT(Dadxl(1,Kk),Dsbcor(1,Kkk))
      if(Velflg) fnfb = fnfb + DOT(Dadxm(1,Kk),Dsbcor(4,Kkk))
c account for dependence of earth and moon coordinates on mass(10)
      if(ikkk.eq.10 .and. Kp(40).gt.0) then
         if(Nplnt.eq.3) then
            fnfb = fnfb - DOT(Dadx(1,Kk),Xpert(1,10))
         else if(Npltg(3).gt.0) then
            fnfb = fnfb - DOT(Dadp(1,Kk,Npltg(3)),Xpert(1,10))
         endif
      endif
      Fn(1) = Fn(1) + fnfb
c*  start=4700
c is this initial condition partial
      icn  = -30 - ikkk
      psum = 0._10
      if(icn.gt.0) then
         if(Kp(100).lt.0) then
c no change for true orbit integration
c        else if(Kp(100).gt.0) then
c subtract mean orbit acceleration partial (use osculating orbit for now)
         else
c subtract elliptic orbit acceleration partial
            psum = psum +
     .             Gamt30*(3._10*Ylpt(Kk)*(Ylpt(1)*Dylpt(1,icn)+Ylpt(2)
     .             *Dylpt(2,icn)+Ylpt(3)*Dylpt(3,icn))
     .             /Rylpt2 - Dylpt(Kk,icn))/Rylpt3
         endif
      endif
c*  start=4800
c check if alpha is time variable gravitational constant
      if(Kp(62).ge.1) then
         if(ikkk.eq.32) then
            sumpar = Gama3*Tvary
            goto 700
         endif
      endif
 
c check if alpha is pertinent mass
      if(ikkk.le.0 .or. ikkk.gt.30) then
         goto 800
      else
         if(Massfl(ikkk).eq.0) goto 800
         if(Massfl(ikkk).lt.1) then
 
c mass of another satellite in system
            sumpar = -Gamat3/Massp1
         else if(Massfl(ikkk).eq.1) then
 
c mass of asteroid (or earth-moon if nplnt=3)
            sumpar = Gamat3/Massp1
         else
 
c mass of central body (or earth-moon if ncentr=10)
            sumpar = Gamat3/Mass(ikkk)
         endif
      endif
  700 psum = psum + sumpar*Sbforc(Kk)
      if(Ncentr.gt.0) psum = psum - sumpar*sumhc(Kk)
  800 Fn(1) = Fn(1) + psum
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=9900
c
c          return right side of differential equation
      return
      end
