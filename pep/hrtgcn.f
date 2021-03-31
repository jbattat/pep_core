      subroutine HRTGCN(ncall)
 
      implicit none

c           subtract effect of target body harmonics on the central body
c           called from appropriate orbiter routine
c        j.f.chandler 1977 april
c        ncall=-2 setup once per iteration of a given step for partials
c        ncall=-1 setup once per iteration of a given step for motion
c        ncall= 0 setup once per step
c        ncall= k evaluate acceleration for equation k.gt.0
c
c array dimensions
      include 'globdefs.inc'
c common
      include 'cnthar.inc'
      include 'emmips.inc'
      include 'fcntrl.inc'
      include 'fmmips.inc'
      include 'intstf.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'sbembr.inc'
      include 'sbrot.inc'
      include 'sbstuf.inc'
      include 'sbthng.inc'
      include 'trghar.inc'
      include 'yvectplp.inc'
 
c local
      logical   dflg
      real*10 rpch(i_mxtrg,5),tlec(4,i_mxtrg),tlec1(4,i_mxtrg),
     . tlec2(4,i_mxtrg),tslac1(3,i_mxtrg),tglec(5,i_mxtrg),
     . tglec1(5,i_mxtrg),tglec2(5,i_mxtrg),tlnc1(3,i_mxtrg),
     . tclnc(i_mxtrg,3),tslnc(i_mxtrg,3),tclac(i_mxtrg),tslac(i_mxtrg),
     . tclacr(i_mxtrg),rpch1(i_mxtrg),tzone(3,i_mxtrg,4),
     . sumht(i_mxtrg,3),tesct1(3,i_mxtrg,5),tesct2(3,i_mxtrg,5),
     . tsnc1(i_mxtrg,5),tsnc2(i_mxtrg,5),rpkh(i_mxtrg,3),
     . tlek(3,i_mxtrg),tlek1(3,i_mxtrg),tlek2(3,i_mxtrg),
     . tslak1(3,i_mxtrg),tclak(i_mxtrg),tslak(i_mxtrg),rpkh1(i_mxtrg),
     . tkzone(3,3),dadct(3,3,i_mxtrg),talnc1(3,i_mxtrg)
      real*10 di2,di3,di4,djh,DOT,dum,sumpar,td,termg,tg,
     .          tgaa,tgar,tgax,tgrr,tgrrp,tgrx,tgxx,tkz,tm1
      real*10 tm2,tq,tqp,zd,zgaa,zgax,zgxx
      integer   i,i1,ib,ih,ii,ikkk,it,itg,itg1,jb,jh,kt,l,ll,
     .          ncall,ntt,ntz
      logical   pntflh
      integer*4 iout/6/
 
      if(ncall .lt. 0) then
c*  start=2000
c
c        set-up once per iteration of a given step for motion
c*  start=3000
c
c        set-up once per iteration for partials
         if(ncall.lt.-1) pntflh = Pntflg
      else if(ncall.eq.0) then
c
c*  start=1000
c set-up once per step
         if(Numtar.gt.0) then
            dflg =.true.
            if(Ntopct.gt.0) then
 
c clear summed effect of central body harmonics
               do i = 1, Ntopct
                  do ii = 1, 3
                     tkzone(ii, i) = 0._10
                  end do
               end do
            endif
 
c determine target body harmonic quantities
            do ll = 1, Numtar
               do ii = 1, 3
                  sumht(ll, ii) = 0._10
               end do
c
c determine target body effect on central body harmonics
               if(Ntopct.gt.0) then
 
c determine latitude
                  tslak(ll) = (Cntrot(3,1)*Tccor(1,ll) + Cntrot(3,2)
     .                        *Tccor(2,ll) + Cntrot(3,3)*Tccor(3,ll))
     .                        /Rtc(ll)
                  tclak(ll) = SQRT(1._10 - tslak(ll)**2)
                  do i = 1, 3
                     tslak1(i, ll) = Cntrot(3, i) - Tccor(i, ll)
     .                               *tslak(ll)/Rtc(ll)
                  end do
 
                  rpkh1(ll)   = Rtc(ll)/Crad
                  rpkh(ll, 1) = rpkh1(ll)**2
                  if(Ntopct.gt.1) then
                     do i = 2, Ntopct
                        rpkh(ll, i) = rpkh(ll, i - 1)*rpkh1(ll)
                     end do
                  endif
 
                  ntz = Ntopct + 1
                  ntt = 0
c
c determine p(n), p'(n)
c
                  call LEGNDR(tslak(ll), tclak(ll), ntz, ntt, tlek(1,ll)
     .                        , tlek1(1,ll), dum, dum)
c
c*  start=1100
c
c effect of target body+cntr hrmncs on equations of motion
                  ntz = Ntopct
 
c if(rtc(ll).gt.hsitb(ll) ) ntz=1
                  do i = 1, ntz
                     di2 = i + 2
                     do ii = 1, 3
                        tkz = (di2*Tccor(ii,ll)*tlek(i,ll)/Rtc(ll)
     .                        - tlek1(i,ll)*tslak1(ii,ll))/rpkh(ll, i)
     .                        /Rtc2(ll)
                        sumht(ll, ii) = sumht(ll, ii) - Czhar(i)*tkz
                        tkzone(ii, i) = tkzone(ii, i) + tkz*Masst(ll)
                     end do
                  end do
               endif
 
               if(Ntopt(ll).gt.0) then
 
c determine latitude
                  tslac(ll) = -(Trgrot(3,1,ll)*Tccor(1,ll) + Trgrot(3,2,
     .                        ll)*Tccor(2,ll) + Trgrot(3,3,ll)
     .                        *Tccor(3,ll))/Rtc(ll)
                  tclac(ll) = SQRT(1._10 - tslac(ll)**2)
                  do i = 1, 3
                     tslac1(i, ll) = Trgrot(3, i, ll) + Tccor(i, ll)
     .                               *tslac(ll)/Rtc(ll)
                  end do
                  tclacr(ll)  = tclac(ll)*Rtc(ll)
                  rpch1(ll)   = Rtc(ll)/Trad(ll)
                  rpch(ll, 1) = rpch1(ll)**2
                  if(Ntopt(ll).gt.1) then
                     ntz = Ntopt(ll)
                     do i = 2, ntz
                        rpch(ll, i) = rpch(ll, i - 1)*rpch1(ll)
                     end do
                  endif
c
c determine longitude
                  ntt = Nttess(ll)
                  if(ntt.gt.1) then
                     tslnc(ll, 1) = -(Trgrot(2,1,ll)*Tccor(1,ll) +
     .                              Trgrot(2,2,ll)*Tccor(2,ll)
     .                              + Trgrot(2,3,ll)*Tccor(3,ll))
     .                              /tclacr(ll)
                     tclnc(ll, 1) = -(Trgrot(1,1,ll)*Tccor(1,ll) +
     .                              Trgrot(1,2,ll)*Tccor(2,ll)
     .                              + Trgrot(1,3,ll)*Tccor(3,ll))
     .                              /tclacr(ll)
                     do i = 2, ntt
                        tslnc(ll, i) = tslnc(ll, i - 1)*tclnc(ll, 1)
     .                                 + tclnc(ll, i - 1)*tslnc(ll, 1)
                        tclnc(ll, i) = tclnc(ll, i - 1)*tclnc(ll, 1)
     .                                 - tslnc(ll, i - 1)*tslnc(ll, 1)
                     end do
                     do i = 1, 3
                        tlnc1(i, ll) = (Trgrot(2,i,ll)*tclnc(ll,1) -
     .                                 Trgrot(1,i,ll)*tslnc(ll,1))
     .                                 /tclac(ll)
                     end do
                  endif
                  ntz = Ntzone(ll)
c
c determine p(n), p'(n), p(n,h), p'(n,h)
c
                  call LEGNDR(tslac(ll), tclac(ll), ntz, ntt, tlec(1,ll)
     .                        , tlec1(1,ll), tglec(1,ll), tglec1(1,ll))
c
c
c        effect of target  body  harmonics on equations of motion
c*  start=1200
c
c        effect of target body  zonal   harmonics on
c        equations of motion
                  if(Ntzone(ll).gt.1) then
                     ntz = Ntzon1(ll)
                     if(Rtc(ll).gt.Hsitb(ll)) ntz = 1
                     do i = 1, ntz
                        di2 = i + 2
                        do ii = 1, 3
                           tzone(ii, ll, i)
     .                        = (-di2*Tccor(ii,ll)*tlec(i,ll)/Rtc(ll)
     .                        - tlec1(i,ll)*tslac1(ii,ll))/rpch(ll, i)
     .                        /Rtc2(ll)
                           sumht(ll, ii) = sumht(ll, ii) + Tzhar(ll, i)
     .                        *tzone(ii, ll, i)
                        end do
                     end do
                  endif
c
c effect of target body tesseral harmonic on eq. of motion
                  if(Nttess(ll).gt.1) then
                     ih  = 0
                     ntt = Nttes1(ll)
                     if(Rtc(ll).le.Hsitb(ll)) then
                        do i = 1, ntt
                           di2 = i + 2
                           i1  = i + 1
                           do jh = 1, i1
                              djh = jh
                              ih  = ih + 1
                              tsnc1(ll, ih) = Tchar(ll, ih)
     .                           *tclnc(ll, jh) + Tshar(ll, ih)
     .                           *tslnc(ll, jh)
                              tsnc2(ll, ih) = -Tchar(ll, ih)
     .                           *tslnc(ll, jh) + Tshar(ll, ih)
     .                           *tclnc(ll, jh)
                              do ii = 1, 3
                                 tesct1(ii, ll, ih)
     .                              = (di2*Tccor(ii,ll)*tglec(ih,ll)
     .                              /Rtc(ll) + tglec1(ih,ll)
     .                              *tslac1(ii,ll))/rpch(ll, i)/Rtc2(ll)
                                 tesct2(ii, ll, ih) = tglec(ih, ll)
     .                              *tlnc1(ii, ll)/rpch(ll, i)/Rtc2(ll)
                                 sumht(ll, ii) = sumht(ll, ii)
     .                              + tsnc1(ll, ih)*tesct1(ii, ll, ih)
     .                              + djh*tsnc2(ll, ih)
     .                              *tesct2(ii, ll, ih)
                              end do
                           end do
                        end do
                     endif
                  endif
               endif
 
c*  start=1300
            end do
         endif
c
c*  start=4000
c
c
c        compute acceleration for motion
      else if(Kkk.gt.0) then
c
c*  start=5000
c
c           compute accelerations for partials
c        effect of target  body harmonics on partial derivatives
         if(Numtar.gt.0) then
 
            ikkk = Icntrl(Kkk)
            itg  = (ikkk - 1)/100
            if(itg.le.0) then
c*  start=6000
c check for other parameters
               do ll = 1, Numtar
                  if(Ntopt(ll).gt.0 .or. Ntopct.gt.0) then
 
c check if alpha is mass of target planet
                     if(ikkk.ne.Ntrg(ll)) then
 
c check if alpha is mass of target planet's center
                        kt = Itgast(ll)
                        if(kt.gt.0) then
                           if(ikkk.eq.Ncnast(kt)) then
                              Fn(1) = Fn(1) - Gamat*Astmab(kt)
     .                                *sumht(ll, Kk)
                              go to 10
                           endif
                        endif
 
c check if alpha is time variable gravitational constant
                        if(ikkk.eq.32) Fn(1) = Fn(1)
     .                      - Gama*Tvary*Masst(ll)*sumht(ll, Kk)
                     else
                        Fn(1) = Fn(1) - Gamat*Masstc(ll)*sumht(ll, Kk)
                     endif
                  endif
   10          end do
            else
               itg1 = ikkk - itg*100
               kt   = -1
 
c see if alpha is center or target parameter
               do l = 1, Numtar
                  if(itg.eq.Ntrg(l)) then
                     kt = l
                     go to 20
                  endif
               end do
               if(itg.ne.Ncentr) return
               kt = 0
 
c check for harmonic coefficient
   20          if(itg1.gt.30) then
 
c*  start=5600
                  sumpar = 0._10
                  if(kt.gt.0) then
 
c alpha is target coefficient
                     do ll = 1, Numtar
                        if(Ntopt(ll).gt.0) then
                           if(Ntzone(ll).gt.1 .and. ikkk.eq.Jtzone(ll))
     .                     then
                              ntz = Ntzon1(ll)
                              if(Rtc(ll).gt.Hsitb(ll)) ntz = 1
                              do i = 1, ntz
 
c check if alpha is zonal harmonic coefficient
                                 if(Itzone(ll,i).eq.Kkk) then
                                    sumpar = tzone(Kk, ll, i)
c*  start=5800
c subtract term for harmonic coefficient
                                    Fn(1) = Fn(1)-Gamat*Masst(ll)*sumpar
                                    return
                                 endif
                              end do
                           endif
 
c effect of target body tesseral harmonics on partials
                           if(Nttess(ll).gt.1) then
                              if(Rtc(ll).le.Hsitb(ll)) then
 
c check if alpha is tess. harmonic coeff.
                                 if( ikkk.eq.Jtcos(ll) .or.
     .                               ikkk.eq.Jtsin(ll) ) then
                                    ih  = 0
                                    ntt = Nttes1(ll)
                                    do i = 1, ntt
                                       i1 = i + 1
                                       do jh = 1, i1
                                         ih = ih + 1
 
c check if alpha is tess. harmonic cosine coefficient
                                         if(ikkk.eq.Jtcos(ll) .and.
     .                                   Itcos(ll,ih).eq.Kkk) then
                                         djh    = jh
                                         sumpar = tclnc(ll, jh)
     .                                      *tesct1(Kk, ll, ih)
     .                                      - tslnc(ll, jh)
     .                                      *tesct2(Kk, ll, ih)*djh
                                         Fn(1)  = Fn(1) -
     .                                            Gamat*Masst(ll)*sumpar
                                         return
                                         endif
 
c check if alpha is tess. harmonic sine coefficient
                                         if(ikkk.eq.Jtsin(ll) .and.
     .                                   Itsin(ll,ih).eq.Kkk) then
                                         djh    = jh
                                         sumpar = tslnc(ll, jh)
     .                                      *tesct1(Kk, ll, ih)
     .                                      + tclnc(ll, jh)
     .                                      *tesct2(Kk, ll, ih)*djh
                                         Fn(1)  = Fn(1) -
     .                                            Gamat*Masst(ll)*sumpar
                                         return
                                         endif
                                       end do
                                    end do
                                 endif
                              endif
                           endif
                        endif
                     end do
 
c alpha is center coefficient
                  else if(ikkk.eq.Jczone) then
                     if(Ntopct.gt.0) then
                        ntz = Ntopct
                        do i = 1, ntz
                           if(Iczone(i).eq.Kkk) then
                              Fn(1) = Fn(1) + Gamat*tkzone(Kk, i)
                              return
                           endif
                        end do
                     endif
                  endif
               else if(itg1.le.6) then
                  if(dflg) then
                     dflg = .false.
c
c effect of target  body harmonics on partial derivatives
c compute grav. grad. matrices
                     do ll = 1, Numtar
 
c*  start=5100
                        if(Ntopt(ll).gt.0 .or. Ntpctp.gt.0) then
                           tg = Masst(ll)*Gamat/Rtc3(ll)
 
c clear out grav. grad. matrix
                           do ib = 1, 9
                              dadct(ib, 1, ll) = 0._10
                           end do
                           if(Ntpctp.gt.0) then
c
c effect of target body+cntr hrmncs on partials
c
c determine p''(n)
                              ntz = Ntpctp + 1
                              ntt = 0
                              call LEGND2(tslak(ll), tclak(ll), ntz,
     .                           ntt, tlek(1,ll), tlek1(1,ll),
     .                           tlek2(1,ll), dum, dum, dum)
 
                              zgaa = 0._10
                              zgax = 0._10
                              zgxx = 0._10
                              zd   = 0._10
                              ntz  = Ntpctp
 
c if(rtc(ll).gt.hsitb(ll) ) ntz=1
                              di3 = 3._10
                              di4 = 4._10
                              do i = 1, ntz
                                 di2  = di3
                                 di3  = di4
                                 di4  = i + 4
                                 tm2  = Czhar(i)/rpkh(ll, i)*tg
                                 zgaa = zgaa - tlek2(i, ll)*tm2
                                 tm1  = tm2*(di2*tlek(i,ll) + tslak(ll)
     .                                  *tlek1(i,ll))
                                 tm2  = tm2*(di3*tlek1(i,ll) + tslak(ll)
     .                                  *tlek2(i,ll))
                                 zgax = zgax - tm2
                                 zgxx = zgxx - di4*tm1 - tslak(ll)*tm2
                                 zd   = zd - tm1
                              end do
 
                              zgax = zgax/Rtc(ll)
                              zgxx = zgxx/Rtc2(ll)
                              zd   = -zd
                              do ib = 1, 3
                                 do jb = 1, 3
                                    termg = zgxx*Tccor(ib, ll)
     .                                 *Tccor(jb, ll)
     .                                 - zgax*(Tccor(ib,ll)*Cntrot(3,jb)
     .                                 + Tccor(jb,ll)*Cntrot(3,ib))
     .                                 + zgaa*Cntrot(3, ib)
     .                                 *Cntrot(3, jb)
                                    dadct(ib, jb, ll) = termg
                                    if(ib.eq.jb) go to 30
                                    dadct(jb, ib, ll) = termg
                                 end do
   30                            dadct(ib, ib, ll) = dadct(ib, ib, ll)
     .                              + zd
                              end do
                           endif
c
c*  start=5200
                           if(Ntopt(ll).gt.0) then
c determine target  body harmonic quantities for partials
c determine p''(n), p''(n,h)
                              ntz = Ntzonp(ll) + 1
                              ntt = Nttesp(ll) + 1
                              call LEGND2(tslac(ll), tclac(ll), ntz,
     .                           ntt, tlec(1,ll), tlec1(1,ll),
     .                           tlec2(1,ll), tglec(1,ll), tglec1(1,ll),
     .                           tglec2(1,ll))
c
c effect of target body  zonal harmonics on partials
                              if(Ntzonp(ll).gt.0) then
                                 zgaa = 0._10
                                 zgax = 0._10
                                 zgxx = 0._10
                                 zd   = 0._10
                                 ntz  = Ntzonp(ll)
                                 if(Rtc(ll).gt.Hsitb(ll)) ntz = 1
                                 di3 = 3._10
                                 di4 = 4._10
                                 do i = 1, ntz
                                    di2  = di3
                                    di3  = di4
                                    di4  = i + 4
                                    tm2  = Tzhar(ll, i)/rpch(ll, i)
                                    zgaa = zgaa + tlec2(i, ll)*tm2
                                    tm1  = tm2*(di2*tlec(i,ll)
     .                                  + tslac(ll)*tlec1(i,ll))
                                    tm2  = tm2*(di3*tlec1(i,ll)
     .                                  + tslac(ll)*tlec2(i,ll))
                                    zgax = zgax + tm2
                                    zgxx = zgxx + di4*tm1 + tslac(ll)
     .                                 *tm2
                                    zd   = zd + tm1
                                 end do
                                 zgaa = zgaa*tg
                                 zgax = zgax*tg/Rtc(ll)
                                 zgxx = zgxx*tg/Rtc2(ll)
                                 zd   = -zd*tg
                                 do ib = 1, 3
                                    do jb = 1, 3
                                       termg = zgxx*Tccor(ib, ll)
     .                                    *Tccor(jb, ll)
     .                                    + zgax*(Tccor(ib,ll)
     .                                    *Trgrot(3,jb,ll)
     .                                    + Tccor(jb,ll)*Trgrot(3,ib,ll)
     .                                    ) + zgaa*Trgrot(3, ib, ll)
     .                                    *Trgrot(3, jb, ll)
                                       dadct(ib, jb, ll)
     .                                    = dadct(ib, jb, ll) + termg
                                       if(ib.eq.jb) go to 40
                                       dadct(jb, ib, ll)
     .                                    = dadct(jb, ib, ll) + termg
                                    end do
   40                               dadct(ib, ib, ll)
     .                                 = dadct(ib, ib, ll) + zd
                                 end do
                              endif
c
c effect of target body tesseral harmonics on partials
                              if(Nttesp(ll).gt.0) then
                                 if(Rtc(ll).le.Hsitb(ll)) then
                                    ntt   = Nttesp(ll)
                                    tgxx  = 0._10
                                    tgax  = 0._10
                                    tgrx  = 0._10
                                    tgrr  = 0._10
                                    tgar  = 0._10
                                    tgaa  = 0._10
                                    tgrrp = 0._10
                                    td    = 0._10
                                    do ib = 1, 3
                                       talnc1(ib, ll)
     .                                    = -Trgrot(2, ib, ll)
     .                                    *tslnc(ll, 1)
     .                                    - Trgrot(1, ib, ll)
     .                                    *tclnc(ll, 1)
                                    end do
                                    ih  = 0
                                    di3 = 3._10
                                    di4 = 4._10
                                    do i = 1, ntt
                                       di2 = di3
                                       di3 = di4
                                       di4 = i + 4
                                       i1  = i + 1
                                       do jh = 1, i1
                                         djh   = jh
                                         ih    = ih + 1
                                         tq    = tsnc1(ll, ih)/rpch(ll,
     .                                           i)
                                         tqp   = djh*tsnc2(ll, ih)
     .                                        /rpch(ll, i)
                                         tm1   = di2*tglec(ih, ll)
     .                                        + tslac(ll)*tglec1(ih, ll)
                                         tm2   = di3*tglec1(ih, ll)
     .                                        + tslac(ll)*tglec2(ih, ll)
                                         tgar  = tgar + tglec1(ih, ll)
     .                                       *tqp
                                         tgrrp = tgrrp + tglec(ih, ll)
     .                                      *tqp
                                         tgaa  = tgaa + tglec2(ih,
     .                                           ll)*tq
                                         tgrx  = tgrx + tqp*tm1
                                         tgax  = tgax + tq*tm2
                                         tgxx  = tgxx +
     .                                       tq*(di4*tm1 +
     .                                       tslac(ll)*tm2)
                                         tgrr = tgrr +
     .                                      djh*djh*tq*tglec(ih, ll)
                                         td   = td + tq*tm1
                                       end do
                                    end do
                                    tgar  = tgar*tg
                                    tgrrp = tgrrp*tg/tclac(ll)
                                    tgaa  = tgaa*tg
                                    tgrx  = tgrx*tg/Rtc(ll)
                                    tgax  = tgax*tg/Rtc(ll)
                                    td    = -td*tg
                                    tgrr  = -tgrr*tg
                                    do ib = 1, 3
                                       do jb = 1, 3
                                         termg = tgxx*Tccor(ib, ll)
     .                                      *Tccor(jb, ll)
     .                                      + tgax*(Tccor(ib,ll)
     .                                      *Trgrot(3,jb,ll)
     .                                      + Tccor(jb,ll)
     .                                      *Trgrot(3,ib,ll))
     .                                      + tgrx*(Tccor(ib,ll)
     .                                      *tlnc1(jb,ll) + Tccor(jb,ll)
     .                                      *tlnc1(ib,ll))
     .                                      + tgar*(Trgrot(3,ib,ll)
     .                                      *tlnc1(jb,ll)
     .                                      + Trgrot(3,jb,ll)
     .                                      *tlnc1(ib,ll))
     .                                      + tgaa*Trgrot(3, ib, ll)
     .                                      *Trgrot(3, jb, ll)
     .                                      + tgrrp*(tlnc1(jb,ll)
     .                                      *talnc1(ib,ll)
     .                                      + tlnc1(ib,ll)*talnc1(jb,ll)
     .                                      ) + tgrr*tlnc1(ib, ll)
     .                                      *tlnc1(jb, ll)
                                         dadct(ib, jb, ll)
     .                                      = dadct(ib, jb, ll) + termg
                                         if(ib.eq.jb) go to 50
                                         dadct(jb, ib, ll)
     .                                      = dadct(jb, ib, ll) + termg
                                       end do
   50                                  dadct(ib, ib, ll)
     .                                    = dadct(ib, ib, ll) + td
                                    end do
                                 endif
                              endif
                           endif
                           if(pntflh) then
                              if(Line.ge.56) then
                                 call NEWPGT(iout, Npage, 24000)
                                 Line = 2
                              endif
                              Line = Line + 4
                              write(iout, 60) ll,
     .                              (dadct(ib,1,ll), ib = 1, 9)
   60                         format('0  IN HRTGCN  DADCT', i1,
     .                               (t26,1p,3D25.15))
                           endif
                        endif
                     end do
                     pntflh = .false.
                  endif
 
c*  start=5500
                  if(itg1.le.6) then
                     if(kt.le.0) then
 
c alpha is i.c. of central planet
                        do ll = 1, Numtar
                           if(Ntopt(ll).gt.0 .or. Ntpctp.gt.0) then
 
c no effect if target is concentric
                              it = Itgast(ll)
                              if(it.le.0 .or. Ncnast(it).ne.Ncentr)
     .                            Fn(1) = Fn(1)
     .                            - DOT(Dccor(1,itg1), dadct(1,Kk,ll))
                           endif
                        end do
 
c alpha is i.c. of  target planet
                     else if(Ntopt(kt).gt.0 .or. Ntpctp.gt.0) then
 
c some other body parameter(not grav.harmonic)
                        Fn(1) = Fn(1)
     .                          + DOT(Dympt(1,itg1,kt), dadct(1,Kk,kt))
                     endif
                  endif
               endif
            endif
         endif
      else if(Numtar.gt.0) then
c
c subtract effect of target body harmonic on central body
         do ll = 1, Numtar
            if(Ntopt(ll).gt.0 .or. Ntopct.gt.0) Fn(1) = Fn(1)
     .          - Gamat*Masst(ll)*sumht(ll, Kk)
         end do
      endif
c*  start=9000
c return with corrections
      return
      end
