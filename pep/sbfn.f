      real*10 function SBFN(k, j, s)

      implicit none

c
c ash/amuchastegui/friedman - june 1969 - real*8 function sbfn
c evaluation of right side of satellite-probe differenial equations
c updated feb. 1980   kcl: tcon added
c

c arguments
      integer k,j
      real*10 s
c           k  =   equation number
c           j  =   iteration number
c           s  =   time (julian date+0.5)
c
c array dimensions
      include 'globdefs.inc'
c           common
      include 'cnthar.inc'
      include 'ellips.inc'
      include 'fcntrl.inc'
      include 'fmmips.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'orblun.inc'
      include 'output.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'prtcod.inc'
      include 'sbembr.inc'
      include 'sbrot.inc'
      include 'sbstuf.inc'
      logical*1 lggfgs(16)
      equivalence (lggfgs, Ggfgs)
      include 'sbthng.inc'
      include 'stint.inc'
      include 'tapdtplp.inc'
      include 'trghar.inc'
      include 'xprcom.inc'
      include 'yvectplp.inc'
c
c local
      logical*1 tlah(8)/8*.false./
      character*8 halt
      equivalence (halt, tlah)
      real*10 julkm, rpbh1(i_mxtrg)
      real*10 dum,rpbkm,rsbh1,tfact,td,tg,rtb5,rtc5,gm
      integer i,i1,icn,jgg,kt,l,l1,ll,nctp1,nczp1,nt,ntt,ntz
c external functions
      real*10 DOT,MSCFR1

c
c nplnt  = 1,...,30 with ncentr.lt.0  integrated body is a planet
c             or asteroid orbiting the sun
c nplnt  =11,...,30 with ncentr.gt.0  integrated body is a natural
c             satellite
c nplnt=31,32,...  integrated body is an artificial space probe
c (this routine would only be called in one of these cases)
c
c ncentr= planet number of central body
c             0=sun     1=mercury   2=venus
c             3=earth   4=mars      5=jupiter
c             6=saturn  7=uranus    8=neptune
c             9=pluto   10=moon     negative=sun
c
c for i=1-6 we have in this routine
c sbcor(i) =  coordinates of integrated body relative to central body
c bcor(i)  =  coordinates of integrated body relative to the sun
c ccor(i)  =  coordinates of central body relative to the sun
c for i=1-3 we have
c pccor(i,l)= coordinates of body l relative to central body for
c             kp(30+l).ge.0 and l.ne.ncentr and l.ne.nplnt (l=1,...,10)
c pbcor(i,l)= coordinates of body l relative to integrated body for
c             kp(30+l).ge.0 and l.ne.ncentr and l.ne.nplnt (l=1,...,10)
c l=3 denotes earth-moon barycenter if kp(40).lt.0 and ncentr.ne.10
c l=3 denotes earth if kp(40).ge.0 or ncentr.eq.10
c mass1(l) =  (mass of body l)/(mass of sun), l=1,...,10
c mass(3)  =  (mass of earth+moon)/(mass of sun)
c mass(10) =  (mass of moon)/(mass of earth+moon)
c masse =     (mass of earth)/(mass of earth+moon)
c massp =     (mass of integrated body)/(mass of central body)
c
c rsb2 = sbcor(1)**2+sbcor(2)**2+sbcor(3)**2
c rsb  = sqrt(rsb2)
c rsb3 = rsb2*rsb
c similar definitions for rb,rc,rpc(l),rpb(l) quantities
c kpb  = logical flags to control whether pbcor,pccor are calculated
c
c gama  =     gravitational constant times mass of sun
c     gamat = gama *(1+prmter(32)*(s-t0))
c
c gama3 =    -gravitational constant times mass of central body times
c             (1+massp)
c     gamat3= gama3*(1+prmter(32)*(s-t0))
c
      Fn(1) = 0._10
      Fn(2) = 0._10
      if(Kkp(60).gt.1)
     . write(6,99333) k,j,s-2440000.5_10
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c computations done only once for a given step
      if(k.eq.4) then
         jgg = j
         if(j.ne.2) then
            Npstep = Npstep + 1

c determine time variable gravitational constant
            Tvary = s - Tvary0
            if(Kp(62).ge.0) then
               tfact  = 1._10 + prmter(32)*Tvary
               Gamat  = Gama*tfact
               Gamat3 = Gama3*tfact
            endif
c
c determine perturbing planet coordinates
            Jd    = s
            Fract = Jd
            Fract = s - Fract
            call PRTCRD(Jd, Fract)
c determine perturbing small body coordinates
            if(Jpert.gt.0) call ASTCRD(Jd, Fract)
c determine coordinates and partials of target and central planet
c position and velocity w.r.t. their orbital elements
            do kt=0,Numtar
               call PLPCRDT(Jd,Fract,kt)
               nt=Nplpt(kt)
               if(Nplnt.eq.3 .and. Kp(40).ge.0 .and. nt.eq.10) then
                  do i=1,6
                     Pcorsav(i,3)=Xpert(i,10)
                     if(Parnum(3).gt.0) Parsav(i,3)=
     .                Dtcor(i,Parnum(3),kt)
                  end do
               endif
            end do
            nt=Nplpt(1)
            do i=1,6
               if(Ncentr.gt.0 .and. Ncentr.lt.10) then
                  Pcorsav(i,3)=Xpert(i,Ncentr)
                  if(Parnum(3).gt.0) Parsav(i,3)=Dccor(i,Parnum(3))
               endif
               if(Numtar.gt.0 .and. nt.lt.10) then
                  Pcorsav(i,2)=Xpert(i,nt)
                  if(Parnum(2).gt.0) Parsav(i,2)=Dtcor(i,Parnum(2),1)
               endif
            end do
            if(Kp(100).lt.0) then
c using true orbit quantities
c
c determine mean orbit coordinates (treat as osculating for now)
c           else if(Kp(100).gt.0) then
            else
c
c determine elliptic orbit coordinates for encke's method
               Tlpt = s - Tlpt0
               call ELIPT(Kp1,Tlpt)
               if(Kkp(60).gt.1) write(6,99332) (Ylpt(i),i=1,3),
     .          (Y(i,j),i=1,3)
99332          format('YLPT',1p3e24.17/'Y',3e12.4)
            endif
c
c determine central body relative to sun
            do i = 1,Jndex1
               Ccor(i) = 0._10
               if(Ncentr.gt.0 .and. Kp(Ncentr+30).ge.0) then
                  if(Ncentr.eq.3) then
                     Ccor(i) = Xpert(i,3) - Mass(10)*Xpert(i,10)
                  else if(Ncentr.lt.10) then
                     Ccor(i) = Xpert(i,Ncentr)
                  else if(Ncentr.eq.10) then
                     Ccor(i) = Xpert(i,3) + Masse*Xpert(i,10)
                  endif
               endif
            end do
            if(Ncentr.gt.0) then
               Rc2 = Ccor(1)**2 + Ccor(2)**2 + Ccor(3)**2
               Rc  = SQRT(Rc2)
               Rc3 = Rc2*Rc
               Rc5 = Rc2*Rc3
               if(Kp(Ncentr+30).ge.0) then
                  do i = 1,Jndex1
                     Ccor3(i) = Ccor(i)/Rc3
                  end do
               endif
            endif
c
c determine perturbing planets (not sun, central body or
c integrated body) relative to central body
            do l = 1,10
               if(Kpb(l)) then
                  do i = 1,3
                     if(l.lt.10) then
                        if(l.eq.3) then
                           if(Ncentr.eq.10) then
                              Pccor(i,l) = -Xpert(i,10)
                              goto 30
                           else if(Kp(40).ge.0) then
                              Pccor(i,l)
     .                           = (Xpert(i,3) - Ccor(i))
     .                           - Mass(10)*Xpert(i,10)
                              goto 30
                           endif
                        endif
                        Pccor(i,l) = Xpert(i,l) - Ccor(i)
                     else if(Ncentr.eq.3) then
                        Pccor(i,l) = Xpert(i,10)
                     else
                        Pccor(i,l) = (Xpert(i,3) - Ccor(i))
     .                     + Masse*Xpert(i,10)
                     endif
   30             end do
                  if(Ncentr.le.0) then
                     if(Kp(40).ge.0 .and. (l.eq.3.or.l.eq.10) ) goto 40
                     Rpc2(l) = Rpert2(l)
                     Rpc(l)  = Rpert(l)
                     Rpc3(l) = Rpert3(l)
                     do i = 1,3
                        Pccor3(i,l) = Xpert3(i,l)
                     end do
                     goto 60
                  endif
   40             Rpc2(l) = DOT(Pccor(1,l),Pccor(1,l))
                  Rpc(l)  = SQRT(Rpc2(l))
                  Rpc3(l) = Rpc2(l)*Rpc(l)
                  Rpc5(l) = Rpc2(l)*Rpc3(l)
                  do i = 1,3
                     Pccor3(i,l) = Pccor(i,l)/Rpc3(l)
                  end do
               endif
   60       end do
c
c determine perturbing asteroids or satellites relative to
c central body of this integration
            do l = 1,Numast

c get elliptic orbit position if not on tape
               if(Kpast(l).le.1) call JLIPT(s-Tlpast(l),Ast999(1,l),0,
     .          Yast(1,l),Rastc(l),Rastc2(l),Rastc3(l),dum)
               Kast = Ncnast(l)
               if(Kast.ne.Ncentr) then
                  do i = 1,3
                     if(Kast.gt.0) Yast(i,l) = Yast(i,l)+Xpert(i,Kast)
                     Yast(i,l) = Yast(i,l) - Ccor(i)
                  end do
                  Rastc2(l) = Yast(1,l)**2 + Yast(2,l)**2 + Yast(3,l)**2
                  Rastc(l)  = SQRT(Rastc2(l))
                  Rastc3(l) = Rastc2(l)*Rastc(l)
               endif
               do i = 1,3
                  Yastc3(i,l) = Yast(i,l)/Rastc3(l)
               end do
            end do
c
c move perturbing planet or elliptic asteroid coordinates
c relative to central body into target planet locations
            do kt = 1,Numtar
               if(Itgast(kt).gt.0) then
                  nt = Itgast(kt)
                  do i = 1,3
                     Tccor(i,kt) = Yast(i,nt)
                  end do
                  Rtc2(kt) = Rastc2(nt)
                  Rtc(kt)  = Rastc(nt)
                  Rtc3(kt) = Rastc3(nt)
               else
                  nt = Ntrg(kt)
                  do i = 1,3
                     Tccor(i,kt) = Pccor(i,nt)
                  end do
                  Rtc2(kt) = Rpc2(nt)
                  Rtc(kt)  = Rpc(nt)
                  Rtc3(kt) = Rpc3(nt)
               endif
            end do
c
c determine rotation matrices for target bodies
            do kt = 1,Numtar
               if(Ntzone(kt).gt.1 .or. Nttess(kt).gt.1)
     .          call ROTMAT(kt,Trgrot(1,1,kt),s)
            end do
c
c set-up once per step for solar probe
            if(lcentr.le.0) then
c determine rotation matrix for integrated body if motion
c affected by its second zonal harmonic
               if(Kp(82).ge.0) call ROTMAT(-1, Bodrot, s)

c need to fill xpert for solgr
               if(Nplnt.le.9) then
                  do i = 1, 6
                     Xpert(i, Nplnt) = Y(i, j)
                     if(Kp(100).ge.0) Xpert(i,Nplnt) = Y(i,j) + Ylpt(i)
                  end do
               endif
               call SOLPRB(0)
            endif
            if(Ncentr.gt.0) then
c
c determine rotation matrix for central bodies
               if(Nczone.gt.0 .or. Nctess.gt.0)
     .            call ROTMAT(0, Cntrot, s)
            endif
c
c set-up once per step for earth satellite
            if(lcentr.gt.0) then
               if(lcentr.eq.3) then
                  call ERTORB(0)
c
c set-up once per step for moon satellite
               else if(lcentr.eq.10) then
                  call MONORB(0)
c
c set-up once per step for planet orbiter
c (artificial space probe)
               else if(Nplnt.gt.30) then
                  call PLNORB(0)
c
c set-up once per step for natural planet satellite
               else
                  call PLNSAT(0)
               endif
            endif
         endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           computations done only once for a given iteration of a
c           given step
c
c
c determine satellite or probe relative to central body
         do i = 1, Index1
            if(Kp(100).lt.0) then
               Sbcor(i) = Y(i, j)
            else
               Sbcor(i) = Y(i, j) + Ylpt(i)
            endif
            Bcor(i) = Sbcor(i)
            Pcorsav(i,1)=Sbcor(i)
         end do
         Rsb2 = Sbcor(1)**2 + Sbcor(2)**2 + Sbcor(3)**2
         Rsb  = SQRT(Rsb2)
c setup for central body mascons at this iteration at this time
         dum = MSCFR1(Sbcor,Cntrot)

c stepsize control
         Rsbsng = Rsb

         Rsb3 = Rsb2*Rsb
         Rsb4 = Rsb2*Rsb2
         Rsb5 = Rsb2*Rsb3
         Rb   = Rsb
         Rb2  = Rsb2
         Rb3  = Rsb3
         Rb5  = Rsb5
c
c determine integrated body relative to sun
         if(Ncentr.gt.0) then
            do i = 1, Jndex1
               Bcor(i) = Ccor(i) + Sbcor(i)
            end do
            Rb2 = Bcor(1)**2 + Bcor(2)**2 + Bcor(3)**2
            Rb  = SQRT(Rb2)
            Rb3 = Rb2*Rb
            Rb5 = Rb2*Rb3
c
c write out distance to central body if kkp(81)>0, dist<con1(12)
c and we are not in starting procedure before first tabular
c point   (ntab = -1)
            if(Kkp(81).gt.0 .and. j.ne.2 .and. Ntab.ge.0 .and.
     .             Rsb.le.Con1(12)) then
               rpbkm = Rsb*Aultsc*Ltvel
               if(Line.ge.58) then
                  call NEWPGT(Iout, Npage, 24000)
                  Line = 2
               endif
               Line  = Line + 1
               julkm = (s - 0.5_10) - 2400000._10
               write(Iout, 70) julkm, Ncentr, rpbkm
   70          format(1x, f16.10, ' DISTANCE TO CENTRAL BODY',
     .                i3, ' IS', f14.5, ' KILOMETERS')
            endif
         endif
c
c determine earth-moon barycenter integration quantities
         if(Nplnt.eq.3 .and. Kp(40).ge.0) then
            do i = 1, 3
               Xes(i) = Sbcor(i) - Mass(10)*Xpert(i, 10)
               Xms(i) = Sbcor(i) + Masse*Xpert(i, 10)
            end do
            Res2 = DOT(Xes, Xes)
            Rms2 = DOT(Xms, Xms)
            Res  = SQRT(Res2)
            Rms  = SQRT(Rms2)
            Res3 = Res2*Res
            Rms3 = Rms2*Rms
            Res5 = Res2*Res3
            Rms5 = Rms2*Rms3
            do l = 1, 9
               if(Kp(l+30).ge.0 .and. l.ne.3) then
                  do i = 1, 3
                     Xpes(i, l) = Pccor(i, l) - Xes(i)
                     Xpms(i, l) = Pccor(i, l) - Xms(i)
                  end do
                  Rpes2(l) = DOT(Xpes(1,l), Xpes(1,l))
                  Rpms2(l) = DOT(Xpms(1,l), Xpms(1,l))
                  Rpes(l)  = SQRT(Rpes2(l))
                  Rpms(l)  = SQRT(Rpms2(l))
                  Rpes3(l) = Rpes2(l)*Rpes(l)
                  Rpms3(l) = Rpms2(l)*Rpms(l)
                  Rpes5(l) = Rpes2(l)*Rpes3(l)
                  Rpms5(l) = Rpms2(l)*Rpms3(l)
               endif
            end do
         endif
c
c determine perturbing planets (not sun, central body or
c integrated body) relative to integrated body
         do l=1,10
            if(Kpb(l)) then
               do i = 1, 3
                  Pbcor(i, l) = Pccor(i, l) - Sbcor(i)
               end do
               Rpb2(l) = Pbcor(1,l)**2 + Pbcor(2,l)**2 + Pbcor(3,l)**2
               Rpb(l)  = SQRT(Rpb2(l))
               Rpb3(l) = Rpb2(l)*Rpb(l)
               Rpb5(l) = Rpb2(l)*Rpb3(l)
               if(Nplnt.eq.3 .and. Kp(40).ge.0) then
                  do i=1,3
                     Pbcor3(i,l)=Masse*Xpes(i,l)/Rpes3(l)
     .                           + Mass(10)*Xpms(i,l)/Rpms3(l)
                  end do
               else
                  do i=1,3
                     Pbcor3(i,l)=Pbcor(i,l)/Rpb3(l)
                  end do
               endif
            endif
         end do
c
c determine perturbing asteroids or satellites relative
c to integrated body
         do l = 1,Numast
            do i = 1,3
               Yastb(i,l) = Yast(i,l) - Sbcor(i)
            end do
            Rastb2(l) = Yastb(1,l)**2 + Yastb(2,l)**2 + Yastb(3,l)**2
            Rastb(l)  = SQRT(Rastb2(l))
            Rastb3(l) = Rastb2(l)*Rastb(l)
         end do
c
c move perturbing planet or elliptic asteroid coordinates
c relative to integrated body into target planet locations
         do kt = 1,Numtar
            if(Itgast(kt).gt.0) then
               nt = Itgast(kt)
               do i = 1,3
                  Tbcor(i,kt) = Yastb(i,nt)
               end do
               Rtb2(kt) = Rastb2(nt)
               Rtb(kt)  = Rastb(nt)
               Rtb3(kt) = Rastb3(nt)
            else
               nt = Ntrg(kt)
               do i = 1,3
                  Tbcor(i,kt) = Pbcor(i,nt)
               end do
               Rtb2(kt) = Rpb2(nt)
               Rtb(kt)  = Rpb(nt)
               Rtb3(kt) = Rpb3(nt)
            endif
         end do
c
c write out distance to target body if dist. < con1(12)
c and we are not in starting procedure before first tabular
c point   (ntab = -1)
         do kt = 1,Numtar

c if(ncentr.gt.0) goto 120
            if(j.eq.2 .or. Ntab.lt.0) goto 150
            if(Rtb(kt).le.Con1(12) .and.
     .       (Nplpt(kt).ne.10 .or. Nplnt.ne.3)) then
               rpbkm = Rtb(kt)*Aultsc*Ltvel
               if(Line.ge.58) then
                  call NEWPGT(Iout,Npage,24000)
                  Line = 2
               endif
               Line  = Line + 1
               julkm = (s - 0.5_10) - 2400000._10
               write(Iout,110) julkm,Ntrg(kt),rpbkm
  110          format(1x,f16.10,' DISTANCE TO TARGET BODY',i3,
     .                ' IS ',f14.5,' KILOMETERS')
            endif
         end do
c
c
c determine target body harmonic quantities
  150    do kt = 1,Numtar
            if(Ntzone(kt).gt.1 .or. Nttess(kt).gt.1) then
c formulas similar to those for central body
c with sbcor replaced by  -tbcor
c
c determine latitude
               Tslat(kt) = -(Trgrot(3,1,kt)*Tbcor(1,kt) + Trgrot(3,2,kt)
     .                     *Tbcor(2,kt) +
     .                     Trgrot(3,3,kt)*Tbcor(3,kt))/Rtb(kt)
               Tclat(kt) = SQRT(1._10 - Tslat(kt)**2)
               do i = 1,3
                  Tslat1(i,kt) = Trgrot(3,i,kt) + Tbcor(i,kt)
     .                           *Tslat(kt)/Rtb(kt)
               end do
               Tclatr(kt)  = Tclat(kt)*Rtb(kt)
               rpbh1(kt)   = Rtb(kt)/Trad(kt)
               Rpbh(kt,1) = rpbh1(kt)**2
               if(Ntopt(kt).gt.1) then
                  ntz = Ntopt(kt)
                  do i = 2,ntz
                     Rpbh(kt,i) = Rpbh(kt,i-1)*rpbh1(kt)
                  end do
               endif
c
c determine longitude
               if(Nttess(kt).gt.1) then
                  Tslng(kt,1) = -(Trgrot(2,1,kt)*Tbcor(1,kt) +
     .                            Trgrot(2,2,kt)*Tbcor(2,kt) +
     .                            Trgrot(2,3,kt)*Tbcor(3,kt))
     .                              /Tclatr(kt)
                  Tclng(kt,1) = -(Trgrot(1,1,kt)*Tbcor(1,kt) +
     .                            Trgrot(1,2,kt)*Tbcor(2,kt) +
     .                            Trgrot(1,3,kt)*Tbcor(3,kt))
     .                              /Tclatr(kt)
                  ntt = Nttess(kt)
                  if(ntt.gt.1) then
                     do i = 2,ntt
                        Tslng(kt,i) = Tslng(kt,i-1)*Tclng(kt,1)
     .                        + Tclng(kt,i-1)*Tslng(kt,1)
                        Tclng(kt,i) = Tclng(kt,i-1)*Tclng(kt,1)
     .                        - Tslng(kt,i-1)*Tslng(kt,1)
                     end do
                  endif
                  do i = 1,3
                     Tlng1(i,kt) = (Trgrot(2,i,kt)*Tclng(kt,1) -
     .                              Trgrot(1,i,kt)*Tslng(kt,1))
     .                                 /Tclat(kt)
                  end do
               endif
c
c determine p(n), p'(n), p(n,h), p'(n,h)
c
               if(Jct(78).gt.0) then
                  call LEGNDS(Tslat(kt),Tclat(kt),Ntzone(kt),
     .                        Nttess(kt),Tleg(1,kt),Tleg1(1,kt),
     .                        Tgleg(1,kt),Tgleg1(1,kt))
               else
                  call LEGNDR(Tslat(kt),Tclat(kt),Ntzone(kt),
     .                        Nttess(kt),Tleg(1,kt),Tleg1(1,kt),
     .                        Tgleg(1,kt),Tgleg1(1,kt))
               endif
            endif

         end do
c
c set-up once per iteration for sun centered probe motion
         if(lcentr.le.0) call SOLPRB(-1)
         if(Ncentr.gt.0) then
c
c determine central body harmonic quantities
            if(Nczone.gt.1 .or. Nctess.gt.1) then
c
c determine latitude
               Cslat = (Cntrot(3,1)*Sbcor(1) + Cntrot(3,2)*Sbcor(2)
     .          + Cntrot(3,3)*Sbcor(3))/Rsb
               Cclat = SQRT(1._10 - Cslat**2)
               do i = 1,3
                  Cslat1(i) = Cntrot(3,i) - Sbcor(i)*Cslat/Rsb
               end do
               Cclatr  = Cclat*Rsb
               rsbh1   = Rsb/Crad
               Rsbh(1) = rsbh1**2
               if(Ntopc.gt.1) then
                  do i = 2,Ntopc
                     Rsbh(i) = Rsbh(i - 1)*rsbh1
                  end do
               endif
c
c determine longitude
               if(Nctess.gt.1) then
                  Cslng(1) =(Cntrot(2,1)*Sbcor(1) + Cntrot(2,2)*Sbcor(2)
     .                     + Cntrot(2,3)*Sbcor(3))/Cclatr
                  Cclng(1) =(Cntrot(1,1)*Sbcor(1) + Cntrot(1,2)*Sbcor(2)
     .                     + Cntrot(1,3)*Sbcor(3))/Cclatr
                  do i = 2,Nctess
                     Cslng(i) = Cslng(i-1)*Cclng(1)+Cclng(i-1)*Cslng(1)
                     Cclng(i) = Cclng(i-1)*Cclng(1)-Cslng(i-1)*Cslng(1)
                  end do
                  do i = 1,3
                     Clng1(i) = (Cntrot(2,i)*Cclng(1) -
     .                Cntrot(1,i)*Cslng(1))/Cclat
                  end do
               endif
c
c determine p(n), p'(n), p(n,h), p'(n,h)
c
               if(Jct(78).gt.0) then
                  call LEGNDS(Cslat,Cclat,Nczone,Nctess,Cleg,Cleg1,
     .             Cgleg,Cgleg1)
               else
                  call LEGNDR(Cslat,Cclat,Nczone,Nctess,Cleg,Cleg1,
     .             Cgleg,Cgleg1)
               endif
            endif
         endif
c
c set-up once per iteration for earth satellite motion
         if(lcentr.gt.0) then
            if(lcentr.eq.3) then
               call ERTORB(-1)
c
c set-up once per iteration for lunar probe motion
            else if(lcentr.eq.10) then
               call MONORB(-1)
c
c set-up once per iteration for planetary probe motion
c (artificial space probe)
            else if(Nplnt.gt.30) then
               call PLNORB(-1)
c
c set-up once per iteration for natural satellite motion
            else
               call PLNSAT(-1)
            endif
         endif
         goto 200

      else
         if(k.eq.10) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           computations for partials done only once for a given
c           iteration of a given step
c
c           determine partial derivatives of motion  relative to central
c           body
c     n1    = n/6-1       called kount in sbset
            jgg = j
            if(jgg.gt.2) jgg = 1

c set grav. grad. flags
            Ggfg = Ggfgs(jgg)

c set flag for grav.grad. matrix print (dadx#)
            if(Mpstep.ne.0) then
               if(j.lt.2 .or. Mpstep.ge.0) then
                  if(mod(Npstep,Mpstep).eq.0) Pntflg = .true.
               endif
            endif
            do l1 = 1, Kount
               l = l1*6
               do i = 1, Index2
                  l = l + 1
                  Dsbcor(i, l1) = Y(l, j)
                  if(Icntrl(l1).le.-31 .and. Kp(100).ge.0) then
                     icn = -30 - Icntrl(l1)
                     Dsbcor(i, l1) = Dsbcor(i, l1) + Dylpt(i, icn)
                  endif
               end do
            end do
            do i = 1, Kount
               Rsbdsb(i) = DOT(Sbcor(1), Dsbcor(1,i))
            end do
c
c determine target  body harmonic quantities for partials
            do kt=1,Numtar
               if(Ntopt(kt).gt.0) then
                  if(Jct(78).gt.0) then
                     call LEGNS2(Tslat(kt),Tclat(kt),Ntzone(kt),
     .                Nttess(kt),Tleg(1,kt),Tleg1(1,kt),Tleg2(1,kt),
     .                Tgleg(1,kt),Tgleg1(1,kt),Tgleg2(1,kt))
                  else
                     call LEGND2(Tslat(kt),Tclat(kt),Ntzone(kt),
     .                Nttess(kt),Tleg(1,kt),Tleg1(1,kt),Tleg2(1,kt),
     .                Tgleg(1,kt),Tgleg1(1,kt),Tgleg2(1,kt))
                  endif
               endif

c determine gradient of force w.r.t. target body coords for indirect terms
               if(Nqt(kt).gt.1) then
                  gm=Gamat*Masst(kt)
                  rtc5=Rtc3(kt)*Rtc2(kt)
                  rtb5=Rtb3(kt)*Rtb2(kt)
                  nt=Ntrg(kt)
                  if(Nplnt.eq.3 .and. Kp(40).ge.0) then
                     if(nt.eq.10) then
                        td=Massp1/Res3 - Massp1/Rms3
                        do l=1,9
                           if(Kpb(l)) td=td+Mass1(l)*
     .                      (1._10/Rpes3(l) - 1._10/Rpms3(l))
                        end do
                        td=Gamat*Masse*Mass(10)*td
                     else
                        td=gm*(Masse/Rpes3(nt)+Mass(10)/Rpms3(nt))
     .                   - gm/Rtc3(kt)
                     endif
                  else if(nt.eq.3 .and. Kp(40).ge.0) then
                     td=Gamat*(Mass1(3)/Rpb3(3)+Mass1(10)/Rpb3(10)
     .                -Mass1(3)/Rpc3(3)-Mass1(10)/Rpc3(10))
                  else if(nt.eq.10) then
                     td=Gamat*Mass1(3)*Mass(10)*(-1._10/Rpb3(3)
     .                +1._10/Rpb3(10)+1._10/Rpc3(3)-1._10/Rpc3(10))
                  else
                     td=gm/Rtb3(kt)-gm/Rtc3(kt)
                  endif
                  do i=1,3
                     do i1=1,i
                        if(Nplnt.eq.3 .and. Kp(40).ge.0) then
                           if(nt.eq.10) then
                              tg=Massp1*(Xms(i)*Xms(i1)/Rms5
     .                         - Xes(i)*Xes(i1)/Res5)
                              do l=1,9
                                 if(Kpb(l)) tg=tg+Mass1(l)*(
     .                            Xpms(i,l)*Xpms(i1,l)/Rpms5(l)-
     .                            Xpes(i,l)*Xpes(i1,l)/Rpes5(l))
                              end do
                              tg=3._10*Gamat*Masse*Mass(10)*tg
                           else
                              tg=3._10*gm*(Tccor(i,kt)*Tccor(i1,kt)/rtc5
     .                         - Masse*Xpes(i,nt)*Xpes(i1,nt)/Rpes5(nt)
     .                         - Mass(10)*Xpms(i,nt)*Xpms(i1,nt)/
     .                         Rpms5(nt))
                           endif
                        else if(nt.eq.3 .and. Kp(40).ge.0) then
                           tg=3._10*Gamat*(Mass1(3)*(
     .                      Pccor(i,3)*Pccor(i1,3)/Rpc5(3)
     .                      -Pbcor(i,3)*Pbcor(i1,3)/Rpb5(3))
     .                      +Mass1(10)*(
     .                      Pccor(i,10)*Pccor(i1,10)/Rpc5(10)
     .                      -Pbcor(i,10)*Pbcor(i1,10)/Rpb5(10)))
                        else if(nt.eq.10) then
                           tg=3._10*Gamat*Mass1(3)*Mass(10)*(
     .                      Pbcor(i,3)*Pbcor(i1,3)/Rpb5(3)
     .                      -Pbcor(i,10)*Pbcor(i1,10)/Rpb5(10)
     .                      -Pccor(i,3)*Pccor(i1,3)/Rpc5(3)
     .                      +Pccor(i,10)*Pccor(i1,10)/Rpc5(10))
                        else
                           tg=3._10*gm*(Tccor(i,kt)*Tccor(i1,kt)/rtc5 -
     .                      Tbcor(i,kt)*Tbcor(i1,kt)/rtb5)
                        endif
                        if(i.eq.i1) then
                           Dadp(i,i,kt)=td+tg
                           if(kt.eq.1) then
                              Dadxsav(i,i,1,2,1)=Dadp(i,i,1)
                           else if(Nplnt.eq.3 .and. Kp(40).ge.0 .and.
     .                         nt.eq.10) then
                              Dadxsav(i,i,1,3,1)=Dadp(i,i,kt)
                           endif
                        else
                           Dadp(i,i1,kt)=tg
                           Dadp(i1,i,kt)=tg
                           if(kt.eq.1) then
                              Dadxsav(i,i1,1,2,1)=tg
                              Dadxsav(i1,i,1,2,1)=tg
                           else if(Nplnt.eq.3 .and. Kp(40).ge.0 .and.
     .                         nt.eq.10) then
                              Dadxsav(i,i1,1,3,1)=tg
                              Dadxsav(i1,i,1,3,1)=tg
                           endif
                        endif
                     end do
                  end do
               endif
            end do
c
c
c set-up once per iteration for sun centered probe for partials
            if(lcentr.le.0) call SOLPRB(-2)
            if(Ncentr.gt.0) then
c
c determine central body harmonic quantities for partials
c this should be in sbset    (next 2 lines)
               nczp1 = min0(Nczone,Nczonp+1)
               nctp1 = min0(Nctess,Nctesp+1)

               if(nczp1.gt.1 .or. nctp1.gt.1) then
                  if(Jct(78).gt.0) then
                     call LEGNS2(Cslat, Cclat, nczp1, nctp1, Cleg,
     .                Cleg1, Cleg2, Cgleg, Cgleg1, Cgleg2)
                  else
                     call LEGND2(Cslat, Cclat, nczp1, nctp1, Cleg,
     .                Cleg1, Cleg2, Cgleg, Cgleg1, Cgleg2)
                  endif
               endif
            endif
c
c
c set-up once per iteration for earth satellite for partials
            if(lcentr.gt.0) then
               if(lcentr.eq.3) then
                  call ERTORB(-2)
c
c set-up once per iteration for lunar satellite for partials
               else if(lcentr.eq.10) then
                  call MONORB(-2)
c
c set-up once per iteration for planet satellite for partials
c (artificial space probe)
               else if(Nplnt.gt.30) then
                  call PLNORB(-2)
c
c set-up once per iteration for natural planet satellite
c for partials
               else
                  call PLNSAT(-2)
               endif
            endif
         endif
      endif
c
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           determine class of equation k
  200 Kkk = (k-1)/6
      Kk  = k - Kkk*6 - 3
      if(Kk.gt.0) then
c
c evaluate right side of equation
         call SBFN1(k)

         Ggfg = halt
      else
         Fn(1) = Y(k + 3, j)
      endif
      SBFN = Fn(1)
      if(Kkp(60).gt.0)
     . write(6,99333) k,j,s-2440000.5_10,Fn(1)
99333 format('k,j,s,fn',2i3,f11.4,1p6d22.15)
      return
      end
