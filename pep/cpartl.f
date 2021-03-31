      subroutine CPARTL(itype, ngo, kick)
 
      implicit none

c m.e.ash    july 1969    subroutine cpartl
c interpolate for partials of coordinates from tape and calculate
c partial of observations
c * * * implicit real*8 * * * *
c           itype indicates which body's partials to compute:
c           1-em, 2-mn, 3-pl, 4-sb, 5-sc, 6-pr, 7-er
c           ngo = 1 calculate partial of coordinates only
c           ngo = 2 calculate partial of observation as well as partial
c                   of coordinates
c           ngo = 3 calculate partial of coordinates only, but combine
c                   the various contributions in DEREM, except when kick=3
c           kick= 1 radar or radio tracking observation
c           kick= 2 optical observation
c           kick= 3 transit or occultation observation
c           kick= 4 interferometer or pulsar observation
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 epitch(3), eroll(3), eyaw(3)
      equivalence (Raddum, epitch), (Raddum(4), eroll),
     .            (Raddum(7), eyaw)
      real*10 r3yaw,r2roll,r2r3,r1r,r1rr,rpitch,ryaw
      equivalence (Angdum, r3yaw), (Angdum(2), r2roll),
     .            (Angdum(3), r2r3), (Angdum(4), r1r),
     .            (Angdum(5), r1rr)
      equivalence (Angdum(6), rpitch), (Angdum(7), ryaw)
      real*10 der3(6), rsc1
      equivalence (Angdum, der3), (Angdum(7), rsc1)
      real*10 rsc2, xsd(3, 2), der4(3)
      equivalence(angdum(8),rsc2),(xsd(1,1),raddum(5)),(der4(1),der3(4))
      real*10 x(6, 3)
      equivalence (Xm(1,1),x(1,1))
      real*10 xo(3), rr1, r, fake, cake
      equivalence (Angdum, xo), (Angdum(4), rr1), (Angdum(5), r),
     .            (Angdum(6), fake), (Angdum(7), cake)
      real*10 xop(3), xnp(3), xpta(3), xptd(3)
      equivalence (Raddum, xop, xpta),
     .            (Raddum(4), xptd),
     .            (Raddum(8), xnp)
      real*10 dvect(6)
      equivalence (Angdum, dvect)
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'ltrapobs.inc'
      integer*2 kind
      equivalence (Numpar, kind)
      include 'mnsprt.inc'
      real*10 dxdisg(6, 2), dxdrho(6, 2), dxdtau(6, 2)
      equivalence (Dxdpsi(1,1), dxdisg(1,1)), (Dxdth(1,1), dxdrho(1,1)),
     .            (Dxdphi(1,1), dxdtau(1,1))
      include 'number.inc'
      include 'obscrd.inc'
      real*10 xfactr, freqtr(2), ct
      equivalence (Dstf(6), xfactr), (Dstf(7), freqtr), (Dstf(9), ct)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'partcm.inc'
      include 'pqind.inc'
      include 'radcrd.inc'
      include 'rotdta.inc'
      include 'rtrdvl.inc'
      include 'sbdta.inc'
      include 'sbdtavtp.inc'
      include 'scdta.inc'
      include 'scdtavtp.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'tabval.inc'
      include 'tapdte.inc'
      include 'tapdtm.inc'
      include 'tapdtp.inc'
      include 'trnocc.inc'
      include 'trpcom.inc'
c
c local
      real*10 drmr1(6)
      real*10 quan2(2), quan3(3), quan4(2), yy(7,3)
      real*10 deyaw(3), deroll(3), deptch(3)
      real*10 vctor1(3), vctor2(3), vctor3(3), vc1vc2
      real*10 quan5(2)
      real*10 deriv1(296), ddfdl1, ddfdl2, dddr1, deriv2(296)
      real*10 dtmda,dtprds,quan23
      integer i,i1,iprdbg,itype,j,jj,jndex,jvl,k,kick,kk,
     .        kspra,lxpdbg,m,n,nc2,ngo,nobj

c indices for arrays derem,dermn,derpl,dersb,dersc
      integer*2 iderst(5)/1, 25, 37, 43, 57/
      character*2 name(6)/'EM','MN','PL','SB','SC','PR'/
 
      logical*4 doboth

c extermal functions
      real*10 DOT, HERMTF, TERPF, TERP14
c
c setup common to all bodies
      iprdbg = mod(Jct(6)/256, 2)
      lxpdbg = Index/3
      if(itype.gt.0 .and. itype.le.8) then
         Ivze(itype) = 1
         Ntab1 = Nb1(itype)
         Ntab2 = 3
      else
         call SUICID('INVALID ITYPE, STOP CPARTL  ',7)
      endif
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c  arrays 'ivze' and 'ivus' indicate status of the partials arrays
c  derem,dermn,derpl,dersb,dersc -- ivze(i)= -1,0,1 if the
c  corresponding der** contains garbage, zeroes, or data.
c  ivus(i)= 0,1 if the corresponding der** is ignored or included
c  in cpartc calculations.
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      if(itype.eq.1) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c partial derivatives of earth-moon barycenter position and velocity
         doboth = (Mouse.eq.2)
         do j = 1, Index
            call YCOFT(yy, Embary(j,Kembry,1), i_mxplprt+1, 1, 1)
            derem(j, 1) = TERPF(Pem(1,1), yy(1,2))
            if(doboth) derem(j,2) = TERPF(Pem(1,2), yy(1,Ntab1))
            if(Ncodf.gt.20) then
               Derpr(j,1,2) = derem(j,1)
               if(doboth) Derpr(j,2,2) = derem(j,2)
            endif
         end do
         call PARTVL(derem,2,kick)
         if(iprdbg.ne.0) then
            do j = 1, Mouse
               if(Line+lxpdbg.gt.58) call OBSPAG
               write(Iout,210) j,0,kind,name(itype),Index,
     .          (Derem(i,j),i=1,Index)
               Line = Line + lxpdbg
            end do
         endif
         if(ngo.eq.1 .or. ngo.eq.3) return
 
      else if(itype.eq.2) then
c
c partial derivatives of moon position and velocity
         doboth = (Mouse.eq.2 .and. Klan.ne.17)
         do j = 1, Index
            call YCOFT(yy, Moon(j,Kmon,1), i_mxplprt+1, 2, 1)
            Dermn(j, 1) = TERP14(Pm(1,1), yy(1,2))
            if(doboth) Dermn(j, 2) = TERP14(Pm(1,2), yy(1,Ntab1))
         end do
         n = 2
         if(Klan.eq.17) n = 1
         call PARTVL(Dermn, n, kick)
         if(iprdbg.ne.0) then
            do j=1,n
               if(Line+lxpdbg.gt.58) call OBSPAG
               write(Iout,210) j,0,kind,name(itype),Index,
     .          (Dermn(i,j),i=1,Index)
               Line = Line + lxpdbg
            end do
         endif
         if(ngo.eq.1) return
         if(kick.ne.3) then
            i1 = 1
            do i = 1, Mouse
               do j = 1, Index
                  if(Ivze(1).le.0) then
                     derem(j,i) = -Dermn(j,i1)*Mnstf
                  else
                     derem(j,i) = derem(j,i)-Dermn(j,i1)*Mnstf
                  endif
                  if(Ncodf.gt.20) Derpr(j,i,2) = derem(j,i)
               end do
               if(Klan.ne.17) i1 = 2
            end do
            Ivze(1) = 1
         endif
         if(ngo.eq.3) return

      else if(itype.eq.3) then
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c partial derivatives of planet position and velocity
         do j = 1, Index
            call YCOFT(yy, Planet(j,Kplnt,1), i_mxplprt+1, 3, 1)
            Derpl(j) = TERPF(Ppl(1), yy(1,2))
         end do
         call PARTVL(Derpl,1,kick)
         if(iprdbg.ne.0) then
            if(Line+lxpdbg.gt.58) call OBSPAG
            write(Iout,210) 1,0,kind,name(itype),Index,
     .       (Derpl(i),i=1,Index)
            Line = Line + lxpdbg
         endif
         if(ngo.eq.1) return
         if(kick.ne.3) then
            do i = 1, Mouse
               do j = 1, Index
                  if(Ivze(1).le.0) then
                     derem(j,i) = -Derpl(j)
                  else
                     derem(j,i) = derem(j,i)-Derpl(j)
                  endif
                  if(Ncodf.gt.20 .and. (Nplnt2.eq.Nplnt0 .or.
     .             (Klanb.gt.0 .and. (Nplnt2.eq.Ncp0.or.Ncs1.eq.Ncp0))))
     .             Derpr(j,i,2) = derem(j,i)
               end do
            end do
            Ivze(1) = 1
         endif
         if(ngo.eq.3) return
 
      else if(itype.eq.4) then
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c partial derivitives of probe position and velocity
         if(Ksb(88).gt.0) then
            do j = 1, Index
               call YCOFT(yy, Satprb(j,Ksprb,1), i_mxeqn, 4, 1)
               Dersb(j) = TERPF(Psb(1), yy(1,2))
            end do
         else if(kisb(1).ge.-1) then
            kspra = Ksprb - 1
            Nc    = Nbtrp(4)
            j     = 0
            do i = 3, Index, 3
               jvl = i/3
               do m = 1, 3
                  j = j + 1
                  call YHERMT(dsprb(1,1,m,kspra), yy, jvl, 1, hc2, Nc)
                  Dersb(j) = HERMTF(yy, jvl, 1, Psb)
               end do
            end do
         endif
         call PARTVL(Dersb, 1, kick)
         if(iprdbg.ne.0) then
            if(Line+lxpdbg.gt.58) call OBSPAG
            write(Iout,210) 1,0,kind,name(itype),Index,
     .       (Dersb(i),i=1,Index)
            Line = Line + lxpdbg
         endif
         if(ngo.eq.1) return
         do j = 1, Index
            Dersb(j) = Dersb(j)*Cmfct
         end do
         if(kick.eq.3) then
            do j = 1, Index
               if(Ivze(3).le.0) then
                  Derpl(j) = Dersb(j)
               else
                  Derpl(j) = Derpl(j)+Dersb(j)
               endif
            end do
            Ivze(3) = 1
         else
            if(Ncodf.eq.13 .or. Ncodf.eq.14) then
               do j = 1, Index
                  derem(j, 1) = 0._10
                  derem(j, 2) = -Dersb(j)
c it is assumed that index=3 (der4(1) equivalenced to der3(4))
                  der3(j)     = -Dersb(j)
                  if(Ncodf.eq.14) then
                     derem(j, 2) = 0._10
                     der4(j)     = -Dersb(j)
                  endif
               end do
               Ivze(1) = 1
               go to 200
            else
               do i = 1, Mouse
                  do j = 1, Index
                     if(Ivze(1).le.0) then
                        derem(j,i) = -Dersb(j)
                     else
                        derem(j,i) = derem(j,i)-Dersb(j)
                     endif
                  end do
               end do
               Ivze(1) = 1
            endif
         endif
         if(ngo.eq.3) return
 
      else if(itype.eq.5) then
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c partial derivatives of observing body position and velocity
         doboth = (Mouse.eq.2 .and. Ncodf.ne.13)
         jndex  = Index
         if(Ncodf.eq.4) jndex = 6
         if(Ksc(88).le.0) then
            kspra = Ksprc - 1
            Nc    = Nb1(5)
            nc2   = Nbtrp(5)
            j     = 0
            do i = 3, jndex, 3
               jvl = i/3
               do m = 1, 3
                  j = j + 1
                  call YHERMT(dsprc(1,1,m,kspra), yy, jvl, 1, gc2, Nc)
                  Dersc(j, 1) = HERMTF(yy, jvl, 1, Psc)
                  if(doboth) then
                     if(nc2.ne.Nc) call
     .                  YHERMT(dsprc(1,1,m,kspra), yy, jvl, 1, gc2, nc2)
                     Dersc(j, 2) = HERMTF(yy, jvl, 1, Psc(1,2))
                  endif
               end do
            end do
         else
            do j = 1, jndex
               call YCOFT(yy, Satprc(j,Ksprc,1), i_mxeqn, 5, 1)
               Dersc(j, 1) = TERPF(Psc, yy(1,2))
               if(doboth) Dersc(j, 2) = TERPF(Psc(1,2), yy(1,Ntab1))
            end do
         endif
         n = 1
         if(Ksite(1).gt.0 .and. Ksite(1).ne.3) n = 2
         call PARTVL(Dersc, n, kick)
         if(iprdbg.ne.0) then
            do j=1,n
               if(Line+lxpdbg.gt.58) call OBSPAG
               write(Iout,210) j,0,kind,name(itype),Index,
     .          (Dersc(i,j),i=1,Index)
               Line = Line + lxpdbg
            end do
         endif
         if(ngo.eq.1) return
         if(kick.ne.3) then
            Ivze(1) = 1
            if(Ncodf.ge.13 .and. Ncodf.le.14) then
               do j = 1, Index
                  derem(j, 1) = -Dersc(j, 1)
                  derem(j, 2) = 0._10
c it is assumed that index=3 (der4(1) equivalenced to der3(4))
                  der3(j)     = Dersc(j, 1)
                  if(Ncodf.eq.14) then
                     derem(j, 2) = -Dersc(j, 2)
                     der4(j)     = Dersc(j, 2)
                  endif
               end do
            else if(Nplnt2.gt.0) then
               do i = 1, Mouse
                  do j = 1, jndex
                     Derpr(j, i, 1) = 0._10
                     Derpr(j, i, 2) = -Dersc(j, i)
                  end do
               end do
            else
               do i = 1, Mouse
                  do j = 1, jndex
                     derem(j, i) = Dersc(j, i)
                  end do
               end do
            endif
         endif
         if(ngo.eq.3) return
 
      else if(itype.eq.6) then
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c partial derivatives of planet or moon rotation
         do j = 1, Index
            call YCOFT(yy, Plnmon(j,Kprt,1), i_mxplprt+1, 6, 1)
            drmr1(j) = TERPF(Ppr(1), yy(1,2))
         end do
 
c remove constant part of dx/dth for eckhardt theory
         if(Kpr(100).eq.1) drmr1(2) = drmr1(2) - Mrcom(Kprt-1)
 
c partials on tape are angles - must get cartesian coordinates
         n = 1
         if(Nspot2.gt.0) n = 2
         do k = 1, n
            if(Kpr(100).ge.0) then
               do j = 1, 3
                  Dermr(j) = drmr1(1)*dxdtau(j,k) + drmr1(2)
     .                       *dxdrho(j,k) + drmr1(3)*dxdisg(j,k)
     .                       + Mrcom(Kprt - 1)*Dxdi(j,k)
               end do
 
c analytic theory needs parital wrt i - numerical integ. has dx/di=0
               if(Index.gt.3) then
                  do j = 4,6
                     jj = j - 3
                     Dermr(j) = drmr1(1)*dxdtau(j,k) + drmr1(2)
     .                          *dxdrho(j,k) + drmr1(3)*dxdisg(j,k)
     .                          + drmr1(4)*dxdtau(jj,k) + drmr1(5)
     .                          *dxdrho(jj,k) + drmr1(6)*dxdisg(jj,k)
                  end do
               endif
            else
 
c euler angles stored in different order
               do j = 1,3
                  Dermr(j) = drmr1(1)*Dxdpsi(j,k) + drmr1(2)
     .                       *Dxdth(j,k) + drmr1(3)*Dxdphi(j,k)
               end do
 
c euler velocities not yet programmed
               if(Index.gt.3) call SUICID(
     .          'EULER ANGLE VELOCITY PARTIALS NOT PROGRAMMED',11)
            endif
            call PARTVL(Dermr,1,kick)
            if(iprdbg.ne.0) then
               if(Line+lxpdbg.gt.58) call OBSPAG
               write(Iout,210) 1,0,kind,'MR',Index,
     .          (Dermr(i),i=1,Index)
               Line = Line + lxpdbg
            endif
c no observables with kick=3 use moon or planet rotation, else there
c would have to be an exception here
            do i = 1, Mouse
               do j = 1, Index
                  if(Ivze(1).le.0) then
                     Derpr(j,i,k) = -Dermr(j)
                  else
                     Derpr(j,i,k) = Derpr(j,i,k)-Dermr(j)
                  endif
               end do
            end do
         end do
         Ivze(1)=1
         if(ngo.eq.1 .or. ngo.eq.3) return
      else
         call SUICID('ILLEGAL ITYPE IN CPARTL ', 6)
      endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=1000
c
c           branch to observable partial
      entry CPARTC(kick)
c
c debug printout (see also label=3000)
  200 iprdbg = mod(Jct(6)/256, 2)
      if(iprdbg.ne.0) then
         lxpdbg = Index/3
         nobj   = 1
         if(Ncodf.gt.20) nobj = 2
         do k = 1, nobj
            do j = 1, Mouse
               if(Line+lxpdbg.gt.58) call OBSPAG
               write(Iout, 210) j, k, kind, name(1), Index,
     .                          (Derpr(i,j,k), i = 1, Index)
  210          format(' CPARTL: SITE', i1, ', OBJ', i1, ', DERIV(', i3,
     .                ')', t43, 'DER', a2, '(1-', i1, ')=',
     .                (t54,1p,3D23.15))
               Line = Line + lxpdbg
            end do
         end do
      endif
 
      if(kick.eq.2) then
c*  start=2000
c
c determine if angular measurement is made on body other
c than the earth
         if(Ksite(1).ne.3 .and. Ksite(1).gt.0) then
c
c
c*  start=2300
c           partial derivatives of pitch and roll angle observations
c           made in spacecraft reference system
            if(Ncodf.le.4) then
               do i = 1, 3
                  quan3(i) = 0._10
               end do
               if(Ivze(5).gt.0) then
 
c partials of unit yaw,pitch,roll vectors
                  do i = 1, 3
                     deyaw(i) = (Xscsun(i,1)*DOT(Xscsun(1,1),Dersc(1,1))
     .                          /ryaw**2 - Dersc(i,1))/ryaw
                  end do
                  call CROSS(Dersc(1,1), Xscsun(4,1), vctor1)
                  call CROSS(Xscsun(1,1), Dersc(4,1), vctor2)
                  do i = 1, 3
                     vctor2(i) = vctor2(i)/Secday
                     vctor3(i) = vctor1(i) + vctor2(i)
                  end do
                  vc1vc2 = DOT(epitch, vctor3(1))
                  do i = 1, 3
                     deptch(i) = (vctor1(i) + vctor2(i) - epitch(i)*
     .                           vc1vc2)/rpitch
                  end do
                  call CROSS(deyaw(1), epitch, vctor1)
                  call CROSS(eyaw(1), deptch(1), vctor2)
                  do i = 1, 3
                     deroll(i) = vctor1(i) + vctor2(i)
                  end do
                  quan3(1) = DOT(Xsitep(1,1), deptch(1))
                  quan3(2) = DOT(Xsitep(1,1), deroll(1))
                  quan3(3) = DOT(Xsitep(1,1), deyaw(1))
               endif
               if(Nice.le.0) Deriv(kind, 1) =
     .                   -(r3yaw*(DOT(derem,eroll)+quan3(2)) -
     .                     r2roll*(DOT(derem,eyaw)+quan3(3)))/r2r3
               if(Nice.ge.0) Deriv(kind, 2) = -(DOT(derem,epitch) +
     .                   quan3(1) + r1r*DOT(derem,Xsitp0))/r1rr
c
c*  start=2500
c partial derivatives of spacecraft look angles relative
c to star background
            else if(Nice.gt.0) then
               Deriv(kind, 2) = (-fake*DOT(Xsitep(1,1),derem(1,1)) -
     .                          derem(3,1))/cake
            else
               Deriv(kind, 1) = -(Xsitep(2,1)*derem(1,1) - Xsitep(1,1)
     .                          *derem(2,1))/rr1
               if(Nice.ge.0) Deriv(kind, 2)
     .             = (-fake*DOT(Xsitep(1,1),derem(1,1)) - derem(3,1))
     .             /cake
            endif
            go to 1400
         else if(Ncodf.gt.20) then
c
c differenced observable
c*  start=2600
            if(Nice.le.0) Deriv(kind, 1) = DOT(xpta, derem)/cake
            if(Nice.ge.0) Deriv(kind, 2) = DOT(xptd, derem)/rr1
            go to 1400
         else
c
c partial derivatives of right ascension and/or declination
            if(Ncodf.eq.20) goto 400
            if(Ncodf.lt.5) then
               if(Klanb.gt.0 .and.(Ncp0.eq.3 .or. Ncp0.eq.10)) goto 400
               call CORCHN(Zp, derem(1,1))
            else if(Jct(39).lt.0) then
c
c*  start=2100
c partial derivatives of right ascension and declination rates
               if(Nice.le. 0) Deriv(kind,1)= -((-derem(1,1)*Xsitep(5,1)+
     .          derem(2,1)*Xsitep(4,1) + (xo(1)*derem(5,1) -
     .          xo(2)*derem(4,1)) /8.64E4_10)/rr1 +
     .          2._10*(xo(1)*Xsitep(5,1)-xo(2)*Xsitep(4,1))/rr1**2 *
     .          (xo(1)*derem(1,1)+xo(2)*derem(2,1))*Convhs )
 
               if(Nice.ge. 0) Deriv(kind,2) = -(derem(6,1)/8.64E4_10 -
     .          ((derem(3,1)-2._10*fake*DOT(xo(1),derem(1,1)) )
     .          *Angdum(9) - xo(3)*(DOT(derem(1,1),Xsitep(4,1)) -
     .          DOT(xo(1),derem(4,1)) /8.64E4_10))/r**2 +
     .          (Xsitep(6,1)+fake*Angdum(9)) *
     .          (derem(1,1)*xo(1)+derem(2,1)*xo(2)) /
     .          Angdum(10)**2) /Angdum(10)/Convds
               goto 1400
            else
               do i = 1, 3
                  Zp(i) = derem(i, 1)
               end do
               go to 300
            endif
         endif
      else if(kick.eq.3) then
c*  start=3000
c
c           partial derivatives of transit or occultation observation
c        input partials of 'logical coordinates' depends on ibtrn
c        first is the observed body relative to logical center
c                  - - - assigned to dermn (1,2,3) or derpl (rest)
c        second is observer relative to logical center
c                  - - - assigned to derem
c        for mutual mid-time events, second is dersc - second body rel.
c                       to logical center
c        zero out garbage
         do i = 1, 5
            i1 = iderst(i) - 1
            if(Ivze(i).lt.0) then
               Ivze(i) = 0
               do j = 1, Index
                  derem(i1+j,1) = 0._10
               end do
            endif
         end do
c get logical position partials
c correct partial for phase correction in mid-time event
         if(Ibtrn.ge.8) call PHSCRP
         if(Ibtrn.eq.1 .or. Ibtrn.eq.4 .or. Ibtrn.eq.9) go to 600
         if(Ibtrn.eq.3) go to 500
         if(Ibtrn.eq.5 .or. Ibtrn.eq.7) go to 1400
         if(Ibtrn.eq.6) then
            if(Ivze(5).ne.1) go to 600
            do i = 1, Index
               Derpl(i) = Derpl(i) - Dersc(i, 1)
            end do
            Ivze(3) = 1
            go to 500
         else if(Ibtrn.eq.8) then
            if(Ivze(1).gt.0) then
               do i = 1, Index
                  Derpl(i)    = Derpl(i) - derem(i, 1)
                  Dersc(i, 1) = Dersc(i, 1) - derem(i, 1)
               end do
               Ivze(3) = 1
               Ivze(5) = 1
            endif
            go to 600
         else
            do i = 1, 3
               Dermn(i, 1) = Derpl(i) - derem(i, 1)
            end do
            go to 600
         endif
      else if(kick.eq.4) then
c*  start=4000
c
c partial derivatives of interferometer observables
         Quan1(2) = 0._10
         quan2(2) = 0._10
         quan3(2) = 0._10
         quan4(2) = 0._10
 
         if(Nplnt0.lt.0) then
c*  start=4200
c first object is a star
            if(Nice.gt.0) go to 1000
            if(Lprspt.eq.1) then
 
c partial wrt star coordinate
               Deriv(kind, 1) = DOT(dvect(1), Derpr(1,1,1))*Aultsc
               if(Nice.ge.0) go to 1000
               go to 1100
            else
               do i = 1, 2
                  quan5(i) = -DOT(Xspcd(1,1), Derpr(1,i,1))
               end do
               Deriv(kind, 1) = (quan5(1) - quan5(2))*Aultsc
               if(Nice.ge.0) go to 1000
               go to 1100
            endif
         else
            do i = 1, Mouse
               quan2(i) = DOT(Xsitp0(1,i), Derpr(1,i,1))
            end do
            if(Mouse.eq.1) quan2(2) = 0._10
            if(Nice.gt.0) go to 900
            ddfdl1 = (quan2(1) - quan2(2))*Aultsc
            Deriv(kind, 1) = ddfdl1
            if(Nice.ge.0) go to 900
            if(nintrf.lt.0) go to 1100
 
c for nice.lt.0 (difnct) deriv1 is partial at start of series
            if(Nk1.lt.0) then
 
c store temporarily ddfdl1 for beginning of first interval
               deriv1(kind) = ddfdl1
               if(nintrf.eq.1) deriv1(kind) = 0._10
               go to 1100
            else if(Nk1.eq.0) then
               deriv1(kind) = deriv1(kind)*freqtr(1)
            endif
            Deriv(kind, 1) = ddfdl1*freqtr(1)
            if(Ncodf.gt.20) go to 1200
            Deriv(kind, 1) = Deriv(kind, 1) - deriv1(kind)
            Deriv(kind, 1) = Deriv(kind, 1)*xfactr
            go to 1400
         endif
      else
c
c*  start=1100
c partial derivatives of time delay and/or doppler shift
         quan2(2) = 0._10
         quan3(2) = 0._10
         quan4(2) = 0._10
         do i = 1, Mouse
            quan2(i) = DOT(Xsitp0(1,i), derem(1,i))
         end do
         if(Ncodf.ge.13 .and. Ncodf.le.14) then
c
c*  start=1200
c partial derivatives of 2 spacecraft time delay and/or
c doppler shift
            quan23 = DOT(xsd(1,1), der3(1))/rsc1
            if(Nice.ge.0) call SUICID(
     .' CANNOT CALCULATE PARTIAL FOR 2 SPACECRAFT DOPPLER, STOP IN CPART
     .L  ', 17)
            Deriv(kind, 1) = (quan2(1) + quan2(2) + quan23)*Aultsc
            if(Ncodf.eq.14) Deriv(kind, 1) = Deriv(kind, 1)
     .          + DOT(xsd(1,2), der4(1))/rsc2*Aultsc
            do j = 1, 6
               der3(j) = 0._10
            end do
         else
c
c correct delay partials for light time
            if(Jct(67).gt.0) then
               do i = 1, 3
                  derem(i, 1) = derem(i, 1) + V2fac(i)*quan2(1)
               end do
               quan2(1) = DOT(Xsitp0, derem)
               if(Mouse.ne.1) then
                  do i = 1, 3
                     derem(i, 2) = derem(i, 2) + quan2(1)*Vxfac(i)
     .                             - quan2(2)*V3fac(i)
                  end do
                  quan2(2) = DOT(Xsitp0(1,2), derem(1,2))
               endif
            endif
            if(Nice.le.0) then
               Deriv(kind, 1) = (quan2(1) + quan2(2))*Aultsc
               if(Ncodf.eq.18) Deriv(kind, 1) = Deriv(kind, 1)*Dpphdt
               if(Nice.lt.0) go to 1400
            endif
            do i = 1, Mouse
               quan3(i) = DOT(Dxsit0(1,i), derem(1,i))*Dopfct(i)
               quan4(i) = DOT(Xsitp0(1,i), derem(4,i))*Dopfct(i)
            end do
            Deriv(kind, 2) = -Freq*((quan3(1)+quan3(2))*Aultsc + (quan4(
     .                       1)+quan4(2))*Aultvl)
         endif
         go to 1400
      endif
      entry CPORTO
      dtmda = (xo(1)*Zp(2) - xo(2)*Zp(1))/Angdum(8)
      do i = 1, 3
         Zp(i) = Zp(i) - Xsitep(i+3,1)*dtmda
      end do
      entry CPARTO
  300 if(Nice.le.0) Deriv(kind, 1) = (xo(2)*Zp(1) - xo(1)*Zp(2))/rr1
      if(Nice.ge.0) Deriv(kind, 2) = (fake*DOT(xo,Zp) - Zp(3))/cake
      go to 1400
c
c*  start=2200
c partial derivatives of azimuth and/or elevation
  400 if(Nice.le.0) Deriv(kind, 1) = DOT(xop, derem)/cake
      if(Nice.ge.0) Deriv(kind, 2) = DOT(xnp, derem)/rr1
      go to 1400
  500 if(Ivze(5).eq.1) then
         do i = 1, 3
            derem(i, 1) = derem(i, 1) - Dersc(i, 1)
         end do
         Ivze(1) = 1
      endif
c*  start=3300
c now calculate partials
  600 kk = 1
      if(Nice.gt.0) go to 800
  700 if(Ibtrn.eq.4 .or. Ibtrn.eq.6) then
 
c transit of planet (limbs) or planet mutual limb event
         Deriv(kind, kk) = -(DOT(Vect(1,kk),Derpl) + DOT(Xpr(1,kk),derem
     .                     ))/Dfdt3(kk)
      else if(Ibtrn.eq.5 .or. Ibtrn.eq.7) then
         go to 1400
      else if(Ibtrn.eq.8 .or. Ibtrn.eq.9) then
 
c mutual event mid-time
         Deriv(kind, 1) = -(DOT(Vect,Derpl) + DOT(Xpr,Derpl(4)) + DOT(
     .                    Vect(1,2),Dersc) + DOT(Xpr(1,2),Dersc(4,1)))
     .                    /Dfdt3(1)
 
c partial of separation
         if(Nice.ge.0) Deriv(kind, 2)
     .       = (Rmrx(2)*(Ftct(1)*DOT(Reh(1,1),Derpl)+Ftct(2)
     .       *DOT(Reh(1,2),Dersc)+Rmxc(2)*Deriv(kind,1))
     .       - (DOT(Xpr(1,1),Derpl)+DOT(Xpr(1,2),Dersc))/Rmxc(1))
     .       /Rmrx(1)
         go to 1400
      else
 
c occultation of star or planet
         if(Re(kk).ne.0._10) dtprds = DOT(Xpr(1,kk), derem)*Rmrx(kk)
         Deriv(kind, kk) = (DOT(Vect(1,kk),Dermn) - Re(kk)*dtprds)
     .                     /Dfdt3(kk)
      endif
  800 if(kk.eq.2) go to 1400
      if(Nice.lt.0) go to 1400
      kk = 2
      go to 700
  900 do i = 1, Mouse
         Quan1(i) = DOT(Xsitep(4,i), Xsitp0(1,i))
         quan2(i) = Quan1(i)*quan2(i)/Rsitp(i)
         quan3(i) = DOT(Xsitep(4,i), Derpr(1,i,1))
         quan4(i) = DOT(Xsitp0(1,i), Derpr(4,i,1))
      end do
      dddr1 = (quan3(1) - quan3(2) + quan2(2) - quan2(1))
     .        *Aultsc + (quan4(1) - quan4(2))*Aultvl
      Deriv(kind, 2) = dddr1
      goto 1100
c delay rate partial
 1000 if(Lprspt.eq.1) then
         Deriv(kind, 2) = DOT(dvect(4), Derpr(1,1,1))*Aultsc
      else
         do i = 1, 2
            quan3(i) = -DOT(Xspcd(1,1), Derpr(4,i,1))
         end do
         Deriv(kind, 2) = (quan3(1) - quan3(2))*Aultvl
      endif
c
c*  start=4400
c a second object is observed - difference partials
 1100 if(Ncodf.lt.30) go to 1400
 1200 if(Nplnt2.lt.0) then
c*  start=4600
c second object is a star
         if(Nice.le.0) then
            if(Lprspt.eq.1) then
               Deriv(kind, 1) = Deriv(kind, 1)
     .                          - DOT(Dvect2(1), Derpr(1,1,2))*Aultsc
            else
               do i = 1, 2
                  quan5(i) = -DOT(Xspcd(1,2), Derpr(1,i,2))
               end do
               Deriv(kind, 1) = Deriv(kind, 1) - (quan5(1) - quan5(2))
     .                          *Aultsc
            endif
         endif
         if(Nice.ge.0) then
 
c delay rate partial for second object
            if(Lprspt.eq.1) then
               Deriv(kind, 2) = Deriv(kind, 2)
     .                          - DOT(Dvect2(4), Derpr(1,1,2))*Aultsc
            else
               do i = 1, 2
                  quan3(i) = -DOT(Xspcd(1,2), Derpr(4,i,2))
               end do
               Deriv(kind, 2) = Deriv(kind, 2) - (quan3(1) - quan3(2))
     .                          *Aultvl
            endif
         endif
      else
         do i = 1, 2
            quan2(i) = DOT(Ysitp0(1,i), Derpr(1,i,2))
         end do
         if(Nice.le.0) then
            ddfdl2 = (quan2(1) - quan2(2))*Aultsc
            if(nintrf.lt.0) Deriv(kind, 1) = ddfdl1 - ddfdl2
         endif
         if(Nice.ge.0) then
            do i = 1, 2
               Quan1(i) = DOT(Ysitep(4,i), Ysitp0(1,i))
               quan2(i) = Quan1(i)*quan2(i)/Rsitp2(i)
               quan3(i) = DOT(Ysitep(4,i), Derpr(1,i,2))
               quan4(i) = DOT(Ysitp0(1,i), Derpr(4,i,2))
            end do
            Deriv(kind, 2) = dddr1 - (quan3(1) - quan3(2) + quan2(2)
     .            - quan2(1))*Aultsc - (quan4(1) - quan4(2))*Aultvl
         else
            if(nintrf.ge.0) then
               if(Nk1.lt.0) then
 
c store temporarily ddfdl2 for beginning of first interval
                  deriv2(kind) = ddfdl2
                  go to 1300
               else if(Nk1.eq.0) then
                  deriv1(kind) = deriv1(kind) - deriv2(kind)*freqtr(2)
                  if(nintrf.eq.1) deriv1(kind) = 0._10
               endif
               Deriv(kind, 1) = Deriv(kind, 1) - ddfdl2*freqtr(2)
     .                          - deriv1(kind)
               Deriv(kind, 1) = Deriv(kind, 1)*xfactr
            endif
         endif
      endif
 1300 do k = 1, 2
         do j = 1, 2
            do i = 1, Index
               Derpr(i, j, k) = 0._10
            end do
         end do
      end do
      Ivze(1) = 0
 
      Lprspt = 0
c*  start=9000
c
c note that used data becomes garbage
 1400 do i = 1, 6
         if(Ivze(i).gt.0) Ivze(i) = -1
      end do
      return
      end
