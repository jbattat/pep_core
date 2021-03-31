      subroutine SBOUT

      implicit none
c
c becker/ash   oct. 1968   subroutine sbout
c output routine for satellite-probe numerical integration
c
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'adams.inc'
      include 'dumdum.inc'
      real*10 sb1(6,i_mxeqn,5)
      equivalence (sb1,A)
      include 'ellips.inc'
      include 'fcntrl.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'nmstrt.inc'
      include 'output.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'rstart.inc'
      include 'sbstuf.inc'
      include 'sscon.inc'
      include 'stint.inc'
c data from satellite - probe numerical integration to be
c written on tape
      common/WRKCOM/ Sb(6,i_mxeqn,5),Tabpt(3,8)
      real*10 Sb,Tabpt
      real*10 dtabpt(3,8,i_mxeqn-1)
      equivalence (dtabpt,Sb(1,1,2))
      include 'zeroes.inc'

c local variables
      real*10 djedz,fsb,fsb1,fsb1st,hc1,hc2,rdv,rdvsv
      integer   i,ifmt,ii,ij,intl,iorder,iordr1,ivl,j,je,
     .          jf,jorder,jordr1,jplnts,jx,k,kk,kk95,kkp2
      integer   kp98,kx,l5sav,lx,mgo,mtabx,n1sav,n2sav,
     .          nexz,ngo,npcont,numprt,numpt1,nwrite
      real*10 juldat,frect,juldas

      character*19 sbfmt/'(F13.6,xxP,6F20.16)'/
      character*8 sbfmth(17)/'( '' JED-', 'xxxxx00''', ',8X,',
     .          ''' X*10**', 'xx'',11X,', ''' Y*10**', 'xx'',11X,',
     .          ''' Z*10**', 'xx'',09X,', ''' DX/DT*', '10**xx'',',
     .          '7X,', ''' DY/DT*', '10**xx'',', '7X,', ''' DZ/DT*',
     .          '10**xx'')'/
      integer*4 jedz/24000/, nexpln/10000/

c external functions
      real*10 DOT
c
c*  start=2000
c insert initial conditions into first tabular point
      if(Ntab.ge.0) goto 800
      Ntab  = 1
      rdvsv = 0._10
      ifmt = Kkp(10)
      ifmt = max0(ifmt, -9)
      ifmt = min0(ifmt, 99)
      nexz = 1
      if(ifmt.lt.0 .or. ifmt.gt.9) nexz = 2
      call MVC('  ',1,2, sbfmt,8)
      call EBCDIX(ifmt, sbfmt, 8, nexz)
      call EBCDIX(jedz, sbfmth, 9, 5)
      call MVC(sbfmt,8,2, sbfmth, 33)
      call MVC(sbfmt,8,2, sbfmth, 49)
      call MVC(sbfmt,8,2, sbfmth, 65)
      call MVC(sbfmt,8,2, sbfmth, 85)
      call MVC(sbfmt,8,2, sbfmth,109)
      call MVC(sbfmt,8,2, sbfmth,133)
      call NEWPGS(Name, Jdp1, Jdp2, sbfmth)
      djedz      = 100._10*jedz
      if(Kout.gt.0) call SBEXPS(djedz, jedz)
      npcont     = 1000000
      kp98   = Kp(98)
      kp98   = iabs(kp98)
      iorder = 8
      jorder = 4
      ivl    = 6
      iordr1 = iorder - 3
      numprt = Iparp - 1
      numpt1 = numprt
      if(numprt.le.0) then
         jorder = 1
         numpt1 = 1
      endif
      jordr1 = jorder - 3
      jx     = 1
      do i = 1, 6
         Sb(i, 1, 1) = X0(i)
      end do
      lx = 6
      if(jx.lt.Iparp) then
         if(Ki(1).ne.0) then
            do j = 1, 6
               if(Ki(j+1).ge.0) then
                  jx = jx + 1
                  do i = 1, 6
                     if(Ki(1).gt.0) lx = lx + 1
                     Sb(i, jx, 1) = Dx0(i, j)
                  end do
               endif
            end do
         endif
         jx = jx + 1
         if(jx.le.Iparp) then
            do j = jx, Iparp
               do i = 1, 6
                  lx = lx + 1
                  Sb(i, j, 1) = V0(lx)
               end do
            end do
         endif
      endif
      Line   = 60
      ngo    = 1
      frect  = 0.5_10 + Hmx
      fsb1st = (T - Hmx) - Jdp0
c
c setup initial time record for output every step of
c nordsieck integration
      if(Kp(99).ge.0 .and. Kp(88).le.0) then
         juldat = T - frect
         je     = juldat + 0.5000001_10
         fsb    = (juldat - je) + 0.5_10
         do i = 1, 3
            ii = i + 3
            Tabpt(i, 1) = Sb(i, 1, 1)
            Tabpt(i, 2) = Sb(ii, 1, 1)
            Tabpt(i, 3) = ff0(ii, 1)
            if(iordr1.gt.0) then
               do j = 1, iordr1
                  Tabpt(i, j + 3) = aa0(ii, 1, j)
               end do
            endif
            if(numprt.gt.0) then
               do k = 1, numprt
                  kk = k + 1
                  dtabpt(i, 1, k) = Sb(i, kk, 1)
                  dtabpt(i, 2, k) = Sb(ii, kk, 1)
                  dtabpt(i, 3, k) = ff0(ii, kk)
                  if(jordr1.gt.0) then
                     do j = 1, jordr1
                        dtabpt(i, j + 3, k) = aa0(ii, kk, j)
                     end do
                  endif
               end do
            endif
         end do
         hc1 = Hmx
         if(Jdstrt.gt.0) hc1 = Hcc1
         hc2  = Hmx
         intl = 0
      endif
c
c print out tabular point data
  100 if(Kp(99).gt.0) goto 700
      mgo    = 1
      npcont = npcont + 1
      if(npcont.lt.kp98) goto 700
      npcont = 0
      if(Line.lt.58) goto 300
  200 call NEWPGT(Iout, Npage, jedz)
      Line = 2
      if(mgo.eq.2) goto 500
  300 juldat = T - frect

cxx   juldas=juldat-2400000
      juldas = juldat - djedz
      write(Iout, sbfmt) juldas, (Sb(i,1,Ntab), i = 1, 6)
      Line = Line + 1
      if(Kp(98).le.0) goto 700
      if(Iparp.le.1) goto 700
      mgo = 2
      jx  = 1
  400 jx  = jx + 1
      if(Line.ge.58) goto 200
  500 write(Iout, 600) jx, (Sb(i,jx,Ntab), i = 1, 6)
  600 format(i11, 1x, 1p, 6D20.13)
      Line = Line + 1
      if(jx.lt.Iparp) goto 400
  700 if(ngo.ne.2) then
         ngo   = 2
         frect = 0.5_10
c
c*  start=6000
c write record of satellite-probe tape
      else if(Iboth.gt.0) then
c
c*  start=5000
c backspace buffer, write probe data set
         n1sav = N1
         n2sav = N2
         l5sav = L5
         if(Kp(99).lt.0) goto 1300
         jf  = Jd
         fsb = Fract
         je  = Jd
         if(Kp(88).gt.0) then
            do j = 1, Iparp
               do i = 1, 6
                  sb1(i, j, 1) = Sb(i, j, 1)
               end do
            end do
            goto 1200
         else
            hc2    = -hc1
            hc1    = -hc1
            nwrite = 2
            goto 1100
         endif
      else
         mtabx = Ntab
         if(Kp(88).le.0) goto 1500
         if(Ntab.lt.5) goto 1600
         goto 1500
      endif
c
c*  start=3000
c insert integration results into tabular point
  800 Ntab = Ntab + 1
      if(Kp(88).le.0) Ntab = 1
      jx = 1
      kx = 6
c
c straight forward integration was employed
      if(Kp(100).lt.0) then
         do i = 1, 6
            Sb(i, 1, Ntab) = Y(i, 3)
         end do
         if(Ki(1).lt.0) then
            if(Kp(100)+2.ne.0) then
            endif
            Tlpt = T - Tlpt0
            goto 900
         else if(Ki(1).ne.0) then
            do j = 1, 6
               if(Ki(j+1).ge.0) then
                  jx = jx + 1
                  do i = 1, 6
                     kx = kx + 1
                     Sb(i, jx, Ntab) = Y(kx, 3)
                  end do
               endif
            end do
         endif
         goto 1000
      else if(Kp(100).eq.0) then
      endif
c
c enckes method of integration was employed
      do i = 1, 6
         Sb(i, 1, Ntab) = Y(i, 3) + Ylpt(i)
      end do
      if(Ki(1).ge.0) then
         if(Ki(1).ne.0) then
            do j = 1, 6
               if(Ki(j+1).ge.0) then
                  jx = jx + 1
                  do i = 1, 6
                     kx = kx + 1
                     Sb(i, jx, Ntab) = Y(kx, 3) + Dylpt(i, j)
                  end do
               endif
            end do
         endif
         goto 1000
      endif
  900 call ELIPT(-1, Tlpt)
      if(Icnd.eq.1) then
         do i = 1, 6
            Dylpt(i, 4) = Dylpt(i, 4) - Dylpt(i, 5)
            Dylpt(i, 5) = Dylpt(i, 5) - Dylpt(i, 6)
         end do
      endif
      do j = 1, 6
         if(Ki(j+1).ge.0) then
            jx = jx + 1
            do i = 1, 6
               Sb(i, jx, Ntab) = Dylpt(i, j)
            end do
         endif
c
c mean orbit integration was employed
      end do
c
c partials w.r.t. non-initial conditions
 1000 jx = jx + 1
      if(jx.le.Iparp) then
         do j = jx, Iparp
            do i = 1, 6
               kx = kx + 1
               Sb(i, j, Ntab) = Y(kx, 3)
            end do
         end do
      endif
c
c move higher derivitives  for kp(88)=0
      if(Kp(99).lt.0) goto 100
      if(Kp(88).gt.0) goto 100
      if(intl.gt.0) hc2 = Hc
      jplnts = Jplnt
      if(Iboth.ge.0) Jplnt = Ibuf
      nwrite = 1
 1100 write(Jplnt) je, fsb, iorder, jorder, numprt, numpt1,
     .             ((Tabpt(i,j),i=1,3), j = 1, iorder),
     .             (((dtabpt(i,j,k),i=1,3),j=1,jorder), k = 1, numpt1),
     .             hc1, hc2, zero9
      if(nwrite.ne.2) then
         Jplnt = jplnts
         je    = T
         fsb   = T - je
         if(intl.gt.0) then
            hc1 = Hc
         else
            intl = 1
         endif
         do i = 1, 3
            ii = i + 3
            Tabpt(i, 1) = Sb(i, 1, 1)
            Tabpt(i, 2) = Sb(ii, 1, 1)
            Tabpt(i, 3) = ff(ii, 1, 3)
            if(iordr1.gt.0) then
               ij = ii
               do j = 1, iordr1
                  if(j.eq.5) ij = ii + 2
                  Tabpt(i, j + 3)    = aa(ij, 1, j)
               end do
            endif
            if(numprt.gt.0) then
               do k = 1, numprt
                  kk = k + 1
                  dtabpt(i, 1, k) = Sb(i, kk, 1)
                  dtabpt(i, 2, k) = Sb(ii, kk, 1)
                  dtabpt(i, 3, k) = ff(ii, kk, 3)
                  if(jordr1.gt.0) then
                     ij = ii
                     do j = 1, jordr1
                        if(j.eq.5) ij  = ii + 2
                        dtabpt(i, j + 3, k) = aa(ij, kk, j)
                     end do
                  endif
               end do
            endif
         end do
         goto 100
      endif
 1200 do while(.true.)
         backspace Ibuf
         if(Kp(88).gt.0) then
            read(Ibuf) je, fsb1, ivl,
     .                 (((Sb(i,j,k),i=1,ivl),j=1,Iparp), k = 1, 5)
            backspace Ibuf
            do k = 2, 5
               kk = 5 - (k - 2)
               do j = 1, Iparp
                  do i = 1, 6
                     sb1(i, j, k) = Sb(i, j, kk)
                  end do
               end do
            end do
            write(Jplnt) jf, fsb, ivl,
     .                   (((sb1(i,j,k),i=1,ivl),j=1,Iparp), k = 1, 5),
     .                   zero9
            do j = 1, Iparp
               do i = 1, 6
                  sb1(i, j, 1) = Sb(i, j, 1)
               end do
            end do
            jf  = je
            fsb = fsb1
            if(fsb.eq.fsb1st) then
               if(je.eq.Jdp0) then
                  rewind Ibuf
                  goto 1300
               endif
            endif
         else
            read(Ibuf) je, fsb, iorder, jorder, numprt, numpt1,
     .                 ((Tabpt(i,j),i=1,3), j = 1, iorder),
     .                 (((dtabpt(i,j,k),i=1,3),j=1,jorder), k = 1,
     .                 numpt1), hc2, hc1
            hc2 = -hc2
            hc1 = -hc1
            backspace Ibuf
            if(fsb.ne.fsb1st) goto 1100
            if(je.ne.Jdp0) goto 1100
            rewind Ibuf
            goto 1300
         endif
      end do
 1300 Iboth = -1
      Nsign = -Nsign
      Hc    = -Hsave
      Hmx   = -Hmx
      Pint5 = 4._10*Hmx
      Tstop = Jdp2
      Tstop = Tstop + Epsp(2) + Hmx
      if(Kp(88).gt.0) Tstop = Tstop + Pint5
      Both = ABS(T - T0)
      T    = T0
      L4   = -10
      L1   = -1
      Ntab = -1
      do i = 1, N
         A(i)  = A0(i)
         B(i)  = B0(i)
         C(i)  = C0(i)
         D(i)  = D0(i)
         Ee(i) = Ee0(i)
         Dydt(i, 1) = Dydt0(i)
      end do
      N1 = n1sav
      N2 = n2sav
      L5 = l5sav
      write(Iout, 1400)
 1400 format(/
     .' INTEGRATION COMPLETED IN ONE DIRECTION FROM EPOCH, STARTED IN OT
     .HER'/)
      Line = Line + 3
      goto 1600
 1500 Ntab = 0
      Nrec = Nrec + 1
      if(Iboth.ge.0) then
         if(Nsign*((Jdp1-T)+Epsp(1)).le.0._10) Iboth = 1
      endif
      if(Kp(99).ge.0) then
         if(Kp(88).gt.0) then
            juldat = T - Pint5
            je     = juldat
            fsb    = juldat - je
            if(Iboth.ge.0) then
               jplnts = Jplnt
               Jplnt  = Ibuf
            endif
            write(Jplnt) je, fsb, ivl,
     .                   (((Sb(i,j,k),i=1,ivl),j=1,Iparp), k = 1, 5),
     .                   zero9
            if(Iboth.ge.0) Jplnt = jplnts
         endif
      endif
 1600 if(Kout.gt.0) then

         if(Kkp(2).lt.0) then
            kkp2 = -Kkp(2)
            rdv  = DOT(Sb(1,1,mtabx), Sb(4,1,mtabx))
            if(rdv*rdvsv.le.0) then
               nexpln = nexpln + 1
               if(nexpln.ge.kkp2) then
                  nexpln = 0
                  call SBEXP(Sb, T, mtabx, sbfmt)
               endif
            endif
            rdvsv = rdv
         else if(Kkp(2).ne.0) then
            nexpln = nexpln + 1
            if(nexpln.ge.Kkp(2)) then
               nexpln = 0
               call SBEXP(Sb, T, mtabx, sbfmt)
            endif
         endif

c printout for close approaches
         if(Kkp(6).gt.0) call SBEXPT(T)
      endif
c
c setup starting proceedure for thrust initiation or
c termination as signaled by subroutine lesthp
      Jthrst = 0
      Lthrst = 0
      T0     = T0sav
      if(Kthrst.gt.0) then
         Kthrst = 0
         Jthrst = 1
         Lthrst = 1
c do not set approximating polynomial coefficients to zero
c in starting proceedure
         T0   = T
         kk95 = Kkp(95)
         Hc   = 2._10**kk95

c station keeping thrusting assumes forward in time integration
         do i = 1, N
            V0(i) = Y(i, 3)
         end do
         if(Kp(100).ge.0) call SUICID(
     .' CANNOT RESTART ENCKE INTEGRATION FOR THRUST, STOP IN SBOUT ',15)
      endif
c
c*  start=9900
      return
      end
