      subroutine MOROUT
 
      implicit none
c
c r.king and r.cappallo  july 1977    subroutine morout
c output routine for simultaneous moon orbit and rotation
c integration.  based on m.e.ash subroutine monout, oct 69.
c
      include 'globdefs.inc'
c
c common
      include 'adams.inc'
      include 'bdydta.inc'
      include 'cassin.inc'
      include 'dumdum.inc'
      include 'ellips.inc'
      include 'empcnd.inc'
      real*10 Meqinc
      equivalence (Mrcond(11),Meqinc)
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'metuna.inc'
      include 'morstf.inc'
      include 'namtim.inc'
      include 'nmstrt.inc'
      include 'orblun.inc'
      include 'output.inc'
      include 'plndta.inc'
      include 'stint.inc'
c data from moon numerical integration to be written on tape
c plus moon rotation and quantities transferred from n-body
      common/WRKCOM/ Moon1(6,i_mxplprt+1,8),Librt1(6,i_mxplprt+1,5),
     . Corrt1(6,i_mxplprt+1,5),Dpsi1(8),Deps1(8),
     . Lib(3,8),Dpsi2(8),Deps2(8),Lib1(3,8)
      real*10 Moon1,Librt1,Corrt1
      real*4 Dpsi1,Deps1,Dpsi2,Deps2,Lib,Lib1
 
c suffix '2' means readback buffer, overlay with coefficients
      real*10 moon2(6,i_mxplprt+1,8), librt2(6,i_mxplprt+1,5)
      equivalence (moon2(1,1,1),librt2(1,1,1),A(1))
      include 'zeroes.inc'
 
      integer*4 i,ix,j,j1,je,jf,jx,k,km98,kx,mgo,ngo,npcont,nut
      integer*4 mvl/6/
      real*10 juldat,frect
      real*10 mfract,tm
      character*1 labl(12),pnt/'.'/,blank/' '/
      character*12 lab4
      equivalence (labl,lab4)
      logical*4 orbprt,rotprt
      integer ibufo,ibufr,ibufc,ilibs,imns,icors
c
c
c-----------------------------------------------------------------------
c           insert initial conditions into first tabular point
c
c           insert initial conditions for orbit
      if(Ntab .ge. 0 .or. Ntabr .ge. 0) goto 1400
 
c first assign buffers. must decode simul. buffers
      ibufo = Ibuf
      ibufr = Ibuf
      if(Kkmr(13).gt.0) ibufr=Kkmr(13)
      ibufc = Kkmr(12)
      if(Orbint .and. Rotint) then
         if(Ibuf.gt.100) then
            ibufo = Ibuf/100
            if(ibufo.eq.ibufr) then
               ibufo = Ibuf - 100*ibufo
            else if(Kkmr(13).le.0) then
               ibufr = Ibuf - 100*ibufo
            endif
         endif
      endif
      if(Orbint) then
         Ntab = 1
         jx   = 1
         do i = 1, 6
            Moon1(i, 1, 1) = X0m(i, 1)
         end do
         if(Kim(1) .ne. 0) then
            do j = 1, 6
               if(Kim(j+1) .ge. 0) then
                  jx = jx + 1
                  do i = 1, 6
                     Moon1(i, jx, 1) = Dx0m(i, j, 1)
                  end do
               endif
            end do
         endif
         jx = jx + 1
         if(jx .le. Iparm) then
            ix = 6*(jx - 1)
            do j = jx, Iparm
               do i = 1, 6
                  ix = ix + 1
                  Moon1(i, j, 1) = V0(ix)
               end do
            end do
         endif
         if(Rotint) then
c insert libration initial conditions as well
            call MONROT(0,X0m(1,2),X0m(2,2),X0m(3,2),
     .       X0m(4,2),X0m(5,2),X0m(6,2))
            call ECLPRC(Jdmn0,0._10,1)
            call MONROT(2,X0m(1,2),X0m(2,2),X0m(3,2),
     .       X0m(4,2),X0m(5,2),X0m(6,2))
            Lib(1,1)=Tauc
            Lib(2,1)=Rhoc
            Lib(3,1)=Meqinc*Tauc-Isig
         endif
      endif
c
c insert initial conditions for rotation
      if(Rotint) then
         Ntabr = 1
         jx    = 1
         j1    = 1
         if(Orbint) j1 = 2
         do i = 1, 6
            Librt1(i,1,1) = X0m(i,j1)
            if(Corint) Corrt1(i,1,1) = X0m(i,j1+1)
         end do
         if(Kir(1) .ne. 0) then
            do j = 1, 6
               if(Kir(j+1) .ge. 0) then
                  jx = jx + 1
                  do i = 1, 6
                     Librt1(i,jx,1) = Dx0m(i,j,j1)
                  end do
               endif
            end do
         endif
         jx = jx + 1
         if(jx .le. Iparmr) then
            ix = 6*(jx - 1)
            do j = jx, Iparmr
               do i = 1, 6
                  ix = ix + 1
                  Librt1(i, j, 1) = V0(ix + Korb)
               end do
            end do
         endif
         if(Corint) then
            ix = 6
            do j = 2,Iparmr
               do i = 1, 6
                  ix = ix + 1
                  Corrt1(i, j, 1) = V0(ix + Krot)
               end do
            end do
         endif
      endif
      Line  = 60
      ngo   = 1
      frect = 0.0_10
      if(Hc .gt. 0.0_10) frect = 1.0_10
 
c set up printing every nth tabular point
      km98 = Km(98)
      if( .not. Orbint) km98 = Kmr(98)
      km98   = iabs(km98)
      npcont = 999999
      orbprt = Orbint .and. (Km(99) .le. 0)
      rotprt = Rotint .and. (Kmr(99) .le. 0)
c
c-----------------------------------------------------------------------
c
c print out tabular point data
  100 if( .not. (orbprt .or. rotprt)) goto 1300
      mgo    = 1
      npcont = npcont + 1
      if(npcont .lt. km98) goto 1300
      npcont = 0
      if(Line .lt. 58) goto 500
  200 call NEWPGT(Iout, Npage, 0)
      write(Iout, 300)
  300 format(4x, 'JED', 14x, 'X', 19x, 'Y', 19x, 'Z', 17x, 'DX/DT',
     .       15x, 'DY/DT', 15x, 'DZ/DT')
      if(rotprt) write(Iout, 400)
  400 format(20x, 'PSI', 16x, 'THETA', 16x, 'PHI', 15x, 'DPSI/DT', 12x,
     .       'DTHETA/DT', 12x, 'DPHI/DT')
      Line = 3
      if(mgo .eq. 2) goto 800
      if(mgo .eq. 3) goto 1200
      if(mgo .eq. 4) goto 1250
  500 juldat = T - frect
      i = juldat
      call EBCDI(i, labl, 8)
      i = 1000.5_10 + 1E3_10*(juldat - i)
      call EBCDI(i, labl(9), 4)
      labl(9) = pnt
      if(orbprt) then
         write(Iout, 550) lab4, (Moon1(i,1,Ntab), i = 1, 6)
  550    format(A12, 6F20.16)
         Line = Line + 1
         if( .not. rotprt) goto 600
         call MOVEBL(blank, 1, labl, 12)
      endif
      tm = Librt1(3, 1, Ntabr)
 
c print phi modulo 2 pi, but keep original for tape
      Librt1(3, 1, Ntabr) = MOD(tm, Twopi)
      write(Iout, 550) lab4, (Librt1(i,1,Ntabr), i = 1, 6)
      Librt1(3, 1, Ntabr) = tm
      Line = Line + 1
      if(Corint) then
         tm = Corrt1(3, 1, Ntabr)
         Corrt1(3, 1, Ntabr) = MOD(tm, Twopi)
         write(Iout, 550) blank, (Corrt1(i,1,Ntabr), i = 1, 6)
         Corrt1(3, 1, Ntabr) = tm
         Line = Line + 1
      endif
 
c print partials, if there and desired
  600 if( .not. orbprt .or. Km(98) .le. 0 .or. Iparm .le. 1)
     .    goto 1000
      mgo = 2
      jx  = 1
  700 jx  = jx + 1
      if(Line .ge. 58) goto 200
  800 write(Iout, 900) jx, (Moon1(i,jx,Ntab), i = 1, 6)
  900 format(i11, 1x, 1p, 6D20.13)
      Line = Line + 1
      if(jx .lt. Iparm) goto 700
 1000 if( .not. rotprt .or. Kmr(98) .le. 0 .or. Iparmr .le. 1)
     .    goto 1300
      mgo = 3
      jx  = 1
 1100 jx  = jx + 1
      if(Line .ge. 58) goto 200
 1200 write(Iout, 900) jx, (Librt1(i,jx,Ntabr), i = 1, 6)
      Line = Line + 1
      if(jx .lt. Iparmr) goto 1100
      if(.not.Corint) goto 1300
      mgo = 4
      jx  = 1
 1210 jx  = jx + 1
      if(Line .ge. 58) goto 200
 1250 write(Iout, 900) jx, (Corrt1(i,jx,Ntabr), i = 1, 6)
      Line = Line + 1
      if(jx .lt. Iparmr) goto 1210
 1300 if(ngo .ne. 2) then
         ngo   = 2
         frect = 0.5_10
c
c-----------------------------------------------------------------------
c
c           write records of output data sets
c
      else if(Iboth .gt. 0) then
c
c-----------------------------------------------------------------------
c
c           if integration in both directions, backspace buffer written
c           with results in one direction, and write output data set in
c           other direction.
         jf = Jd
         k  = 2*iabs(jf - Jdbd(1)) + 1
         if( .not. Orbint) goto 1900
         if(Km(99) .lt. 0) goto 1900
         Dpsi2(1) = Psid(k)
         Deps2(1) = Epsd(k)
         if(.not. Rotint) then
            do i = 1, 3
               Lib1(i, 1) = Librt(k, i)
            end do
         else
            do i = 1, 3
               Lib1(i,1) = Lib(i,1)
            end do
         endif
         do j = 1, Iparm
            do i = 1, 6
               moon2(i, j, 1) = Moon1(i, j, 1)
            end do
         end do
         do while( .true. )
            backspace ibufo
            read(ibufo) je, mfract, mvl,
     .                  (((Moon1(i,j,k),i=1,mvl),j=1,Iparm), k = 1, 8),
     .                  (Dpsi1(i), Deps1(i), i = 1, 8), Lib
            backspace ibufo
            do k = 2, 8
               Kk = 8 - (k - 2)
               Dpsi2(k) = Dpsi1(Kk)
               Deps2(k) = Deps1(Kk)
               do i = 1, 3
                  Lib1(i, k) = Lib(i, Kk)
               end do
               do j = 1, Iparm
                  do i = 1, 6
                     moon2(i, j, k) = Moon1(i, j, Kk)
                  end do
               end do
            end do
            write(Imn) jf, mfract, mvl,
     .                 (((moon2(i,j,k),i=1,mvl),j=1,Iparm), k = 1, 8),
     .                 (Dpsi2(i), Deps2(i), i = 1, 8), Lib1, zero9
            do j = 1, Iparm
               do i = 1, 6
                  moon2(i, j, 1) = Moon1(i, j, 1)
               end do
            end do
            Dpsi2(1) = Dpsi1(1)
            Deps2(1) = Deps1(1)
            do i = 1, 3
               Lib1(i, 1) = Lib(i, 1)
            end do
            jf = je
            if(je .eq. Jdmn0) then
               rewind ibufo
               goto 1900
            endif
         end do
      else
c this code reached from fortran statement label = 99
c
c write record of moon tape
         if(Orbint) then
            if(Ntab .ge. 8) then
               Ntab = 0
               Nrec = Nrec + 1
               if(Iboth .ge. 0) then
                  if(Nsign*(Jdm1-Jd) .le. 0) Iboth = 1
               endif
               if(Km(99) .ge. 0) then
                  juldat = T - Mint
                  je     = juldat + 1E-5_10
                  mfract = juldat - je
 
c tabular interval must be 0.5 (intm=-1) to transfer nutation
                  nut = 2*iabs(je - Jdbd(1))
                  if(Msign .lt. 0) nut = nut + 2
                  do j = 1, 8
                     k = nut + ISIGN(j, Msign)
                     Dpsi1(j) = Psid(k)
                     Deps1(j) = Epsd(k)
                     if(.not. Rotint) then
                        do i = 1, 3
                           Lib(i, j) = Librt(k, i)
                        end do
                     endif
                  end do
                  if(Iboth .ge. 0) then
                     imns = Imn
                     Imn  = ibufo
                  endif
                  write(Imn) je, mfract, mvl,
     .                       (((Moon1(i,j,k),i=1,mvl),j=1,Iparm),
     .                       k = 1, 8), (Dpsi1(i), Deps1(i), i = 1, 8),
     .                       Lib, zero9
                  if(Iboth .ge. 0) Imn = imns
               endif
            endif
         endif
c
c write record of moon rotation tape and moon core rotation tape
         if(Rotint) then
            if(Ntabr .ge. 5) then
               Ntabr = 0
               Nrec  = Nrec + 1
               if(Iboth .ge. 0) then
                  if(Nsign*((Jdm1-T)+Epsm(1)).le.0._10) Iboth = 1
               endif
               if(Kmr(99) .ge. 0) then
                  juldat = T - Mint5
                  je     = juldat
                  mfract = juldat - je
                  ilibs  = Ilib
                  if(Iboth .ge. 0) Ilib = ibufr
                  write(Ilib) je, mfract, mvl,
     .             (((Librt1(i,j,k),i=1,mvl),j=1,Iparmr),k=1,5),zero9
                  Ilib = ilibs
                  if(Corint .and. Icor.gt.0) then
                     icors=Icor
                     if(Iboth.ge.0) Icor=ibufc
                     write(Icor) je, mfract, mvl,
     .                (((Corrt1(i,j,k),i=1,mvl),j=1,Iparmr),k=1,5),zero9
                     Icor=icors
                  endif
               endif
            endif
         endif
         return
      endif
c
c-----------------------------------------------------------------------
c
c           insert integration results into tabular point
c
c           results for orbit
 1400 if( .not. Orbint) goto 1800
      Ntab = Ntab + 1
      jx   = 1
      kx   = 6
c
c straight forward integration was employed
      if(Km(100) .lt. 0) then
         do i = 1, 6
            Moon1(i, 1, Ntab) = Y(i, 3)
         end do
         if(Kim(1) .lt. 0) then
            if(Km(100) + 2 .lt. 0) then
               Tlpt = T - Tlpt0
               goto 1500
            else if(Km(100) + 2 .eq. 0) then
               goto 1600
            endif
         else if(Kim(1) .eq. 0) then
            goto 1700
         endif
         do j = 1, 6
            if(Kim(j+1) .ge. 0) then
               jx = jx + 1
               do i = 1, 6
                  kx = kx + 1
                  Moon1(i, jx, Ntab) = Y(kx, 3)
               end do
            endif
         end do
      else if(Km(100) .eq. 0) then
c
c enckes method of integration was employed
         do i = 1, 6
            Moon1(i, 1, Ntab) = Y(i, 3) + Ylpt(i)
         end do
         if(Kim(1) .lt. 0) goto 1500
         if(Kim(1) .ne. 0) then
            do j = 1, 6
               if(Kim(j+1) .ge. 0) then
                  jx = jx + 1
                  do i = 1, 6
                     kx = kx + 1
                     Moon1(i, jx, Ntab) = Y(kx, 3) + Dylpt(i, j)
                  end do
               endif
            end do
         endif
      else
c
c mean lunar orbit integration was employed
         do i = 1, 6
            Moon1(i, 1, Ntab) = Y(i, 3) + Ylun(i)
         end do
         if(Kim(1) .lt. 0) goto 1600
         if(Kim(1) .ne. 0) then
            do j = 1, 6
               if(Kim(j+1) .ge. 0) then
                  jx = jx + 1
                  do i = 1, 6
                     kx = kx + 1
                     Moon1(i, jx, Ntab) = Y(kx, 3) + Dylun(i, j)
                  end do
               endif
            end do
         endif
      endif
      goto 1700
 1500 call ELIPT(-1,Tlpt)
      do j = 1, 6
         if(Kim(j+1) .gt. 0) then
            jx = jx + 1
            do i = 1, 6
               Moon1(i, jx, Ntab) = Dylpt(i, j)
            end do
         endif
      end do
      goto 1700
 1600 call LUNORB(Jd, Fract, -1)
      do j = 1, 6
         if(Kim(j+1) .gt. 0) then
            jx = jx + 1
            do i = 1, 6
               Moon1(i, jx, Ntab) = Dylun(i, j)
            end do
         endif
      end do
c
c partials w.r.t. non-initial conditions
 1700 jx = jx + 1
      if(jx .le. Iparm) then
         do j = jx, Iparm
            do i = 1, 6
               kx = kx + 1
               Moon1(i, j, Ntab) = Y(kx, 3)
            end do
         end do
      endif
c
c insert libration for moon tape if integrated
      if(Rotint) then
         if(Kout.eq.0 .or. Fract.ne.0._10) then
            call ECLPRC(Jd,Fract,1)
            call MONROT(2,Y(Korb+1,3),Y(Korb+2,3),Y(Korb+3,3),
     .       Y(Korb+4,3),Y(Korb+5,3),Y(Korb+6,3))
         endif
         Lib(1,Ntab)=Tauc
         Lib(2,Ntab)=Rhoc
         Lib(3,Ntab)=Meqinc*Tauc-Isig
      endif
c
c insert results for rotation and core rotation
 1800 if(Rotint) then
         Ntabr = Ntabr + 1
         do i = 1, 6
            Librt1(i,1,Ntabr) = Y(i+Korb,3)
            if(Corint) Corrt1(i,1,Ntabr) = Y(i+Krot,3)
         end do
         kx    = 6
c
c partials
         do j = 2,Iparmr
            do i = 1, 6
               kx = kx + 1
               Librt1(i,j,Ntabr) = Y(kx+Korb,3)
               if(Corint) Corrt1(i,j,Ntabr) = Y(kx+Krot,3)
            end do
         end do
      endif
      goto 100
 
 1900 if(Kmr(99) .ge. 0) then
         if(Rotint) then
c copy and reverse buffered output that was written backwards from the
c starting point, first: libration
            jf = Jd
            tm = Fract
            do j = 1, Iparmr
               do i = 1, 6
                  librt2(i, j, 1) = Librt1(i, j, 1)
               end do
            end do
            do while( .true. )
               backspace ibufr
               read(ibufr) je, mfract, mvl,
     .                     (((Librt1(i,j,k),i=1,mvl),j=1,Iparmr),
     .                     k = 1, 5)
               backspace ibufr
               do k = 2, 5
                  Kk = 5 - (k - 2)
                  do j = 1, Iparmr
                     do i = 1, 6
                        librt2(i, j, k) = Librt1(i, j, Kk)
                     end do
                  end do
               end do
               write(Ilib) jf, tm, mvl,
     .                     (((librt2(i,j,k),i=1,mvl),j=1,Iparmr),
     .                     k = 1, 5), zero9
               do j = 1, Iparmr
                  do i = 1, 6
                     librt2(i, j, 1) = Librt1(i, j, 1)
                  end do
               end do
               jf = je
               tm = mfract
               if(je .eq. Jdmr0) then
                  rewind ibufr
                  goto 1950
               endif
            end do

c next: core rotation
 1950       if(Icor.le.0 .or. .not.Corint) goto 2000
            jf = Jd
            tm = Fract
            do j = 1, Iparmr
               do i = 1, 6
                  librt2(i, j, 1) = Corrt1(i, j, 1)
               end do
            end do
            do while( .true. )
               backspace ibufc
               read(ibufc) je, mfract, mvl,
     .          (((Corrt1(i,j,k),i=1,mvl),j=1,Iparmr),k=1,5)
               backspace ibufc
               do k = 2, 5
                  Kk = 5 - (k - 2)
                  do j = 1, Iparmr
                     do i = 1, 6
                        librt2(i, j, k) = Corrt1(i, j, Kk)
                     end do
                  end do
               end do
               write(Icor) jf, tm, mvl,
     .          (((librt2(i,j,k),i=1,mvl),j=1,Iparmr),k=1,5), zero9
               do j = 1, Iparmr
                  do i = 1, 6
                     librt2(i, j, 1) = Corrt1(i, j, 1)
                  end do
               end do
               jf = je
               tm = mfract
               if(je .eq. Jdmr0) then
                  rewind ibufc
                  goto 2000
               endif
            end do
         endif
      endif
 
 2000 Iboth = -1
      Nsign = -Nsign
      Hc    = -Hsave
      Hmx   = -Hmx
      Mint  = -Mint
      Mint5 = -Mint5
      Msign = -Msign
      Tstop = Jdm2 + ISIGN(4, Nsign)
      if( .not. Orbint) Tstop = Jdm2 + Epsm(2) + Mint5 + Hmx
      Both  = ABS(T - T0)
      T     = T0
      L4    = -10
      L1    = -1
      Ntab  = -1
      Ntabr = -1
      do i = 1, N
         A(i)  = A0(i)
         B(i)  = B0(i)
         C(i)  = C0(i)
         D(i)  = D0(i)
         Ee(i) = Ee0(i)
         Dydt(i, 1) = Dydt0(i)
      end do
      write(Iout, 2100)
 2100 format(/
     .' INTEGRATION COMPLETED IN ONE DIRECTION FROM EPOCH, STARTED IN OT
     .HER'/)
      Line = Line + 3
 
      return
      end
