      subroutine BODOUT
 
      implicit none
c
c ash/smith/connolly    may 1968    subroutine bodout
c output routine for n-body integration
c
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'adams.inc'
      include 'bdctrl.inc'
      include 'bddtaint.inc'
      include 'bdydta.inc'
      include 'bodstf.inc'
      include 'cassin.inc'
      include 'dumdum.inc'
      include 'empcnd.inc'
      real*10 Meqinc
      equivalence (Mrcond(11),Meqinc)
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'nmstrt.inc'
      include 'output.inc'
      include 'stint.inc'
 
c data from n-body numerical integration to be written on tape
      common/WRKCOM/ Merc1(6,10),Venpl1(6,5,8),Moon1(6,40),Mrot1(6,40),
     .        Merc2(6,10),Venpl2(6,5,8),Moon2(6,40),Mcor1(6,40),
     .        Nutat1(2,40),Librt1(3,40),Nutat2(2,40),
     .        Librt2(3,40)
      real*10 Merc1,Merc2,Venpl1,Venpl2,Moon1,Moon2,Mrot1,Mcor1
      real*4  Nutat1,Nutat2,Librt1,Librt2

c local
      integer*4 i,ibodys,igo,ihmx,ivl,j,jdir,je,jf,k,kk,
     .          kprint,l,mntb,mvl,nb,ngo,np,np1,np2
      real*10 juldat,frect
      real*10 fsb,fsb1
      data ivl,mvl/6,3/
 
      if(Ntab.ge.0) goto 300
      if(Kout.gt.0 .and. MOD(Intxpr/64,2).eq.1) write(Kout,50)
   50 format(4x,'JED',17x,'X',25x,'Y',25x,'Z',24x,
     . 'DX/DT',21x, 'DY/DT',21x,'DZ/DT')
      Ntab = 1
      Mtab = 1
      mvl  = 6
      nb   = Nbodyp - 1
      l    = 0
      do i = 1, 6
         l = l + 1
         Merc1(i,1) = V0(l)
      end do
      do j = 1, nb
         do i = 1, 6
            l = l + 1
            Venpl1(i,1,j) = V0(l)
         end do
      end do
      if(Jmoon.gt.0) then
         Ntabr= 1
         do i = 1, 6
            Moon1(i,1) = X0m(i,1)
            Mrot1(i,1) = X0m(i,2)
            Mcor1(i,1) = X0m(i,3)
         end do
         if(ABS(Mrot1(3,1)).gt.Twopi) Mrot1(3,1)=MOD(Mrot1(3,1),Twopi)
c insert libration initial conditions as well
         call MONROT(0,X0m(1,2),X0m(2,2),X0m(3,2),
     .    X0m(4,2),X0m(5,2),X0m(6,2))
         call ECLPRC(Jdbdy0,0._10,1)
         call MONROT(2,X0m(1,2),X0m(2,2),X0m(3,2),
     .    X0m(4,2),X0m(5,2),X0m(6,2))
         Librt1(1,1)=Tauc
         Librt1(2,1)=Rhoc
         Librt1(3,1)=Meqinc*Tauc-Isig
         if(Jct(21).le.0) then
            call PEPNUT((Jdbdy0-2451545.5_10)/36525._10,Nutat1(1,1),
     .       Nutat1(2,1))
         else
            call IAU2000A(Jdbdy0-0.5_10,Nutat1(1,1),Nutat1(2,1))
         endif
      else
c if not integrating moon, indicate moon is always ready (from other)
         Ntabr=40
      endif         
      Line   = 60
      kprint = Kbdy(38)
      ngo    = 1
      frect  = 0.5_10 + Hmx
c
c print out tabular point data
  100 if(Kbdy(39).le.0) then
         kprint = kprint + 1
         if(kprint.ge.Kbdy(38)) then
            kprint = 0
            do igo=0,Nbodyp
               if(Line.ge.58) then
                  call NEWPG
                  write(Iout,110)
  110             format(4x,'JED', 14x, 'X', 19x, 'Y', 19x, 'Z', 15x,
     .                   'DX/DT*10**2', 9x, 'DY/DT*10**2', 9x,
     .                   'DZ/DT*10**2')
                  Line = 2
               endif
               if(igo.le.0) then
                  juldat = T - frect
                  write(Iout,120) juldat,(Merc1(i,Mtab),i = 1,6)
  120             format(f12.3,3F20.16,2p,3F20.16)
                  Line = Line+1
                  if(Kout.gt.0 .and. MOD(Intxpr/64,2).eq.1)
     .             write(Kout,125) juldat,(Merc1(i,Mtab),i = 1,6)
  125             format(f12.3,1p6e26.19)
               else if(igo.le.nb) then
                  write(Iout,130) (Venpl1(i,Ntab,igo),i = 1,6)
  130             format(12x,0p3F20.16,2p3F20.16)
                  Line = Line+1
                  if(Kout.gt.0 .and. MOD(Intxpr/64,2).eq.1)
     .             write(Kout,135) (Venpl1(i,Ntab,igo),i = 1,6)
  135             format(12x,1p6e26.19)
               else if(Jmoon.gt.0) then
                  write(Iout,130) (Moon1(i,Ntabr),i=1,6),
     .             (Mrot1(i,Ntabr),i=1,6)
                  Line = Line+2
                  if(Kout.gt.0 .and. MOD(Intxpr/64,2).eq.1) then
                     write(Kout,135) (Moon1(i,Ntabr),i=1,6),
     .                (Mrot1(i,Ntabr),i=1,6)
                     if(Corint) write(Kout,135) (Mcor1(i,Ntabr),i=1,6)
                  endif
               endif
            end do
         endif
      endif
      if(ngo.eq.2) goto 600
      ngo   = 2
      frect = 0.5_10
c
c insert integration results into tabular point
  300 if(Jmoon.gt.0) then
         Ntabr=Ntabr+1
         l=(Jmoon-1)*6
         do i=1,6
            l=l+1
            Moon1(i,Ntabr)=Y(l,3)
            Mrot1(i,Ntabr)=Y(l+6,3)
            Mcor1(i,Ntabr)=Y(l+12,3)
         end do
         if(ABS(Mrot1(3,Ntabr)).gt.Twopi) Mrot1(3,Ntabr)=
     .    MOD(Mrot1(3,Ntabr),Twopi)
         call ECLPRC(Jd,Fract,1)
         call MONROT(2,Y(l+1,3),Y(l+2,3),Y(l+3,3),
     .       Y(l+4,3),Y(l+5,3),Y(l+6,3))
         Librt1(1,Ntabr)=Tauc
         Librt1(2,Ntabr)=Rhoc
         Librt1(3,Ntabr)=Meqinc*Tauc-Isig
         if(Jct(21).le.0) then
            call PEPNUT((T-2451545.5_10)/36525._10,Nutat1(1,Ntabr),
     .       Nutat1(2,Ntabr))
         else
            call IAU2000A(T-0.5_10,Nutat1(1,Ntabr),Nutat1(2,Ntabr))
         endif
         if(MOD(Ntabr,4).ne.1) goto 600
      endif
      Mtab = Mtab + 1
      l    = 0
      do i = 1, 6
         l = l + 1
         Merc1(i,Mtab) = Y(l,3)
      end do
      if(Mtab.ne.((Mtab/2)*2)) then
         Ntab = Ntab + 1
         do j = 1, nb
            do i = 1, 6
               l = l + 1
               Venpl1(i,Ntab,j) = Y(l,3)
            end do
         end do
         goto 100
      endif
c we now have a pair of mercury tabular points corresponding to a single
c set of tabular points for the other planets.  move the corresponding
c octet of moon, nutation, library from /bdydta/ to complete the set
c if we are writing an output tape
      if(Kbdy(39).ge.0 .and. Jmoon.le.0) then
         mntb = (Mtab-2)*4
         jf   = T
         je   = T-Hmx
         jdir = 1
         np   = IABS(je - Jdbd(1))*2
         if((jf-je)*Ibdsgn.le.0) then
            jdir = -1
            np   = np + 2
         endif
         do i = 1, Nbody
            if(Nplbdy(i).eq.10) goto 750
         end do
         np1 = np + jdir
         np2 = np + jdir*8
         do i = 1, 3
            if(mvl.gt.Mvel(i) .and.
     .       np1.le.Recend(i) .and. np2.ge.Recbeg(i)) then
               mvl = 3
               goto 700
            endif
         end do
  700    k=np
         do j = 1, 8
            k = k + jdir
            do i = 1, mvl
               Moon1(i,mntb+j) = Mon(i,k)
            end do
         end do
  750    k=np
         do j = 1, 8
            k = k + jdir
            Nutat1(1,mntb+j) = Psid(k)
            Nutat1(2,mntb+j) = Epsd(k)
            do i = 1, 3
               Librt1(i,mntb+j) = Librt(k,i)
            end do
         end do
      endif
      goto 600
c
c setup integration in other direction from epoch
  400 Iboth = -1
      Nsign = -Nsign
      Hc    = -Hsave
      Hmx   = -Hmx
      Bint  = -Bint
      ihmx  = Intbdy
      if(ihmx.le.0) ihmx=1
      Tstop = Jdbdy2 + ISIGN(10*ihmx,Nsign)
      Both  = ABS(T - T0)
      T     = T0
      L4    = -10
      L1    = -1
      Ntab  = -1
      do i = 1, N
         A(i)  = A0(i)
         B(i)  = B0(i)
         C(i)  = C0(i)
         D(i)  = D0(i)
         Ee(i) = Ee0(i)
         Dydt(i,1) = Dydt0(i)
      end do
      write(Iout,500)
  500 format(
     .'0INTEGRATION COMPLETED IN ONE DIRECTION FROM EPOCH, STARTED IN OT
     .HER'/)
      Line = Line + 3
      return
c
c decide whether to write record of n-body tape
  600 if(Iboth.gt.0) then
c
c backspace buffer, write n-body data set
         if(Kbdy(39).lt.0) goto 400
 
c move integration result tabular point
         jf  = T
         fsb = T - jf
         do i = 1, 6
            Merc2(i,1) = Merc1(i,1)
            do k = 1, nb
               Venpl2(i,1,k) = Venpl1(i,1,k)
            end do
         end do
         if(Jmoon.gt.0) then
            do i=1,6
               Moon2(i,1)=Moon1(i,1)
            end do
            do i=1,3
               Librt2(i,1)=Librt1(i,1)
            end do
            if(Jct(21).le.0) then
               call PEPNUT((T-2451545.5_10)/36525._10,Nutat2(1,1),
     .          Nutat2(2,1))
            else
               call IAU2000A(T-0.5_10,Nutat2(1,1),Nutat2(2,1))
            endif
         else

c move one tabular point of moon,nutation,libration from /bdydta/
            np = iabs(jf - Jdbd(1))*2 + 1
            do i = 1, 6
               Moon2(i,1) = Mon(i,np)
            end do
            do i = 1, 3
               Librt2(i,1) = Librt(np,i)
            end do
            Nutat2(1,1) = Psid(np)
            Nutat2(2,1) = Epsd(np)
         endif
         do while( .true. )
 
c backspace and read buffer
            backspace Ibuf
            read(Ibuf) je,fsb1,ivl,mvl,
     .                 ((Merc1(i,j),i=1,ivl),j = 1,10),
     .                 (((Venpl1(i,j,k),i=1,ivl),j=1,5),k = 1,nb),
     .                 ((Moon1(i,j),i=1,mvl),j = 1,40),Nutat1,
     .                 Librt1
            backspace Ibuf
 
c move coordinates
            do i = 1, ivl
               do k = 2, 10
                  kk = 10 - (k - 2)
                  Merc2(i,k) = Merc1(i,kk)
               end do
               do k = 2, 5
                  kk = 5 - (k - 2)
                  do j = 1, nb
                     Venpl2(i,k,j) = Venpl1(i,kk,j)
                  end do
               end do
            end do
            do k = 2, 40
               kk = 40 - (k - 2)
               do i = 1, mvl
                  Moon2(i,k) = Moon1(i,kk)
               end do
 
c move nutation
               do i = 1, 2
                  Nutat2(i,k) = Nutat1(i,kk)
               end do
 
c move libration
               do i = 1, 3
                  Librt2(i,k) = Librt1(i,kk)
               end do
            end do
 
c write output n-body tape
            write(Ibody) jf,fsb,ivl,mvl,
     .                   ((Merc2(i,j),i=1,ivl),j = 1,10),
     .                   (((Venpl2(i,j,k),i=1,ivl),j=1,5),k = 1,nb),
     .                   ((Moon2(i,j),i=1,mvl),j = 1,40),Nutat2,
     .                   Librt2
 
c move one tabular point
            do i = 1, ivl
               Merc2(i,1) = Merc1(i,1)
               do k = 1, nb
                  Venpl2(i,1,k) = Venpl1(i,1,k)
               end do
            end do
            do i = 1, mvl
               Moon2(i,1) = Moon1(i,1)
            end do
            do i = 1, 2
               Nutat2(i,1) = Nutat1(i,1)
            end do
            do i = 1, 3
               Librt2(i,1) = Librt1(i,1)
            end do
 
c test for end
            jf  = je
            fsb = fsb1
            if(fsb.eq.0.0_10) then
               if(je.eq.Jdbdy0) then
                  rewind Ibuf
                  goto 400
               endif
            endif
         end do
      else if(Mtab.ge.10 .and. Ntabr.ge.40) then
         Ntab = 0
         Mtab = 0
         if(Jmoon.gt.0) Ntabr=0
         jf   = T
         if(Iboth.ge.0) then
            if(Nsign*(Jdbdy1-jf).le.0) Iboth = 1
         endif
         if(Kbdy(39).lt.0) return
         juldat = T - Bint
         je     = juldat
         fsb    = juldat - je
c
c write n-body data set
         ibodys = Ibody
         if(Iboth.ge.0) ibodys = Ibuf
         write(ibodys) je,fsb,ivl,mvl,
     .                 ((Merc1(i,j),i=1,ivl),j = 1,10),
     .                 (((Venpl1(i,j,k),i=1,ivl),j=1,5),k = 1,nb),
     .                 ((Moon1(i,j),i=1,mvl),j = 1,40),Nutat1,
     .                 Librt1
      endif
 
      return
      end
