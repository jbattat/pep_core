      subroutine ASTRD1(lice, nplnt, kp)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ints, j, lice
 
c*** end of declarations inserted by spag
 
 
c  J.F.Chandler - 1976 July
c  Read 1st 5 records of asteroid/satellite tape,
c  2 header, 3 data records.
c  If LICE = 0 print out header stuff.
c  NPLNT is integrated body, ignore it, even if on the tape.
c  KP is array of control integers for  ast/sat group.
c
      integer*2 nplnt, kp(30)
c
c common
      include 'b2dtaint.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
 
c temporary storage
      common/YVECT/ Bmas(12),Relft(12),Heada(32),Name(6,12),Jd0(12),
     . Jvlb2,Itrt,Npg,Intr,Epsb2,Ncp(12),Kb2(40),Kkb2(80)
      character*4 Name,Heada
      real*10 Bmas,Relft
      real*4    Epsb2
      integer*4 Jd0,Itrt,Npg,Jvlb2,Intr
      integer*2 Ncp,Kb2,Kkb2
 
      Nast = 0
      if(Jpert.le.0) return
      Itrwnd(Jpert) = 1
      read(Jpert,end=500,err=100) Heada
      goto 300
  100 read(Jpert) Heada
      write(Iout,200) Jpert
  200 format(
     .'  *** ERROR IGNORED ON FIRST RECORD OF PERTURBING SATELLITE DATAS
     .ET', i3, ' ***')
      Line = Line + 2
  300 read(Jpert,end=500,err=600) Nast, (Np2(i),i=1,Nast),
     .     (Ncp(i),i=1,Nast), (Inb2(i),i=1,Nast), Jdb21, Jdb22,
     .     (Jd0(i),i=1,Nast), ((B2ta(j,i),j=1,6),i=1,Nast),
     .     ((Name(j,i),j=1,6),i=1,Nast), Na1, Na4, Na8, Intr,
     .     Jvlb2, Epsb2, Kb2, Itrt, Npg, (Bmas(i), i=1, Nast),
     .     (Relft(i), i=1, Nast), Kkb2
      Intb2 = Intr
 
c set up bit pattern to select non-empty classes (1, 4, or 8)
      Naf = 0
      if(Na1.gt.0) Naf = 1
      if(Na4.gt.0) Naf = Naf + 2
      if(Na8.gt.0) Naf = Naf + 4
      Ib2sgn = ISIGN(1, Jdb22 - Jdb21)
      Ddir   = Ib2sgn
      Fstep  = Intb2
      if(Intb2.lt.1) Fstep = 2._10**Intb2
      if(lice.le.0) then
         call NEWPG
         write(Iout, 350) Jpert, Heada
  350    format('0  HEADER INFO FROM PERTURBING SATELLITE DATA SET', i3/
     .          1x, 32A4)
         do i = 1, Nast
            Jd0(i) = Jd0(i) - 1
            end do
         write(Iout, 400) Ib2sgn, Jdb21, Jdb22, Na1, Na4, Na8, Intb2
  400    format(' IDIR=', i2, '  JDT1=', i8, '  JDT2=', i8, '   N1=',
     .          i2, '   N4=', i2, '   N8=', i2, '   INT=', i3,
     .          '.  INITIAL CONDITIONS ARE')
         write(Iout, 450) (Np2(j), (Name(i,j),i=1,6), Jd0(j),
     .                    (B2ta(i,j),i=1,6), j=1, Nast)
  450    format(/i4, 1x, 6A4, ' JD0=', i8, '.5  A=', 1pd22.15, 4x,
     .          'E=', 1pd22.15, 5x, 'INC=', 1pd22.15/44x, 'ASC=',
     .          1pd22.15, '  PER=', 1pd22.15, 3x, 'ANOM0=', 1pd22.15)
         Line = Line + 3 + 3*Nast
      endif
      if(Kkb2(70).ne.Jct(13)) call SUICID('S-BODY HAS WRONG REFERENCE FR
     .AME, STOP IN ASTRD1',12)
      goto 700
 
c error stops
  500 call SUICID(
     .' END ENCOUNTERED IN HEADERS OF PERTURBING ASTEROID DATA SET, STOP
     . IN ASTRD1 ', 19)
  600 call SUICID(
     .' REDUNDANCY ENCOUNTERED IN SECOND RECORD OF PERTURBING ASTEROID D
     .ATA SET, STOP IN ASTRD1', 22)
 
c setup computations
  700 ints   = Fstep*32._10 - 1._10
      Jd0(1) = max0(Jdb21, Jdb22) - ints
      Jdb21  = min0(Jdb21, Jdb22) + ints
      Jdb22  = Jd0(1)
      D1(1)  = Fstep*Ddir
      D4(1)  = D1(1)*4._10
      D8(1)  = D1(1)*8._10
      D1(2)  = D1(1)**2
      D4(2)  = D4(1)**2
      D8(2)  = D8(1)**2
      do i = 1, Nast
         Nas(i) = 0
         j = Np2(i)
         Kpa(i) = kp(j)
         if(nplnt.eq.j) Kpa(i) = -1
         end do
      call RDATAP(0)
      return
      end
