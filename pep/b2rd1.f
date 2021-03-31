      subroutine B2RD1(lice)

      implicit none

c j.f.chandler  1977 june        subroutine b2rd1
c read first 5 records of small-body tape

c arguments
      integer*2 lice

c lice = 0 print out data on first two records
c lice = 1 do not print out

c array dimensions
      include 'globdefs.inc'

c common
      include 'b2dta.inc'
      include 'comdat.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'

c temporary storage
      common/YVECT/ Tpl(12),Head(16),Name(6,12),Kb2(40),
     . Bmas(12),Cmas(12),Relft(12),Fr0(12),Itrt,Npg,Kkb2(80),
     . Jvlb2,Epsb2,Jdb20(12),Nc2(12)
      character*8 Head
      character*4 Name
      real*10 Tpl,Bmas,Cmas,Relft,Fr0
      real*4    Epsb2
      integer*4 Itrt,Npg,Jdb20,Jvlb2
      integer*2 Kb2,Kkb2,Nc2

c local
      real*10 amarg,b2int
      integer*4 i,imarg,intr,j,m1,n2rec

c for now . . . there is no second n-body integration
      integer*4 i2bod
      equivalence (i2bod,Jpert)

      character*8 eend/' * * END'/, eerr/' * ERROR'/

c
c read first two records
      n2rec = 1
      read(i2bod,err=100,end=700) Head
      goto 300
  100 write(Iout,200) i2bod
  200 format('0 * * * ERROR ON RECORD 1 OF S-BODY DATA SET',i3,
     .       ' * * *')
      Line = Line + 2
      read(i2bod) Head
  300 n2rec = 2
      read(i2bod,err=600,end=700) Nast,(Np2(i),i = 1,Nast),
     .     (Nc2(i),i = 1,Nast),(Inb2(i),i = 1,Nast),Jdb21,Jdb22,
     .     (Jdb20(i),i = 1,Nast),((B2ta(i,j),i=1,6),j = 1,Nast),
     .     ((Name(i,j),i=1,6),j = 1,Nast),Na1,Na4,Na8,intr,
     .     Jvlb2,Epsb2,Kb2,Itrt,Npg,(Bmas(i),i = 1,Nast),
     .     (Relft(i),i = 1,Nast),Kkb2,(Icn2(i),i = 1,Nast),
     .     (Cmas(i),i = 1,Nast),(Fr0(i),i = 1,Nast)
      Intb2 = intr

      if(Nast.le.0) goto 390
c
c write out data on first two records
      if(lice.le.0) then
         call NEWPG
         write(Iout,320) Nast,i2bod,Itrt,Npg,Head
  320    format('0INFORMATION FROM FIRST TWO RECORDS OF',i3,
     .    '-BODY DATA SET (S-BODY)',i3,' WITH ITERAT=',i3,
     .    ' PAGE=',i5,' AND TITLE:'/ 1x,16A8)
         write(Iout,340) Jdb21,Jdb22,Na1,Na4,Na8,Intb2
  340    format(' JD1=',i7,'  JD2=',i7,' (NA1=',i2,' NA4=',i2,
     .    ' NA8=',i2,') INT=',i3,'  INITIAL CONDITIONS ARE')
         Line = Line + 4
         do i = 1,Nast
            Tb20(i) = Jdb20(i) + Fr0(i)
            Tpl(i)  = Tb20(i) - 0.5_10
            end do
         write(Iout,360) (Np2(j),(Name(i,j),i=1,6),Tpl(j),
     .       (B2ta(i,j),i=1,3),Nc2(j),Icn2(j),
     .       (B2ta(i,j),i=4,6),Bmas(j),Cmas(j),j = 1,Nast)
  360    format('0',i3,1x,6A4,' JD0=',0pf11.3,' A=',1pd22.15,
     .    4x,'E=',d22.15,5x,'INC=',d22.15/ 5x,'NCP=',i2,
     .    ' ICND=',i2,25x,'ASC=',d22.15,'  PER=',d22.15,
     .    3x,'ANOM0=',d22.15/ 9x,'MASS=',d22.15,
     .    '  CMASS=',d22.15)
         write(Iout,380) (i,Kb2(i),i = 1,40)
  380    format(10('   K(',i2,')=',i4))
         Line = Line + 4*Nast + 4
      endif

c check for valid numbers of bodies of each type
      if(Na1.gt.6 .or. Na4.gt.8 .or. Na8.gt.8 .or. Nast.gt.12) goto 900
      if(Na8.gt.0 .and. (Na4.gt.4 .or. Na1.gt.5)) goto 900
      if(Na4.gt.0 .and. Na1.gt.4) goto 900

      if(Kkb2(70).ne.Jct(13)) call SUICID('S-BODY HAS WRONG REFERENCE FR
     .AME, STOP IN B2RD1 ',12)

c calculate interval quantities
      b2int = Intb2
      if(Intb2.le.0) b2int = 2._10**Intb2
      Naf = 0
      if(Na1.gt.0) Naf = 1
      if(Na4.gt.0) Naf = Naf + 2
      if(Na8.gt.0) Naf = Naf + 4
      amarg = b2int
      if(Na4.gt.0) amarg = b2int*4._10
      if(Na8.gt.0) amarg = b2int*8._10
      imarg  = amarg*5._10
      imarg  = imarg + 1
      Ib2sgn = ISIGN(1,Jdb22 - Jdb21)
      m1     = Jdb21
      Jdb21  = min0(Jdb21,Jdb22) + imarg
      Jdb22  = max0(m1,Jdb22) - imarg
      D(1)   = b2int*Ib2sgn
      D(2)   = D(1)*4._10
      D(3)   = D(1)*8._10
      B2int5 = D(3)*5._10

c compare bodies+centers with input list
c
c read first three data records
      call B2RED1(0)
      return

c
c error messages
  390 write(Iout,400) i2bod,Nast
  400 format(' S-BODY DATA SET',i3,' HAS NAST=',i3)
      goto 1100
  600 write(Iout,800) eerr,n2rec,i2bod
      goto 1100
  700 write(Iout,800) eend,n2rec,i2bod
  800 format('0* *',a8,' ON RECORD',i2,' OF S-BODY DATA SET',i3,' * *')
      goto 1100
  900 write(Iout,1000) i2bod
 1000 format('0S-BODY DATA SET',i3,' HAS TOO MANY BODIES')

 1100 call SUICID('ERRORS ON S-BODY DATA SET, STOP IN B2RD1',10)

      end
