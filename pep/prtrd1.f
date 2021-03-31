      subroutine PRTRD1(lice,nplnt,ncentr,kp,lvel,lnut,llib)

      implicit none
c
c M.E.Ash   Aug 1968   Subroutine PRTRD1 (n-body format)
c First five records of perturbing planet tape are read.
c
c   LICE =0 printout of first two records of pert.planet tape
c   LICE =1 no printout of first two records
c   NPLNT is the planet to be integrated about central body NCENTR.
c   Position always calculated for body j if KP(j).ge.0.
c   In addition, velocity also calculated if LVEL=1 for all such
c   j.  If LVEL=-1, this velocity only calculated for central
c   body j=NCENTR.
c
c arguments
      integer*4 lice,lvel,lnut,llib
      integer*2 kp(10),nplnt,ncentr

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'bdctrl.inc'
      include 'bddtaint.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'param.inc'
      include 'prtpin.inc'

c temporary storage
      common/YVECT/ Mass1(10),Relftb(10),Head(32),Name(6,10),Je0(10),
     .  Kbd(40),Epsbd,Jvlbd,Itrt,Npg,Intbd,Kkbd(80)
      character*4 Head,Name
      real*10 Mass1,Relftb
      real*4    Epsbd
      integer*4 Je0,Jvlbd,Itrt,Npg
      integer*2 Intbd,Kbd,Kkbd

      integer   i,j,rec
      character*24 nbrown/'  MOON  (BROWN MEAN THR)'/

      Ntab1 = 1
      Mtab1 = 2

      if(Ipert.gt.0) then

c read first record
         Itrwnd(Ipert) = 1
         read(Ipert,end=500,err=100) Head
         goto 300
      endif
      if(Ict(50).eq.0) call SUICID(
     .       'NO PERTURBING BODY TAPE, STOP AT PRTRD1 ',10)

c set up arrays as if n-body tape were being read
      Nbdy     = 1
      Nplbd(1) = 10
      Ncpbd(1) = 3
      Intb(1)  = -1
      Jdbd0(1) = 0
      Jdbd1    = 0
      Jdbd2    = 9999999
      call ZFILL(Betabd(1,1),16*6)
      Nmoon  = -1
      Ibdsgn = 1
      Dir    = 1._10
      call MVC(nbrown,1,24,Nammon,1)
      goto 800

  100 read(Ipert) Head
      write(Iout,200) Ipert
  200 format(
     .' **** ERROR IGNORED ON FIRST RECORD OF PERTURBING PLANET DATA SET
     .',i3,' ****')
      Line = Line + 1
  300 read(Ipert,end=500,err=600) Nbdy,(Nplbd(i),i=1,Nbdy),
     .     (Ncpbd(i),i=1,Nbdy),(Intb(i),i=1,Nbdy),Jdbd1,Jdbd2,
     .     (Jdbd0(i),i=1,Nbdy),((Betabd(i,j),i=1,6),j=1,Nbdy),
     .     ((Name(i,j),i=1,6),j=1,Nbdy),Nmoon,Nbdy1,Intbd,
     .     Jvlbd,Epsbd,Kbd,Itrt,Npg,(Mass1(i),i=1,Nbdy),
     .     (Relftb(i),i=1,Nbdy),Kkbd
      Ibdsgn = ISIGN(1,Jdbd2 - Jdbd1)
      call MVC(Name(1,Nbdy),1,24,Nammon,1)
      if(lice.le.0) then
c
c printout first two id records if lice=0
         call NEWPG
         write(Iout,350) Ipert,Head
  350    format(
     .'0   INFORMATION FROM FIRST TWO RECORDS OF PERTURBING PLANET PERIP
     .HERAL DATA SET',i3/ 1x,32A4)
         do i = 1,10
            Je0(i) = Jdbd0(i) - 1
         end do
         write(Iout,400) Ibdsgn,Jdbd1,Jdbd2,Nmoon
  400    format(' IDIR=',i2,' JDT1=',i8,' JDT2=',i8,' NMOON=',
     .          i2,' INITIAL CONDITIONS ARE')
         write(Iout,450) ((Name(i,j),i=1,6),Je0(j),(Betabd(i,j),i=1,
     .                    6),j = 1,Nbdy)
  450    format('0',6A4,' JD0=',i7,'.5  A =',1pd22.15,'  E =',
     .          1pd22.15,'   INC=',1pd22.15/ 40x,'ASC=',1pd22.15,
     .          ' PER=',1pd22.15,' ANOM0=',1pd22.15)
         Line = Line + 33
      endif
      if(Kkbd(70).ne.Jct(13)) call SUICID('N-BODY HAS WRONG REFERENCE FR
     .AME, STOP IN PRTRD1',12)
      goto 700
c
c error stops
  500 call SUICID(
     .' END ENCOUNTERED READING FIRST TWO RECORDS OF PERTURBING PLANET D
     .ATA SET, STOP IN PRTRD1',22)
  600 call SUICID(
     .' REDUNDENCY ENCOUNTERED READING SECOND RECORD OF PERTURBING PLANE
     .T DATA SET, STOP IN PRTRD1 ',23)
c
c setup computations
  700 if(Jdbdy0.lt.0 .and. Ipert.eq.Ibody) Jdbd2 = Jdbdy2
      Je0(1) = max0(Jdbd1,Jdbd2) - 20
      Jdbd1  = min0(Jdbd1,Jdbd2) + 20
      Jdbd2  = Je0(1)
      Dir    = Ibdsgn
      if(Ibdsgn.gt.0) then
         Jdbd1 = Jdbd1 - 1
      else
         Jdbd2 = Jdbd2 + 1
      endif
  800 Dmerc(1) = 2.0_10*Dir
      Dbody(1) = 4.0_10*Dir
      Dmoon(1) = 0.5_10*Dir
      Dmerc(2) = Dmerc(1)**2
      Dbody(2) = Dbody(1)**2
      Dmoon(2) = Dmoon(1)**2
      if(Nmoon.gt.0) then
         if(Mdstau.le.0._10)
     .        call SUICID(' NMOON=1, MDSTAU=0._10, STOP IN PRTRD1   ',
     .       10)
      endif
      rec=0
      do i=1,3
         Recbeg(i)=rec+1
         rec=rec+40
         Recend(i)=rec
      end do
      do i = 1,10
         Kiss(i) = kp(i)
         if(nplnt.eq.i) Kiss(i) = -1
      end do
      Kiss(11) = lnut
      Kiss(12) = llib
      Ler    = lvel
      Lps    = ncentr
      Limvel = 3
      if(Ler.gt.0) Limvel = 3*(Ler + 1)
      Limnut = lnut
      Limlib = llib
      Mascnt = 1._10
      do i = 1,9
         if(Kp(i).ge.0 .or. i.eq.nplnt .or. i.eq.ncentr)
     .    Mascnt = Mascnt + Mass(i)
      end do
c
c read first three data records
      call RDYCAL(0)
      return
      end
