      subroutine BDRD1(lice)

      implicit none

c
c M.E.Ash   Sept 1968    Subroutine BDRD1
c First five records of n-body tape are read.
c Also, asteroid center-of-mass data set, if any.
c
      integer*2 lice
c LICE =0 printout of data on first two records of n-body tape
c LICE =1 no such printout
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'bdctrl.inc'
      include 'bddta.inc'
      include 'comdat.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'namtim.inc'
      include 'param.inc'
c
c temporary storage
      common/YVECT/ Mass1(10), Relftb(10), Name(6,10), Kbd(40),
     .      Title(32), Epsbd, Jvlbd, Itrt, Npg, Nmoon, Intbd, Kkbd(80)
      character*4 Name, Title
      real*4    Epsbd
      integer*4 Jvlbd, Itrt, Npg, Nmoon
      integer*2 Intbd, Kbd, Kkbd
      real*10 Mass1, Relftb

c local variables
      real*10 dum
      integer*4 i,icmsgn,idatin,j,k,kk,levl,m1,nbcom,nbdrec
      character*56 messag
      character*1 blank/' '/
      integer*2 kcom70
c
c initialize conglomerate asteroid mass
      Sumcom = 0._10

c read first two records of n-body peripheral data set
      nbdrec = 1
      idatin = Ibody
      read(Ibody, err=100, end=900) Title
      goto 300
  100 write(Iout, 200) Ibody
  200 format('0**** ERROR ON RECORD 1 OF N-BODY DATA SET', i3, ' ****')
      Line = Line + 2
      read(Ibody) Title
c if system changes this second read of first record might not be
c needed to get data into storage in case of error
  300 nbdrec = 2
      read(Ibody, err=700, end=900) Nbdy, (Npl(i), i=1,Nbdy),
     .     (Ncp(i), i=1,Nbdy), (Intb(i), i=1,Nbdy), Jdbd1, Jdbd2,
     .     (Jdbd0(i), i=1,Nbdy), ((Beta(i,j),i=1,6), j=1,Nbdy),
     .     ((Name(i,j),i=1,6), j=1,Nbdy), Nmoon, Nbdy1, Intbd,
     .     Jvlbd, Epsbd, Kbd, Itrt, Npg, (Mass1(i), i=1,Nbdy),
     .     (Relftb(i), i=1,Nbdy), Kkbd

      if(Nbdy.le.0) then
         write(messag, 550) Ibody, Nbdy
  550    format(' N-BODY DATA SET', i3, ' HAS NBDY=', i3,
     .          ', STOP IN BDRD1 ')
         call SUICID(messag, 12)
      endif
c
c print out first two records of n-body peripheral data set
      if(lice.le.0) then
         call NEWPG
         write(Iout, 320) Nbdy, Ibody, Itrt, Npg, Title
  320    format('0INFORMATION FROM FIRST TWO RECORDS OF', i3,
     .          '-BODY DATA SET', i3, ' PRODUCED ON ITERATION', i3,
     .          ' PAGE', i5, ' OF RUN WITH TITLE'/1x, 32A4)
         write(Iout, 340) Jdbd1, Jdbd2, Nmoon, Nbdy1, Intbd, Jvlbd,
     .                    Epsbd
  340    format(' JDBD1=', i7, '  JDBD2=', i7, '  NMOON=', i2,
     .          ' NBDY1=', i2, '   INTBD=', i2, '  JVLBD=', i4,
     .          '  EPSBD=', 1pe12.5)
         Line = Line + 4
         do j = 1, Nbdy
            write(Iout, 350) j, (Name(i,j), i=1,6), Npl(j),
     .                       Ncp(j), Intb(j), Jdbd0(j),
     .                       (Beta(i,j), i=1,6)
  350       format(i3, '. ', 6A4, 3x, 'NPL=', i2, 3x, 'NCP=', i2,
     .             3x, 'INT=', i2, 3x, 'JD0=', i7, 6x, 'A =',
     .             1pd22.15, 5x, ' E =', 1pd22.15/13x, 'INC=',
     .             1pd22.15, 5x, 'ASC=', 1pd22.15, 5x, 'PER=',
     .             1pd22.15, 4x, 'ANOM=', 1pd22.15)
            Line = Line + 2
         end do
         write(Iout, 360) (i, Kbd(i), i=1,40)
  360    format(10('   K(',i2,')=',i4))
         write(Iout, 380) (blank, i, Mass1(i), i=1,Nbdy)
  380    format(4(a1,'MASS1(',i2,')=',1pd22.15))
         write(Iout, 400) (blank, i, Relftb(i), i=1,Nbdy)
  400    format(4(a1,'RELFT(',i2,')=',1pd22.15))
         Line = Line + 4 + 2*(1 + (Nbdy-1)/4)
      endif
      if(Kkbd(70).ne.Jct(13)) call SUICID('N-BODY HAS WRONG REFERENCE FR
     .AME, STOP IN BDRD1 ',12)
c
c read asteroid center-of-mass data set headers
      if(Kpert.gt.0) then
         nbdrec = 1
         idatin = Kpert
         read(Kpert, end=900) (Title(i), i=1,22), levl
         nbdrec = 2
         read(Kpert, end=900) Jdcom1, Jdcom2, Intcom, nbcom, Sumcom,
     .        kcom70
         if(lice.le.0) then
            write(Iout, 410) Kpert, (Title(i), i=1,22), levl,
     .                       Jdcom1, Jdcom2, Intcom, nbcom, Sumcom
  410       format('-INFORMATION FROM CENTER-OF-MASS DATA SET', i3,
     .             ' WITH'/' TITLE= ', 2A4, 1x, 18A4, 1x, 2A4,
     .             ' LEVEL=', a4/' JDCOM1=', i7, '  JDCOM2=', i7,
     .             '  INTCOM=', i2, '  NBCOM=', i4, '  SUMCOM=',
     .             1pd15.7)
            Line = Line + 5
         endif
         if(kcom70.ne.Jct(13)) call SUICID(
     .    'AST.COM HAS WRONG REFERENCE FRAME, STOP IN BDRD1',12)
      endif
c
c calculate interval quantities
      if(Jdbdy9.lt.0) Jdbd2 = Jdbdy2
      Ibdsgn = ISIGN(1, Jdbd2 - Jdbd1)
      m1     = Jdbd1
      Jdbd1  = min0(Jdbd1, Jdbd2) + 20
      Jdbd2  = max0(m1, Jdbd2) - 20
      Nbdy2  = Nbdy - 2
      Intbd5 = 5*Intb(2)
c
c setup for center of mass of solar system computations
c
c set up individual elliptic approximations for partials
      do i = 1, Nbdy1
         call JNITL(Gauss*SQRT(1._10+Mass(i)), Beta(1,i), Setpss(1,i),
     .              0, dum)
      end do
      Jdcnts = 0
      Jdcnt1 = 0
      Jdcnt2 = 0
      Jdent1 = 0
      Jdent2 = 0
      Nmmm   = 0
c total mass computation moved to PLCMS

c extra setup if asteroid center-of-mass included
      if(Kpert.gt.0) then
         icmsgn = ISIGN(1, Jdcom2 - Jdcom1)
         m1     = Jdcom1
         Jdcom1 = min0(Jdcom1, Jdcom2) + Intbd5
         Jdcom2 = max0(m1, Jdcom2) - Intbd5

c check consistency
         if(Intcom.ne.Intb(2) .or. icmsgn.ne.Ibdsgn) call SUICID(
     .'INCONSISTENCY BETWEEN N-BODY AND CENTER-OF-MASS DATA SETS, STOP I
     .N BDRD1', 18)
      endif
c
c see if there is reintegration of n-bodies (if possible)
      if(Nplbdy(1).eq.0) goto 500
      do k = 1, Nbdy1
         kk = -3
         if(Npl(k).ne.3) then
            kk = -2
            if(Npl(k).ne.10) then
               do kk = 1, Numpln
                  if(Npl(k).eq.Nplnt(kk)) goto 420
               end do
               goto 450
            endif
         endif
  420    do i = 1, 6
            if(Lpl(i,kk).gt.0) then
               Jdbdy0 = Jdbd0(k)
               goto 500
            endif
         end do
  450 end do
c
c read first three data records of n-body data set
  500 call BDRED1(0)
c
c check for consistency of moon distance unit
      if(Nmoon.gt.0) then
         if(prmter(99).le.0.0_10) call SUICID(
     .' NMOON=1 ON N-BODY DATA SET, PRMTER(99)=0 INSTEAD OF 0.021275213,
     . STOP IN BDRD1 ', 20)
      endif

      return

c
c error messages


  700 write(messag, 800) nbdrec, Ibody
  800 format(' ERROR ON RECORD', i2, ' OF N-BODY DATA SET', i3,
     .       ', STOP IN BDRD1 ')
      call SUICID(messag, 14)

  900 write(messag, 1000) nbdrec, idatin
 1000 format(' END ON RECORD', i2, ' OF DATA SET', i3,
     .       ', STOP IN BDRD1  ')
      call SUICID(messag, 12)

      end
