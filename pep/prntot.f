      subroutine PRNTOT(noprnt,iseq,nstop)
 
      implicit none

c ash/forni  october 1966    subroutine prntot
c revised oct 1969   march 1970  may 1973
c reno'd 1978 april
c printout of input data

c arguments
      integer*4 noprnt,iseq,nstop

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'aprtbf.inc'
      integer*4 ibufa(5),ibufb(5)
      equivalence (ibufa,Ibuf1),(ibufb,Ibuf6)
      include 'bdctrl.inc'
      include 'dtparm.inc'
      include 'empcnd.inc'
      include 'ethhar.inc'
      include 'fcntrl.inc'
      include 'france.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'monhar.inc'
      include 'namtim.inc'
      include 'nrmwgt.inc'
      include 'obsdta.inc'
      include 'param.inc'
      include 'plndta.inc'
      include 'plnhar.inc'
 
c locally nshape called nshp which is documented external name
      integer*4 nshp(4)
      equivalence (Nshape,nshp)
      include 'scoef4.inc'
      include 'smlbdy.inc'
      include 'timstf.inc'
 
      real*10 ovlyfs(4,122)
      real*4    transf(u_stdsz,4,1000/u_stdsz)
      integer*4 latdim(4),londim(4),ngdpts(4)
      integer*2 lshp4(4,122)
      equivalence (Pzhar,transf,ovlyfs),
     .            (Lpzhar,lshp4),(Scontl(1,7),latdim),
     .            (Scontl(1,8),londim),(Scontl(1,9),ngdpts)
c equivalencing for fourier series shape coefficients
c comments in bdyred,bodred,bodntl
      real*10 pa(4,20),pb(4,20),pc(4,20),pd(4,20),pap(4,20),
     .          pcp(4,20),papp(4),pbpp(4)
      integer*2 lpa(4,20),lpb(4,20),lpc(4,20),lpd(4,20),
     .          lpap(4,20),lpcp(4,20),lpapp(4),lpbpp(4)
      equivalence (ovlyfs(1,1),pa),(ovlyfs(1,21),pb),
     .            (ovlyfs(1,41),pc),(ovlyfs(1,61),pd),
     .            (ovlyfs(1,81),pap),(ovlyfs(1,101),pcp),
     .            (ovlyfs(1,121),papp),(ovlyfs(1,122),pbpp)
      equivalence (lshp4(1,1),lpa),(lshp4(1,21),lpb),
     .            (lshp4(1,41),lpc),(lshp4(1,61),lpd),
     .            (lshp4(1,81),lpap),(lshp4(1,101),lpcp),
     .            (lshp4(1,121),lpapp),(lshp4(1,122),lpbpp)
 
c shared external work area
      common /WRKCOM/ Grid(1000),Imass(30),
     . Idt((u_mxtes*(u_mxtes+1))/2-1),Idi((u_mxtes*(u_mxtes+1))/2-1)
      real*4 Grid
      real*10 Imass
      integer*4 Idi,Idt

c local
      integer*2 idex(9)/9*0/,jdex(9)/9*0/,nelem(2)
      character*8 sinlit/'  SINE  '/,coslit/' COSINE '/
      character*8 teslit/'TESSERAL'/,zonlit/'  ZONAL '/,bl8/'        '/
      character*2 bl,blank
      equivalence (bl8,bl,blank)
      real*4    sec
      integer*2 ncsp(4)/0,3,3,10/
      character*8 earth/' EARTH  '/
      character*2 symsp(5)/'EM','MN','ER','MR','PL'/
 
c literal data for planet
      character*1 xx1/'J'/,xx2/'C'/,xx3/'S'/
      character*4 xshp(8)/'FA( ', 'FB( ', 'FC( ', 'FD( ', 'FAP(',
     .                    'FCP(', 'FAPP', 'FBPP'/
      character*96 harmsg/'                GRAVITATIONAL POTENTIAL HARMO
     .NIC COEFFICIENTS -          (NPLNT=..., NTESS=...) '/
      character*8 pname
      character*133 lineout,linealt

      real*10 diff,frp,frq
      integer*4 i,icount,igrid,ihr,ilat,ilon,imin,
     . intdy,intp22,ithing,j,k,k5,klam,l,lgrid,linc,ll,lll,
     . ltest,mstop,n,n1,n1a,n1b,n1bs,n2,ncentr,ncount,ndt1,
     . nki,nlat,nlon,noprtp,npall,npt,npln,npz
      character*40 erms/'****** INPUT DATA ERRORS, STOP IN PRNTOT'/
c*       start=1000
c
      call PLOG9
c
c write title page for input link
      call PAGSET('PLANETARY EPHEMERIS PROGRAM INPUT DATA  ', 10)
      call NEWPG
      call CONVRT(ihr,imin,sec,Itimst)
      write(Iout,100) Idayr,ihr,imin,sec
  100 format('0DAY OF YEAR=', 1x, a3, '    CLOCK TIME AT EXECUTION=',
     .       i3, 'H', i3, 'M', f6.2, 'S'/'-'/'-'/'-'/)
      Line   = 13
      noprtp = noprnt
      if(Ict(1).eq.-99) noprnt = 0
      if(noprnt.le.0) then
         write(Iout,150)
  150    format(
     .34x,'**   **           **   **********    **       **   **********
     .**'/
     .34x,'**   ***          **   ***********   **       **   **********
     .**'/
     .34x,'**   ****         **   **       **   **       **        ** '/
     .34x,'**   ** **        **   **       **   **       **        ** '/
     .34x,'**   **  **       **   **       **   **       **        ** '/
     .34x,'**   **   **      **   **       **   **       **        ** '/
     .34x,'**   **    **     **   ***********   **       **        ** '/
     .34x,'**   **     **    **   **********    **       **        ** '/
     .34x,'**   **      **   **   **            **       **        ** '/
     .34x,'**   **       **  **   **            **       **        ** '/
     .34x,'**   **        ** **   **            **       **        ** '/
     .34x,'**   **         ****   **            **       **        ** '/
     .34x,'**   **          ***   **            ***********        ** '/
     .34x,'**   **           **   **             *********         ** ')
         call PLINK
         write(Iout,200)
  200    format('-'/
     .          '0SPECIFICATION OF THE UNITS OF MASS, LENGTH AND TIME'/
     .          '0MASS OF SUN', 35x, '= 1.0000000000000')
         if(Gauss.eq.0.017202098950000_10) then
            write(Iout,250) Gauss
  250       format(' SQRT(GRAVITATIONAL CONSTANT TIMES MASS OF SUN)=',
     .       f16.13, ' ((AU)**3/2)/(DAY)')
         else
            write(Iout,260) Gauss
  260       format(' SQRT(GRAVITATIONAL CONSTANT TIMES MASS OF SUN)=',
     .       f24.21, ' ((AU)**3/2)/(DAY)')
         endif
         write(Iout,300)
  300    format(
     .  ' INDEPENDENT TIME VARIABLE OF EPHEMERIDES      = A.1 + 32.15S')
c
c check consistency of control constants
         call NEWPG
      endif
      mstop = 0
      call CHECK(mstop)
      call PERCHK(mstop)
 
c check for inconsistent old-style dt setup
      if(Ict(1).ge.0 .and. Numdt.gt.0 .and. Numdt.le.80 .and.
     .    Jddt0.le.0) then
         ndt1 = Numdt + 1
         do i = ndt1, 100
            if(Ldt(i).ne.0 .or. Dt(i).ne.0.0) goto 400
         end do
         do i = 1, Numdt
            if(Ldt(i+100).ne.0 .or. Dt(i+100).ne.0.0) then
               write(Iout,310)
  310          format(
     .            ' ***** OBSOLETE DT SETUP: USE 1-200,201-400,401-600',
     .            ' NOT 1-100,101-200,201-300')
               mstop = mstop + 1
               if(Mout.gt.0) write(Mout,310)
               goto 400
            endif
         end do
      endif
 
  400 nstop = nstop + mstop
      if(noprnt.gt.0) goto 1990
c
c printout general control constants
      if(mstop.gt.0) call NEWPG
      ithing = 3
      write(Iout,450) In,Iout,Jout,Kout,Ipunch,Igraph,Intern,
     . Ibuf,Imat,Ipert,Jpert,Kpert,Ivcnt,Iplcon,Iobcon,Iobs,
     . ithing,Imat1,Imat2,Jmat,Imat3,Imat4,Ictat,Libhoc
  450 format('0PERIPHERAL DATA SET NUMBERS NOT BELONGING TO A SPECIFIC',
     . ' BODY OR OBSERVATION LIBRARY'/
     . '      IN =', i3, '    IOUT =', i3,
     .'    JOUT =', i3, '    KOUT =', i3, '   IPUNCH=', i3,
     .'   IGRAPH=', i3, '   INTERN=', i3, '    IBUF =', i3,
     .'    IMAT =', i3, '    IPERT=', i3/'    JPERT=', i3, '    KPERT=',
     .i3, '    IVCNT=', i3, '   IPLCON=', i3, '   IOBCON=', i3,
     .'    IOBS =', i3, '   ITHING=', i3, '    IMAT1=', i3,
     .'    IMAT2=', i3, '    JMAT =', i3/'    IMAT3=', i3, '    IMAT4=',
     .i3, '    ICTAT=', i3,'   LIBHOC=', i3)
      write(Iout,500) Nummt0,(Imat0(i),i = 1,Nummt0)
  500 format(' NUMMT0=', i3, '   IMAT0=', 20I4)
      Line = Line + 6
      write(Iout,550) Nparam,(i,Ict(i),i = 1,80)
  550 format(
     . '0CONTROL CONSTANTS FOR PROGRAM FLOW AND LEAST-SQUARES ANALYSIS',
     . 5x, 'NPARAM=', i4/(10('  ICT(',i2,')=',i3)))
      write(Iout,600) (i,Eps(i),i = 1,30)
  600 format(6('  EPS(',i2,')=',1pe12.5))
      Line = Line + 15
c
c printout parameters and corresponding least-squares control
c constants
      write(Iout,650) (i,Lprm(i),i = 1,u_nmprm)
  650 format(
     . '0CONTROL CONSTANTS FOR LEAST SQUARES PARAMETER ADJUSTMENT AND TH
     .E CORRESPONDING PARAMETER VALUES'/
     . (10(' LPRM(',i2,')=',i3)))
      Line = Line + 12
      write(Iout,700) (i,Mass(i),i = 1,28)
  700 format((4('  MASS(',i2,')=',1pd22.15)))
      write(Iout,750) (prmter(i),i = 29,44)
  750 format('  MASS(29)=', 1pd22.15, '  MASS(30)=', d22.15,
     .       ' RELFCT31 =', d22.15, ' GMVARY32 =',
     .       d22.15/' SUNHAR33 =', d22.15, ' PRMTR(34)=', d22.15,
     .       ' PRMTR(35)=', d22.15, ' PRMTR(36)=',
     .       d22.15/' PRMTR(37)=', d22.15, ' PRMTR(38)=', d22.15,
     .       ' SUNPOE39 =', d22.15, ' PRMTR(40)=',
     .       d22.15/' BETA  41 =', d22.15, ' GAMMA 42 =', d22.15,
     .       ' BETA'' 43 =', d22.15, ' GAMMA''44 =', d22.15)
      write(Iout,800) (i,prmter(i),i = 45,48)
  800 format(4(' PRMTR(',i2,')=',1pd22.15))
      write(Iout,850) (prmter(i),i = 49,56)
  850 format('  ASBA 49 =', 1pd22.15, ' ASBMAS50 =', d22.15,
     .       ' AULTSC51 =', d22.15, ' LTVARY52 =',
     .       d22.15/' RELDEL53 =', 1pd22.15, ' RELDOP54 =', d22.15,
     .       ' PRMTR(55)=', d22.15, ' PRMTR(56)=', d22.15)
      write(Iout,800) (i,prmter(i),i = 57,68)
      write(Iout,900) (i,prmter(i),i = 69,71),Ctvary
  900 format(3(' PRMTR(',i2,')=',1pd22.15), ' CTVARY72 =', d22.15)
      write(Iout,800) (i,prmter(i),i = 73,88)
      write(Iout,940) (prmter(i),i = 89,100)
  940 format(' PRTMR(89)=', 1pd22.15, ' PRMTR(90)=', 1pd22.15,
     .       ' ECINC 91 =', 1pd22.15, ' SEQINC92 =', 1pd22.15 /
     .       ' SEQASC93 =', 1pd22.15, ' SUNRAD94 =', 1pd22.15,
     .       ' PRMTR(95)=', 1pd22.15, ' PRMTR(96)=', 1pd22.15 /
     .       ' PRMTR(97)=', 1pd22.15, ' MDSTAU98 =', 1pd22.15,
     .       ' MDSTSC99 =', 1pd22.15, ' LTVEL100 =', 1pd22.15)
      Line = Line + 25
c
c printout inverse masses, constants in /aprtbf/, obs lib data
c set nos., nbody integration control constants
      call NEWPG
c
c printout inverse masses
      call PAGCHK(60,10,0)
      write(Iout,950)
  950 format('0INVERSE MASSES CALCULATED FROM PREVIOUS PAGE')
      do i = 1, 30
         Imass(i) = 0.0_10
         if(Mass(i).gt.0._10) Imass(i) = 1._10/Mass(i)
      end do
      write(Iout,1000) (blank,i,Imass(i),i = 1,30)
 1000 format(4(1A1,'IMASS(',i2,')=',1pd22.15))
c
c printout constants in /aprtbf/ labeled common
      call PAGCHK(60,9,0)
      write(Iout,1050) Ipert0,Ipert1,Ipert2,Jpert0,Jpert1,
     .                  Jpert2, Kpert0, Kpert1, Kpert2,
     .                  (i,ibufa(i),i = 1,5),
     .                  (i,ibufb(i-5),i = 6,10)
 1050 format(
     .'0ALTERNATE PERTURBING PLANET DATA SETS, BUFFERS, ACCURACY CONSTAN
     .TS'/3x, 'IPERT0=', i3, 3x, 'IPERT1=', i3, 3x, 'IPERT2=', i3, 3x,
     . 'JPERT0=', i3, 3x, 'JPERT1=', i3, 3x, 'JPERT2=', i3, 3x,
     . 'KPERT0=', i3, 3x, 'KPERT1=', i3, 3x, 'KPERT2=', i3/
     . 9(4x,'IBUF',i1,'=',i3), 3x, 'IBUF', i2, '=', i3)
      write(Iout,1100) (i,Epsa(i),i = 1,30)
 1100 format(6(' EPSA(',i2,')=',1pe12.5))
c
c printout additional control constants
      call PAGCHK(60,4,0)
      write(Iout,1140) Jpunch,Kpunch,Lout,Mout,Nout,Ieng,
     .                  Jeng, Keng, Typout, Extprc, noprtp, iseq
 1140 format(
     .'0ADDITIONAL PERIPHERAL DATA SET NUMBERS NOT BELONGING TO A SPECIF
     .IC BODY OR OBSERVATION LIBRARY + SWITCHES'/3x, 'JPUNCH=', i3, 3x,
     .'KPUNCH=', i3, 3x, ' LOUT =', i3, 3x, ' MOUT =', i3, 3x,
     .' NOUT =', i3, 3x, ' IENG =', i3, 3x, ' JENG =', i3, 3x,
     .' KENG =', i3/3x, 'TYPOUT=', i3, 3x, 'EXTPRC=', i3, 3x, 'NOPRNT=',
     .i3, 3x, ' ISEQ =', i3)
      call PAGCHK(60,12,0)
      write(Iout,1150) (i,Jct(i),i = 1,100)
 1150 format('0ADDITIONAL CONTROL CONSTANTS FOR PROGRAM FLOW AND LEAST-S
     .QUARES ANALYSIS'/(10('  JCT(',i2,')=',i3)))
c
c printout observation library data set numbers
      linc = 3 + (Numobt - 1)/3
      call PAGCHK(60,linc,0)
      write(Iout,1200) Numobt,
     .                  (blank,i,Iobs0(i),i,Iobs1(i),i,Iobs2
     .                  (i),i = 1,Numobt)
 1200 format(
     . '0OBSERVATION LIBRARY PERIPHERAL DATA SET NUMBERS (NUMOBT=', i2,
     . ')'/(3(a1,1x,'IOBS0(',i2,')=',i2,2x,'IOBS1(',i2,')=',i2,2x,
     . 'IOBS2(',i2,')=',i2)))
 
c printout normal equation weights
      linc = 3 + (Numobt - 1)/4 + (Nummt0 + 3)/4
      call PAGCHK(60,linc,0)
      write(Iout,1250) (bl,i,Wgtobs(i),i = 1,Numobt)
 1250 format('0WEIGHTS APPLIED TO NORMAL EQUATIONS FORMED FROM OBSE',
     .  'RVATION DATA SETS OR RESTORED FROM SAVED NORMAL EQUATIONS'/
     .  4(a1,'WGTOBS(',i1,')=',1pd22.15))
      if(Nummt0.gt.0) write(Iout,1300) (bl,i,Wgtmt0(i),i = 1,Nummt0)
 1300 format(4(a1,'WGTMT0(',i1,')=',1pd22.15))
c
c*  start=1500
c printout n-body integration control constants
      write(Iout,1350) Nbody,Ibody,Jdbdy1,Jdbdy0,Jdbdy2,
     .                  Jvlbdy, Epsbdy, Intbdy
      linc = 8 + (Nbody - 1)/9
      call PAGCHK(60,linc,0)
 1350 format('0INPUT CONTROL CONSTANTS FOR N-BODY INTEGRATION'/
     .       3x, 'NBODY =', i3, 3x, 'IBODY =', i3, 3x, 'JDBDY1=', i8,
     .       3x, 'JDBDY0=', i8, 3x, 'JDBDY2=', i8, 3x, 'JVLBDY=',
     .       i6, 3x, 'EPSBDY=', 1pe12.5, 3x, 'INTBDY=', i3)
      write(Iout,1400) (i,Kbdy(i),i = 1,40)
 1400 format(10(' KBDY(',i2,')=',i3))
      if(Nbody.le.0) then
         write(Iout,1420)
 1420    format(' THERE IS NO N-BODY INTEGRATION OR N-BODY TAPE')
      else
         write(Iout,1440) (blank,i,Nplbdy(i),i = 1,Nbody)
 1440    format(9(1A1,'NPLBDY(',i2,')=',i2))
      endif
c
c initialize index table for tesseral coefficients
      j = 0
      do n = 2, 50
         do i = 1, n
            Idt(j + i) = n
            Idi(j + i) = i
         end do
         j = j + n
      end do
 
      Line  = 60
      npall = 4 + Numpln
c
c printout planet constants
      klam = -1
      do 1800 k = 1, npall
         j = k - 4
         pname  = Aplnt(j)
         npln   = Nplnt(j)
         ncentr = Npcent(j)
         if(k.le.4) then
 
c special body or rotation (earth or moon)
            k5     = k
         else
 
c ordinary planet
            k5     = 5
         endif
 
         call PAGCHK(60,29,0)
         write(Iout,1460) pname,j,npln,Iplnt(j),symsp(k5),
     .                     Jd1(k),symsp(k5),Jdpl0(j),symsp(k5),
     .                     Jd2(k),Int(k),ncentr,
     .                     (symsp(k5),i,Lpl(i,j),i = 1,u_nmbod)
 1460    format('0', a8, 4x, 'NPLNT(', i2, ')=', i3, 4x, 'IPLNT=',
     .          i3, 9x, 'JD', a2, '1=', i8, 7x, 'JD', a2, '0=', i8,
     .          7x, 'JD', a2, '2=', i8, 3x, '  INT =', i3, 3x,
     .          'NCENTR=', i3, /(10('  L',a2,'(',i2,')=',i3)))
         intp22 = Intp2(k)
         frp    = Intp1(k)
         frp    = frp*2._10**intp22
         frq    = Jdpl0(j) + (frp - 0.5_10)
         write(Iout,1470) frq,Ihrp(k),Iminp(k),Secp(k),
     .    Intp1(k),Intp2(k),frp,Icnd(j)
 1470    format(' INITIAL EPOCH (COORD.TIME) ', f16.8, ' IHR=',
     .    i2, ' IMIN=', i2, ' SEC=', f7.4, ' INT1=', i10,
     .    ' INT2=', i3, ' FRACT=', f16.13, '  ICND=', i3)
         write(Iout,1480) (Pcond(i,j),i = 1,8),
     .                     (i,Pcond(i+6,j),i = 3,22)
 1480    format(9x,'A=', 1pd22.15, 9x, 'E=', d22.15, 7x, 'INC=',
     .          d22.15, 7x, 'ASC=', d22.15/7x, 'PER=', d22.15, 6x,
     .          'ANOM=', d22.15, 4x, 'RADIUS=', d22.15, 6x, 'FLAT=',
     .          d22.15/4('   CON(',i2,')=',d22.15))
         write(Iout,1500) (bl,i,Pcond(i+6,j),i = 23,24),
     .                     (bl,i,Dumcon(i,k),i = 1,12)
 1500    format(2(2x,a1,'CON(',i2,')=',1pd22.15),
     .          2(1x,a1,'CON1(',i2,')=',1pd22.15)/
     .          (4(1x,a1,'CON1(',i2,')=',1pd22.15)))
         write(Iout,1520) (i,Kkk(i,k),i = 1,u_nmprm)
 1520    format(10(3x,'K(',i2,')=',i4))
         nki = Ndumki(k)
         if(nki.gt.33) nki = 33
         write(lineout,1540) Ndumki(k),(Kdumi(i,k),i = 1,nki)
 1540    format(' NUMKI=', i3, '  KI=', 7I2, 26I4)
         write(linealt,1545) (Kdumi(i,k),i=8,nki)
 1545    format(26i5)
         do i=8,nki
            if(linealt((i-8)*5+1:(i-8)*5+3).eq.'-10')
     .       lineout(i*4-2:i*4+1)='-M'//linealt((i-8)*5+4:(i-8)*5+5)
         end do
         write(Iout,1547) lineout
 1547    format(a)
         write(Iout,600) (i,Dumeps(i,k),i = 1,6)
c test for gravitational potential or shape coefs.
         if(klam.gt.Nmphar) goto 1800
         if(k.gt.4 .and. Nplnt(j).ne.Nplhar(klam)) goto 1800
 
c check for non-spherical harmonic shape model (nshp.gt.0)
         if(k.gt.4 .and. Nplnt(j).lt.0 .and. nshp(klam).gt.0) then
 
c planet shape
            if(nshp(klam).eq.1) then
 
c two dimensional fourier series shape model
               call PAGCHK(60,33,0)
 
               write(Iout,1550) Aplnt(j),Nplhar(klam),nshp(klam)
 1550 format('0TWO DIMENSIONAL FOURIER SERIES SHAPE MODEL FOR ', a8,
     .       ' (NPLNT=', i3, ', NSHP=', i3, ')')
               write(Iout,1560) Scontl(klam,1),Scontl(klam,2)
 1560 format('    LATITUDE BOUNDARIES ON PLANET FOR SHAPE MODEL ARE -',
     .    4x, 'LOWER LATITUDE=', f7.2, ', UPPER LATITUDE=', f7.2)
 
c coefficients dependent upon both latitude and longitude
               write(Iout,1570)
               write(Iout,1580) (xshp(1),i,pa(klam,i),i = 1,20)
               write(Iout,1590) (xshp(1),i,lpa(klam,i),i = 1,20)
               write(Iout,1580) (xshp(2),i,pb(klam,i),i = 1,20)
               write(Iout,1590) (xshp(2),i,lpb(klam,i),i = 1,20)
               write(Iout,1580) (xshp(3),i,pc(klam,i),i = 1,20)
               write(Iout,1590) (xshp(3),i,lpc(klam,i),i = 1,20)
               write(Iout,1580) (xshp(4),i,pd(klam,i),i = 1,20)
               write(Iout,1590) (xshp(4),i,lpd(klam,i),i = 1,20)
 1570          format(
     .    '0COEFFICIENTS DEPENDENT UPON BOTH LONGITUDE AND LATITUDE')
 1580          format(4(3x,a4,i2,')=',1pd22.15))
 1590          format(10('  L',a4,i2,')=',i2))
               call PAGCHK(60,19,0)
 
c coefficients dependent upon longitude only
               write(Iout,1600)
 1600          format(' COEFFICIENTS DEPENDENT UPON LONGITUDE ONLY')
               write(Iout,1580) (xshp(5),i,pap(klam,i),i = 1,20)
               write(Iout,1590) (xshp(5),i,lpap(klam,i),i = 1,20)
               write(Iout,1580) (xshp(6),i,pcp(klam,i),i = 1,20)
               write(Iout,1590) (xshp(6),i,lpcp(klam,i),i = 1,20)
 
c coefficients dependent upon latitude only
               write(Iout,1610)
 1610          format(' COEFFICIENTS DEPENDENT UPON LATITUDE ONLY')
               write(Iout,1620) xshp(7),papp(klam),xshp(8),pbpp(klam)
               write(Iout,1630) xshp(7),lpapp(klam),xshp(8),
     .                           lpbpp(klam)
 1620          format(3x,2(a4,'=',1pd20.13,5x))
 1630          format(3x,2('L',a4,'=',i3,5x), /)
c
c*  start=2500
c
            else if(nshp(klam).eq.2) then
c        altitude grid-local shape model-nshp=2
c
c        array transf(2,4,500) contains the 1000 real*4
c        elements of the altitude grid for each of
c        4 possible planets which are equivilenced to
c        the real*8 harmonic coefficients in /scoef4/
c
c        check for non-zero elements in grid and transfer
c        to array grid for easier handling
               nelem(1) = 0
               nelem(2) = 0
               igrid    = 0
               do i = 1,1000/u_stdsz
                  do l = 1,u_stdsz
                     igrid = igrid + 1
                     if(igrid.gt.ngdpts(klam)) goto 1640
                     Grid(igrid) = transf(l,klam,i)
                     if(Grid(igrid).ne.0.0) nelem(1)= nelem(1) + 1
                     if(Lpzhar(klam,igrid).gt.0) nelem(2)= nelem(2) + 1
                  end do
               end do
 
c check for new page
 1640          ltest = 59 - (ngdpts(klam)/6 + 8)
               if(nelem(1).eq.0) ltest = 50
               if(Line.ge.ltest) call NEWPG
 
               write(Iout,1650) Aplnt(j),Nplhar(klam),nshp(klam)
 1650          format('0ALTITUDE GRID - LOCAL SHAPE MODEL FOR ',
     .                a8, ' (NPLNT=', i3, ', NSHP=', i3, ')')
               write(Iout,1660) (Scontl(klam,i),i = 1,6)
 1660          format(3x,'RANGE ON PLANET(IN DEGREES):  TLAT= ',
     .                f7.2, ' TO ', f7.2, 4x, 'TLON= ', f7.2,
     .                ' TO ', f7.2, /, 3x,
     .                'GRID SPACING(IN DEGREES):  TLATIN= ',
     .                f7.2, 4x, 'TLONIN= ', f7.2)
               write(Iout,1670) latdim(klam),londim(klam),ngdpts(klam)
 1670          format(3x,'GRID DIMENSIONS ARE:   (', i2, ',', i3, ')'//
     .           3x, 'ALTITUDE IN KILOMETERS AT GRID POINTS (NGDPTS=',
     .           i4, ')')
 1700          format('0  L-VECTOR FOR GRID')
 
               Line   = Line + 7
               nlat   = latdim(klam)
               nlon   = londim(klam)
               ncount = 6
 
c print out grid and l-vector (unless all zero)
               do lgrid = 1, 2
                  if(nelem(lgrid).eq.0) then
 
c all grid points zero
                     write(Iout,1690)
 1690                format(4x,'*** ALL VALUES = 0 ***')
                     Line = Line + 1
                  else
                     igrid  = 0
                     icount = 0
 
c lgrid=1 for grid print;lgrid=2 for l-vector print
                     do ilat = 1, nlat
                        do ilon = 1, nlon
                           icount = icount + 1
                           if(ilat.eq.nlat .and.
     .                        ilon.eq.nlon) ncount = icount
                           idex(icount) = ilat
                           jdex(icount) = ilon
                           if(icount.ge.ncount) then
 
                              call PAGCHK(58,1,0)
                              if(lgrid.eq.1)
     .                          write (Iout,1680) (idex(i),jdex(i),
     .                          Grid(igrid+i),i = 1,ncount)
 1680 format(4x,6('(',i2,',',i3,')=',f9.4,', '))
                              if(lgrid.eq.2)
     .                          write (Iout,1710) (idex(i),jdex(i),
     .                          Lpzhar(klam,igrid+i),i = 1,ncount)
 1710 format(4x,9('L(',i2,',',i3,')=',i2,', '))
 
                              igrid  = igrid + ncount
                              icount = 0
                           endif
                        end do
                     end do
                  endif
 
                  if(lgrid.ne.2) then
 
c print out l-vector for grid
                     ltest = 59 - (ngdpts(klam)/9 + 3)
                     if(Line.ge.ltest .and. Line.ge.30) call NEWPG
                     write(Iout,1700)
                     Line   = Line + 2
                     ncount = 9
                  endif
 
               end do
            else if(nshp(klam).eq.3) then
               if(line.gt.58) call NEWPG
               write(Iout,1715) Aplnt(j),Nplhar(klam),nshp(klam)
 1715          format('0EXTERNAL SHAPE MODEL FOR ',
     .                a8, ' (NPLNT=', i3, ', NSHP=', i3, ')')
               Line=Line+2
            endif
            klam = klam + 1
            goto 1800
         endif
         if(k.gt.4) then
            npz = Npzone(klam)
            npt = Nptess(klam)
         else if(npln.eq.3) then
            pname = earth
            npz   = Nezone
            npt   = Netess
         else if(npln.eq.10) then
            npz = Nmzone
            npt = Nmtess
         else
            goto 1800
         endif
 
         call MVC(pname,1,8,harmsg,65)
         call MVC(bl8,1,8,harmsg,1)
         call MVC('  ZONAL ', 1, 8, harmsg, 9)
         call MOVEBL('ZONE', 4, harmsg(87:90), 4)
         call EBCDIX(npln,harmsg,81,3)
         call EBCDIX(npz,harmsg,92,3)
         call PAGSET(harmsg,-24)
         n1 = npz - 1
         ll = 2
         if(n1.gt.0) ll = 4 + n1/4 + n1/10
c
c planets gravitational zonal harmonic coefficients
         call PAGHED(ll)
         if(n1.gt.0) then
            n2 = n1 + 1
            if(npln.eq.3) then
 
c special case: earth zonal harmonics
               write(Iout,1720) (xx1,i,Ezhar(i-1),i = 2,n2)
               write(Iout,1730) (blank,xx1,i,Lezhar(i-1),i = 2,n2)
            else if(npln.eq.10) then
 
c special case: moon zonal harmonics
               write(Iout,1720) (xx1,i,Mzhar(i-1),i = 2,n2)
               write(Iout,1730) (blank,xx1,i,Lmzhar(i-1),i = 2,n2)
            else
               write(Iout,1720) (xx1,i,Pzhar(klam,i-1),i = 2,n2)
 1720          format(2(4(7x,a1,i1,'=',1pd22.15)/),
     .                (4(6x,a1,i2,'=',1pd22.15)))
               write(Iout,1730) (blank,xx1,i,Lpzhar(klam,i-1),i =
     .                           2, n2)
 1730          format(8(a2,4x,'L',a1,i1,'=',i3),
     .                2(a2,'   L',a1,i2,'=',i3)/
     .                (10(a2,'   L',a1,i2,'=',i3)))
            endif
         endif
         n1 = (npt*(npt+1))/2 - 1
         n1a = 1
         n1b = min0(n1,220)
         n1bs = n1b
         ll = 0
         if(n1.gt.0) ll = 1 + n1b/4
         lll = 1 + n1/10
         call MVC(teslit,1,8,harmsg,1)
         call MVC(coslit,1,8,harmsg,9)
         call MOVEBL(teslit,4,harmsg(87:90), 4)
         call EBCDIX(npt,harmsg,92,3)
         call PAGSET(harmsg,-24)
c
c planets gravitational potential tesseral harmonic
c cosine coefficients
         call PAGHED(ll)
         if(n1.gt.0) then
            if(npln.eq.3) then
 
c special case: earth tesseral harmonics
               write(Iout,1740) (xx2,Idt(i),Idi(i),Echar(i),i=1,n1b)
               do while (n1b.lt.n1)
                  n1a = n1b + 1
                  n1b = min0(n1,n1b + 220)
                  call PAGE(1+(n1b-n1a+1)/4,1)
                  write(Iout,1745) (xx2,Idt(i),Idi(i),Echar(i),
     .                              i=n1a,n1b)
               end do
               call PAGCHK(60,lll,1)
               write(Iout,1750) (bl,xx2,Idt(i),Idi(i),Lechar(i),
     .                           i = 1, n1)
            else if(npln.eq.10) then
 
c special case: moon tesseral harmonics
               write(Iout,1740) (xx2,Idt(i),Idi(i),Mchar(i),i=1,n1b)
               do while (n1b.lt.n1)
                  n1a = n1b + 1
                  n1b = min0(n1,n1b + 220)
                  call PAGE(1+(n1b-n1a+1)/4,1)
                  write(Iout,1745) (xx2,Idt(i),Idi(i),Mchar(i),
     .                              i=n1a,n1b)
               end do
               call PAGCHK(60,lll,1)
               write(Iout,1750) (bl,xx2,Idt(i),Idi(i),Lmchar(i),
     .                           i = 1, n1)
            else
               write(Iout,1740) (xx2,Idt(i),Idi(i),Pchar(klam,i),
     .                           i = 1, n1b)
 1740          format(11(4(4x,a1,i1,'(',i1,')=',1pd22.15)/),
     .                2(4(3x,a1,i2,'(',i1,')=',1pd22.15)/), 3x, a1,
     .                i2, '(', i1, ')=', 1pd22.15, 2x, a1, i2, '(',
     .                i2, ')=', 1pd22.15,
     .                2(3x,a1,i2,'(',i1,')=',1pd22.15)/
     .                (4(2x,a1,i2,'(',i2,')=',1pd22.15)))
               do while (n1b.lt.n1)
                  n1a = n1b + 1
                  n1b = min0(n1,n1b + 220)
                  call PAGE(1+(n1b-n1a+1)/4,1)
                  write(Iout,1745) (xx2,Idt(i),Idi(i),Pchar(klam,i),
     .                              i=n1a,n1b)
 1745             format((4(2x,a1,i2,'(',i2,')=',1pd22.15)))
               end do
               call PAGCHK(60,lll,1)
               write(Iout,1750) (blank,xx2,Idt(i),Idi(i),Lpchar
     .                           (klam,i),i = 1,n1)
 1750          format(4(10(a2,'  L',a1,i1,'(',i1,')=',i2)/),
     .                4(a2,'  L',a1,i1,'(',i1,')=',i2),
     .                6(a2,' L',a1,i2,'(',i1,')=',i2)/
     .                3(a2,' L',a1,i2,'(',i1,')=',i2),
     .                7(a2,'L',a1,i2,'(',i2,')=',i2)/
     .                (10(a2,'L',a1,i2,'(',i2,')=',i2)))
            endif
         endif
c
c planets gravitational potential tesseral harmonic
c sine coefficients
         call MVC(sinlit,1,8,harmsg,9)
         call PAGSET(harmsg,-24)
         call PAGHED(ll)
         if(n1.gt.0) then
            n1b = n1bs
            if(npln.eq.3) then
 
c special case: earth tesseral harmonics
               write(Iout,1740) (xx3,Idt(i),Idi(i),Eshar(i),i=1,n1b)
               do while (n1b.lt.n1)
                  n1a = n1b + 1
                  n1b = min0(n1,n1b + 220)
                  call PAGE(1+(n1b-n1a+1)/4,1)
                  write(Iout,1745) (xx3,Idt(i),Idi(i),Eshar(i),
     .                              i=n1a,n1b)
               end do
               call PAGCHK(60,lll,1)
               write(Iout,1750) (bl,xx3,Idt(i),Idi(i),Leshar(i),
     .                           i = 1, n1)
            else if(npln.eq.10) then
 
c special case: moon tesseral harmonics
               write(Iout,1740) (xx3,Idt(i),Idi(i),Mshar(i),i=1,n1b)
               do while (n1b.lt.n1)
                  n1a = n1b + 1
                  n1b = min0(n1,n1b + 220)
                  call PAGE(1+(n1b-n1a+1)/4,1)
                  write(Iout,1745) (xx3,Idt(i),Idi(i),Mshar(i),
     .                              i=n1a,n1b)
               end do
               call PAGCHK(60,lll,1)
               write(Iout,1750) (bl,xx3,Idt(i),Idi(i),Lmshar(i),
     .                           i = 1, n1)
            else
               write(Iout,1740) (xx3,Idt(i),Idi(i),Pshar(klam,i),
     .                           i = 1, n1b)
               do while (n1b.lt.n1)
                  n1a = n1b + 1
                  n1b = min0(n1,n1b + 220)
                  call PAGE(1+(n1b-n1a+1)/4,1)
                  write(Iout,1745) (xx3,Idt(i),Idi(i),Pshar(klam,i),
     .                              i=n1a,n1b)
               end do
               call PAGCHK(60,lll,1)
               write(Iout,1750) (blank,xx3,Idt(i),Idi(i),Lpshar
     .                           (klam,i),i = 1,n1)
            endif
         endif
         klam = klam + 1
 
c no printout coded for nshp gt 3
 1800    continue
c
c print limited asteroids
c*  start=2800
      do i = 1, Numsml
         call PAGCHK(60,4,0)
         write(Iout,1810) Nbsml(i),Scond(7,i),Denpts(i),
     .                     Jdsml0(i),(Scond(j,i),j = 1,6)
 1810    format('-SMALL BODY', i5, '  RADIUS=', 1pe14.7,
     .          ' KM, DEN=PRMTER(', i2, ')  JD0=', i7/' A=',
     .          e14.7, '  E=', e14.7, '  I=', e14.7, '  ASC=',
     .          e14.7, '  PER=', e14.7, '  ANOM=', e14.7)
      end do
c
c printout et-ut2, au-ut1 or wobble table
c*  start=3000
      if(Numdt.gt.0) then
         intdy = 0
         if(Line.gt.50) call NEWPG
         do i = 1, Numdt
            if(i.ne.1) then
               if(Line.le.57) goto 1860
               call NEWPG
            endif
            if(Jddt0.gt.1) then
               write(Iout,1820)
 1820          format(
     .  '0EPHEMERIS TIME MINUS UNIVERSAL TIME TABLE (ET-UT2 OR CT-UT2)'/
     .  '  N     DT', 7x, 'JDDT  INTDY INTYR  LDT')
               Line = Line + 3
            else if(Jddt0.ge.1) then
               write(Iout,1830)
 1830          format('0A1-UT1 TABLE'/'  N     DT', 7x,
     .                'JDDT  INTDY INTYR  LDT')
               Line = Line + 3
            else if(Jddt0.lt.0) then
               write(Iout,1840)
 1840          format('0WOBBLE TABLE'/'  N   XWOB', 12x, 'YWOB', 12x,
     .                ' JDDT   INTDY INTYR LDT(N+200) LDT(N+400)')
               Line = Line + 3
            else
               write(Iout,1850)
 1850          format('0A1-UT1 AND WOBBLE TABLE'/'  N   A1-UT1',
     .                11x, 'XWOB', 12x, 'YWOB', 12x,
     .            ' JDDT   INTDY INTYR LDT(N) LDT(N+200) LDT(N+400)')
               Line = Line + 3
            endif
 1860       if(i.ne.1) intdy = Jddt(i) - Jddt(i - 1)
            diff = intdy
            diff = diff/365.25_10
            if(Jddt0.lt.0) then
               write(Iout,1870) i,Dt(i + 200),Dt(i + 400),
     .                           Jddt(i),intdy,diff,Ldt(i + 100),
     .                           Ldt(i + 200)
 1870          format(i4,'.', 2x, 1pe14.7, 2x, 1pe14.7, i8, i6,
     .                f7.3, 2I5)
            else if(Jddt0.eq.0) then
               write(Iout,1880) i,Dt(i),Dt(i + 200),Dt(i + 400),
     .                           Jddt(i),intdy,diff,Ldt(i),
     .                           Ldt(i + 200),Ldt(i + 400)
 1880          format(i4,'.', 1p, 3E16.7, i8, i6, f7.3, 3I7)
            else
               write(Iout,1890) i,Dt(i),Jddt(i),intdy,diff,
     .                           Ldt(i)
 1890          format(i4,'.', f9.4, i8, i6, f7.3, i3)
            endif
            Line = Line + 1
         end do
      else
         call PAGCHK(60,2,0)
         write(Iout,1900)
 1900    format('0 THERE IS NO INPUT ET-UT2, A1-UT1, OR WOBBLE TABLE')
      endif
c
c printout for observing sites   and
c printout for spots on observed bodies
      call PRNCRD(nstop)
c
c printout for pulsars
      call PRNPSR(nstop)
c
c printout for sky corrections
      call PRNSTR
c
c printout for observation constants
      call PRNOBS(nstop,iseq)
 1990 continue
 
      call PRNACM(nstop)
      call PRNFIL(noprnt,nstop)
      if(nstop.gt.0) then
c
c input data contains errors, program terminated
         call EBCDI(nstop,erms,6)
         call SUICID(erms,10)
      endif
 
      if(Ict(1).eq.-99) call SUICID(' NORMAL STOP - JUST PRINTOUT', 7)
c
c input data inserted into storage without errors
c write quantities onto disk so they do not clutter up
c core storage
      call INWRIT
c
c printout time required in input link
      call TIMRIT('  SETTING UP, READING IN AND WRITING OUT INPUT DATA '
     .            , 13)
      return
c
      end
