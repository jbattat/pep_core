      subroutine ADJUST(keep)
 
      implicit none

c     ash/forni    june 1969     subroutine adjust
c     rdr/kcl jan. 1980 add effect of the a priori to the norm
c     printout error analysis for all observations used in forming
c     normal equations.  make adjustments to site coordinates, equator-
c     equinox corrections, solar system parameters in this routine
c     and then to other parameters by calling appropriate routines in
c     the appropriate order. printout adjustments to inverse masses.
c     decide if convergence has been obtained in least squares iteration
c
c arguments
c if keep true then keep this solution
      logical*4 keep
 
c array dimensions
      include 'globdefs.inc'

c commons
      include 'adjstf.inc'
      include 'aprtbf.inc'
      include 'cureph.inc'
      include 'dtparm.inc'
      include 'empcnd.inc'
      include 'eqenox.inc'
      include 'ethhar.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'monhar.inc'
      include 'namtim.inc'
      include 'param.inc'
      include 'plnhar.inc'
      integer*4 ngdpts(4)
      equivalence (Scontl(1,9),ngdpts)
      include 'psrstf.inc'
      include 'rtsidesl.inc'
      include 'scail.inc'
      include 'scoef4.inc'
      include 'stcord.inc'
c
c external work area
      common/WRKCOM/ Apsol(1000),Vect(1000),Pointr(1000),
     . Fcthr1((u_mxtes*(u_mxtes+1))/2-1),
     . Fcthrt((u_mxtes*(u_mxtes+1))/2-1),
     . Fctbod(u_nmbod),Wdrs(u_nmbod),Fctprm(u_nmprm-30),
     . Fctpln(u_nmbod),Fctprb(u_nmbod),Fctsit(6),Fcteqn(3),
     . Ppr16(2,10),Tfr8(5,30),Rjr8(5,2),Tfc8(30),Rjc8(2,2),Ppr4(10),
     . Tfc2(30),Tfi2(4,30),Rjc2(2),Rji2(3,2)
      character*8 Wdrs,Rjc8,Tfc8
      character*4 Ppr4
      character*2 Rjc2,Tfc2
      real*10   Ppr16
      real*10 Apsol,Vect,Fcthr1,Fcthrt,Fctbod,Fctprm,Fctpln,
     .       Fctprb,Fctsit,Fcteqn,Tfr8,Rjr8
      integer*2 Tfi2,Rji2,Pointr
      real*10 fctemb(u_nmbod),fcter(u_nmbod),fctmon(u_nmbod),
     . fctmr(u_nmbod),fctplb(u_nmbod),fctpsr(16)
      equivalence (Fctbod,fctemb,fcter,fctmon,fctmr,fctplb,fctpsr)

c external functions
      real*10 DOTN,LEGSCL
c
c internal to subroutine adjust
      real*10 ams,convd1,ddt,double,dum,fnobs,pstmen,pstsum,t1,t2
      integer   i,i4,iaprio,ibfap,ibr,igmvry,imtf,ip,ippr,
     .          ipsrpr,j,jj,k,kk,l,ll,ll1,ll1p,ll2,m
      integer   msave,ne,nphase,nrbias,nsave,nspot
      character*8 blank/'        '/
      character*8 peryr/'(PER YR)'/
      real*10 a(3),mass2
      real*4    snwv
      character*32 tmesg/'   ADJUSTING*******  PARAMETERS '/
      character*2 astrik(3)/'* ','& ','  '/
      character*8 dtwrd(5)/' ET-UT2 ',' A1-UT1 ','  XWOB  ',
     .    '  YWOB  ',' WOBBLE '/
 
      character*12 wrdc(6)/' RADIUS   KM',' LONGITUDEDG',' LATITUDE DG',
     .                     ' VERT.V MM/Y',' WEST.V MM/Y',' NRTH.V MM/Y'/
 
      character*12 nameqn(3)/'DEQUINOX  S ','DEQUATOR  " ',
     .                       'D-DECL    " '/
 
      character*8 wrdp(70)/
     .     'RELFCT  ', 'GMVARY  ', 'SUNHAR  ', 'PRMTER34', 'PRMTER35',
     .     'PRMTER36', 'PRMTER37', 'PRMTER38', 'SUNPOE  ', 'PRMTER40',
     .     'BETA    ', 'GAMMA   ', 'PRMTER43', 'PRMTER44', 'PRMTER45',
     .     'PRMTER46', 'ASCABELT', 'INCABELT', 'DSTABELT', 'MASABELT',
     .     'AULTSC  ', 'LTVARY  ', 'RELDEL  ', 'RELDOP  ', 'PRMTER55',
     .     'PRMTER56', 'PRMTER57', 'PRMTER58', 'PRMTER59', 'PLASMA C',
     .     'PLASMA V', 'ATM FACT', 'ION FACT', 'PRMTER64', 'PRMTER65',
     .     'PRMTER66', 'PRMTER67', 'PRMTER68', 'PRMTER69', 'PRMTER70',
     .     'PRMTER71', 'CTVARY  ', 'PRMTER73', 'PRMTER74', 'PRMTER75',
     .     'PRMTER76', 'PRMTER77', 'PRMTER78', 'PRMTER79', 'PRMTER80',
     .     'SFATTOCT', 'PHATTOCT', 'PRMTER83', 'PRMTER84', 'PRMTER85',
     .     'PRMTER86', 'PRMTER87', 'PRMTER88', 'PRMTER89', 'PRMTER90',
     .     'ECINC   ', 'SEQINC  ', 'SEQASC  ', 'SUNRAD  ', 'SUNRTR  ',
     .     'PRMTER96', 'PRMTER97', 'MDSTAU  ', 'MDSTSC  ', 'LTVEL   '/
 
      character*8 wrdem(u_nmbod)/'A   (AU)', 'E       ', 'INC (DG)',
     .          'ASC (DG)', 'PER (DG)', 'ANOM(DG)', 'ERAD(KM)',
     .          'EFLAT   ', 'CON( 3) ', 'CON( 4) ', 'CON( 5) ',
     .          'CON( 6) ', 'CON( 7) ', 'CON( 8) ', 'CON( 9) ',
     .          'CON(10) ', 'CON(11) ', 'CON(12) ', 'CON(13) ',
     .          'CON(14) ', 'CON(15) ', 'CON(16) ', 'CON(17) ',
     .          'CON(18) ', 'CON(19) ', 'CON(20) ', 'CON(21) ',
     .          'CON(22) ', 'CON(23) ', 'RELFCT  '/
      integer*2 npln3/3/, npln0/0/
 
      character*8 erotat/' EROTAT '/
      character*8 wrder(u_nmbod)/'COND(1) ', 'COND(2) ', 'COND(3) ',
     .          'COND(4) ', 'COND(5) ', 'COND(6) ', 'PSID(1) ',
     .          'PSID(2) ', 'PSID(3) ', 'PSI(1)  ', 'PSI (2) ',
     .          'PSI (3) ', 'CON( 7) ', 'CON( 8) ', 'CON( 9) ',
     .          'CON(10) ', 'R0 A1UT1', 'R1 A1UT1', 'R2 A1UT1',
     .          'B0 A1UT1', 'B1 A1UT1', 'D0 A1UT1', 'D1 A1UT1',
     .          'A0 A1UT1', 'A1 A1UT1', 'C0 A1UT1', 'C1 A1UT1',
     .          'PRCN    ', 'OBLQ    ', 'NUTN    '/
      integer*2 mpln3/-3/
      character*8 wrder1(20)/'FNUT COS', 'FNUT SIN', 'DEPS  14',
     .          'DPSI  14', 'DEPS 183', 'DPSI 183', 'DEPS 365',
     .          'DPSI 365', 'DEPS 19Y', 'DPSI 19Y', 'DEPSX 14',
     .          'DPSIX 14', 'DEPSX183', 'DPSIX183', 'DEPSX365',
     .          'DPSIX365', 'DEPSX19Y', 'DPSIX19Y', 'K/C   14',
     .          'K/C   30'/
      character*5 ec/'EC00('/, es/'ES00('/
 
      character*8 wrdmn(u_nmbod)/'A   (AU)', 'E       ', 'INC (DG)',
     .          'ASC (DG)', 'PER (DG)', 'ANOM(DG)', 'MRAD(KM)',
     .          'CON( 2) ', 'CON( 3) ', 'CON( 4) ', 'CON( 5) ',
     .          'CON( 6) ', 'CON( 7) ', 'CON( 8) ', 'CON( 9) ',
     .          'CON(10) ', 'CON(11) ', 'CON(12) ', 'CON(13) ',
     .          'CON(14) ', 'CON(15) ', 'CON(16) ', 'CON(17) ',
     .          'CON(18) ', 'CON(19) ', 'ETADELTA', 'CON(21) ',
     .          'CON(22) ', 'CON(23) ', 'RELFCT  '/
      integer*2 npln10/10/
 
      character*8 mrotat/' MROTAT '/
      character*8 wrdmr(u_nmbod)/'COND(1) ', 'COND(2) ', 'COND(3) ',
     .          'COND(4) ', 'COND(5) ', 'COND(6) ', 'CON( 1) ',
     .          'CON( 2) ', 'BETA    ', 'GAMMA   ', 'CON( 5) ',
     .          'CON( 6) ', 'CON( 7) ', 'CON( 8) ', 'CON( 9) ',
     .          'CON(10) ', 'CON(11) ', 'CON(12) ', 'CON(13) ',
     .          'CON(14) ', 'CON(15) ', 'CON(16) ', 'CON(17) ',
     .          'CON(18) ', 'CON(19) ', 'CON(20) ', 'CON(21) ',
     .          'CON(22) ', 'CON(23) ', 'CON(24) '/
      integer*2 mpln10/-10/
 
      character*5 mc/'MC00('/, ms/'MS00('/
 
      character*8 wrdr(6)/'R   (AU)', 'RA  (DG)', 'DECL(DG)',
     .    'V* AU/DY', 'FLAZ(DG)', 'FLANG(D)'/
      character*8 wrds(u_nmbod)/'A   (AU)', 'E       ', 'INC (DG)',
     .    'ASC (DG)',
     .          'PER (DG)', 'ANOM(DG)', 'PRAD(KM)', 'E1      ',
     .          'E2      ', 'E3      ', 'NU  (DG)', 'PHI0(DG)',
     .          'PERIOD  ', 'DEC (DG)', 'R.A.(DG)', 'I0  (DG)',
     .          'PSI0(DG)', 'MU      ', 'CON(13) ', 'CON(14) ',
     .          'CON(15) ', 'CON(16) ', 'CON(17) ', 'CON(18) ',
     .          'CON(19) ', 'CON(20) ', 'CON(21) ', 'CON(22) ',
     .          'CON(23) ', 'RELFCT  '/
      character*4 wrdl2(16)/', 1)', ', 2)', ', 3)', ', 4)', ', 5)',
     .          ', 6)', ', 7)', ', 8)', ', 9)', ',10)', ',11)', ',12)',
     .          ',13)', ',14)', ',15)', ',16)'/
      character*8 wrdb(u_nmbod)/'X   (AU)', 'Y   (AU)', 'Z   (AU)',
     .          'X'' AU/DY', 'Y'' AU/DY', 'Z'' AU/DY', 'CON( 1) ',
     .          'CON( 2) ', 'CON( 3) ', 'CON( 4) ', 'CON( 5) ',
     .          'CON( 6) ', 'CON( 7) ', 'CON( 8) ', 'CON( 9) ',
     .          'CON(10) ', 'CON(11) ', 'CON(12) ', 'CON(13) ',
     .          'CON(14) ', 'CON(15) ', 'CON(16) ', 'CON(17) ',
     .          'CON(18) ', 'CON(19) ', 'CON(20) ', 'CON(21) ',
     .          'CON(22) ', 'CON(23) ', 'RELFCT  '/
 
      character*4 charz(4)/'(1) ', '(2) ', '(3) ', '(4) '/
      character*4 chart(4)/',1) ', ',2) ', ',3) ', ',4) '/
      character*5 pc/'PC00('/, ps/'PS00('/
 
      character*8 strnam/'  STAR  '/, bodynm
      integer*2 onemns/-1/
      real*10 rhs(1)
      equivalence (rhs, Sigma)
      real*10 norm,normsq,pstrms,pstssq,prerms,normsv,totssq
      character*8 embary/' EMBARY '/,earth/' EARTH  '/
      character*8 plsr/'PLSR....'/
c
c write out error analysis
      call PAGSET('ADJUSTMENT TO PARAMETERS',6)
      call PAGSET('        L VECTOR       PARAMETER   BODY  NPLNT'//
     .            '  JD0          OLD VALUE        ADJUSTMENT    '//
     .            '     NEW VALUE           SIGMA    FRACT ',-33)
      call PAGSET(-1,Jout)
      call NEWPG
      ams = Measmt
      if(Measmt.le.0) ams = 1._10
      a(1)  = Ermeas(1)/ams
      a(2)  = Ermeas(2)/ams
      fnobs = Ermeas(3)/ams
      a(3)  = SQRT(fnobs)
      if(Ict(42).gt.0) then
         write(Iout,50) Ephnum,Datea,Dateb
   50    format('0EPOCH NUMBER',i4,5x,'FROM',f13.4,'  TO',f13.4)
         Line = Line + 2
         if(Jout.gt.0) write(Jout, 50) Ephnum, Datea, Dateb
      endif
      write(Iout, 100) Measmt, a, fnobs, Ermeas(3)
  100 format('0ERROR ANALYSIS FOR THE', i8,
     .       ' MEASUREMENTS USED TO FORM THE NORMAL EQUATIONS'/
     .       '          AVERAGE (OBS-TH)/ERROR', 1pe14.5/
     .       '       AVERAGE ABS(OBS-TH)/ERROR', 1pe14.5/
     .       ' ROOT MEAN SQUARE (OBS-TH)/ERROR', 1pe14.5/
     .       '     AVERAGE ((OBS-TH)/ERROR)**2', 1pe14.5/
     .       '         SUM ((OBS-TH)/ERROR)**2', 1pe14.5)
      if(Jout.gt.0) write(Jout, 100) Measmt, a, fnobs, Ermeas(3)
      write(Iout, 200) Nparam, Measmt
  200 format('0ADJUSTMENT TO', i5, ' PARAMETERS BY SOLUTIONS OF THE ',
     . 'NORMAL EQUATIONS FORMED FROM THE RESULTS OF', i8,
     . ' MEASUREMENTS')
      if(Jout.gt.0) write(Jout, 200) Nparam, Measmt
      Line = Line + 10
      call PAGHED(0)
c
c initialization
      Keepit = keep
      Fradj  = Eps(14)
      if(Fradj.eq.0._10) Fradj = 1._10
      igmvry = 0
      Iboda  = 0
      ipsrpr = 0
      Aukm   = Aultsc*Ltvel
      convd1 = 1._10/Convd
      do i = 1, 3
         Statf(i) = 0._10
         Nkind(i) = 0
      end do
      Il = 1
      Nf = 0
      N  = 0
      call FRSTAT(0, 3, ' INITIAL')
c
c set scale factors for converting from unnormalized
c gravitational potential harmonic adjustments to normalized
c grav.pot.harmonic output (tesserals)
c note: Fcthr1 must precede Fcthrt and be the same size, so that shape
c models can use it for unit scaling.
      k = 0
      do j = 2, u_mxtes
         do i = 1, j
            k = k + 1
            Fcthr1(k) = 1._10
            Fcthrt(k) = 1._10/LEGSCL(j, i)
         end do
      end do
c
c adjustment to observing site coordinates
c*  start=100
      Fctsit(1) = Aukm
      Fctsit(2) = Aultsc/Convd
      do i=4,6
         Fctsit(i)=1._10
      end do
      do k = 1, Numsit
         Fctsit(3) = Fctsit(2)
         if(Kscrd(k).lt.0) Fctsit(3) = Fctsit(1)
         call LINCHK
         do i = 1,6
            if(Lscrd(i,k).gt.0) then
               call ADJAST(Scord(i,k), Fctsit(i))
               write(Iout, 210) astrik(Ntype), N, i, k, Lscrd(i, k),
     .                          wrdc(i), sitd(k),
     .                          Scord(i, k), Adj, Nwv, Sig, Fract
  210          format(1x, a1, i4, '. LSCRD(', i1, ',', i2, ')=', i2,
     .                a12, ' OF ', a8, 8x, 1pd22.15, 1pd16.8, 1pd22.15,
     .                1pd10.3, 0pf8.3)
               if(Jout.gt.0) write(Jout, 210) astrik(Ntype), N, i,
     .            k, Lscrd(i, k), wrdc(i), sitd(k),
     .            Scord(i, k), Adj, Nwv, Sig, Fract
               if(Keepit) Scord(i, k) = Nwv
            endif
         end do
         call FRSTAT(1, 2, sitd(k))
      end do
      call FRSTAT(2, 2, 'OBS SITE')
c
c adjustment to equator-equinox corrections for optical
c observation series
c*  start=200
      Fcteqn(1) = 1._10
      Fcteqn(2) = 1._10
      Fcteqn(3) = 1._10
      do k = 1, Numeqn
         call LINCHK
         do i = 1, 3
            if(Leqn(i,k).gt.0) then
               double = deqnx(k, i)
               call ADJAST(double, Fcteqn(i))
               snwv = Nwv
               write(Iout, 220) astrik(Ntype), N, i, k, Leqn(i, k),
     .                          nameqn(i), Eqnsit(k),
     .                          Eqnser(k), deqnx(k, i), Adj, snwv, Sig,
     .                          Fract
  220          format(1x, a1, i4, '. LEQN(', i1, ',', i3, ')= ', i1,
     .                1x, a12, ' FOR ', a4, 1x, a4, 4x, 6x, 1pe12.5,
     .                5x, 1pd16.8, 5x, 1pe12.5, 5x, 1pd10.3, 0pf8.3)
               if(Jout.gt.0) write(Jout, 220) astrik(Ntype), N, i,
     .            k, Leqn(i, k), nameqn(i), Eqnsit(k),
     .            Eqnser(k), deqnx(k, i), Adj, snwv, Sig, Fract
               if(Keepit) deqnx(k, i) = snwv
            endif
         end do
         call FRSTAT(1, 2, 'E-E SITE')
      end do
      call FRSTAT(2, 2, 'EQ - EQ ')
c
c adjustments to sky corrections (star catalog errors)
      call ADJSKY
c
c adjustment to solar system parameters
c*  start=300
      do i = 1, 70
         Fctprm(i) = 1._10
      end do
      Fctprm(1)  = (8.64E4_10/(Gauss*Aultsc))**2
      Fctprm(2)  = 1._10/365.25_10
      Fctprm(3)  = (Aukm/Sunrad)**2
      Fctprm(11) = Fctprm(1)/Aultsc
      Fctprm(12) = Fctprm(1)/Aultsc
      Fctprm(17) = convd1
      Fctprm(18) = convd1
      Fctprm(21) = Aultsc
      Fctprm(23) = Fctprm(1)/Aultsc
 
c imtf=  inverse mass counter
      imtf = 0
      call LINCHK
      do k = 1, u_nmprm
         if(Lprm(k).le.0) goto 300
         l = Lprm(k)
         if(l.gt.30) then
            m = l - 30
            call ADJAST(prmter(l), Fctprm(m))
            write(Iout, 240) astrik(Ntype), N, k, Lprm(k), wrdp(m),
     .                       blank, prmter(l), Adj, Nwv, Sig, Fract
  240       format(1x, a1, i4, '. LPRM(', i3, ')  =', i2, 1x, 2A8, 15x,
     .             1pd22.15, 1pd16.8, 1pd22.15, 1pd10.3, 0pf8.3)
            if(Jout.gt.0) write(Jout, 240) astrik(Ntype), N, k,
     .                              Lprm(k), wrdp(m), blank, prmter(l)
     .                              , Adj, Nwv, Sig, Fract
 
            if(l.eq.32 .or. l.eq.72) then
               igmvry = igmvry + 1
               Rjc2(igmvry) = astrik(ntype)
               Rji2(1,igmvry) = N
               Rji2(2,igmvry) = k
               Rji2(3,igmvry) = Lprm(k)
               Rjc8(1,igmvry) = wrdp(m)
               Rjc8(2,igmvry) = peryr
               Rjr8(1,igmvry) = prmter(l)*365.25_10
               Rjr8(2,igmvry) = Adj*365.25_10
               Rjr8(3,igmvry) = Nwv*365.25_10
               Rjr8(4,igmvry) = Sig*365.25_10
               Rjr8(5,igmvry) = Fract
            endif
         else
            call ADJAST(prmter(l), 1._10)
            bodynm = blank
            do i = 1, Numpln
               if(l.eq.Nplnt(i)) then
                  bodynm = Aplnt(i)
                  goto 260
               endif
            end do
            if(l.eq.10) bodynm = Aplnt(17)
            if(l.eq.3) bodynm  = embary
  260       write(Iout, 280) astrik(Ntype), N, k, Lprm(k), l, bodynm,
     .                       Mass(l), Adj, Nwv, Sig, Fract
  280       format(1x, a1, i4, '. LPRM(', i3, ')  =', i2, ' MASS(', i2,
     .             ') OF ', 1A8, 11x, 1pd22.15, 1pd16.8, 1pd22.15,
     .             1pd10.3, 0pf8.3)
            if(Jout.gt.0) write(Jout, 280) astrik(Ntype), N, k,
     .                              Lprm(k), l, bodynm, Mass(l), Adj,
     .                              Nwv, Sig, Fract
c
c store information for inverse mass print out
            imtf = imtf + 1
            Tfc2(imtf)=  astrik(Ntype)
            Tfi2(1,imtf) = N
            Tfi2(2,imtf) = k
            Tfi2(3,imtf) = Lprm(k)
            Tfi2(4,imtf) = l
            Tfc8(imtf) =  bodynm
            Tfr8(5,imtf) = -Fract
            mass2 = Mass(l)
            if(mass2.ne.0._10) then
               Tfr8(1,imtf) = 1._10/mass2
               Tfr8(4,imtf) = Sig/(mass2*mass2)
            else
               Tfr8(1,imtf) = 0._10
               Tfr8(4,imtf) = 0._10
            endif
            if(Nwv.ne.0._10) then
               Tfr8(3,imtf) = 1._10/Nwv
            else
               Tfr8(3,imtf) = 0._10
            endif
             Tfr8(2,imtf) = Tfr8(3,imtf) - Tfr8(1,imtf)
         endif
         if(Keepit) prmter(l) = Nwv
      end do
  300 call FRSTAT(1, 1, 'SOL SYS ')
c
c adjustment to earth-moon barycenter
c and earth parameters
c*  start=500
      do i = 1, u_nmbod
         fctemb(i) = 1._10
      end do
      do i = 3, 6
         fctemb(i) = convd1
      end do
      fctemb(30) = Fctprm(1)
      call ADJBDY(Lem, Econd, fctemb, wrdem, 'LEM(', ')   ', embary,
     .            npln3, Jdem0, u_nmbod, npln0, npln0)
      call FRSTAT(1, 3, embary)
c
c adjustment to earth rotation initial conditions
c and parameters
      if(Jct(29).ne.0) then
         do i = 1, 20
            wrder(i + 6) = wrder1(i)
         end do
      endif
      do i = 1, u_nmbod
         fcter(i) = 1._10
      end do
      call ADJBDY(Ler, Ercond, fcter, wrder, 'LER(', ')   ', erotat,
     .            mpln3, Jder0, u_nmbod, npln0, npln0)
      call FRSTAT(1, 3, erotat)
c
c adjustment to et-ut2, a1-ut1, or wobble parameters
c*  start=600
      if(Numdt.gt.0) then
         call LINCHK
         jj = 0
         j  = 1
         if(Jddt0.le.1) j = 2
         do while( .true. )
            do kk = 1, Numdt
               k = kk + jj*200
               if(Ldt(k).gt.0) then
                  ddt = Dt(k)
                  call ADJAST(ddt, 1._10)
                  snwv = Nwv
                  write(Iout, 310) astrik(Ntype), N, k, Ldt(k),
     .                             dtwrd(j), Jddt(kk), Dt(k), Adj,
     .                             snwv, Sig, Fract
  310             format(1x, a1, i4, '. LDT(', i3, ')   =', i2, a8,
     .                   'VALUE AT JULDAT ', i8, 5x, 1pe12.5, 5x,
     .                   1pd16.8, 5x, 1pe12.5, 5x, 1pd10.3, 0pf8.3)
                  if(Jout.gt.0) write(Jout, 310) astrik(Ntype), N,
     .               k, Ldt(k), dtwrd(j), Jddt(kk), Dt(k), Adj, snwv,
     .               Sig, Fract
                  if(Keepit) Dt(k) = Nwv
               endif
            end do
            call FRSTAT(1, 1, dtwrd(j))
            if(Jddt0.gt.0) goto 400
            if(jj.ge.2) goto 400
            jj = jj + 1
            j  = j + 1
         end do
      endif
c
c adjustment to earth gravitational potential harmonic coeff.
c*  start=700
  400 if(Nezone.gt.1) then
         call LINCHK
         ll = Nezone - 1
         call ADJHAR(Lezhar, Ezhar, 1, 1, ll, 1, 'EZ', '    ', earth,
     .               npln3, Fcthr1)
         call FRSTAT(1, 2, 'ZON HARM')
      endif
      if(Netess.gt.1) then
         call LINCHK
         ll1 = 1
         do ll = 2, Netess
            call EBCDIX(ll,ec,3,2)
            call ADJHAR(Lechar(ll1), Echar(ll1), 1, 1, ll, 0, ec,
     .                  ')   ', earth, npln3, Fcthrt(ll1))
            ll1 = ll1 + ll
         end do
         call FRSTAT(1, 2, 'COS HARM')
         call LINCHK
         ll1 = 1
         do ll = 2, Netess
            call EBCDIX(ll,es,3,2)
            call ADJHAR(Leshar(ll1), Eshar(ll1), 1, 1, ll, 0, es,
     .                  ')   ', earth, npln3, Fcthrt(ll1))
            ll1 = ll1 + ll
         end do
         call FRSTAT(1, 2, 'SIN HARM')
      endif
c
c adjustment to moon initial conditions and parameters
c*  start=800
      do i = 1, u_nmbod
         fctmon(i) = 1._10
      end do
      do i = 3, 6
         fctmon(i) = convd1
      end do
      fctmon(7)  = Aukm
      fctmon(30) = Fctprm(1)
      call FRSTAT(2, 3, 'ALL HARM')
      call FRSTAT(3, 3, earth)
      call ADJBDY(Lmn, Mcond, fctmon, wrdmn, 'LMN(', ')   ', Aplnt(17),
     .            npln10, Jdmn0, u_nmbod, npln3, npln0)
      call FRSTAT(1, 3, 'MOON IC ')
c
c adjustment to moon rotation initial conditions and
c parameters
      do i = 1, u_nmbod
         fctmr(i) = 1._10
      end do
      call ADJBDY(Lmr, Mrcond, fctmr, wrdmr, 'LMR(', ')   ', mrotat,
     .            mpln10, Jdmr0, u_nmbod, npln0, npln0)
      call FRSTAT(1, 3, mrotat)
c
c adjustment to moon gravitational potential harmonic coeff.
      if(Nmzone.gt.1) then
         call LINCHK
         ll = Nmzone - 1
         call ADJHAR(Lmzhar, Mzhar, 1, 1, ll, 1, 'MZ', '    ', Aplnt(17)
     .               , npln10, Fcthr1)
         call FRSTAT(1, 2, 'ZON HARM')
      endif
      if(Nmtess.gt.1) then
         call LINCHK
         ll1 = 1
         do ll = 2, Nmtess
            call EBCDIX(ll,mc,3,2)
            call ADJHAR(Lmchar(ll1), Mchar(ll1), 1, 1, ll, 0, mc,
     .                  ')   ', Aplnt(17), npln10, Fcthrt(ll1))
            ll1 = ll1 + ll
         end do
         call FRSTAT(1, 2, 'COS HARM')
         call LINCHK
         ll1 = 1
         do ll = 2, Nmtess
            call EBCDIX(ll,ms,3,2)
            call ADJHAR(Lmshar(ll1), Mshar(ll1), 1, 1, ll, 0, ms,
     .                  ')   ', Aplnt(17), npln10, Fcthrt(ll1))
            ll1 = ll1 + ll
         end do
         call FRSTAT(1, 2, 'SIN HARM')
      endif
c
c adjustment to moon spot coordinates
      nspot = 0
      call FRSTAT(2, 3, 'ALL HARM')
      call ADJSPT(nspot, npln10, Aplnt(17))
      call FRSTAT(2, 3, 'MOONSPOT')
c
c adjustment to moon radar observation biases
      nrbias = 0
      call ADJRBS(nrbias, npln10, Aplnt(17))
      call FRSTAT(2, 3, 'MOONBIAS')
c
c adjustment to moon optical observations phase corrections
      nphase = 0
      call ADJPHS(nphase, npln10, Aplnt(17))
      call FRSTAT(2, 3, 'MOONFAZE')
      call FRSTAT(3, 3, Aplnt(17))
c
c adjustment to sun spot coordinates
      call ADJSPT(nspot, npln0, Aplnt(18))
      call FRSTAT(2, 3, 'SUNSPOT ')
c
c adjustment to sun radar observation biases
      call ADJRBS(nrbias, npln0, Aplnt(18))
      call FRSTAT(2, 3, 'SUNBIAS ')
c
c adjustment to sun optical observation phase corrections
      call ADJPHS(nphase, npln0, Aplnt(18))
      call FRSTAT(2, 3, 'SUNFAZE ')
      call FRSTAT(3, 3, Aplnt(18))
c
c adjustment to planet initial conditions, parameters
c and optical observation phase corrections
c*  start=1000
      do i = 1, u_nmbod
         Fctprb(i) = 1._10
         Fctpln(i) = 1._10
      end do
      Fctprb(29) = Aultsc
      Fctprb(30) = Fctprm(1)
      do i = 3, 6
         Fctpln(i) = convd1
      end do
      Fctpln(7)  = Aukm
      Fctpln(8)  = Aultsc
      Fctpln(9)  = Aultsc
      Fctpln(10) = Aultsc
      Fctpln(11) = Aultsc*convd1
      Fctpln(12) = Fctpln(11)
      Fctpln(13) = Aultsc
      Fctpln(14) = -Fctpln(11)
      Fctpln(15) = Fctpln(11)
      Fctpln(30) = Fctprm(1)
 
      do k = 1, Numpln
         if(Nplnt(k).eq.0) goto 600
 
         if(Nplnt(k).ge.0) then
            if(Nplnt(k).le.30) then
               do i = 1, u_nmbod
                  fctplb(i) = Fctpln(i)
                  Wdrs(i)   = wrds(i)
               end do
               goto 450
            endif
         endif
         do i = 1, u_nmbod
            fctplb(i) = Fctprb(i)
            Wdrs(i)   = wrdb(i)
         end do
  450    if(Nplnt(k).ge.0 .and. Icnd(k).ge.0) then
            if(Icnd(k).ge.2) then
               fctplb(1) = 1._10
               fctplb(2) = convd1
               fctplb(3) = convd1
               fctplb(4) = 1._10
               fctplb(5) = convd1
               fctplb(6) = convd1
               do i = 1, 6
                  Wdrs(i) = wrdr(i)
               end do
            else
               do i = 1, 6
                  fctplb(i) = Fctpln(i)
                  Wdrs(i)   = wrds(i)
               end do
            endif
            goto 500
         endif
         do i = 1, 6
            fctplb(i) = Fctprb(i)
            Wdrs(i)   = wrdb(i)
         end do
 
  500    call ADJBDY(Lpl(1,k), Pcond(1,k), fctplb, Wdrs, 'LPL(',
     .               wrdl2(k), Aplnt(k), Nplnt(k), Jdpl0(k), u_nmbod,
     .               Npcent(k), Icnd(k))
         call FRSTAT(1, 3, 'PLNT IC ')
c
c adjustment to planet gravitational potential harmonic coeff.
c and planet shape models
c*  start=1200
         do i = 1, 4
            if(Nplhar(i).eq.Nplnt(k)) then
               if(Nplnt(k).le.0) then
                  ibr = Nshape(i) + 1
                  if(ibr.eq.2) then
c
c fourier shape
                     ll = 122
                     call ADJFOR(i, ll, Aplnt(k), Nplhar(i))
                  else if(ibr.eq.3) then
c
c grid shape
                     ll = ngdpts(i)
                     call ADJGRD(i, ll, Aplnt(k), Nplhar(i))
                  else
                     goto 510
                  endif
                  goto 550
               endif
  510          if(Npzone(i).gt.1) then
                  call LINCHK
                  ll = Npzone(i) - 1
                  call ADJHAR(Lpzhar, Pzhar, 4, i, ll, 1, 'PZ', charz(i)
     .                        , Aplnt(k), Nplhar(i), Fcthr1)
                  call FRSTAT(1, 2, 'ZON HARM')
               endif
               if(Nptess(i).gt.1) then
                  ll2 = Nptess(i)
                  call LINCHK
                  ll1  = 1
                  ll1p = 0
c -(u_mxtes*(u_mxtes+1))/2-1 effectively uses Fcthr1 in place of Fcthrt
                  if(Nplnt(k).lt.0) ll1p = -(u_mxtes*(u_mxtes+1))/2-1
                  do ll = 2, ll2
                     call EBCDIX(ll,pc,3,2)
                     call ADJHAR(Lpchar(1,ll1), Pchar(1,ll1), 4, i, ll,
     .                           0, pc, chart(i), Aplnt(k), Nplhar(i),
     .                           Fcthrt(ll1p+ll1))
                     ll1 = ll1 + ll
                  end do
                  call FRSTAT(1, 2, 'COS HARM')
                  call LINCHK
                  ll1 = 1
                  do ll = 2, ll2
                     call EBCDIX(ll,ps,3,2)
                     call ADJHAR(Lpshar(1,ll1), Pshar(1,ll1), 4, i, ll,
     .                           0, ps, chart(i), Aplnt(k), Nplhar(i),
     .                           Fcthrt(ll1p+ll1))
                     ll1 = ll1 + ll
                  end do
                  call FRSTAT(1, 2, 'SIN HARM')
               endif
               goto 550
            endif
 
         end do
  550    call FRSTAT(2, 3, 'ALL SHAP')
c
c do not call spot and bias routines for planet rotation
         if(Nplnt(k).gt.0) then
c
c adjustment to planet spot coordinates
            call ADJSPT(nspot, Nplnt(k), Aplnt(k))
            call FRSTAT(2, 3, 'PL SPOT ')
c
c adjustment to planet radar observation biases
            call ADJRBS(nrbias, Nplnt(k), Aplnt(k))
            call FRSTAT(2, 3, 'PL BIAS ')
c
c adjustment to planet optical observation phase correction
            call ADJPHS(nphase, Nplnt(k), Aplnt(k))
            call FRSTAT(2, 3, 'PL FAZE ')
 
            call FRSTAT(3, 3, Aplnt(k))
         endif
      end do
c
c
c adjustment to star coordinates
  600 call ADJSPT(nspot, onemns, strnam)
c
c adjustment to star quantities in radar bias common
      call ADJRBS(nrbias, onemns, strnam)
c
c adjustment to star quantities in optical phase corr.common
      call ADJPHS(nphase, onemns, strnam)
c
c pulsar parameters
      if(Numpsr.gt.0) then
         do i = 1, 16
            fctpsr(i) = 1._10
         end do
         do k = 1, Numpsr
            plsr(5:8) = Sptpsr(k)
            call LINCHK
            do i = 1, 16
               l = Lpsrcn(i, k)
               if(l.le.0) goto 650
               call ADJAST(Psrcn(l,k), fctpsr(l))
               write(Iout, 610) astrik(Ntype), N, i, k, l, wrdb(l + 6),
     .                          Sptpsr(k), Psrcn(l, k), Adj, Nwv, Sig,
     .                          Fract
  610          format(1x, a1, i4, '. LPSRCN(', i2, ',', i2, ')=', i2,
     .                1x, a8, 'OF ', a4, 14x, 1pd22.15, d16.8, d22.15,
     .                d10.3, 0pf8.3)
               if(Jout.gt.0) write(Jout, 610) astrik(Ntype), N, i,
     .            k, l, wrdb(l + 6), Sptpsr(k), Psrcn(l, k), Adj, Nwv,
     .            Sig, Fract
 
c save extended-precision values for later print
               if(l.eq.6 .and. Plspr(k).ne.0._10) then
                  ipsrpr = ipsrpr + 1
                  if(ipsrpr.gt.u_mxpsr)
     .                 call SUICID('TOO MANY PULSARS IN ADJUST  ', 7)
                  Ppr4(ipsrpr)     = Sptpsr(k)
                  Ppr16(1, ipsrpr) = Plspr(k)
                  Ppr16(2, ipsrpr) = Ppr16(1, ipsrpr)
                  Ppr16(1, ipsrpr) = Ppr16(1, ipsrpr) + Psrcn(l, k)
                  Ppr16(2, ipsrpr) = Ppr16(2, ipsrpr) + Nwv
               endif
               if(Keepit) Psrcn(l, k) = Nwv
            end do
            call FRSTAT(1, 2, plsr)
  650    end do
         call FRSTAT(2, 2, 'PSR PRMS')
      endif
c
c shall we start a new page for completion message
      call PAGCHK(60, 12, 0)
c
c write completion message
c*  start=2000
      write(Iout, 700) astrik(1), Nkind(1), Eps(9), astrik(2), Nkind(2)
     .                 , Eps(10), Nkind(3)
  700 format(//1x, a1, i4, ' PARAMETERS NOT CONVERGED WITH EPS(9) =',
     .       1pe12.5, 2x, 'CRITERION'/ 1x, a1, i4,
     .       ' PARAMETERS NOT CONVERGED WITH EPS(10)=', 1pe12.5,
     .       '  CRITERION'/ i6, ' PARAMETERS CONVERGED')
      if(Jout.gt.0) write(Jout, 700) astrik(1), Nkind(1), Eps(9),
     .                        astrik(2), Nkind(2), Eps(10), Nkind(3)
c
c write statistics for fractional adjustments
      ne    = N - Nf
      fnobs = Nf
      if(Nf.eq.0) fnobs = 1._10
      a(1)  = Statf(1)/fnobs
      a(2)  = Statf(2)/fnobs
      fnobs = Statf(3)/fnobs
      a(3)  = SQRT(fnobs)
      write(Iout, 800) ne, Nf, N, a, fnobs, Statf(3)
  800 format(/i6, ' STANDARD DEVIATIONS WERE ZERO. THE REMAINING', i5,
     .       ' OF THE TOTAL OF', i5,
     .       ' ADJUSTMENTS HAD THE FOLLOWING FRACTIONAL STATISTICS'/
     .       10x, 'AVERAGE FRACT', 1pe14.5/ 5x, 'AVERAGE ABS(FRACT)',
     .       1pe14.5/1x, 'ROOT MEAN SQUARE FRACT', 1pe14.5/7x,
     .       'AVERAGE FRACT**2', 1pe14.5/11x, 'SUM FRACT**2', 1pe14.5)
      if(Jout.gt.0) write(Jout, 800) ne, Nf, N, a, fnobs, Statf(3)
      if(Ict(15).eq.0) then
 
c skip this  when using saved solution tape  jmat
         if(Ict(5).gt.1 .and. Iterat.eq.1) goto 1100
c
c calculate convergence norm and predicted residual sum of
c squares:  must read rhs back from temporary storage
         do i = 1, Nparam
            read(Ibuf)
         end do
         read(Ibuf) j, (rhs(i), i = 1, Nparam)
         rewind Ibuf
 
c normsq = delta x * rhs
         normsq = Sumzns
         do i = 1, Nparam
            normsq = normsq + rhs(i)*Solut(i)/Scale(i)
         end do
         norm   = SQRT(normsq)
         normsv = normsq
c
c new addition to correct the effect of the a priori facility
         if(Imat2.le.0 .or. Jct(60).gt.0) then
            if(Ibuf2.eq.0 .or. Ict(44).eq.0) goto 900
 
c skip title
            read(Ibuf2)

c skip priori sumsq
            read(Ibuf2)
 
c read pointer
            read(Ibuf2) nsave,(Pointr(i),i=1,nsave)
            do i = 1, nsave
               ip = Pointr(i)
               Vect(i) = Solut(ip)
            end do
            ibfap = Ibuf2
         else
 
c first, check if sne dataset includes a priori information
            ibfap = Imat2
            call FRMHED(ibfap, dum, dum, 0, ippr, 0)
            read(ibfap) (i4,i=1,5),(dum,i=1,5),iaprio
            if(iaprio.le.0) goto 850
            nsave = Nparam
            call BSKIP(ibfap, nsave)
            do i = 1, nsave
               Vect(i) = Solut(i)
            end do
         endif
 
c read v: the contribution of a priori estimate to rhs
         call QREAD(ibfap, msave, Apsol, nsave)
 
c form the dot product
         t1 = -2._10*DOTN(Vect, Apsol, nsave)
 
c read b one row at a time and form the dot product
         t2 = 0._10
         do while( .true. )
            call QREAD(ibfap, i, Apsol, nsave)
            t2 = t2 + Vect(i)*DOTN(Vect, Apsol, nsave)
            if(i.ge.nsave) then
 
c the two extra terms are now added to the norm
               normsq = normsq + t1 + t2
               goto 850
            endif
         end do
  850    rewind ibfap
         Itrwnd(ibfap) = 0
c
c pstssq = pressq - normsq
  900    pstssq = Ermeas(3) - normsq
         if(pstssq.lt.0._10) then
            call PAGCHK(60, 3, 0)
            write(Iout, 920)
            pstssq = -pstssq
  920       format(' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'/
     .             ' CHANGING SIGN OF PRDCT SUM TO AVOID NEG. SQRT.'/
     .             ' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX')
         endif
         pstrms = SQRT(ABS(pstssq)/ams)
         if(pstssq.lt.0._10) pstrms = -pstrms
         prerms = SQRT(Ermeas(3)/ams)
c
c room on this page?
         call PAGCHK(60,13,0)
c
c norm message
         totssq = Ermeas(3)+Sumaps
         write(Iout,950) Nparam,norm,normsv,Ermeas(3),prerms,
     .    pstssq,pstrms,totssq,totssq-normsv
  950    format('-CONVERGENCE NORM FOR THE ADJUSTMENT OF',i4,
     .          ' PARAMETERS'/
     .          '   NORM   ',1pd14.5/,
     .          '   NORM**2',1pd14.5/
     .          '   PREFIT SUM (O-C/ERROR)**2',1pd14.5/
     .          '   PREFIT RMS (O-C/ERROR)   ',1pd14.5/
     .          '   PRDICT SUM (O-C/ERROR)**2',1pd14.5/
     .          '   PRDICT RMS (O-C/ERROR)   ',1pd14.5/
     .          '   PREFIT TOTAL SUM-SQUARE  ',1pd14.5/
     .          '   PRDICT TOTAL SUM-SQUARE  ',1pd14.5)
         if(Jout.gt.0) write(Jout,950) Nparam,norm,normsv,Ermeas(3),
     .    prerms,pstssq,pstrms,totssq,totssq-normsv
         write(Iout,1000) Eps(15)
         if(Jout.gt.0) write(Jout,1000) Eps(15)
 1000    format('0CONVERGENCE CRITERION FOR NORM**2 IS',1pe14.5)
 
         if(Jct(60).le.0) then
            pstsum = Ermeas(1) - DOTN(Solut,Vectk,Nparam) - Sumzsm
            pstmen = pstsum/ams
            call PAGCHK(60,3,0)
            write(Iout,1020) pstsum,pstmen
 1020    format('0  PRDICT SUM  (O-C/ERROR)',1pd16.5/
     .          '   PRDICT MEAN (O-C/ERROR)',1pd16.5)
            if(Jout.gt.0) write(Jout,1020) pstsum,pstmen
         endif
      endif
c print out adjustment to inverse masses
c*  start=2500
 1100 if(imtf.ne.0) then
         if(Line.gt.49) call NEWPG
         if(Line.ne.1) then
            write(Iout, 1120)
 1120       format(' ')
            if(Jout.gt.0) write(Jout, 1120)
            Line = Line + 1
         endif
         do i = 1, imtf
            call PAGCHK(60, 1, 0)
            if(Line.le.2 .or. i.le.1) then
               write(Iout, 1130)
               if(Jout.gt.0) write(Jout, 1130)
 1130          format('0ADJUSTMENT TO INVERSE MASSES')
               Line = Line + 2
               call PAGHED(0)
            endif
            write(Iout,280) Tfc2(i),(Tfi2(j,i),j=1,4),
     .       Tfc8(i),(Tfr8(j,i),j=1,5)
            if(Jout.gt.0) write(Jout,280) Tfc2(i),(Tfi2(j,i),j=1,4),
     .       Tfc8(i),(Tfr8(j,i),j=1,5)
         end do
      endif
c
c print out adjustment to gmvary or ctvary per year
      if(igmvry.gt.0) then
         call PAGCHK(58, 2 + igmvry, 0)
         write(Iout, 1150)
 1150    format('0ADJUSTMENT TO PARAMETERS IN UNITS SPECIFIED')
         if(Jout.gt.0) write(Jout, 1150)
         call PAGHED(0)
         do i = 1, igmvry
            write(Iout,240) Rjc2(i),(Rji2(j,i),j=1,3),
     .       (Rjc8(j,i),j=1,2),(Rjr8(j,i),j=1,5)
            if(Jout.gt.0) write(Jout,240) Rjc2(i),(Rji2(j,i),j=1,3),
     .       (Rjc8(j,i),j=1,2),(Rjr8(j,i),j=1,5)
         end do
      endif
c
c printout adjustment to satellite init.cond.in kilometers
      if(Iboda.gt.0) then
         if(Line.gt.49) call NEWPG
         if(Line.ne.1) then
            write(Iout, 1120)
            if(Jout.gt.0) write(Jout, 1120)
            Line = Line + 1
         endif
         do i = 1, Iboda
            call PAGCHK(60,1,0)
            if(Line.le.2 .or. i.le.1) then
               write(Iout,1160)
               if(Jout.gt.0) write(Jout,1160)
 1160          format(
     .              '0ADJUSTMENTS TO SATELLITE INITIAL CONDITIONS (KM)')
               Line = Line + 2
               call PAGHED(0)
            endif
            write(Iout,1180) (cbod4(j,i),Jboda(j,i),j=1,3),cbod8(1,i),
     .       cbod8(2,i),Jboda(4,i),Jboda(5,i),(Boda(j,i),j = 1,5)
            if(Jout.gt.0) write(Jout,1180) (cbod4(j,i),Jboda(j,i),
     .       j=1,3),cbod8(1,i),cbod8(2,i),Jboda(4,i),Jboda(5,i),
     .       (Boda(j,i),j = 1,5)
 1180       format(1x,a1,i4,'. ',a4,i2,a4,' =',i2,1x,a8,
     .             ' OF ',a8,i3,i8,1pd22.15,1pd16.8,1pd22.15,
     .             1pd10.3,0pf8.3)
         end do
      endif
 
      if(ipsrpr.gt.0) then
         call PAGCHK(58,2+ipsrpr,0)
         write(Iout, 1200)
 1200    format('0EXTENDED-PRECISION PULSAR PERIODS (SEC)')
         if(Jout.gt.0) write(Jout, 1200)
         call PAGHED(0)
         write(Iout, 1250) (Ppr4(i), (Ppr16(j,i),j=1,2), i = 1, ipsrpr)
 1250    format(36x, a4, 1p, 2d40.31)
         if(Jout.gt.0) write(Jout, 1250)
     .                           (Ppr4(i), (Ppr16(j,i),j=1,2), i = 1,
     .                           ipsrpr)
      endif
 
      if(Jout.gt.0) then
         write(Jout, 1300) N
 1300    format(/i6, ' PARAMETERS ADJUSTED')
      endif
c     at this point file jout used to be rewound (if no prdict),
c     now jout is left open to the end of the job, allowing all
c     iterations to be written on jout and allowing jout to share
c     the same fortran unit number with mout .
c
c           print out timer information
      call EBCDIX(N, tmesg, 13, 7)
      call TIMRIT(tmesg, 8)
c
c decide if least squares iteration is to continue
      if(Iterat.le.Ict(1)) then
         if(Jct(1).ge.1) then
            if(normsq.lt.Eps(15)) Ict(1) = Iterat
         else
            if(Nkind(1)+Nkind(2).le.0) Ict(1) = Iterat
         endif
      endif
c
c*  start=9990
      return
      end
