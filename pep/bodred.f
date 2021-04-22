      subroutine BODRED(in0,nstop,jdpad)


      implicit none

c
c     m.e.ash   june 1966    subroutine bodred
c     control and data constants are initialized and read for
c     earth, moon, planets, asteroids, satellites, artificial space
c     probes, earth rotation, moon rotation and planet rotation.
c
c arguments
      integer in0,nstop,jdpad
c in0  - fortran unit for editing input stream
c nstop- cumulative count of errors
c jdpad- increment for modifying jd1-jd2 ranges

c array dimensions
      include 'globdefs.inc'

c common
      include 'fcntrlx.inc'
      include 'inodta.inc'
      include 'mnrtlb.inc'
      include 'param.inc'

c body parameters
      common/WRKCOM/ Cond(6),Con(u_nmbod-6),Con1(12),Tcon(30),Sec,
     . Int1,Int2,Eps(6),Name,Jd1,Jd0,Jd2,Itape,L(u_nmbod),I4fill,
     . Nshp,Scntrl(9),J2,J3,J4,J5,J6,J7,J8,J9,J10,J11,J12,J13,J14,
     . J15,J16,J17,J18,J19,J20,J21,J22,J23,J24,J25,J26,J27,J28,J29,
     . J30,J31,J32,J33,J34,J35,J36,J37,J38,J39,J40,J41,J42,J43,J44,
     . J45,J46,J47,J48,J49,J50,C2(2),C3(3),C4(4),C5(5),C6(6),C7(7),
     . C8(8),C9(9),C10(10),C11(11),C12(12),C13(13),C14(14),C15(15),
     . C16(16),C17(17),C18(18),C19(19),C20(20),C21(21),C22(22),
     . C23(23),C24(24),C25(25),C26(26),C27(27),C28(28),C29(29),
     . C30(30),C31(31),C32(32),C33(33),C34(34),C35(35),C36(36),
     . C37(37),C38(38),C39(39),C40(40),C41(41),C42(42),C43(43),
     . C44(44),C45(45),C46(46),C47(47),C48(48),C49(49),C50(50),
     . S2(2),S3(3),S4(4),S5(5),S6(6),S7(7),S8(8),S9(9),S10(10),
     . S11(11),S12(12),S13(13),S14(14),S15(15),S16(16),S17(17)

      common/WRKCOM/ S18(18),S19(19),S20(20),S21(21),S22(22),S23(23),
     . S24(24),S25(25),S26(26),S27(27),S28(28),S29(29),S30(30),
     . S31(31),S32(32),S33(33),S34(34),S35(35),S36(36),S37(37),
     . S38(38),S39(39),S40(40),S41(41),S42(42),S43(43),S44(44),
     . S45(45),S46(46),S47(47),S48(48),S49(49),S50(50),
     . Lj2,Lj3,Lj4,Lj5,Lj6,Lj7,Lj8,Lj9,Lj10,Lj11,Lj12,Lj13,Lj14,
     . Lj15,Lj16,Lj17,Lj18,Lj19,Lj20,Lj21,Lj22,Lj23,Lj24,Lj25,Lj26,
     . Lj27,Lj28,Lj29,Lj30,Lj31,Lj32,Lj33,Lj34,Lj35,Lj36,Lj37,Lj38,
     . Lj39,Lj40,Lj41,Lj42,Lj43,Lj44,Lj45,Lj46,Lj47,Lj48,Lj49,Lj50,
     . Lc2(2),Lc3(3),Lc4(4),Lc5(5),Lc6(6),Lc7(7),Lc8(8),Lc9(9),
     . Lc10(10),Lc11(11),Lc12(12),Lc13(13),Lc14(14),Lc15(15),
     . Lc16(16),Lc17(17),Lc18(18),Lc19(19),Lc20(20),Lc21(21),
     . Lc22(22),Lc23(23),Lc24(24),Lc25(25),Lc26(26),Lc27(27),
     . Lc28(28),Lc29(29),Lc30(30),Lc31(31),Lc32(32),Lc33(33),
     . Lc34(34),Lc35(35),Lc36(36),Lc37(37),Lc38(38),Lc39(39)

      common/WRKCOM/ Lc40(40),Lc41(41),Lc42(42),Lc43(43),Lc44(44),
     . Lc45(45),Lc46(46),Lc47(47),Lc48(48),Lc49(49),Lc50(50),Ls2(2),
     . Ls3(3),Ls4(4),Ls5(5),Ls6(6),Ls7(7),Ls8(8),Ls9(9),Ls10(10),
     . Ls11(11),Ls12(12),Ls13(13),Ls14(14),Ls15(15),Ls16(16),
     . Ls17(17),Ls18(18),Ls19(19),Ls20(20),Ls21(21),Ls22(22),
     . Ls23(23),Ls24(24),Ls25(25),Ls26(26),Ls27(27),Ls28(28),
     . Ls29(29),Ls30(30),Ls31(31),Ls32(32),Ls33(33),Ls34(34),
     . Ls35(35),Ls36(36),Ls37(37),Ls38(38),Ls39(39),Ls40(40),
     . Ls41(41),Ls42(42),Ls43(43),Ls44(44),Ls45(45),Ls46(46),
     . Ls47(47),Ls48(48),Ls49(49),Ls50(50),
     . Denptr,K(u_nmprm),Int,Nplnt,Ncentr,Nzone,Ntess,Ihr,Imin,
     . Kk(100),Icnd,Numki,Ki(99),Small,Pulsar

      character*8 Name
      real*10 Cond,Con,Con1,Tcon,Sec,J2,J3,J4,J5,J6,J7,J8,J9,
     .        J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,J20,
     .        J21,J22,J23,J24,J25,J26,J27,J28,J29,J30,J31,
     .        J32,J33,J34,J35,J36,J37,J38,J39,J40,J41,J42,
     .        J43,J44,J45,J46,J47,J48,J49,J50,C2,C3,C4,C5,
     .        C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,
     .        C18,C19,C20,C21,C22,C23,C24,C25,C26,C27,C28,
     .        C29,C30,C31,C32,C33,C34,C35,C36,C37,C38,C39,
     .        C40,C41,C42,C43,C44,C45,C46,C47,C48,C49,C50,
     .        S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,S13,S14,
     .        S15,S16,S17,S18,S19,S20,S21,S22,S23,S24,S25,
     .        S26,S27,S28,S29,S30,S31,S32,S33,S34,S35,S36,
     .        S37,S38,S39,S40,S41,S42,S43,S44,S45,S46,S47,
     .        S48,S49,S50
      real*4    Eps,Scntrl
      logical*4 Small,Pulsar
      integer*4 Int1,Int2,Jd1,Jd0,Jd2,Itape,I4fill,Nshp
      integer*2 L,Lj2,Lj3,Lj4,Lj5,Lj6,Lj7,Lj8,Lj9,Lj10,Lj11,
     .          Lj12,Lj13,Lj14,Lj15,Lj16,Lj17,Lj18,Lj19,Lj20,
     .          Lj21,Lj22,Lj23,Lj24,Lj25,Lj26,Lj27,Lj28,Lj29,
     .          Lj30,Lj31,Lj32,Lj33,Lj34,Lj35,Lj36,Lj37,Lj38,
     .          Lj39,Lj40,Lj41,Lj42,Lj43,Lj44,Lj45,Lj46,Lj47,
     .          Lj48,Lj49,Lj50,Lc2,Lc3,Lc4,Lc5,Lc6,Lc7,Lc8,
     .          Lc9,Lc10,Lc11,Lc12,Lc13,Lc14,Lc15,Lc16,Lc17,
     .          Lc18,Lc19,Lc20,Lc21,Lc22,Lc23,Lc24,Lc25,Lc26,
     .          Lc27,Lc28,Lc29,Lc30,Lc31,Lc32,Lc33,Lc34,Lc35,
     .          Lc36,Lc37,Lc38,Lc39,Lc40,Lc41,Lc42,Lc43,Lc44,
     .          Lc45,Lc46,Lc47,Lc48,Lc49,Lc50,Ls2,Ls3,Ls4,Ls5,
     .          Ls6,Ls7,Ls8,Ls9,Ls10,Ls11,Ls12,Ls13,Ls14,Ls15,
     .          Ls16,Ls17,Ls18,Ls19,Ls20,Ls21,Ls22,Ls23,Ls24,
     .          Ls25,Ls26,Ls27,Ls28,Ls29,Ls30,Ls31,Ls32,Ls33,
     .          Ls34,Ls35,Ls36,Ls37,Ls38,Ls39,Ls40,Ls41,Ls42,
     .          Ls43,Ls44,Ls45,Ls46,Ls47,Ls48,Ls49,Ls50,
     .          Denptr,K,Int,Nplnt,Ncentr,Nzone,Ntess,
     .          Ihr,Imin,Kk,Icnd,Numki,Ki

      real*10 zh(u_mxzon-1),ch((u_mxtes*(u_mxtes+1))/2-1),
     . sh((u_mxtes*(u_mxtes+1))/2-1)
      equivalence (zh,J2),(ch,C2),(sh,S2)
      integer*2 lzh(u_mxzon-1),lch((u_mxtes*(u_mxtes+1))/2-1),
     . lsh((u_mxtes*(u_mxtes+1))/2-1)
      equivalence (lzh,Lj2),(lch,Lc2),(lsh,Ls2)

      real*10 a,e,inc,asc,per,anom,radius
      equivalence (Cond,a),(Cond(2),e),(Cond(3),inc),
     .            (Cond(4),asc),(Cond(5),per),(Cond(6),anom)
      equivalence (Con,radius)

      integer*2 jtype

      real*10 alpha,delta,beta,angle,r,v
      equivalence (alpha,Cond(2)),(delta,Cond(3)),(beta,Cond(6)),
     .            (angle,Cond(5)),(r,Cond),(v,Cond(4))
      equivalence (Nshp,ntype)

c           --- harmonic+shape storage ---
c        for nmlst2 with positive nplnt, gravitational
c        potential harmonics is assumed.
c        if a -nplnt nmlst2 is input, this is assumed
c        to be planet shape. model depends on nshp.
c        nshp=0,(default) spherical harmonic expansion
c        nshp=1,          fourier series expansion
c        nshp=2,          altitude grid (local model)
c        nshp=3,          external model expected to be supplied on obslibs
c        (nothing for nshp .gt.3 yet)
c        note: see later comments on shape for full
c        explanation.
c
      real*10 j(50)
      real*4    shape(1000), grid(1000)
      integer*2 lj(50),lshape(1000),lgrid(1000)
      equivalence (J2,j(2),shape,grid), (Lj2,lj(2),lshape,lgrid)

c        overlay equivalencing for fourier series coefs
c          s.brody april 1975-coding for two dimensional fourier series
c
c          ovely(lovely)equivalenced starting with j2(lj2)
c          for the next 252 real*8(integer*2)elements.
c *** historical reasons for 252 elements. doesn't hurt. ***
c          fa thru fbpp variables are the corresponding
c          fourier coefficients for the two dimensional fourier
c          series representation (note:  read p as primed and
c          pp as double primed). lfa thru lfbpp are the
c          corresponding l-vectors. see later comments on
c          shape for full explanation.
c
      integer*2 lfa(20),lfb(20),lfc(20),lfd(20),lfap(20),lfcp(20),
     .          lfapp,lfbpp
      integer*2 lovely(252)
      real*10 fa(20),fb(20),fc(20),fd(20),fap(20),fcp(20),fapp,
     .          fbpp
      real*10 ovely(252)
      equivalence (ovely,J2),(lovely,Lj2)
      equivalence (ovely,fa),(ovely(21),fb),
     .            (ovely(41),fc),(ovely(61),fd),
     .            (ovely(81),fap),(ovely(101),fcp),
     .            (ovely(121),fapp),(ovely(122),fbpp)
      equivalence (lovely,lfa),(lovely(21),lfb),
     .            (lovely(41),lfc),(lovely(61),lfd),
     .            (lovely(81),lfap),(lovely(101),lfcp),
     .            (lovely(121),lfapp),(lovely(122),lfbpp)

      integer*4 zwrkcm/48700/   !r8=27340,r10=48700


      integer*2 imn0, idy0, iyr0, imn1, idy1, iyr1, imn2, idy2, iyr2
c     imnk,idyk,iyrk = month,day,year converted to jdk (k=0,1,2)
c                      if iyrk.gt.0
c                     year is 0 to 99 in 1900s
c                     year is greater than or equal to 100 in 2000s
c
      real*4    tlat(2), tlon(2), tlatin, tlonin
      equivalence (tlat, Scntrl), (tlon, Scntrl(3)),
     .            (tlatin, Scntrl(5)), (tlonin, Scntrl(6)),
     .            (latdim, Scntrl(7)), (londim, Scntrl(8)),
     .            (ngdpts, Scntrl(9))

      integer*4 sclzon, scltes

c     sclzon = 0 unscaled zonal harmonics are input (default)
c     sclzon = 1 scaled zonal harmonics are input (change to unscaled
c                before return to calling program)
c     scltes = 0 unscaled tesseral harmonics are input (change to scaled
c                before return to calling program)
c     scltes = 1 scaled tesseral harmonics are input (default)


      namelist/NMLST2/Cond,Con,Con1,Eps,Name,Jd1,Jd0,Jd2,L,K,Int,Nplnt,
     .Ncentr,Itape,a,e,inc,asc,per,anom,radius,j,J2,J3,J4,J5,J6,J7,J8,
     .J9,J10,J11,J12,J13,J14,J15,J16,J17,J18,J19,J20,J21,J22,J23,J24,
     .J25,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,C18,
     .C19,C20,C21,C22,C23,C24,C25,S2,S3,S4,S5,S6,S7,S8,S9,S10,S11,S12,
     .S13,S14,S15,S16,S17,S18,S19,S20,S21,S22,S23,S24,S25,lj,Lj2,Lj3,
     .Lj4,Lj5,Lj6,Lj7,Lj8,Lj9,Lj10,Lj11,Lj12,Lj13,Lj14,Lj15,Lj16,Lj17,
     .Lj18,Lj19,Lj20,Lj21,Lj22,Lj23,Lj24,Lj25,Lc2,Lc3,Lc4,Lc5,Lc6,Lc7,
     .Lc8,Lc9,Lc10,Lc11,Lc12,Lc13,Lc14,Lc15,Lc16,Lc17,Lc18,Lc19,Lc20,
     .Lc21,Lc22,Lc23,Lc24,Lc25,Ls2,Ls3,Ls4,Ls5,Ls6,Ls7,Ls8,Ls9,Ls10,
     .Ls11,Ls12,Ls13,Ls14,Ls15,Ls16,Ls17,Ls18,Ls19,Ls20,Ls21,Ls22,Ls23,
     .Ls24,Ls25,Nzone,Ntess,Int1,Int2,Ihr,Imin,Sec,Kk,Icnd,sclzon,
     .scltes,jcnd,kcnd,lcnd,jtype,alpha,delta,beta,angle,r,v,
     .mean,mcentr,mcnd,fract0,iftkm,lfa,lfb,
     .lfc,lfd,lfap,lfcp,lfapp,lfbpp,fa,fb,fc,fd,fap,fcp,fapp,fbpp,grid,
     .lgrid,tlat,tlon,tlatin,tlonin,Nshp,imn0,idy0,iyr0,imn1,idy1,iyr1,
     .imn2,idy2,iyr2,jdtype,incnd,inunit,Tcon,Numki,Ki,Small,Denptr,
     .Pulsar,ntype,zh,ch,sh,lzh,lch,lsh


c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c NPLNT = body number
c            1 = Mercury
c            2 = Venus
c            3 = Earth-Moon barycenter
c            4 = Mars
c            5 = Jupiter
c            6 = Saturn
c            7 = Uranus
c            8 = Neptune
c            9 = Pluto
c           10 = Moon
c    11,...,30 = other possible natural planets,asteroids or satellites
c    31,32,... = artificial space probes
c    -1,-2,-3,-4,...,-9,-10,-11,... = rotation for planet iabs(NPLNT)
c            0 = indicates end of sequence of &nmlst2 namelists,
c                subroutine input gets out of loop calling subroutine
c                bodred. in other parts of program, planet number 0 is
c                associated with the sun.
c --- However, if SMALL=t, then NPLNT is an asteroid number,
c     starting with 1 for Ceres, 2 for Juno, etc.  Asteroids may be
c     input either as ordinary minor planets with NPLNT=11 to 30 or
c     as limited asteroids with any NPLNT number.
c
c NCENTR = central body of numerical integration
c         -1,0 = Sun
c            1 = Mercury
c            2 = Venus
c            3 = Earth (not Earth-Moon barycenter)
c            4,5,...same as above
c
c if NPLNT between 1 and 30 and NCENTR=0, individual body integration
c            with sun as central body is done in subroutine FN
c if NPLNT between 1 and 30 and NCENTR.lt.0, individual body integration
c            with sun as central body is done in subroutine SBFN
c if NPLNT.gt.30 or NCENTR.gt.0, individual body integration with
c            central body NCENTR is done in subroutine SBFN
c if NPLNT=10, moon integration is always done in MONFN
c if NPLNT= 3, earth-moon barycenter integration is done in FN or SBFN
c            depending on whether ICT(40)=0 or 1
c There is no indivdual body integration if JD0=0 or if the specific
c            body is one of the NPLBDY(.) bodies in the n-body /BDCTRL/
c            labeled common
c
c (mass of probe NPLNT)/(mass of central body NCENTR)=0 if
c              NPLNT.gt.30
c if NCENTR.le.0, (mass of body NPLNT)/(mass of sun)=MASS(NPLNT) if
c              NPLNT between 1 and 30 (MASS vector is in /PARAM/
c              labeled common)
c if NCENTR.gt.0, (mass of body NPLNT)/(mass of planet-satellite system
c                 NCENTR) =MASS(NPLNT) if NPLNT between 1 and 30
c Note that when we integrate the motion of Jupiter, for instance, we
c are actually integrating for the motion of the center of mass of the
c Jovian planet-satellite system relative to the sun and that mass(5) is
c the total mass of the Jovian planet-satellite system divided by the
c mass of the sun
c
c
c NAME = 8 characters giving name of body NPLNT
c
c
c   JD0 =initial epoch of NPLNT integration is midnight beginning of
c        day of day with julian day number   JD0  (if negative there is
c        checkpoint restart at the epoch   -JD0  unless
c          JD0 .le.-3000000. if JD0=-1 there is checkpoint restart
c          at record just before end of file on ephemeris tape.)
c        (if JD0=0 or JD0=-3000000, there is no
c        numerical integration, comparison of theory and observation
c        assumes motion is already on tape, but if initial conditions
c        are adjusted,   JD0  is set to positive value on second record
c        of tape so there will be reintegration on next least squares
c        iteration, whereas if no initial conditions are adjusted,
c          JD0  is set to zero in the compar link so there will be no
c        reintegration on subsequent least squares iterations)
c        (  JD0 .le.-3000000 is signal to compar link on first
c        iteration to override end time on tape by input   JD2 )
c
c        Initial time is midnight ephemeris time on day with Julian day
c        number JD0 (Julian date at noon on given day), except that
c        integrations in SBFN (NPLNT.gt.30 or NCENTR.ne.0) can have
c        initial epoch at day JD0 and nonzero fraction of a day
c        INT1*(2**INT2). INT2 negative and less than 30 in absolute
c        value. INT1 positive and less than 2**iabs(INT2).  Initial
c        epoch time can be input as IHR,IMIN,SEC ephemeris time or as
c        FRACT0 fraction of day, rather than INT1,INT2, but it will be
c        converted to INT1 with INT2=-24 (accurate to .01 second)
c
c   JD1 =julian day number of first data record on output tape
c   JD2  julian day number of last data record on output tape
c   (see also eps(1) and eps(2) below)
c   see also description of jdpad in PRMRED
c
c must have JD0 lying between JD1 and JD2, unless JD0 is
c        negative indicating checkpoint restart, in which case   -JD0
c        must be between JD1 and JD2, or unless JD0 is zero
c        indicating no integration or less than or equal to -3000000
c        indicating no integration and input JD2 to override that on
c        second record of integration data set in COMPAR link (needed
c        because such data set might have been produced by check point
c        restart going beyound point indicated on second record)
c        (if JD0=-1, consistency check for -JD0 between JD1,JD2 not made
c
c time direction of output tape is from JD1 to JD2, forward in
c        time if JD1 .lt. JD2, backward in time if JD1 .gt. JD2
c
c if JD0 .lt.0 & .gt.-3000000, output tape from a previous integration
c        is read to the time -JD0 (or to the record just before end
c        of file if JD0=-1 ) and the coordinates on the tape at
c        that epoch are used to restart that integration, proceeding
c        from -JD0 to JD2. Integration cannot go in both
c        directions from epoch -JD0 in restart mode and a
c        previous integration cannot be restarted at a time in the
c        first direction of the previous integration if it was in the
c        two direction mode, since output tape is not written until it
c        gets to second direction. In the checkpoint restart mode with
c        JD0 negative, the only data which can be input for the
c        integration are JD1,JD0,JD2,K(99),EPS(2),EPS(6) because rest
c        of the integration data is taken from the second record of the
c        tape of the previous integration.  If K(99) .lt.0, the tape
c        of the previous integration is rewound, whereas if
c        K(99) .ge.0, the restarted integration continues writing
c        the tape of the previous integration.
c
c if JD0 is positive between JD1 and JD2, we have
c       (a) if JD0=JD1, numerical integration starts at JD0
c           and proceeds to JD2, writing output tape as integration
c           proceeds if K(99).ge.0
c       (b) if JD0 is strictly between JD1 and JD2, integration
c           starts at JD0 and goes in direction of JD1 (forward
c           in time if JD0 .lt. JD1, backward in time if JD0 .gt.
c           JD1). The integration polynomial coefficients determined
c           at JD0 by the starting proceedure are saved and the
c           numerical integration proceeds towards JD1, writing onto
c           disk buffer IBUF if K(99).ge.0.  When JD1 is reached,
c           IBUF is backspaced and read to write output tape ITAPE in
c           time direction from JD1 to JD0 if K(99).ge.0.  Then
c           using the saved starting polynomial coefficients, the
c           integration goes from JD0 to JD2, writing output
c           tape ITAPE as the integration proceeds if K(99).ge.0.
c           When the integration is completed, we have an output tape
c           from JD1 to JD2, even though integration went from
c           JD0 to JD1 and then from JD0 to JD2. These
c           contortions are gone through in order to simplify the logic
c           needed to read the tape in the COMPAR link.
c
c           Note: in a simultaneous integration of lunar orbit,
c           libration, and core rotation, three separate buffers are
c           needed for the two-direction mode.  IBUF is used for the
c           orbit, while KK(13) and KK(12) of the Moon rotation are
c           used for the libration and core rotation, respectively.
c
c
c INT  gives interval between tabular points.
c if  INT.gt.0, tabular interval is  INT days (only value allowed
c               if NPLNT=1,...,9,11,...,30 and NCENTR=0)
c if  INT.le.0, tabular interval is 2**INT days (for the Moon, NPLNT=10,
c               it is expected that INT=-1 for half day tabular intervl)
c These constant output tabular intervals apply throughout the
c integration for integration methods K(88)=1,2,3. for K(88)=0 method,
c output interval from initial epoch to first tabular point is
c controlled by INT. Thereafter, tabular output is every step of
c variable step Nordsieck integration. Output every step immediately
c after the initial epoch is supressed because the starting interval for
c the Nordsieck method of integration is very small, but it grows to the
c size warranted by the stability requirements of the differential
c equations integrated and the input accuracy constant EPS(3) as the
c integration proceeds away from the initial epoch.
c
c
c ITAPE =   data set number for output tape
c           first record is 88 characters describing computer run,second
c           record is data about body integration, including
c           IPAR=number of partials+1 (calculated in setup for integ.)
c           third and subsequent records have integration results.
c
c           (1) if NPLNT=1,...,9,11,...,30 with NCENTR=0 there are five
c           tabular points per record with format
c               JD,FRACT,IVL,(((CRD(i,j,k),i=1,IVL),j=1,IPAR),k=1,5)
c           where JD is Julian day number of the data at the first
c           tabular point,FRACT is double precision fraction of day from
c           midnight beginning of day for this JD (0 if INT.gt.0)
c           and where for k=1,5 the double precision CRD array is
c           CRD(i,1,k) = position,velocity at tabular point k, i=1,6
c           CRD(i,j,k) = partial derivative of position,velocity with
c                        respect parameter j at tabular point k,
c                        i=1,6 and j=2,IPAR
c           if  IPAR=1, no partials on tape
c           IPAR  calculated in integration setup routine from KI(1-n)
c           the meaning of the partials j=2,ipar  is given by  KI(1-n).
c           IVL=6 in integration output routine so velocity will be on
c           tape for checkpoint restart purposes. To save space on tape
c           or disk output can be copied with IVL=3 before recent times
c           since optical observations do not need accurate velocity
c           and the velocity they need is gotten by numerical
c           differentiation from position tabulation.
c           Everett eighth difference interpolation is used to
c           interpolate from this tape in the COMPAR link. There are
c           three records in storage in the COMPAR link and there are
c           five tabular points per record so that can interpolate for
c           any time in the middle record (need more than four tabular
c           points per record because of the existence of receive,
c           reflection and send times of an observation)
c
c           (2) if NPLNT=10 for the Moon there are eight tabular points
c           per record with format
c               JD,FRACT,IVL,(((CRD(i,j,k),i=1,IVL),j=1,IPAR),k=1,8),
c                  nutation,libration,k=1,8)
c           where the meaning of JD,FRACT,IVL,CRD is the same as in
c           case(1) with FRACT =0 if INT=-1,-2,-3, where the
c           nutation of the Earth and the physical libration of the
c           Moon are copied from the perturbing planet data set for
c           use in the COMPAR link.
c
c           (3) if NPLNT=1,...,9,11,...,30 with NCENTR=-1 or
c           NCENTR.gt.0, or if NPLNT.gt.30 there are two cases.
c           if K(88).gt.0, it is the same as case(1) above.
c           if K(88)=0 (Nordsieck variable output interval) there is
c           one tabular point per record with format
c               JD,FRACT,IORD,JORD,NUMPRT,NUMPT1,((SPRB(i,j),i=1,3),
c                 j=1,IORD),(((DSPRB(i,j,k),i=1,3),j=1,JORD),
c                 k=1,NUMPT1),HC1,HC2
c           where meaning of JD,FRACT is the same as above and where
c           for the given tabular point at time JD,FRACT
c           IORD  = order of time derivative - 1 for motion on tape
c           JORD  = order of time derivative - 1 for partial on tape
c           NUMPRT= number of partials on tape
c           NUMPT1=NUMPRT if NUMPRT.gt.0, 1 otherwise
c           SPRB(i,1)= postion coordinate i, i=1,2,3
c           SPRB(i,2)= velocity coordinate
c           SPRB(i,3)= acceleration coordinate
c           SPRB(i,4)= jerk coordinate
c                   etc. up to j=IORD
c           DSPRB(i,j,k)= similar meaning for partial of motion with
c                         respect to parameter k
c           HC1   = time distance from previous tabular point to
c                   this tabular point
c           HC2   = time distance from this tabular point to next
c                   next tabular point
c           There are six records in storage in compar link for Hermite
c           interpolation for a time point in the midst of the middle
c           four tabular points so that either a two or four point
c           Hermite interpolation method can be used of order limited
c           only by the number of derivatives on the tape and by our
c           being able to derive the formulas. The derivatives after
c           acceleration come from the Nordsieck intgration coefficients
c           which are not the real higher derivatives in so far as
c           they tend to lag behind the current tabular point, the jerk
c           lagging least of all.
c
c
c  COND(i), i=1,6,= initial conditions for equations of motion
c    (a) Normally, COND(1-6) are the
c        initial osculating eliptic orbital elements with non-time
c        variable gravitational constant and with augmented mass of
c        central body (so that equations for partial of motion with
c        respect to given body mass will have non-zero init.cond.)
c            COND(1)=A    =semi-major axis (astronomical units)
c            COND(2)=E    =eccentricity
c            COND(3)=INC  =inclination (deg)
c            COND(4)=ASC  =right ascension of ascending node (deg)
c            COND(5)=PER  =argument of perihelion   (deg)
c            COND(6)=ANOM =initial mean anomaly (deg)
c    (b) If NPLNT.lt.0 then COND(1-6) are initial Euler angles and
c            rates (rad and rad/day).
c    (c) Input controls may specify other types of initial conditions.
c The equations of motion are integrated in a coordinate system
c referred to the mean equinox and equator of the reference epoch.
c (See JCT(13).)
c
c The above conventions for initial conditions are superceded by the
c following (default value of icnd is given by above)
c ICND =-1 initial conditions are cartesian position and velocity
c          (units are astronomical units and astronomical units per day)
c          (referred to mean equinox and equator of reference epoch)
c ICND = 0 initial conditions are elliptic orbital elements
c ICND = 1 initial conditions are elliptic orbital elements with sums of
c          angles   ASC, ASC+PER, ASC+PER+ANOM
c ICND = 2 initial conditions are
c          COND(1) = distance from central body (astronomical units)
c          COND(2) = right ascension of position (epoch coordinates)
c                    (between 0 and 360 deg)
c          COND(3) = declination of body (epoch coordinates)
c                    (between -90 and 90 deg)
c          COND(4) = magnitude of velocity minus the velocity in
c                    circular orbit at radius COND(1) (au/day)
c          COND(5) = azimuth of velocity vector on plane normal to
c                    radius vector (between 0 and 360 deg)
c          COND(6) = flight path angle measured from radius vector to
c                    velocity vector (between 0 and 180 deg)
c          ADBARV
c                a = right ascension = COND(2)
c                d = declination     = COND(3)
c                b = beta            = COND(5)
c                a = angle           = COND(6)
c                r = radius          = COND(1)
c                v = velocity        = COND(4)+circular orbit velocity
c ADBARV initial conditions can be read in &nmlst2 namelist by the names
c ALPHA,DELTA,BETA,ANGLE,R,V with appropriate values for ICND,KCND,LCND
c and also JTYPE
c
c The initial conditions must be in one of the above forms in the input
c stream, but several variants are allowed, namely, distance and vel-
c ocity may be in other units, the reference frame may be other than
c 1950 mean equatorial, and the velocity 'v' of adbarv notation may be
c the total velocity.  The variants, if any, are all resolved to the
c standard forms in the input procedure.  For convenience, the elements
c may be input in one form (say, Cartesian) and converted to another
c (say, elliptic).  The desired final form is always specified via
c 'ICND' (as above), and the input form, if different, via 'INCND' with
c (with the same logic as 'ICND' above).
c
c INCND=-1 input elements are Cartesian
c INCND= 0 input elements are elliptic
c INCND= 1 input elements are elliptic with angle sums
c INCND= 2 input elements are ADBARV
c
c INUNIT=-1 input distances, velocities are au, au/40 d
c INUNIT= 0 input distances, velocities are au, au/d
c INUNIT= 1 input distances, velocities are km, km/s
c INUNIT= 2 input distances, velocities are ft, ft/s
c
c LCND=0 normal ADBARV notation
c LCND=1 input ADBARV (if any) has v=total velocity, not excess
c
c JTYPE specifies the reference frame of input elements.  The default
c is -1 (1950 mean equinox and equator).  See below for other options.
c
c For Moon rotation there is special logic that applies only to the
c core's angle rates:
c INCND=0 input rates are Euler angle rates
c INCND=3 input rates are angular velocity components of the core in
c         the frame of the mantle (to be converted to "0" if ICND=0)
c         If ICND=3, then retain these values and integrate only the
c         core angular velocity, ignoring what would otherwise be the
c         Euler angles
c
c * * * * * obsolete * * * * * see 'INUNIT' and 'INCND' above
c ICND switches internal to BODRED are
c ICND =-8 change input position and velocity in km and km/sec to
c          the ICND=1 initial conditions
c ICND =-7 change input position and velocity in km and km/sec to
c          the ICND=0 initial conditions
c ICND =-6 change input position and velocity in au and au/day to
c          the ICND=1 initial conditions
c ICND =-5 change input position and velocity in au and au/day to
c          the ICND=0 initial conditions
c ICND =-4 same as ICND=1 but unit of COND(1) is kilometers instead of
c          au.  COND(1) is changed to astronomical units (au) and icnd
c          is set equal to 1
c ICND =-3 same as ICND=0 but unit of COND(1) is kilometers instead of
c          au.  COND(1) is changed to astronomical units (au) and icnd
c          is set equal to 0
c ICND =-2 same as INCD=-1 but units of position and velocity are
c          kilometers and kilometers per second instead of au and au/day
c          COND(1-6) are transformed to au and au/day and ICND is set
c          equal to -1
c ICND = 3 same as ICND=2 but units for COND(1) & COND(4) are kilometers
c          and km/sec instead of au and au/day, respectively.
c          COND(1) is transformed to au and COND(4) is transformed to
c          au/day and icnd is set equal to 2
c ICND = 4 change input position and velocity in au and au/day to
c          the ICND=2 initial conditions in au and au/day and set
c          ICND equal to 2
c ICND = 5 change input position and velocity in km and km/sec to
c          the ICND=2 initial conditions in au and au/day and set
c          ICND equal to 2
c * * * * * obsolete * * * * *
c
c * * * * * obsolete * * * * * see 'INUNIT' and 'INCND' above
c JCND switches internal to BODRED for preprocessing of initial
c conditions before ICND preprocessing is done
c JCND = 0 no element preprocessing in subroutine CHNCNE
c JCND = 1 change initial osculating orbital elements (ICND=0 type)
c          to Cartesian coordinates (semi-major axis unit is
c          astronomical units)
c JCND = 2 same as JCND=1 but semi-major axis unit is kilometers
c          (change to astronomical units)
c * * * * * obsolete * * * * *
c
c * * * * * obsolete * * * * * see 'INUNIT' and 'INCND' above and JTYPE
c KCND switches internal to BODRED for preprocessing of initial
c conditions before JCND,ICND preprocessing is done
c KCND = 0 no coordinate system preprocessing in CHNCNE
c KCND = 1 change initial osculating orbital elements (ICND=0 type)
c          from Euler angles inc,asc,per referred to true
c          equinox and equator of date to mean equinox and equator
c          of reference epoch
c KCND = 2 change initial cartesian coordinates referred to true equinox
c          and equator of date to mean reference equinox and equator
c KCND = 3 change initial ADBARV coordinates referred to true equinox
c          and equator of date to initial cartesian coordinates referred
c          to the mean reference equinox and equator (ICND.gt.2 can
c          change back to ADBARV)
c          (units are au & au/day)
c KCND = 4 units changed from km & km/sec and then KCND=3 is done
c * * * * * obsolete * * * * *
c
c JTYPE overrides above change to mean equinox and equator of epoch.
c Instead of true equinox and equator of date, original coordinates are
c JTYPE =-1 mean equinox and equator of 1950.0
c JTYPE = 0 mean equinox and equator of date
c JTYPE = 1 true equinox and equator of date (default assumed if
c           obsolete KCND is given)
c JTYPE = 2 mean equinox and ecliptic of date
c JTYPE = 3 mean equinox and ecliptic of 1950.0
c JTYPE = 4 mean lunar plane of date (x axis along intersection of mean
c           lunar and ecliptic planes of date)
c JTYPE = 5 mean equinox and true equator of date JDTYPE (0 hr)
c JTYPE = 6 mean equinox and equator of date JDTYPE (12 hr)
c JTYPE = 7 mean equinox and ecliptic of date JDTYPE (12 hr)
c
c JDTYPE defaults to the appropriate epoch according to JCT(13)
c
c LCND = 0 nothing
c LCND = 1 input ADBARV COND(4) is velocity and not difference between
c          velocity and circular orbit velocity  (units au & au/day)
c          change to difference velocity
c LCND = 2 change from km & km/sec to au & au/day and then do LCND=1
c          (obsolete, see 'INUNIT' above)
c
c * * * * * obsolete * * * * * see 'INUNIT' and 'INCND' above
c Before any of the above are done, the following is done.
c IFTKM = 0 nothing
c IFTKM = 1 for elliptic orbital element initial condtions change semi-
c           major axis from feet to kilometers
c IFTKM = 2 for ADBARV initial conditions change radius COND(1) from
c           feet to kilometers and velocity COND(4) from ft/sec to
c           km/sec
c IFTKM = 3 for cartesian initial conditions change position and
c           velocity COND(1-6) from ft and ft/sec to km and km/sec
c * * * * * obsolete * * * * *
c
c MEAN = 0 just the usual &nmlst2 namelist read for this body
c MEAN = 1 read the NORAD 2 card mean element set after this &nmlst2.
c          The mean elements are changed to Cartesian coordinates, and
c          set JTYPE=0, INUNIT=0, INCND=-1.  ICND should have been set
c          to the desired output type of orbital elements (in this case
c          the default is -1).
c
c MCENTR =-2 NCENTR gives the central body, no changes made
c MCENTR.ge.-1 NCENTR gives the central body that the initial conditions
c             are for. Change initial conditions to the MCENTR central
c             body and set NCENTR=MCENTR. Before this is done, by use of
c             ICND,JCND,KCND,LCND,JTYPE we must have ended up at the
c             ICND=-1 type of initial conditions (Cartesian coordinates
c             referred to the mean reference equinox and equator in au
c             and au/day). After the change of central body, the initial
c             conditions can be changed to the ICND=0,1 or 2 type by
c             MCND switch, which does not apply unless MCENTR.ge.0
c             the preferred values of MCND are (-1,0,1,2), but the
c             following (obsolete) values from the (obsolete) ICND
c             conventions can also be used:
c MCND =-6 change ICND=-1 initial conditions to ICND=1 i.c.
c MCND =-5 change ICND=-1 initial conditions to ICND=0 i.c.
c MCND = 4 change ICND=-1 initial conditions to ICND=2 i.c.
c
c
c CON(1) to CON(24) are adjustable parameters belong to given body.
c CON1(1) to CON1(12) are non-adjustable parameters belonging to
c given body.
c
c For all bodies (NPLNT positive)
c     CON(24) = adjustable relativity motion factor belonging to body
c
c For a planet, asteroid or natural satellite (NPLNT.le.30 and
c NPLNT not equal 3 or 10 nor negative)
c     CON(1)=radius of body
c     CON(i), i=1,12,  adjustable parameters giving shape and rotation
c                     of planet
c     CON(2),CON(3),CON(4) = j1,c11,s11 for spherical harmonic shape
c             model
c     CON(2)               = flattening for fourier series shape model
c     one of CON(2),CON(3),CON(4) must be non-zero in order to
c     include planetary shape in processing radar observations of
c     planets.  shape coefficients are specified as harmonic input
c     for body -NPLNT.
c     CON(6)  = phase angle of rotation at epoch CON1(1) (this has
c               different meanings in ellipsoidal shape routines in
c               compar link and in the calculation of the effect of this
c               body's gravitational potential harmonics on another
c               body in subroutines sbfn,sbfn1)
c               for special mars rot. code, this angle is referred to
c               the inertial frame defined by CON1(7) and CON1(8)
c     CON(7)  = rotation period in days at epoch CON1(1)
c     CON(8)  = declination of rotation pole at epoch CON1(1) referred
c               to the mean equinox and equator of reference epoch
c     CON(9)  = right ascension of rotation pole at epoch CON1(1)
c               referred to the mean reference equinox and equator
c     CON(10) = i0, inclination at epoch CON1(1) referred to
c               intermediate inertial frame defined by CON1(7) and
c               CON1(8)
c     CON(11) = psi0, precession phase angle at epoch CON1(1)
c               referred to inertial frame defined by CON1(7) and
c               CON1(8)
c     CON(12) = mu, precession constant
c     CON(13) = time derivative of CON(7) in days per days past CON1(1)
c     CON(14) = a: CON(14) through CON(17)
c     CON(15) = b: are coeficients in the seasonal change of
c     CON(16) = c: spin period of mars. see memo by rdr dated 2/14/78
c     CON(17) = d: units are in radians
c     CON(14) = A1 nongrav. force coefficient for comets (see K(83))
c     CON(15) = A2 nongrav. force coefficient for comets
c     CON(18) = atmospheric scale height (km) of this planet (in patctl)
c     CON(19) = atmospheric zenith delay (s) at this planet's surface,
c               assuming a constant scale height
c     CON(20) = deltap = principle of equivalence violation parameter
c               eta*delta(body) where eta = 4*beta-gamma-3,
c               delta = ratio of gravitational to total energy (neg.)
c               see CON1(10).
c     CON(21) = atmospheric scale height (km) of this planet (in plnorb)
c               [Yes, this is conceptually the same as CON(18), but the
c                two are measured from completely different phenomena
c                and may need to be estimated independently.]
c     CON(22) = ref. density (gm/cm**3) of atmosphere: see CON1(12)
c     CON(23) = radius of body in kilometers for transit or occultation
c                     observations
c     CON1(1) = epoch from which rotation is measured
c               for asteroids or satellites 11,...,30 the epoch is
c               also used for initial conditions for elliptic orbits
c               which perturb motion in integration in subroutines
c               SBFN,SBFN1.  Epoch is Julian ephemeris date.
c     CON1(2-4)= components of unit vector pointing from earth to given
c                     body at epoch defining zero longitude on body.
c     CON1(5)  =constant (degrees) to be added to longitude calculated
c               with CON1(1-4) to get IAU system of longitude (applied
c               in subroutine deldop, affects shape model of planet)
c     CON1(6)  =constant (degrees) to be added in subroutine RADAR to
c               IAU longitude on data card for time delay observation
c               away from the subradar point to get longitude that would
c               be calculated with CON1(1-5)
c               (only one of CON1(5) and CON1(6) should be nonzero)
c               (we assume the external world is on the IAU east (not
c               west) longitude system and we either convert PEP to the
c               external world with CON1(5) or convert the external
c               world to PEP with CON1(6), presuming we were not able
c               to get the IAU sytsem with CON1(1-4) alone).
c     CON1(7) = n, an angle defining intermediate inertial frame.
c               Relates inertial frame defined by a mean orbit of date
c               to the reference inertial frame. See Reasenber and King,
c               "The Rotation of Mars" for more precise definition.
c     CON1(8) = j, same comments as for CON1(7)
c     CON1(9)=drgeps,if exponent of density function .gt.drgeps
c          then print exponent with warning: see CON(21-22)
c          CON1(9) could be converted to a TCON if the space is
c          needed for another purpose
c     CON1(10) = ratio of gravitational to total energy (negative,
c                for principle of equivalence calculations): see CON(20)
c     CON1(12)  ref. altitude for CON(22)  (km)
c
c For all probes (NPLNT.gt.30) in the COMPAR link
c CON(16)= transmitter frequency (hz)
c CON(17)= transmitter frequency rate (hz/sec)
c CON(18)= transmitter frequency accelerating (hz/sec**2)
c CON1(3)= initial epoch for CON(17) and CON(18) (PEP jd.fract)
c
c For a probe (NPLNT.gt.30) with logical center the sun
c     CON(i),i=1,11,  adjustable parameters giving low thrust and other
c                     forces which affect motion of probe
c      gas leak force parameters
c     CON(1) = alpha 1
c     CON(2) = alpha 2
c     CON(3) = f 1
c     CON(4) = f 2
c     CON(5) = f 3
c      radiation pressure force parameters
c     CON(6) = g 1
c     CON(7) = g 2
c     CON(8) = g 3
c
c For a probe (NPLNT.gt.30) with logical center a planet (low thrust
c forces computed in plnorb)
c     CON(1)- CON(4) magnitudes of gas leak terms (dynes)
c     CON(5)- CON(8) inverse time constants for gas leaks (/day)
c     CON1(1) rt. asc. of canopus or other ref. star (deg.)
c     CON1(2) dec. of canopus or ref. star (deg.)
c     CON1(3) beta:rotation of gas leak coor. system from sun-star plane
c     CON1(4) epoch for time varying gas leak forces (j.e.d.)
c
c     radiation pressure parameters
c     CON(9)-CON(15) radiation pressure constants
c     CON1(5) beta2:rotation of rad. pres. system from sun-star plane
c     CON1(6) epoch for time varying radiation pressure quan. (j.e.d.)
c
c     atmospheric drag parameters
c     CON1(10)= cd1 drag coefficient of probe
c     CON1(11)=cd2 lift coefficient of probe
c     in COMPAR link we also have:
c     CON(23) transponder delay for time delay and counted doppler
c          calculation ((seconds)
c     CON(1) for instantaneous Doppler (NCODF=1), assumed radius
c          of observed body for time delay (km)
c
c For a probe (NPLNT.gt.30) with logical center the Earth (low-thrust
c forces computed in ERTORB)
c     radiation pressure parameters (see also KK(78))
c     CON(9) = scale factor for direct acceleration
c     CON(10)= scale factor (x total direct accel.) for acceleration
c              in y-axis direction (for GPS satellites with isotropic
c              or nds model)
c     CON(11)= angle between solar panels and y-axis (for GPS NDS model
c              only) (radians)
c
c     CON(18)= ad hoc along-track acceleration (au/day**2)
c
c for LES-8 and LES-9 Earth satellites in subroutine LESTHP
c     CON1(1) = weight of satellite in pounds (if zero, default launch
c               values used in in computing radiation pressure and
c               thrust accelerations)
c     CON(17) = solvable parameter multiplying LES-8/9 thrust
c               (if zero, assumed one)
c     CON(18) = solvable parameter multiplying LES-8/9 radiation
c               pressure (if zero, assumed one)
c
c for LES-8/9 station keeping simulation (setup at end of subroutine
c SBSET, used in subroutines LESTHP and STATON)
c     CON1(2) = satellite orbit plane inclination to ecliptic (deg)
c     CON1(3) = angle between ascending node of satellite orbit
c               plane on ecliptic and sun perigee location (deg)
c     CON1(4) = desired station east longitude (deg)
c     CON1(5) = damping constant (dimensionless)
c     CON1(6) = coasting constant (deg/day)
c     CON1(7) = saturating drift constant (deg/day)
c     CON1(8) = station keeping thrust level (millipounds)
c     CON1(9) = bias to be added to east longitude sought by station
c               keeping system (to correct for discrepency between
c               longitude derived from sun transits and actual longitude
c               of center of LES-8/9 figure eight ground track)
c     CON1(10)= Julian ephemeris date of enabling of sun transit
c               measurements (half day less than jd,fract)
c
c For an exo-planet or other star companion (NPLNT.GT.30, NCENTR.EQ.-4)
c
c     NAME must begin with the four characters of the primary's spot
c     name and continue with any desired additional identification.
c
c     The elliptic orbit elements COND(1-6) are referred to the plane of
c     the sky, with the three coordinate axes pointing north, east, and
c     toward the observer.  For observations of the primary, the orbit
c     position (suitably scaled by the mass factor) must be subtracted.
c     The initial epoch of the orbit is given by JD0 of the primary.
c
c     CON(3) = orbital period in days, explicitly imposed by scaling the
c              "mass" used in the elliptic orbit formulae as the cube of
c              the input semimajor axis.
c     CON(4) = mass of body divided by the sum of masses of the body and
c              its primary.  (Not yet implemented.  The mass factor for
c              the body is assumed to be incorporated in COND(1).)
c
c For any integration which uses SBFN (probe, satellite or even embary
c if forced there by ICT(40)=1 or planet if forced their by NCENTR=-1)
c     CON1(12)= distance below which distance to target body and
c               if KK(81).gt.0 distance to central body is printed out
c               (unit is astronomical units)
c
c For the Earth (NPLNT=3)
c     CON(1) = mean equatorial radius in kilometers
c     CON(2) = degree of flattening
c     CON(3) = Love number h for solid Earth tides
c     CON(4) = Love number l for solid Earth tides
c     CON(5) = lag time (sec) for solid Earth tides
c
c For the Moon (NPLNT=10)
c     CON(1) = mean radius pointing towards mean subradar point
c     CON(3) = Love number h for solid Moon tides
c     CON(4) = Love number l for solid Moon tides
c     CON(5) = lag time (sec) for solid Moon tides
c     CON(6) = ad hoc node precession rate (rad/day about ecliptic)
c     CON(16-18) constants for effect of tidal friction on moon motion
c                sin(2*lag angle)=CON(16)+CON(17)*t+CON(18)*t**2
c     CON(19) = dels = principle of equivalence violation parameter
c               eta*(delta(Earth)+delta(Moon)) where
c               eta=4*beta-gamma-3 and delta (negative) is ratio of
c               gravitational to total energy
c     CON(20) = deld = eta*(delta(Earth)-delta(Moon))
c
c For the Earth rotation (NPLNT=-3)
c   CON(1-21) are shared by two mutually exclusive sets of parameters
c   if JCT(29) = 0:
c     CON(1) = psidot(1) ad hoc small rotation angular rate
c     CON(2) = psidot(2)
c     CON(3) = psidot(3)
c     CON(4) = psi(1)    ad hoc small rotation angle
c     CON(5) = psi(2)
c     CON(6) = psi(3)
c     CON(11-13) polynomial terms in CT-UT (before 1956 Jan 17.0,
c                applied if JCT(34).gt.0) or in A1-UT1 (after 1956 Jan
c                17.0, applied if erotat CON1(1).gt.0).  See JCT(34)
c     CON(11) = r0  constant term for polynomial variation in A1-UT1
c                   (CT-UT before 1956) (sec).  See A1UT1F or CTUTF
c     CON(12) = r1  linear term of this polynomial variation (sec/cy)
c                   epoch for A1-UT1 is CON1(1), epoch for CT-UT is
c                   1956 Jan 17.0.
c     CON(13) = r2  quadratic term of this polynomial variation
c     CON(14-21) seasonal terms in ut2-ut1 (before 1956 Jan 17.0) or
c                a1-ut1 (after 1956 Jan 17.0), applied if ICT(34).gt.0
c     CON(14)= coefficient of cosine annual term
c     CON(15)= time variation of CON(14)
c     CON(16)= coefficient of cosine semi-annual term
c     CON(17)= time variation of CON(16)
c     CON(18)= coefficient of sine annual term
c     CON(19)= time variation of CON(18)
c     CON(20)= coefficient of sine semi-annual term
c     CON(21)= time variation of CON(20)
c   if JCT(29) > 0:
c     CON(1-2)= free core nutation, cos & sin components.  See JCT(28).
c               period (approx. 460 days) is given by CON1(2)
c     CON(3-10) coefficients of adjustable terms in nutation series
c     CON(3)  = coefficient of cos(2F +2 omega) term in deps  (13.7 days
c     CON(4)  = coefficient of sin(2F +2 omega) term in dpsi
c     CON(5)  = coefficient of cos(2L -2D +2 omega) term in deps (183 da
c     CON(6)  = coefficient of sin(2L -2D +2 omega) term in dpsi
c     CON(7)  = coefficient of cos(Lp) term in deps  (365 days)
c     CON(8)  = coefficient of sin(Lp) term in dpsi
c     CON(9)  = coefficient of cos(omega) term in deps  (18.6 yrs)
c     CON(10) = coefficient of sin(omega) term in dpsi
c     CON(11-18) coefficients of analogous out-of-phase terms--
c                 sin terms in deps, and cos  terms in dpsi
c     CON(19) = k/c for fortnightly terms in ut1
c     CON(20) = k/c for monthly terms in ut1
c   For any value of JCT(29):
c     CON(22)= precession constant (if zero nominal value used)
c     CON(23)= obliquity constant (if zero, nominal value used)
c     CON(24)= nutation constant (if zero, nominal value used)
c     CON1(1) = epoch for polynomial variation, CON(11-13).
c               and for seasonal parameters, CON(14-21). for before 1956
c               hardwired in ctutf as 1956 jan 17.0
c               should be jan 17.0 of whatever year is chosen, so that
c               the pre- and post-1956 seasonal terms are in phase
c             = 0 corrections not applied(see a1ut1f & wobblf)
c     CON1(2) = period(days) of free mode nutation, CON(1-2).
c
c For the Moon rotation (NPLNT=-10)
c     CON(1) = mean equatorial radius
c
c     moment of inertia ratios
c     CON(2) = alpha
c     CON(3) = beta
c     CON(4) = gamma
c
c     CON(5) = meqinc, mean inclination of equator to orbit
c     CON(6) = k2, lunar love number
c     CON(7) = k2*t or k2/q, lunar dissipation parameter.  see k(83).
c     CON(8-10) = coefficients of analytic dissipation terms
c     CON(11-16) = initial conditions of lunar core rotation
c     CON(17)= core-mantle rotation coupling constant
c              (solar-mass-AU-squared-per-day)
c     CON(18)= core moment of inertia (solar-mass-AU-squared)
c     CON(19)= core flattening (ignored in computing mantle moments)
c
c For limited asteroids (SMALL=t)
c     initial conditions are specified as for other bodies
c     NPLNT = asteroid number
c     CON(1) = radius (km)
c     CON1(1)= epoch of initial conditions (must be midnight)
c     denptr = index into prmter array for density parameter.
c              logic elsewhere limits denptr to 34-38.
c
c For pulsars (PULSAR=t)
c     NAME - 1st four characters are pulsar spot name
c     JD0   = reference epoch for pulse phase model
c     NTYPE = code for kind of phase model (not used)
c     CON1(2)= approximate pulse period (sec)
c     CON(1)= parallax (arcsec)
c     CON(2)= ra proper motion (arcsec/yr)
c     CON(3)= dec proper motion (arcsec/yr)
c     CON(4)= dispersion (sec-hz**2)
c     CON(5)= pulse phase at reference epoch jd0 (cycles)
c     CON(6)= pulse period - CON1(2) (sec)
c     CON(7)= pulse period drift rate (sec/sec)
c     CON(8)= pulse period acceleration (sec/sec**2)
c     CON(9)= time derivative of dispersion (sec-Hz**2/d)
c     CON(10)= amplitude of the planet/companion's signature (sec)
c     CON(11)= planet orbit eccentricity
c     CON(12)= planet argument of periapse (deg)
c     CON(13)= planet mean anomaly at JD0 (deg)
c     CON(14)= planet orbital period (days)
c     CON(15)= planet orbital inclination (deg)
c     CON(16)= planet mass (Msun)
c     CON(17)= time derivative of planet's signature (sec/d)
c     CON(18)= time derivative of planet's periapse (deg/d)
c
c
c L(1) to L(30) control whether initial conditions COND(1-6) and
c parameters CON(1-24) belonging to given body are adjusted in the
c iterative least squares process.
c           (a) for i=1,6, l(i)=0 or 1 according to whether the
c               corresponding initial condition  COND(i) is not adjusted
c               or is adjusted in the least-squares analysis
c           (b) for i=7,30,l(i)=k means that parameter CON(k) is
c               adjusted,k=1,24.if l(i)=0,no such parameter is adjusted.
c               L(7-30) must satisfy same logical conditions as
c               LPRM(1-100), namely
c               (1) if j.gt.i and L(i)=0, then L(j)=0
c               (2) if j.gt.i and L(i),L(j) not 0, then L(j).gt.L(i)
c
c
c N.B. EPS read in &NMLST2 is not the same as EPS read in &NMLST1.
c
c EPS(1) = start of integration output tape is Julian day number JD1
c          and fraction of day EPS(1) (EPS(1) assumed zero in
c          non satellite-probe case)
c EPS(2) = end of integration output tape is Julian day number JD2
c          and fraction of day EPS(2) (EPS(2) assumed zero in
c          non satellite-probe case)
c EPS(3) = accuracy constant to control numerical integration step size
c          for Nordsieck method (K(88)=0,1 or K(88)=2,3 in starting
c          proceedure)
c EPS(6) = if JD0 negative for checkpoint restart, numerical integration
c          restarted at Julian day number -JD0 and fraction of day
c          EPS(6) or as close after this time as is practicable.
c
c
c  K(1) to K(100) control integration of equations of motion and
c equations for partials of motion for body with planet number NPLNT.
c K(1-30) were formerly used to select partial derivatives for
c integration.  This function is now filled by KI(1-NUMKI).  If NUMKI
c is not specified, then the program logic sets NUMKI to the index of
c the highest-numbered non-zero KI.  Old-format input streams will be
c recognized by the lack of KI's coupled with values .le.1 for all of
c K(1-7).  Old-format input streams will be translated automatically
c (and internally) to the new format.
c

      include 'maxkidat.inc'

c
c  KI(1)
c  KI(1) =-1 equations of motion are integrated and partial derivatives
c              w.r.t. initial conditions are determined from elliptic
c              or mean orbit formulas (depending on value of  K(100))
c  KI(1) = 0 equations of motion are integrated
c  KI(1) = 1 equations of motion and equations for partial derivatives
c              w.r.t. initial conditions are integrated
c
c  KI(2-7)
c  KI(1+j)=-1   partial derivatives with respect to initial orbital
c                 element j (or initial position or velocity coordinate
c                 j if NPLNT.gt.30) are not determined, j=1,...,6
c  KI(1+j)= 0,1 partial derivatives with respect to initial orbital
c                 element j (or initial position or velocity coordinate
c                 j if NPLNT.gt.30) are determined if KI(1) is not zero,
c                 j=1,...,6
c
c  KI(8-NUMKI)
c  KI(j), j=8,...,NUMKI, determine parameters for which equations for
c              partial derivatives are integrated. we must have
c              (1) if j greater than i, then KI(j) greater than KI(i)
c                  unless (a) KI(j)=0 or unless (b) KI(j) is negative,
c                  in which case KI(i) must also be negative with KI(i)
c                  greater than KI(j)  except that KI(k+1),KI(k+2)
c                  are not considered in this logic if KI(k)=100*n+31
c                  nor are KI(k+1),KI(k+1),KI(k+3),KI(k+4) considered
c                  in this logic if KI(k)=100*n+41 or KI(k)=100*n+51
c                  (k=8,...,NUMKI, n=1,2,...)
c              (2) if j greater than i and KI(i)=0, then KI(j)=0
c            with these conventions, we then have for j=8,...,NUMKI
c  KI(j) = 0 no equations integrated
c  KI(j) =-m (m positive) equations for partial derivatives with respect
c                     to parameters  CON(m) are integrated (m=1,24)
c  KI(j) = m equations for partial derivatives w.r.t. PRMTER(m) are
c              integrated, m=1,...,100 (really m=1,...,50 because
c              prmter(1-50) are supposed to affect motion and
c              prmter(51-100) are supposed not to affect motion but
c              only observations)
c  KI(j) = 100*n+i    equations for partial derivatives with respect to
c                     initial condition i (i=1,6) or parameter i-6 (i=
c                     7,30) for perturbing target planet or central
c                     planet n (n=1,10) are integrated
c  KI(j) = 100*n+nn
c          nn=31      equations for partial derivatives w.r.t. zonal
c                     harmonic coefficients between orders KI(j+1) and
c                     KI(j+2) inclusive for target or central planet n
c                     are integrated
c          nn=40-44   equations for partial derivatives w.r.t. tesseral
c                     harmonic cosine coefficients between orders
c                     KI(j+1),KI(j+2) and KI(j+3),KI(j+4) inclusive for
c                     target or central planet n are integrated.
c                     Warning: only 41 working in COMPAR (PARTL1/HPARTL)
c          nn=45-49   equations for partial derivatives w.r.t. resonant
c                     tesseral harmonic cosine coefficients of order
c                     KI(j+1) and of degree KI(j+2) thru KI(j+3).
c                     Not working in COMPAR.
c          nn=50-54   equations for partial derivatives w.r.t. tesseral
c                     harmonic  sine  coefficients between orders
c                     KI(j+1),KI(j+2) and KI(j+3),KI(j+4) inclusive for
c                     target or central planet n are integrated.
c                     Warning: only 51 working in COMPAR (PARTL1/HPARTL)
c          nn=55-59   equations for partial derivatives w.r.t. resonant
c                     tesseral harmonic sine coefficients of order
c                     KI(j+1) and of degree KI(j+2) thru KI(j+3).
c                     Not working in COMPAR.
c
c  K(1-11)          planet numbers of target bodies for probes
c
c
c  K(30)=-1  effect of any input limited asteroids not included
c  K(30)= 0  effect of all input limited asteroids (if any) included in
c              equations of motion but not partial derivatives
c  K(30)= 1  effect of all input limited asteroids (if any) included in
c              equations of motion and partial derivatives (not implem.)
c  K(31-60)
c  K( 30+j)=-1 effect of mass(j) not included in equations of motion or
c                in equations for partial derivatives, j=1,...,30
c  K( 30+j)= 0 effect of mass(j) included in equations of motion but not
c                in equations for partial derivatives, j=1,...,30
c  K( 30+j)= 1 effect of mass(j) included in equations of motion and
c                in equations for partial derivatives, j=1,...,30
c if NPLNT between 1 and 30, K(30+NPLNT) is ignored
c if ncentr is between 1 and 30 (ncentr assumed 3 for NPLNT=10) we have
c  K(30+ncentr)=-1 effect of sun not included in equations of motion nor
c                  in equations for partials of satellite
c  K(30+ncentr)= 0 effect of sun included in equations of motion but
c                  not in equations for partials of satellite
c  K(30+ncentr)= 1 effect of sun included in equations of motion and
c                  in equations for partials of satellite
c
c  K(61-80)
c  K( 30+j)=-1 effect of PRMTER(j) not included in equations of motion
c                or in equations for partial derivatives, j=31,...,50
c  K( 30+j)= 0 effect of PRMTER(j) included in equations of motion but
c                not in equations for partial derivatives, j=31,...,50
c  K( 30+j)= 1 effect of PRMTER(j) included in equations of motion and
c                in equations for partial derivatives, j=31,...,50
c Note: whether  K(61) controls effect of PRMTER(31) or controls effect
c of  CON(24) depends on value of  K(97).
c Also, for sun-centered bodies, K(61)=2 selects a simplified form for
c the general relativity correction to the force which only includes
c the sun and omits terms proportional to the velocity of the sun
c relative to the solar system center of mass.  This simplified
c formulation is the only one available for non-sun-centered bodies.
c K(61)=2 includes relativity in partials.
c Note: K(80) controls whether asteroid ring is included in equations
c of motion in solprb (PRMTER(50) is mass of asteroid ring).
c Also coded in PLNORB.
c
c  K(81) to  K(86) specific to given body
c
c For a planet (NPLNT.gt.0 and NPLNT.le.30)
c      K(82)=-1,0,1   according to whether the interaction of the
c                     planet's second zonal harmonic with the sun
c                     is not included, is included in the motion but
c                     not the partials, or is included in motion and
c                     partials (this last option not yet coded)
c
c For probe (NPLNT.gt.30) with logical central body the sun (SOLPRB)
c (also for planet, asteroid, or comet)
c      K(81) =-1   no low thrust forces
c      K(81) = 0   low thrust forces included in motion, not in partials
c      K(81) = 1   low thrust forces included in motion and partials
c
c      K(83) =-1   no cometary nongravitational forces
c      K(83) = 0   cometary nongrav. forces for motion, not partials
c      K(83) = 1   cometary nongrav. forces for motion and partials
c
c For a probe with logical central body a planet (ERTROB or PLNORB)
c     K(81)- K(84) control for various low thrust forces
c     K(81)= 10*kreflt+ kdirct
c          kdirct.le.0  no direct radiation pressure
c          kdirct= 1 direct radiation pressure in motion
c
c          reflected radiation not yet coded
c          kreflt.le.0 no radiation from central body
c          kreflt=1 reflected radiation from planet in motion
c          kreflt=2 emitted infared radiation from planet in motion
c          kreflt=3 reflected and emitted radiations included in motion
c     K(82).le.0 no atmospheric drag
c     K(82)=+1   effect of drag on motion and partials
c     K(82).gt.1 reserved for table lookup of density in atmden
c                (not coded)
c     K(83).le.0 no gas leak forces included
c     K(83)=1  effect of gas leaks on motion included
c
c     K(84).le.0  no ad hoc along-track accelerations
c     K(84).gt.0  effect of ad hoc along-track acceleration included
c                (ertorb only)
c
c K(86).lt.0 default number of central body harmonics in sbset included
c            in integration of equations for partials
c            default=j2 only (+c22,s22 if moon=central body)
c K(86).ge.0, K(86)=ntess*100+nzone   ntess tesserals and nzone zonals
c            of central body harmonics included in partials integration
c
c K(85) is exactly the same for target bodies as K(86) is for central
c            bodies. no distinction made in K(85) between target bodies,
c            the decision as to target body harmonic order included
c            in partials integration is the same for all of them.
c            default=no harmonics
c
c For Moon (NPLNT=10)
c  K(81)=n, where n is the degree of the highest zonal harmonic
c           included in the Earth's potential.
c  K(82)=m, where m is the maximum degree of the lunar gravity field
c           used. both zonals and tesserals are carried through the
c           mth order.
c  K(84)=-1,0,1 according to whether tidal friction does not affect
c               motion or partials, affects motion but not partials, or
c               affects motion and partials.
c  K(85)=-1,0,1 according to whether a violation of the principle of
c               equivalence does not affect motion or partials, affects
c               motion but not partials, or affects motion and partials
c               (this switch controls the contribution from CON(19)
c               and CON(20), plus the contribution from the metric
c               parameters beta and gamma, PRMTER(43-44) )
c
c  For Moon rotation (NPLNT=-10)
c  K(81)=-1,0,1 according to whether the sun does not affect the
c         rotation, affects the rotation e.o.m. only, or affects
c         the rotation eom and the variational equations.
c  K(82)=-1,0 according to whether or not the Earth-figure torque
c         is included in the equations of motion.
c  K(83)=-1,0,2 implying, respectively, no lunar elasticity &
c         dissipation (rigid-body), elas. & constant-time-lag diss.,
c         elas. & constant-q dissipation model.
c         1: same as 0, but calculate lunar mean motion from masses
c  K(84)=-1,0,1 same as K(81), but for planets Venus through Saturn
c  K(86)  controls lunar harmonics for partials, same as in probes
c         (but use K(86) of Moon if higher)
c
c  note-- if orbit and rotation are simultaneously integrated
c         (both JD0s non-zero) then the integration control info.
c         comes from NPLNT=10 stream. this includes JD0, JD1, JD2,
c         EPS(3), K(87-92), and also K(31-40).
c
c  K(87)   = positive, step size for Adams-Moulton or second sum
c            numerical integration is  K(87) days
c  K(87)  .le.zero, interval size for Adams-Moulton or second sum
c            numerical integration is 2** K(87) days
c
c  K(88) = Indicator of the type of integration, if any, to be done.
c          Also saved in NAMTIM common as INTTYP for each body.
c          The default value is 3.
c  K(88) =-9 no integration: orbit is approximated by a quasi-elliptic
c            model specified in COND(1-6) and TCON(1-10). JD0 had
c            better be 0!
c  K(88) =-8 no integration: orbit is approximated by a purely elliptic
c            model with elements COND(1-6) and epoch CON1(1).  JD0 had
c            better be 0!
c  K(88) = 0 Nordsieck variable step size method used for numerical
c            integration with variable tabular interval output
c  K(88) = 1 Nordsieck variable step size method used for numerical
c            integration with constant tabular interval output
c  K(88) = 2 Adams-Moulton method used for numerical integration
c  K(88) = 3 Royal Road (second sum) method used for numerical integ.
c Both Adams-Moulton and Royal Road are constant step integration
c methods with constant tabular interval output, and they both call
c Nordsieck method in their starting proceedures (Nordsieck is self
c starting).
c
c  K(89)   = number of predictor and corrector terms in Adams-Moulton
c            or Royal Road (second sum) numerical integration
c
c  K(90) = number of equations controlling integration step size
c  K(91) = interval in starting proceedure for integration is 2** K(91)
c  K(92) = minimum interval of integration is 2** K(92)
c
c  K(93-96)  (old format only) planet numbers for probe targets
c
c  K(97) = 0 relativity motion factor for force controlled by  K(61)
c            is PRMTER(31)
c  K(97) = 1 relativity motion factor for force controlled by  K(61)
c            is  CON(24) inserted into PRMTER(31) in sub.PLANET for fn
c
c K(98).le.0 no printout of partial derivatives
c K(98).gt.0 printout of partial derivatives
c  K(98) = also controls how often data is printed out, namely
c            every iabs(K(98)) tabular points if K(98).ne.0
c            and every tabular point of integration if K(98)=0
c
c  K(99) =-1 printout, no tape for given body integration results
c  K(99) = 0 printout and tape for given body integration results
c  K(99) = 1 tape, no printout for given body integration results
c
c  K(100)=-1 numerical integration of usual equations of motion
c  K(100)= 0 numerical integration of equations for difference between
c                true orbit and initial osculating elliptic orbit
c  K(100)= 1 numerical integration of equations for difference between
c                true orbit and mean orbit
c (For planets, K(100) assumed zero for fn integration)
c (For Moon rotation, K(100)=-1,0,1 signifies Euler angles,
c   osculating Cassini angles, or mean Cassini angles)
c
c KK(1-100)  additional control integers
c
c KK(1)  = logical central body number which controls what subroutines
c          are called in SBFN and SBFN1 (default is NCENTR)
c
c KK(2)  = 0 never any output on KOUT.ne.0 from SBEXP extended print
c            called from SBOUT during SBFN probe or planet integration
c      .lt.0 SBEXP called every -KK(2)th apsis
c      .gt.0 SBEXP called every KK(2)th integration step
c
c KK(3),KK(4),KK(5) define the form of the extended print block for
c                   SBFN integration of a probe or planet
c KK(6) defines the form of an alternate extended print block
c
c KK(3) is a packed bits number
c          1 cartesian coordinates are output
c          2 all partials are output
c          4 apsidal time is printed if near
c          8 elliptic elements referred to the mean equinox and equator
c            (of the Earth) of the reference epoch are output
c not programmed as yet
c         16 print elliptic elements  central body of date
c         32 print elliptic elements  plane of the sky
c
c KK(4) is a packed bits number
c          1 print  Earth - s/c vector
c          2 print  Earth - central body vector
c          4 print   sun  - s/c vecyor
c          8 print   sun  - central body vector
c         16 print target - s/c vector
c
c KK(5) is a packed bits number
c          1 print  latitude & longitude of s/c w.r.t. central body
c          2 print projected motion of s/c on 'stationary' central body
c            print projected motion of s/c on   rotating   central body
c
c KK(6) is a packed bits number
c          1 print vector of s/c w.r.t. Arecibo at highest signal
c          2 print vector of s/c w.r.t. Earth at close approach
c          4 print vector of s/c w.r.t. 1st target at close approach
c          8 ditto 2nd target
c            ...
c       8192 ditto 12th target
c
c KK(7) is a packed bits number for Kout extra print
c         1 print accelerations at every tabular point (coords+partials)
c         2 print relativistic accelerations and (1st only) partials
c         4 print coordinates and gravity gradients w.r.t. them for
c           three bodies: integrated, 1st target, center
c           when integrating embary, 3rd body is moon
c         8 print partials of coordinates of same three bodies w.r.t.
c           1st parameter
c        16 print newtonian indirect partials contributions due to
c           positions and velocities of same three bodies, and then
c           relativistic likewise (note newtonian velocity indirect
c           contribution is always zero)
c        32 change interval of printout to 100 days
c        64 reserved for nbody extra print
c       128 suppress printout of moon libration angles
c
c KK(10) =   scale factor for f-format coordinate output in sbout
c KK(10) = n print (x,y,etc.)*10**n in sbout and sbexp
c          (partials are not scaled)
c
c KK(11) = moon core integration output dataset (for moon rotation)
c KK(12) = moon core integration forward/backward buffer dataset
c KK(13) = moon libration integration forward/backward buffer dataset
c          (supersedes ibuf if nonzero)
c
c KK(50) =   number of integration steps per printing of
c            the dadx# matrices.  if .lt. zero, print for
c            iteration 1 only.  control is bl pntflg, set
c            in sbfn at "once per iteration for partials"
c
c        extra control for Nordsieck variable stepsize integration
c        see TCON(1-4) below
c KK(51) =  0 stepsize control not nused (default)
c          -1 make transition to stepsize control
c           1 lock onto stepsize control
c
c KK(70) = copy of JCT(13), not independently settable
c
c KK(71) =   UTC hour for LES-8/9 station keeping firing time
c KK(72) =   UTC minute for LES-8/9 station keeping firing time
c KK(73) =   UTC second for LES-8/9 station keeping firing time
c
c KK(78) = 0 no data set with radiation pressure model on it
c KK(78).gt.0 data set number for cards from which mvmsrf interpolates
c            solar radiation force
c KK(78).gt.0 in ERTORB, use NDS radiation pressure model
c
c KK(81) = 0 no printout of distance to central body
c KK(81) = 1 printout of distance to central body in kilometers
c            if this distance is less than CON1(12) au (central
c            body not sun)
c
c KK(84) = 0  read no center/target body tape during integration
c KK(84) = 1  read tape for center/target body position, velocity only
c KK(84) = 2  read tape for partials as well
c             (option 1 above not implemented: treated as option 2)
c             default is 2 for the Moon, 0 otherwise
c
c KK(85)     hsitb=distance from target body at which the
c            effect of the harmonics(not j2)on the e. of m. and partials
c            is disabled.  it is nominally set so that the
c            maximum possible harmonic-acceleration is less than
c            the solar-acceleration by hsftr=1.0e15
c            hsitb=nominal+abs(KK(85))/100.    (a.u.)
c
c KK(88) =   intm in EVAL
c          0 correct 2 times - old form - default
c         -1 Royal Road mode - 1 correction
c          1 like -1 plus pseudo correction with dadx matrix
c                  not tested and evidently not working: KK(88)=1
c
c KK(91) = 0 no special LES-8/9 forces in equations of motion
c KK(91) =-1 LES-8/9 radiation pressure, no thrusts, standard attitude
c KK(91).gt.0 LES-8/9 radiation pressure with thrust and attitude
c             history read from data set KK(91), which is written in
c             subroutine lesin from reading cards following this &nmlst2
c
c KK(92).gt.0 print output during thrust in subroutine lesthp
c
c KK(93) = 0 no change of numerical integration method in checkpoint
c            restart
c KK(93) = 1 such a change in subroutine sbset with input values of
c            EPS(3),K(87-92) overriding values on second record of
c            restart tape
c
c KK(94) = 0 no LES-8/9 station keeping thrusts
c KK(94) = 1 LES-8/9 station keeping emulation is turned on
c            in numerical integration of equations of motion
c KK(94) = 2 emulation turned on plus write out thrust events
c
c KK(95) =   starting interval is 2**KK(95) (for KK(95) negative)
c            for variable step size integration starting proceedure
c            at thrust initiation and termination  for LES-8/9
c            station keeping emulation. K(91) is still starting
c            interval for starting proceedure at initial epoch.
c
c KK(96) =   number of days from epoch in first time direction
c            for reintegration of equations of motion after orbit fit
c            convergence
c
c KK(97) = 0 no partials integrated in reintegration of equations of
c            motion after orbit fit convergence
c KK(97) = 1 partials are integrated in reintegration of equations of
c            motion after orbit fit convergence
c
c KK(98) =   number of days from epoch in second time direction
c            for reintegration of equations motion after orbit fit
c            convergence
c
c KK(99) =   output tape number for reintegration of equations of motion
c            after orbit fit convergence if JCT(79)=1
c
c KK(100)= 0 partial derivatives are integrated each iteration if
c            K(1-30) so indicate
c KK(100).gt.0 partial derivatives are not integrated if iterat.gt.
c            KK(100) (iterat is least squares interation conter)
c            (integration will do what K(1-30) say to do when iterat=1
c            no matter what the value of KK(100))
c
c
c        the following four TCONs are used for Nordsieck integration
c        in conjunction with stepsize control option - see KK(51)
c    TCON(1)= scaling variable, approx eq semi-major axis of probe in au
c    TCON(2)= period of probe in days
c    TCON(3)= approximate number of integration steps per period
c    TCON(4)= exponential index of (r/a)
c
c     for non-integrated bodies TCON(1-10) are parameters of
c     a precessing quasi-elliptic orbit (all deg and day)
c TCON(1)= i0, inclination of orbit to laplacian plane
c TCON(2)= ldot, mean motion in longitude
c TCON(3)= pdot, rate of longitude of periapse
c TCON(4)= kdot, rate of argument of node on laplacian plane
c TCON(5)= ndot, rate of node of laplacian plane on earth equator
c TCON(6)= jdot, rate of inclination of laplacian plane
c TCON(7)= ldotdot, 1/2 rate of change of mean motion
c TCON(8)= theta, amplitude of term sin(k-phi)
c TCON(9)= phi, phase of term sin(k-phi)
c TCON(10)=k0, initial argument of node on laplacian plane
c
c In the following gravitational potential is for NPLNT positive.
c For NPLNT negative, this is to be read planet shape.
c
c NZONE =order of zonal harmonic input for gravitational potential
c            of given body
c Ji    =zonal harmonic coefficient for given body, i=2,3,4,...
c LJi   =0,1 according to whether corresponding zonal harmonic is not
c            or is adjusted in least squares analysis
c
c NTESS =order of tesseral harmonic input for gravitiational potential
c            of given body
c Ci(k) =tesseral cosine harmoic coefficient for given body,
c            i=2,3,4,... and k=1,...,i
c Si(k) =tesseral sine harmononic coefficient for given body,
c            i=2,3,4,... and k=1,...,i
c LCi(k)=0,1 according to whether corresponding tesseral cosine
c            harmonic is not or is adjusted in least squares analysis
c LSi(k)=0,1 according to whether corresponding tesseral sine
c            harmonic is not or is adjusted in least squares analysis
c
c
c There are 4 planet shape models programmed in the radar link of PEP.
c CON(2,3 or 4) for NPLNT.gt.0 turns on shape for radar observations of
c a planet.  The model is given as harmonic input for body -NPLNT in
c one of the following manners:
c        1)shperical harmonic expansion, if nshp=0 (default)
c        2)fourier series model, if nshp=1
c        3)altitude grid (local model), if nshp=2
c        4)external model if nshp=3
c
c In the spherical harmonic expansion over the whole planet surface,
c CON(2-4) for NPLNT are J1,C11,S11 with Jn,Cmj,Smj (n.le.NZONE, m.le.
c NTESS, j.le.m) for -NPLNT being the remaining harmonic coefficients.
c normalized coefficients are used for shape, whereas unnormalized are
c used for the gravitational potential. for the spherical harmonic model
c planet radius = CON(1) + spherical harmonic series
c (units are seconds in solution for harmonic coef,  km for CON(1)  )
c
c ***** Fourier series model
c
c If NSHP=1 for -NPLNT the Fourier series model is used with
c TLAT(1) and TLAT(2) equal to the lower and upper latitude
c values for the planet NPLNT being analyzed.
c Let phi=east long on planet,
c     theta=(latitude-tlat(1))/(tlat(2)-tlat(1))*twopi
c planet radius = function of CON(1) mean radius and CON(2) degree of
c flattening + fourier series coefficients for -NPLNT
c the shape model, as it is currently, is a 2-dimensional fourier series
c with 20 longitudinal terms (smallest period = 360/20=18 deg) and
c 1 latitudinal term.  it is:  h(theta,phi)= summation for m=1,20 of
c        (am*cos(m*phi)*cos(theta)  +  bm*cos(m*phi)*sin(theta)
c      +  cm*sin(m*phi)*cos(theta)  +  dm*sin(m*phi)*sin(theta)
c      +  a'm*cos(m*phi)  +  c'm*sin(m*phi) )
c      +  a''*cos(theta)  +  b''*sin(theta)
c notice that the a' and c' coefficients are just longitude dependent,
c and the a'' and b'' are the only two latitude-only terms.
c all 122 fourier series and radius & flattening partial derivatives
c should be calculated when computing the observed minus theoretical
c values of observations.  suggested parameters to be solved for
c in addition to radius and perhaps flattening are:
c Mercury - only low order longitudinal terms, maybe a'1 thru a'5 and
c           c'1 thru c'5.
c Venus   - how about a'1 thru 20 and c'1 thru 20?
c Mars    - probably all 122 fourier shape parameters
c
c Summary: to use the fourier shape section of coding you set NSHP=1,
c TLAT(1&2) equal to lower and upper latitude values for the rectangle
c (almost) on the planet being analyzed in &NMLST2 for -NPLNT.
c Default values presently set are:
c     Mercury   TLAT(1)=-15., TLAT(2)=15.,
c     Venus     TLAT(1)=-15., TLAT(2)=15.,
c     Mars      TLAT(1)=-30., TLAT(2)=30.,
c In the corresponding &NMLST2 with positive NPLNT, CON(2)=1.0E-20_10
c sets npshp=1 so HARSHP is called, and L(7)=1, as usual, adjusts the
c radius. In addition L(8)=2 adjust flattening CON(2), which could
c have a real value instead of the dummy 1.E-20_10 used as a switch.
c Units of radius are km, flattening is dimensionless, and
c units of Fourier seriers coefficients for radius are light-seconds.
c
c      The Fourier series coefficients and L-vectors may be input
c directly as the following &NMLST2 variables:
c    coefficients                                       l-vectors
c   fa(1)  to fa(20)         for a coefficients            lfa
c   fb(1)  to fb(20)         for b coefficients            lfb
c   fc(1)  to fc(20)         for c coefficients            lfc
c   fd(1)  to fd(20)         for d coefficients            lfd
c   fap(1) to fap(20)        for a' coefficients           lfap
c   fcp(1) to fcp(20)        for c' coefficients           lfcp
c   fapp                     for a'' coefficient           lfapp
c   fbpp                     for b'' coefficient           lfbpp
c      to adjust coefficients, set appropriate l-vector to 1
c ***** example *****
c   mercury fourier series shape
c   solve for a'1 thru a'5 and c'1 thru c'5
c      required input is:
c       &nmlst2
c       NPLNT= -1,name= 'shpmrcry',
c       nshp= 1,
c       lfap= 1,1,1,1,1,     (or lfap= 5*1,)
c       lfcp= 1,1,1,1,1,     (or lfcp= 5*1,)
c       &end
c
c ***** altitude grid-local shape model (s.brody 5/76)
c
c if nshp=2 for -NPLNT the shape model is a grid of
c altitudes (in kilometers) relative to an adjusted spheroid
c for a appropriate region of the planet NPLNT. the region
c of applicability for the model is input as the real*4
c nmlst2 variables tlat(1),tlat(2),tlon(1),tlon(2) representing
c minimum and maximum values for latitude and longitude, respectively
c (in degrees). tlatin and tlonin give the increment in
c degrees for lat and long (ie the grid spacing). default
c values for the lat. range are as given above for the fourier
c series model. default for long is 0 to 360 degrees (east)
c grid spacing defaults are 10 degrees. (latdim,londim) the dimensions
c of the grid and the total number of grid points (ngdpts) are computed
c from tlat,tlon,tlatin,tlonin.
c to input a priori estimates for grid points (in km),
c determine the corresponding index into the (1-dimensional)'grid'
c array where the 2-dimensional grid points are stored by
c row (ie., in latitude strips). grid(1) corresponds to the
c grid point of min. long and max lat (ie., the upper left
c point). lgrid is the l-vector for the grid points with
c     lgrid(i)=0, do not adjust grid(i)
c             =1,        adjust grid(i)
c as in the fourier shape model, the flattening coefficient CON(2)
c for the body |NPLNT| should be set to some small non-zero value,
c such as 1E-20_10 to enable the shape computations in the COMPAR link.
c note that the flattening may also be treated as an adjustable
c parameter as well as a switch.
c
c ***** additional shape models may be implemented using nshp gt 3 *****
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c local variables
      integer*2 nplch
      character*8 blank/'        '/
      real*10 dum,fct,fract,fract0,fract9,wct(3)
      integer i,ict66,ictl,iftkm,incnd,incr,indnd,
     . intx,inunit,jx,jcnd,jd09,jdtype,kcnd,
     . latdim,lcnd,londim,m,mcentr,mcnd,mcndi,mean,ngdpts,ntype
c external functions
      real*10 LEGSCL
      integer JULDAY

c new values of incnd, inunit, icnd for old-style icnd:
c input       icnd= ...  -8 -7 -6 -5 -4 -3 -2 -1  0  1  2  3  4  5
      integer*2 ivin(14)/-1,-1,-1,-1, 1, 0,-1,-1, 0, 1, 2, 2,-1,-1/,
     .          ivun(14)/ 1, 1, 0, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 1/,
     .          ivic(14)/ 1, 0, 1, 0, 1, 0,-1,-1, 0, 1, 2, 2, 2, 2/


c
c initialize wrkcom labeled common
      call ZFILL(Cond, zwrkcm)

      Eps(3) = 1.0E-16
      Name   = blank
      Icnd   = 99
      incnd  = 99
      inunit = 0
      jcnd   = 0
      kcnd   = -1
      lcnd   = 0
      jtype  = 99
      iftkm  = 0
      mean   = 0
      mcentr = -2
      mcnd   = -1
      do i = 31, 39
         K(i) = 1
      end do
      do i = 40, 86
         K(i) = -1
      end do
      K(61)  = 1
      K(88)  = 3
      K(89)  = 11
      K(90)  = 6
      K(91)  = -12
      K(92)  = -14
      fract0 = 0._10
      Itape  = -1
      Pulsar = .false.
      Small  = .false.
      Nzone  = -1
      Ntess  = -1
      sclzon = 0
      scltes = 1

      imn0   = 0
      idy0   = 0
      iyr0   = 0
      imn1   = 0
      idy1   = 0
      iyr1   = 0
      imn2   = 0
      idy2   = 0
      iyr2   = 0
      jdtype = 0

c
c spool &nmlst2 from in to in0 with a-format printout
      call PEPTIC(In, Iout, in0, 9,
     .            'INDIVIDUAL BODY PARAMETERS &NMLST2  ', nstop, 0)
c
c read &nmlst2 namelist from in0
      read(in0, NMLST2)
      rewind in0
c
c initialize for specific body and read data again
      if(jtype.eq.99) then
         jtype = -1
         if(kcnd.gt.0) jtype = 1
      else if(jtype.ge.6 .and. Jct(13).eq.1) then
         jdtype = 2451545
      endif
      if(Small) goto 100
      if(Pulsar) goto 300
      if(Nplnt.eq.0) goto 200
      call BODNTL
      read(in0, NMLST2)
      rewind in0
      Kk(70) = Jct(13)
c
c detect explicit new style partials controls
      if(Numki.le.0) then

c detect implicit new style (targets without partials)
         do i = 1, 7
            if(K(i).gt.1) Numki = 8
         end do

c detect implicit new style and count partial controls
         do i = 1, maxki
            if(Ki(i).ne.0) Numki = i
         end do
      endif
      if(Numki.gt.0 .and. Numki.lt.8) Numki = 8

c convert old format controls to new (if necessary)
      call KP2KI(K, Numki, Ki)
c
c convert from calendar date to Julian day number
      if(iyr0.gt.0) Jd0 = JULDAY(imn0, idy0, iyr0)
      if(iyr1.gt.0) Jd1 = JULDAY(imn1, idy1, iyr1)
      if(iyr2.gt.0) Jd2 = JULDAY(imn2, idy2, iyr2)
c
c are norad 2 card mean elements to be read
      if(mean.gt.0) then
         call NORAD(Jd0, Ihr, Imin, Sec, Cond, in0)
         jtype  = 0
         inunit = 0
         incnd  = -1
         if(Icnd.eq.99) Icnd = -1
      endif
c
c fix starting and ending times of integration
      if(Jd0.ne.0 .and. Jd0.ne.-1 .and. Jd0.gt.-3000000) then
         if(iabs(Jd1).lt.100000) Jd1 = Jd1 + iabs(Jd0)
         if(iabs(Jd2).lt.100000) Jd2 = Jd2 + iabs(Jd0)
      endif
      if(Jd1.gt.0 .and. Jd2.gt.0 .and. jdpad.gt.0) then
         if(Jd2.ge.Jd1) then
            Jd2=Jd2+jdpad
            Jd1=Jd1-jdpad
         else
            Jd1=Jd1+jdpad
            Jd2=Jd2-jdpad
         endif
      endif
c
c setup nonzero initial fraction of day
      if(fract0.ne.0._10) then
         fract = fract0*2._10**30
      else if(Int1.ne.0) then
         if(Int2.gt.0) Int2 = -Int2
         goto 50
      else
         fract = 60*(Imin + 60*Ihr)
         fract = (fract + Sec)/8.64E4_10
         if(fract.eq.0._10) goto 100
         fract = fract*2._10**30
      endif
      Int1  = fract
      fract = fract - Int1
      if(fract.gt.0.5_10) Int1 = Int1 + 1
      Int2 = -30
      do i = 1, 29
         intx = Int1/2
         if(Int1.ne.intx*2) goto 50
         Int1 = intx
         Int2 = Int2 + 1
      end do
   50 fract = Int1
      fract = fract*2._10**Int2
      Sec   = fract*8.64E4_10
      Ihr   = Sec/3600._10
      Sec   = Sec - Ihr*3600
      Imin  = Sec/60._10
      Sec   = Sec - Imin*60
c
c change units and type of initial conditions to the
c allowed icnd=-1,0,1,2 and au & au/day
  100 nplch = Nplnt
      if(Small) nplch = 999
      if(Icnd.lt.-8 .or. Icnd.gt.5) Icnd = 0
      if(nplch.gt.0) then

c convert old-style icnd
         indnd = Icnd + 9
         if(incnd.lt.-1 .or. incnd.gt.2) then
            incnd  = ivin(indnd)
            if(inunit.eq.0) inunit = ivun(indnd)
         endif
         Icnd = ivic(indnd)
         if(jcnd.gt.0 .and. incnd.eq.-1) incnd = 0
         if(kcnd.ge.3 .and. incnd.eq.-1) incnd = 2
         if(jcnd.gt.1 .or. kcnd.eq.4 .or. lcnd.eq.2) inunit = 1
         if(kcnd.eq.0) jtype = -1
         if(mcentr.ge.-1 .and. mcentr.ne.Ncentr .and. Jd0.gt.0)Icnd = -1
         if(lcnd.ne.0 .and. incnd.ne.2) then
            write(Iout, 120) incnd, lcnd
  120       format(' *** INCND=', i3, ' LCND=', i3,
     .             ' INCONSISTENT, BODRED ERROR')
            if(Mout.gt.0) write(Mout, 120) incnd, lcnd
            nstop = nstop + 1
         endif
         if(iftkm.gt.0) inunit = 2
         if(inunit.ne.0) then
            fct = 1._10
            if(inunit.eq.1) fct = Aultsc*Ltvel
            if(inunit.eq.2) fct = Aultsc*Ltvel/.3048E-3_10
            incr = 6
            if(incnd.eq.2) incr = 3
            if(incnd.lt.0) incr = 1
            do i = 1, 6, incr
               Cond(i) = Cond(i)/fct
            end do
            if(incr.ne.6) then
               fct = 8.64E4_10
               if(inunit.lt.0) fct = 1._10/40._10
               do i = 4, 6, incr
                  Cond(i) = Cond(i)*fct
               end do
            endif
         endif
      endif
      jd09   = Jd0
      fract9 = fract

c coordinate system epoch is different than initial condition epoch
      if(jtype.eq.-1) then
         jtype = 0
         jd09 = 2433282
         fract9 = 0.923_10
      else if(jtype.eq.3) then
         jtype = 2
         jd09 = 2433282
         fract9 = 0.923_10
      else if(jtype.eq.5) then
         jd09 = jdtype
         fract9 = 0._10
      else if(jtype.eq.6) then
         jtype = 0
         jd09 = jdtype
         fract9 = 0.5_10
      else if(jtype.eq.7) then
         jtype = 2
         jd09 = jdtype
         fract9 = 0.5_10
      endif
      if(jd09.le.0) then
         write(Iout, 130) jd09
  130    format(' **** INVALID DATE FOR TRANSFORMATION',I7,
     .          ', ERROR IN BODRED')
         if(Mout.gt.0) write(Mout, 130) jd09
         nstop = nstop + 1
      endif
      call CHNCNE(nplch, Ncentr, Cond, jd09, fract9, lcnd, jtype,
     .            incnd, Icnd, Jct(13)+1)
c
c special logic for lunar core initial conditions
      if(nplch.eq.-10 .and. incnd.eq.3) then
         if(Icnd.eq.0) then
            call MONROT(0,Cond(1),Cond(2),Cond(3),dum,dum,dum)
            call MNCROT(Con(11),Con(12),Con(13),dum,dum,dum)

c get rotation matrix relating core-fixed to body-fixed
c multiply a body-fixed vector by mcmrot to get a core-fixed vector
            call PRODCT(Mcrtlb,Mrotlb,Mcmrot,3,-3,3)
            call PRODCT(Mcmrot,Con(14),wct,3,3,1)

            Con(14) = (wct(1)*Sphic + wct(2)*Cphic)/Sthetac
            Con(15) = wct(1)*Cphic - wct(2)*Sphic
            Con(16) = wct(3) - Con(14)*Cthetac
            incnd   = 0
         else if(Icnd.eq.3) then
         else
            write(Iout,135) incnd,Icnd
  135       format(' **** INVALID INCND,ICND',2I3,', ERROR IN BODRED')
            if(Mout.gt.0) write(Mout,135) incnd,Icnd
            nstop = nstop + 1
         endif
      else if(nplch.eq.-10 .and. incnd.eq.0 .and. Icnd.eq.3) then
         call MONROT(0,Cond(1),Cond(2),Cond(3),dum,dum,dum)
         call MNCROT(Con(11),Con(12),Con(13),dum,dum,dum)
         call PRODCT(Mcrtlb,Mrotlb,Mcmrot,3,-3,3)
c determine angular velocities in core-fixed frame
         wct(1) = Con(15)*Cphic + Con(14)*Sphic*Sthetac
         wct(2) = -Con(15)*Sphic + Con(14)*Cphic*Sthetac
         wct(3) = Con(14)*Cthetac + Con(16)
c convert to mantle frame
         call PRODCT(Mcmrot,wct,Con(14),-3,3,1)
         incnd=3
      endif

      if(nplch.gt.0) then
c
c calculate i0,psi0 from alpha0, delta0 or vice-versa
         if(.not.Small) then
            ict66 = Ict(66)
            if((mod(ict66,2).eq.1) .and. (Nplnt.eq.4) ) then
               ictl = 0
               if(Con(8).eq.0._10) ictl = 1
               call ROCHNG(Con(10), Con(11), Con(8), Con(9), Con1(7),
     .                     Con1(8), ictl)
            endif
c
c check for sb integration of embary
            if(Nplnt.eq.3 .and. Ncentr.eq.-1) Ict(40) = 1
c
c change central body
            if(mcentr.ge.-1 .and. mcentr.ne.Ncentr .and. Jd0.gt.0) then
               mcndi = -1

c set for conversion after changing center
               indnd = mcnd + 9
               if(indnd.ge.1 .and. indnd.le.14) mcndi = ivic(indnd)
               call CHNCNT(Jd0, fract, Nplnt, Ncentr, mcentr,
     .                     Cond, Icnd, mcndi)
            endif
         else
            if(Icnd.ne.0) then
               write(Iout, 140) Icnd
  140          format(' * * * INVALID ICND=',i3,' FOR LIMITED ASTEROID')
               if(Mout.gt.0) write(Mout, 140) Icnd
               nstop = nstop + 1
            endif
            goto 300
         endif
      endif
c
c change scale of input gravitational potential harmonics
c to output values of input link
      if(Nplnt.gt.0) then
         if(Nzone.gt.1 .and. sclzon.ge.1) then
            do jx = 2, Nzone
               i = jx - 1
               zh(i) = -zh(i)*LEGSCL(jx,0)
            end do
         endif
         if(Ntess.gt.1 .and. scltes.le.0) then
            m = 0
            do jx = 2, Ntess
               do i = 1, jx
                  m = m + 1
                  ch(m) = ch(m)/LEGSCL(jx,i)
                  sh(m) = sh(m)/LEGSCL(jx,i)
               end do
            end do
         endif
      endif
c
c read LES-8/9 thrust and attitude history
  200 if(Kk(91).gt.0) call LESIN(nstop, Name, in0, Kk(91))
c
c additional setup for planet shape
      if(Nplnt.lt.0) then
         if(Nshp.eq.2) then
c for altitude grid model of shape(nshp=2),
c determine grid dimensions and number of grid points
            call SHPDIM(tlat, tlatin, ' LAT', latdim)
            call SHPDIM(tlon, tlonin, ' LON', londim)
            if((tlon(2)-tlon(1)).eq.360.) londim = londim - 1
            ngdpts = latdim*londim
         endif
      endif

  300 continue
      return
      end
