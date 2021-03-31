      subroutine PLNHED
 
      implicit none
c
c           ash/forni  november 1966    subroutine plnhed
c     printout of data page for planet numerical integration
c     write first two records of output ephemeris tape
c     updated feb. 1980   kcl: ptcon added
c
c array dimensions
      include 'globdefs.inc'
c           common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'loadid.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'rstart.inc'
      include 'sbstuf.inc'
      include 'stint.inc'
      include 'tapdtplp.inc'
      include 'zeroes.inc'
c
c see subroutine bodred of input link for explanation of
c                  cond(1-6)       called cond(1-6)   in bodred
c                  con(1-24)       called con(1-24)   in bodred
c                  con1(1-12)      called con1(1-12)  in bodred
c                  tcon(1-30)      called tcon(1-30)  in bodred
c                  epsp(1-6)       called eps(1-6)    in bodred
c                  name            called name        in bodred
c                  jdp1,jdp0,jdp2  called jd1,jd0,jd2 in bodred
c                  jplnt           called itape       in bodred
c                  kp(1-100)       called k(1-100)    in bodred
c                  kkp(1-100)      called kk(1-100)   in bodred
c                  ki(1-n)         called ki(1-n)     in bodred
c                  intp            called int         in bodred
c                  intp1           called int1        in bodred
c                  intp2           called int2        in bodred
c                  nplnt           called nplnt       in bodred
c                  ncentr          called ncentr      in bodred
c                  ihrp,iminp,secp called ihr,imin,secin bodred
c                  icnd            called icnd        in bodred
c
c iparp= number of partial derivatives on output tape + 1
c           if iparp=1, no partials on tape
c           iparp calculated in subroutine setup from values of kp(1-30)
c           the meaning of the partials j=2,iparp is given by kp(1-30).
c
c klan = 0  for earth-moon barycenter integration
c klan = 1,...u_mxpl  for planet integration,indicating that record klan+4
c                   of iplcon data set and arrays with second dimension
c                   klan from /namtim/ and /empcnd/ were used in putting
c                   data into /petuna/
c

c local variables
      real*10 frp,frq
      integer*4 i,ifiltr,itf,itf1,itf2,ll,nstop,ntyp1,ntyp2
      integer*2 nparitr,nparit1,nparpl(20)
      character*44 mesg/'  SETUP FOR ******** NUMERICAL INTEGRATION  '/
      character*8 encowm(3)/' COWELL ',' ENCKE  ','MEAN ORB'/
      character*20 numeth(4)/'NORDSIECK (VAR.INT.)',
     .                       'NORDSIECK (CON.INT.)',
     .                       'ADAMS MOULTON       ',
     .                       'ROYAL ROAD (2ND SUM)'/
 
      ifiltr = 0
      if(Ict(42).gt.0 .and. Jct(56).gt.0) ifiltr = 1
c
c write first two records of planet tape
      if(Kp(99).ge.0 .and. Jdstrt.le.0) then
         nparitr=0
         do i=0,Numtar
            if(Jplntg(i).gt.0) then
               nparitr=nparitr+1
               nparpl(nparitr)=Nplpt(i)
            endif
         end do
         nparit1=nparitr
         if(nparitr.eq.0) then
            nparit1=1
            nparpl(1)=0
         endif
         write(Jplnt) Name,Heding,Date,Lnklvl
         write(Jplnt) Nplnt,Ncentr,Iparp,Intp,Jdp1,Jdp0,Jdp2,
     .    u_nmprm,u_nmbod,
     .    Cond,Con,Con1,prmter,Epsp,Kp,Npage,Iterat,Icnd,Intp1,Intp2,
     .    Ihrp,Iminp,Secp,Kkp,ifiltr,Lnklvl,Numki,(Ki(i),i=1,Numki),
     .    nparitr,nparit1,(nparpl(i),i=1,nparit1),
     .    Izero,Izero,Izero,Izero,Izero,Izero
      endif
c
c printout  planet constants
      call NEWPGS(Name,Jdp1,Jdp2,0)
      call NEWPG
      write(Iout, 100) Name, Klan, Nplnt, Jplnt, Jdp1, Jdp0, Jdp2,
     .                 Intp, Ncentr
  100 format('0', A8, 4x, 'NPLNT(', i2, ')=', i3, 4x, 'JPLNT=', i3,
     .       9x, ' JDP1=', i8, 7x, ' JDP0=', i8, 7x, ' JDP2=', i8, 3x,
     .       '  INT =', i3, 3x, 'NCENTR=', i3)
      frp = Intp1
      frp = frp*2._10**Intp2
      frq = Jdp0
      frq = frq + (frp - 0.5_10)
      write(Iout, 200) frq, Ihrp, Iminp, Secp, Intp1, Intp2, frp
  200 format(' INITIAL EPOCH (COORD.TIME) JED=', f16.8, ' IHR=', i2,
     .       ' IMIN=', i2, ' SEC=', f7.4, ' INT1=', i10, ' INT2=', i3,
     .       ' FRACT=', 1pd21.14)
      itf = 22
      write(Iout, 300) (Cond(i), i = 1, 6), Con(1), Con(2),
     .                 (i, Con(i), i = 3, itf)
  300 format(9x, 'A=', 1pd22.15, 9x, 'E=', d22.15, 7x, 'INC=', d22.15,
     .       7x, 'ASC=', d22.15/7x, 'PER=', d22.15, 6x, 'ANOM=',
     .       d22.15, 4x, 'RADIUS=', d22.15, 6x, 'FLAT=',
     .       d22.15/4('   CON(',i2,')=',d22.15))
      itf1 = itf + 1
      itf2 = itf1 + 1
      write(Iout, 400) (i, Con(i), i = itf1, itf2),
     .                 (i, Con1(i), i = 1, 10)
  400 format(2('   CON(',i2,')=',1pd22.15),
     .       2('  CON1(',i2,')=',1pd22.15), /,
     .       (4('  CON1(',i2,')=',1pd22.15)))
      write(Iout, 500) (i, Con1(i), i = 11, 12), (i, Tcon(i), i = 1, 2)
  500 format(2('  CON1(',i2,')=',1pd22.15),
     .       2('  TCON(',i2,')=',1pd22.15))
      write(Iout, 600) (i, Tcon(i), i = 3, 30)
  600 format((4('  TCON(',i2,')=',1pd22.15)))
      write(Iout, 700) (i, Kp(i), i = 1, 100)
  700 format(10(2x,'KP(',i2,')=',i4))
      write(Iout, 800) (i, Kkp(i), i = 1, 100)
  800 format(10(2x,'KKP(',i2,')=',i3))
      ll = (Numki - 8)/25 + 1
      write(Iout, 900) Numki, (Ki(i), i = 1, Numki)
  900 format(' NUMKI=', i3, '  KI=', 7I2, (t30,25I4))
      write(Iout, 1000) (i, Epsp(i), i = 1, 6)
 1000 format(6(' EPSP(',i2,')=',1pe12.5))
      Line  = 43 + ll
      ntyp1 = Kp(100) + 2
      if(Nplnt.le.30 .and. Ncentr.eq.0) ntyp1 = 2
      ntyp2 = Kp(88) + 1
      call PAGCHK(60, 1, 0)
      write(Iout, 1100) Icnd, Iparp, encowm(ntyp1), numeth(ntyp2)
 1100 format(3x, ' ICND =', i3, 3x, 'IPARP =', i3, 3x,
     .       'NUMERICAL INTEGRATION METHOD IS ', a8, 1x, a20)
      if(Jdstrt.gt.0) then
         call PAGCHK(60, 2, 0)
         write(Iout, 1150) Jdstrt
 1150    format('0CHECKPOINT RESTART AT JDSTRT=', i8)
      endif
c
c check consistency of jd1,jd0,jd2
      nstop = 0
      call TCHECK(Name, T1, T0, T2, Intp, nstop)
      if(nstop.gt.0) call SUICID(
     .       ' JD1,JD0,JD2 NOT CONSISTENT, STOP IN PLNHED ',11)
c
c determine real and task time required for setup
      call MVC(Name, 1, 8, mesg, 13)
      call TIMRIT(mesg, 11)
      return
      end
