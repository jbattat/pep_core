      subroutine MORHED
 
      implicit none
c
c king/cappallo   august 1977   subroutine morhed
c printout of data page for moon orbit and rotation numerical
c integration.  based on subroutine monhed (ash/forni jun 67).
c
c
c array dimensions
      include 'globdefs.inc'
c common
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'loadid.inc'
      include 'metuna.inc'
      include 'morstf.inc'
      include 'namtim.inc'
      include 'param.inc'
      include 'plndta.inc'
      include 'rstart.inc'
      include 'stint.inc'
      include 'tapdtplp.inc'
      include 'zeroes.inc'
c
c see subroutine bodred of input link for explanation of
c   mcond(1-30) or mrcond(1-30)   called cond(1-6),con(1-24) in bodred
c   mcon1(1-12) or mrcon1(1-12)   called con1(1-12)          in bodred
c   epsm(1-6) or epsmr(1-6)       called eps(1-6)            in bodred
c   jdm1,jdmn0 (or jdmr0),jdm2    called jd1,jd0,jd2         in bodred
c   imn or ilib                   called itape               in bodred
c   km(1-100) or kmr(1-100)       called k(1-100)            in bodred
c   kim(1-n) or kir(1-n)          called ki(1-n)             in bodred
c   intm or intmr                 called int                 in bodred
c   mplnt                         called nplnt               in bodred
c   mcentr                        called ncentr              in bodred
c
c ipar= number of partial derivatives on output tape + 1
c           iparm for orbit, iparmr for rotation
c           if ipar=1, no partials on tape
c           iparm and iparmr calculated in subroutine morset from
c           kim(1-n) and kir(1-n), respectively
c
c local variables
      integer*4 i,j,jd0,ll,nstop,ntyp1,ntyp2
      integer*2 nparitr,nparit1,nparpl(20)
      character*8 encowm(3)/' COWELL ',' ENCKE  ','MEAN ORB'/
      character*20 numeth(4)/'NORDSIECK (VAR.INT.)',
     .                       'NORDSIECK (CON.INT.)',
     .                       'ADAMS MOULTON       ',
     .                       'ROYAL ROAD (2ND SUM)'/
 
      jd0 = Jdmr0
      if(Orbint) jd0 = Jdmn0
      call NEWPGS(Aplnt(-2),Jdm1,Jdm2,0)
      nparitr=0
      do i=0,i_mxtrg
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
c
c write first two records of moon orbit tape
      if(Orbint) then
         if(Km(99).ge.0 .and. Jdstrt.le.0) then
            write(Imn) Aplnt(-2),Heding,Date,Lnklvl
            write(Imn) Mplnt,Mcentr,Iparm,Intmn,Jdm1,Jdmn0,Jdm2,
     .       u_nmprm,u_nmbod,
     .       Mcond,Mcon1,prmter,Epsm,Km,Npage,Iterat,
     .       izr2,izr4,izr4,izr2,izr2,Zero(1),
     .       Kkm,izr4,Lnklvl,Numkim,(Kim(i),i=1,Numkim),
     .       nparitr,nparit1,(nparpl(i),i=1,nparit1),
     .       Izero,Izero,Izero,Izero,Izero,Izero
         endif
c
c printout moon constants
         call NEWPG
         write(Iout,50) Mplnt,Imn,Jdm1,Jdmn0,Jdm2,Intmn,Mcentr
   50    format('0MOON    ', 4x, 'MPLNT    =', i3, 4x, 'IMN  =', i3,
     .          9x, ' JDM1=', i8, 7x, 'JDMN0=', i8, 7x, ' JDM2=', i8,
     .          3x, ' INTMN=', i3, 3x, 'MCENTR=', i3)
         write(Iout,100) (Mcond(i),i = 1,8),
     .                    (i,Mcond(i+6),i = 3,u_nmbod-6)
  100    format(7x,'A  =', 1pd22.15, 7x, 'E  =', 1pd22.15, 5x,
     .          'INC  =', 1pd22.15, 5x, 'ASC  =', 1pd22.15/5x,
     .          'PER  =', 1pd22.15, 4x, 'ANOM  =', 1pd22.15, 5x,
     .          'MRAD =', 1pd22.15, '  MCON( 2)=', 1pd22.15/,
     .          (4('  MCON(',i2,')=',1pd22.15)))
         write(Iout,150) (i,Mcon1(i),i = 1,12)
  150    format(4(' MCON1(',i2,')=',1pd22.15))
         write(Iout,200) (i,Km(i),i = 1,100)
  200    format(10(2x,'KM(',i2,')=',i4))
         write(Iout,250) (i,Kkm(i),i = 1,20)
  250    format(10('  KKM(',i2,')=',i3))
         ll = (Numkim - 8)/20 + 1
         write(Iout,300) Numkim,(Kim(i),i = 1,Numkim)
  300    format(' NUMKI=', i3, '  KI=', 7I2, (t30,20I5))
         write(Iout,350) (i,Epsm(i),i = 1,6)
  350    format(6(' EPSM(',i2,')=',1pe12.5))
         Line = 37 + ll
         call PAGCHK(60,3,0)
         ntyp1 = Km(100) + 2
         ntyp2 = Km(88) + 1
         write(Iout,400) Iparm,encowm(ntyp1),numeth(ntyp2)
  400    format(3x,'IPARM =', i3, 3x,
     .          'NUMERICAL INTEGRATION METHOD IS ', a8, 1x, a20)
         if(Jdstrt.gt.0) write(Iout,450) Jdstrt
  450    format('0CHECKPOINT RESTART AT JDSTRT=', i8)
      endif
c
c
c write first two records of moon rotation tape
      if(Rotint) then
         if(Kmr(99).ge.0 .and. Jdstrt.le.0) then
            write(Ilib) Aplnt(0),Heding,Date,Lnklvl
            write(Ilib) Mplntr,Mcentr,Iparmr,Intmr,Jdm1,jd0,Jdm2,
     .       u_nmprm,u_nmbod,
     .       Mrcond,Mrcon1,prmter,Epsm,Kmr,Npage,Iterat,
     .       izr2,izr4,izr4,izr2,izr2,Zero(1),
     .       Kkmr,izr4,Lnklvl,Numkir,(Kir(i),i=1,Numkir),
     .       nparitr,nparit1,(nparpl(i),i=1,nparit1),
     .       Izero,Izero,Izero,Izero,Izero,Izero
c write first two records of core rotation tape
            if(Corint .and. Icor.gt.0) then
               write(Icor) ' MONCOR ',Heding,Date,Lnklvl
               write(Icor) Mplntr,Mcentr,Iparmr,Intmr,Jdm1,jd0,Jdm2,
     .          u_nmprm,u_nmbod,
     .          Mrcond,Mrcon1,prmter,Epsm,Kmr,Npage,Iterat,
     .          izr2,izr4,izr4,izr2,izr2,Zero(1),
     .          Kkmr,izr4,Lnklvl,Numkir,(Kir(i),i=1,Numkir),
     .          nparitr,nparit1,(nparpl(i),i=1,nparit1),
     .          Izero,Izero,Izero,Izero,Izero,Izero
            endif
         endif
c
c printout moon rotation constants
         call NEWPG
         write(Iout,500) Mplntr,Ilib,Jdm1,jd0,Jdm2,Intmr
  500    format('0MOON ROTATION    ', 4x, 'MPLNT    =', i3, 4x,
     .          'ILIB=', i3, 9x, ' JDM1=', i8, 7x, 'JDMR0=', i8, 7x,
     .          ' JDM2=', i8, 3x, '  INTMR=', i3)
         write(Iout,550) (Mrcond(i),i = 1,11),
     .                    (j,Mrcond(j+6),j = 6,u_nmbod-6)
  550    format(2x,' PSI =', 1pd22.15, 5x, 'THETA =', 1pd22.15, 7x,
     .          'PHI =', 1pd22.15/2x, 'DPSI =', 1pd22.15, 4x,
     .          'DTHETA =', 1pd22.15, 6x, 'DPHI =', 1pd22.15, /5x,
     .          'MRAD =', 0pf18.12, 9x, 'ALPHA=', 1pd22.15, 5x,
     .          'BETA =', 1pd22.15, 5x, 'GAMMA=', 1pd22.15/4x,
     .          'MEQINC=', 1pd22.15, 3('   CON(',i2,')=',1pd22.15), /,
     .          (4('   CON(',i2,')=',1pd22.15)))
         write(Iout,600) (i,Mrcon1(i),i = 1,12)
  600    format(4('  CON1(',i2,')=',1pd22.15))
         write(Iout,650) (i,Kmr(i),i = 1,100)
  650    format(10(2x,'KMR(',i2,')=',i3))
         write(Iout,670) (i,Kkmr(i),i = 1,20)
  670    format(10(' KKMR(',i2,')=',i3))
         ll = (Numkir - 8)/20 + 1
         write(Iout,300) Numkir,(Kir(i),i = 1,Numkir)
         write(Iout,350) (i,Epsm(i),i = 1,6)
         Line = 27 + ll
         call PAGCHK(60,3,0)
         ntyp1 = Km(100) + 2
         ntyp2 = Km(88) + 1
         write(Iout,700) Iparmr,encowm(ntyp1),numeth(ntyp2)
  700    format(3x,'IPARMR =', i3, 3x,
     .          'NUMERICAL INTEGRATION METHOD IS ', a8, 1x, a20)
         if(Jdstrt.gt.0) then
            write(Iout,720) Jdstrt
  720       format('0CHECKPOINT RESTART AT JDSTRT=', i8)
         endif
c
c check consistency of jd1,jd0,jd2
         nstop = 0
         call TCHECK(Aplnt(-2),T1,T0,T2,Intmn,nstop)
         if(nstop.gt.0) call SUICID(
     .       ' JD1,JD0,JD2 NOT CONSISTENT, STOP IN MORHED ', 11)
      endif
 
      return
      end
