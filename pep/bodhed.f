      subroutine BODHED
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, nstop, ntyp2
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash        march 1968   subroutine bodhed
c printout title page and write first two records of output tape
c for n-body integration
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'bdctrl.inc'
      include 'bodstf.inc'
      include 'bddtaint.inc'
      character*24 namm24
      equivalence (Nammon,namm24)
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'loadid.inc'
      include 'rstart.inc'
      include 'stint.inc'
      include 'zeroes.inc'
c
c quantities for first two records of n-body tape to give it
c the format of a perturbing planet tape
      character*8 nambod/' N-BODY '/, nambh/'XX-BODY '/
      character*4 filbod(9)/' PER', 'TURB', 'ING ', 'PLNT', ' TAP',
     .          'E FR', 'OM N', 'BODY', ' INT'/
      integer*2 nbdyh
      integer*2 intbh(10)
 
      character*4 numeth(5,4)/'NORD','SIEC','K (V','AR.I','NT.)',
     .          'NORD', 'SIEC', 'K (C', 'ON.I', 'NT.)', 'ADAM', 'S MO',
     .          'ULTO', 'N   ', '    ', 'ROYA', 'L RO', 'AD (', '2ND ',
     .          'SUM)'/
      character*4 blank/'    '/
c
c write first two records of n-body output tape
      if(Kbdy(39).ge.0 .and. Jdstrt.le.0) then
         write(Ibody) nambod,Heding,Date,filbod,Lnklvl

c 2nd record describes moon, regardless of whether it's integrated,
c and doesn't describe lunar libration, regardless
         nbdyh = Nbodyp+1
c set tabular intervals for all bodies
         intbh(1) = Intbdy
         if(Intbdy.lt.0) then
            do i = 2, Nbodyp
               intbh(i) = intbh(1) + 1
            end do
         else
            if(Intbdy.eq.0) intbh(1) = 1
            do i = 2, Nbodyp
               intbh(i) = intbh(1)*2
            end do
         endif
         intbh(nbdyh)  = -1

c 2nd record for case in which moon is not one of integrated bodies
c in that case, we have already read the headers for the moon input
         if(Jmoon.le.0) then
            do i=1,6
               Beta(i,nbdyh)=Betabd(i,Nbdy)
            end do
            Relftb(nbdyh)=0._10
c if moon is integrated, then its name, mass, and relfct are all set
         else
            namm24=Name(Jmoon)//Heding(1)//Heding(2)//Date(1)//Date(2)
         endif
         write(Ibody) nbdyh,(Nplbdy(i),i=1,nbdyh),(Ncbdy(i),i=1,nbdyh),
     .    (intbh(i),i=1,nbdyh),Jdbdy1,Jdbdy2,(Jdbdy0,i=1,nbdyh),
     .    ((Beta(i,j),i=1,6),j=1,nbdyh),
     .    (Name(i),Heding(1),Heding(2),Date,i=1,Nbodyp),namm24,
     .    Nmoon,Nbodyp,Intbdy,Jvlbdy,Epsbdy,(Kbdy(i),i=1,40),
     .    Iterat,Npage,(Mass1(i),i=1,nbdyh),(Relftb(i),i=1,nbdyh),
     .    (Kkbdy(i),i=1,80),Lnklvl,(Lnklvl,i=1,nbdyh),
     .    (izr4,i=1,240)
      endif
c we must have nbody=9 with moon being transfered from perturbing
c planet data set or integrated now
c
c print title page for n-body integration
      i = Nbody
      call EBCDI(i,nambh,2)
      call NEWPGS(nambh,Jdbdy1,Jdbdy2,0)
      call NEWPG
      write(Iout,100) Nbody,Ibody,Jdbdy1,Jdbdy0,Jdbdy2,Jvlbdy,
     .                 Epsbdy, Intbdy
  100 format('0  NBODY =', i3, 3x, 'IBODY =', i3, 3x, 'JDBDY1=', i8,
     .       3x, 'JDBDY0=', i8, 3x, 'JDBDY2=', i8, 3x, 'JVLBDY=', i6,
     .       3x, 'EPSBDY=', 1pe12.5, 3x, 'INTBDY=', i3)
      write(Iout,200) (i,Kbdy(i),i = 1,40)
  200 format(10(1x,' KBDY(',i2,')=',i2))
      write(Iout,300) (i,Kkbdy(i),i = 1,80)
  300 format(10(1x,'KKBDY(',i2,')=',i2))
      write(Iout,400) (blank,i,Mass1(i),i = 1,Nbody)
  400 format(4(1A1,'MASS1(',i2,')=',1pd22.15))
      write(Iout,500) (blank,i,Relftb(i),i = 1,Nbody)
  500 format(4(1A1,'RELFT(',i2,')=',1pd22.15))
      do j = 1, Nbody
         write(Iout,550) j,Name(j),Nplbdy(j),Ncbdy(j),
     .                    (Beta(i,j),i = 1,6)
  550    format(i3,'. ', 1A8, '  NPLBDY=', i2, '  NCBDY=', i2, 5x,
     .          ' A =', 1pd22.15, 5x, ' E =', 1pd22.15, 5x, 'INC=',
     .          1pd22.15/39x, 'ASC=', 1pd22.15, 5x, 'PER=', 1pd22.15,
     .          4x, 'ANOM=', 1pd22.15)
      end do
      ntyp2 = Kbdy(28) + 1
      write(Iout,600) (numeth(i,ntyp2),i = 1,5)
  600 format(3x,' NUMERICAL INTEGRATION METHOD IS COWELL ', 5A4)
      if(Jdstrt.gt.0) then
         write(Iout,650) Jdstrt
  650    format('0CHECKPOINT RESTART AT JDSTRT=', i8)
      endif
c
c check consistency of jd1,jd0,jd2
      nstop = 0
      call TCHECK(nambod,T1,T0,T2,Intbdy,nstop)
      if(nstop.gt.0) call SUICID(
     .               ' JD1,JD0,JD2 NOT CONSISTENT, STOP IN BODHED 91  '
     .               , 12)
c
c determine real and task time required for setup
      call TIMRIT('SETUP FOR  N-BODY NUMERICAL INTEGRATION ', 10)
      return
      end
