      subroutine BODSET(lice)

      implicit none

c
c ash/smith/connolly    oct. 1968   subroutine bodset
c setup routine for n-body integration
c m.e.ash   may 1974   asteroid perturbations added
c
c arguments
      integer*4 lice
c lice=0 perturbing planet data set not used in n-body integration
c lice=1 perturbing planet data set used in n-body integration to
c calculate perturbation on motion of earth-moon barycenter

c array dimensions
      include 'globdefs.inc'

c commons
      include 'adams.inc'
      include 'bdctrl.inc'
      include 'bddtaint.inc'
      include 'bodstf.inc'
      include 'ellips.inc'
      include 'empcnd.inc'
      include 'funcon.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'metuna.inc'
      include 'namtim.inc'
      include 'nmstrt.inc'
      include 'param.inc'
      include 'petina.inc'
      include 'plndta.inc'
      include 'rstart.inc'
      include 'stint.inc'
 
c v0    = initial conditions
c n     = number of differential equations
c m     = number of controling equations

c local
      real*10 eta,hmrc,rd8,xsign
      real*4 epsmr(6)
      integer i,inx,j,k,klanb,klref(18),np
      integer*2 id2
      character*96 messag
      logical*1 usepl(9)
c
c definitions:
c nbody = number of bodies being integrated. may include the moon.
c nbodyp= number of planets being integrated, excluding the moon.
c nbodyx= number of bodies included in GR calculations, aside from the
c         sun. may include the moon.
c jembry= pointer to embary in the list of integrated bodies, if present
c jmoon = pointer to moon in the list of integrated bodies, if present
c jmoonx= pointer to moon in extended list for handling relativity
c
c           see if this is checkpoint restart
      Jdstrt = 0
      if(Jdbdy0.le.0) then
         Jdstrt = -Jdbdy0
         call BDREST
      endif
c
c determine initial conditions and body quantities
      Gama  = Gauss**2
      Mmoon = Mass(3)*Mass(10)
c     Mrm   = Mass(10)
      Merth = Mass(3) - Mmoon
      Masse = 1._10 - Mass(10)
      i     = prmter(97) + 0.5001_10
      Tvary0= i
      N     = 0
      Gamat = Gama
      Jembry = -1
      Jmoon = -1
      Ksepem=0
      Orbint=.false.
      Rotint=.false.
      Corint=.false.
      do j=1,9
         usepl(j)=.false.
      end do
      do j = 1, Nbody
         np = Nplbdy(j)
         if(Kbdy(j).ge.0) Ksepem=1
         Mass1(j) = 0._10
         Relftb(j) = 0._10
         if(Kbdy(21).eq.0) Relftb(j) = Relfct
         if(np.gt.0 .and. np.le.9) then
            usepl(np)=.true.
         endif
 
         if(np.eq.3) then
            Jembry   = j
            klanb=-3
         else if(np.eq.10) then
            Jmoon    = j
            klanb=-2
         else if(np.eq.-10) then
            klanb=0
         else
            do k = 1, Numpln
               if(np.eq.Nplnt(k)) then
                  klanb=k
                  goto 50
               endif
            end do
            write(messag,20) Nbody, j, Nplbdy(j)
   20       format(' NO INPUT INITIAL CONDITIONS FOR', i3,
     .             '-BODY INTEGRATION FOR PLANET NPLBDY(', i2, ')=',
     .             i3, ', STOP IN BODSET  ')
            call SUICID(messag, 24)
         endif
 
   50    Name(j)  = Aplnt(klanb)
         Ncbdy(j) = Npcent(klanb)
         klref(j) = klanb
         if(Kbdy(21).gt.0) Relftb(j) = Pcond(30,klanb)
         do i=1,6
            Beta(i,j) = Pcond(i,klanb)
         end do
         if(((Icnd(klanb).ne.0.and.Kbdy(37).le.0) .or.
     .    (Icnd(klanb).eq.0.and.Kbdy(37).gt.0)) .and. np.gt.0) then
            write(messag,55) j,Nplbdy(j),Icnd(klanb),Kbdy(37)
   55       format(' PLANET NPLBDY(',i2,')=',i3,' HAS ICND=',i2,
     .       ', BUT KBDY(37)=',i2,', STOP IN BODSET ')
            call SUICID(messag,17)
         endif
c epoch of orbital elements
         if(Jdpl0(klanb).gt.0 .and. Jdbdy0.gt.0 .and.
     .    Jdpl0(klanb).ne.Jdbdy0) then
            write(messag,56) j,Nplbdy(j),Jdpl0(klanb),Jdbdy0
   56       format(' PLANET NPLBDY(',i2,')=',i3,' HAS JD0=',i7,
     .       ', BUT JDBDY0=',i7,', STOP BODSET ')
            call SUICID(messag,18)
         endif
         if(np.eq.10) then
            Ksepem = 1
            Orbint = .true.
            Mass1(j) = Mmoon
            if(Jdstrt.le.0 .and. Kbdy(37).le.0) then
               call INITL(Gauss, Mass(3), Beta(1,j), 0)
               call ELIPT(0, 0._10)
            endif
         else if(np.eq.-10) then
c if np=-10 appears twice, the second one is the lunar core
            if(.not.Orbint) call SUICID(
     .'CANNOT INTEGRATE LIBRATION WITHOUT LUNAR ORBIT, STOP BODSET ',15)
            if(Rotint) then
               Corint = .true.
               if(Mrcond(23).eq.0._10) call SUICID(
     .'CANNOT INTEGRATE LUNAR CORE WITH MRCON(17)=0, STOP BODSET   ',15)
               do i=1,6
                  Beta(i,j) = Pcond(i+16,klanb)
               end do
            else
               Rotint = .true.
            endif
         else
            if(Orbint) call SUICID('MOON MUST BE LAST, STOP BODSET  ',8)
            if(np.le.30) Mass1(j) = Mass(np)
            if(Jdstrt.le.0 .and. Kbdy(37).le.0) then
               call INITL(Gauss, Mass1(j), Beta(1,j), 1)
               call ELIPT(0, 0._10)
            endif
         endif
         do i = 1,6
            if(Jdstrt.gt.0) then
               X0(i,j) = Rstrt(i,j)
            else if(Kbdy(37).gt.0 .or. np.eq.-10) then
               X0(i,j) = Beta(i,j)
            else
               X0(i,j) = Ylpt(i)
            endif
            N     = N + 1
            V0(N) = X0(i,j)
         end do
      end do
      if(Jembry.le.0) Ksepem=0
      if(Orbint) then
         if(.not.Rotint) call SUICID(
     .'CANNOT INTEGRATE LUNAR ORBIT WITHOUT LIBRATION, STOP BODSET ',15)
         if(Mrcond(23).ne.0._10 .and. .not.Corint) call SUICID(
     .'LUNAR CORE NOT INCLUDED IN N-BODY INTEGRATION   ',-12)
         if(Jembry.le.0) call SUICID(
     .'CANNOT INTEGRATE LUNAR ORBIT WITHOUT EMBARY, STOP BODSET',14)
      endif
c
c setup constants
      Lpert = lice
      Intxpr= Kkbdy(71)
c
c setup for second harmonic of sun
      if(Kbdy(23).ge.0) call SJ2SET

c account for classes of body, depending on whether the moon is
c included in the integration.  if the moon is included, then it is
c already set up as body "jmoon"
      if(Ksepem.le.0) then
c moon not used and not included in newtonian or relativistic forces
         Nbodyp=Nbody
         Nbodyx=Nbody
         Mass1(Nbody+1)=0._10
      else
         Mass1(Jembry)=Merth
         Jmoonx=Jmoon
         if(Jmoonx.le.0) then
            Nbodyp=Nbody
            Nbodyx=Nbody+1
            Jmoonx=Nbodyx
            Nplbdy(Jmoonx)=10
            Ncbdy(Jmoonx)=3
            Mass1(Jmoonx)=Mmoon
            klref(Jmoonx)=-2
         else
            Nbodyp=Jmoon-1
            Nbodyx=Jmoon
         endif
      endif

c setup for violation of principle of equivalence
c Metric parameters beta or gamma not equal to one give principle of
c equivalence violation (unless ratio of gravitational to total energy
c is zero), as does non-zero Deltas.  In this case, the array xsp3 is
c no longer simply the vector divided by distance cubed, but includes
c the violation factor.  Similarly, mfac is no longer 1 + the planet
c mass, but also includes the violation.
      eta    = 4._10*Betapr - Gamapr - 3._10
      Epfacs = 1._10 + eta*Sunpoe + Deltas
      do j=1,Nbodyx
         klanb=klref(j)
            Epfacp(j)= 1._10 + eta*Con11(10,klanb+4)
         if(Nplbdy(j).eq.10) then
            Epfacp(j)= Epfacp(j) + Pcond(25,klanb)-Pcond(26,klanb)
         else
            Epfacp(j)= Epfacp(j) + Pcond(26,klanb)
         endif
         Mfac(j)= Epfacp(j) + Mass1(j)*Epfacs
      end do
      if(Ksepem.gt.0) then
         Mface = Masse*Epfacp(Jembry)+Merth*Epfacs
         Mfacm = Mass(10)*Epfacp(Jmoonx)+Mmoon*Epfacs
      endif

c
c general relativity setup
      Cvel2 = (864E2_10/Aultsc)**2
      B2g2  = 2._10*(Betapm+Gamapm)
      B2g21 = B2g2+1._10
      B2m1  = 2._10*Betapm-1._10
      G21   = 2._10*Gamapm+1._10
      G22   = G21+1._10
      G11   = Gamapm+1._10
      G7    = (4._10*Gamapm+3._10)/2._10
      A44   = (Gauss*Aultsc/864E2_10)**2
      Nsun  =Nbodyx+1
      Mass1(Nsun)=1._10
c supersede solar-system mass possibly set up in prtrd1
      Mascnt=1._10
      do j=1,Nbodyx
         Mascnt=Mascnt+Mass1(j)
      end do
c
c extra setup for moon
      if(Orbint) then
         read(Iplcon) Econ1
         read(Iplcon) Mcon1,Epsm,Km,Jdm1,Jdm2,Intmn,i,i,id2,id2,rd8,
     .    Kkm,id2,Tconm,Numkim,(Kim(i),i=1,Numkim)
         read(Iplcon) Ercon1
         read(Iplcon) Mrcon1,epsmr,Kmr,Jdmr1,Jdmr2,Intmr,i,i,id2,id2,
     .    rd8,Kkmr,id2,Tconmr,Numkir,(Kir(i),i=1,Numkir)
         rewind Iplcon
c enforce conditions of n-body integration or halt if impossible
c no partials
         Numkim=0
         Kim(1)=0
         Numkir=0
         Kir(1)=0
c type of integration (cowell vs encke)
         Km(100)=Kbdy(19)
         if(Jdstrt.gt.0) call SUICID(
     .    'CANNOT RESTART MOON-NBODY INTEGRATION, STOP BODSET  ',13)
c relativity?
         Km(61)=Kbdy(21)
         if(Kbdy(21).ge.0 .and. Kbdy(36).lt.0) Km(61)=2
         Km(97)=Kbdy(21)
c include G-dot?
         Km(62)=Kbdy(22)
c sun J2?
         Km(63)=Kbdy(23)
c perturbing planets
         do np=1,9
c if moon and earth effects on a planet are done separately, then that
c planet is forcibly included as a perturber of the moon
            if(usepl(np)) Km(30+np)=0
c halt if the moon expects a perturber which is not part of the n-body
c (this should never happen)
            if(Km(30+np).ge.0 .and..not.usepl(np)) then
               write(messag,60) np
   60          format('PLANET ',i1,' NEEDED FOR MOON INTEGRATION,',
     .          ' BUT NOT IN N-BODY LIST, STOP IN BODSET')
               call SUICID(messag,19)
            endif
         end do
c extra print on Kout?
         Kkm(7)=Kkbdy(71)
c
         call MORSET(Nbodyp*6)
      endif
c
c setup numerical integration constants
      L1    = 0
      T1    = Jdbdy1
      T2    = Jdbdy2
      Nsign = Jdbdy2 - Jdbdy1
      Nsign = ISIGN(1,Nsign)
      xsign = Nsign
 
c     m=kbdy(30)
      M     = N
      j = Kbdy(31)
      if(j.le.0) then
         Hc = (2._10**j)*xsign
      else
         Hc = j*Nsign
      endif
      j = Kbdy(32)
      if(j.le.0) then
         Hmn = 2._10**j
      else
         Hmn = j
      endif
      inx=Intbdy
      if(Jmoon.gt.0) inx=-1
      if(inx.le.0) then
         Hmx = (2._10**inx)*xsign
      else
         Hmx = inx*Nsign
      endif
      if(Intbdy.le.0) then
         hmrc = (2._10**Intbdy)*xsign
      else
         hmrc = Intbdy*Nsign
      endif
      Epsi  = Epsbdy
      T0    = Jdbdy0
      Tstop = T2 + 10._10*hmrc
      if(Jmoon.gt.0) then
         Bint  = 39._10*Hmx
      else
         Bint  = 9._10*hmrc
      endif

c extra setup for adams-moulton integration
      Npredt = Kbdy(29)
      Ncoret = Kbdy(29)
      Neq    = N
      Nh     = Kbdy(27)
      Inthmx = inx
      Itype  = Kbdy(28)
      if(Jmoon.gt.0 .and. Nh.ge.-1) call SUICID(
     . 'STEP SIZE KBDY(27) TOO COARSE FOR MOON, STOP BODSET ',13)
c
c change initial epoch if this is checkpoint restart
      Both  = 0._10
      Iboth = -1
      if(Jdstrt.gt.0) then
         T0 = Jdstrt
         T1 = T0
c
c alter numerical integration setup for integrating forward
c and backward from initial epoch
      else if(Jdbdy0.ne.Jdbdy1) then
         Iboth = 0
         Nsign = -Nsign
         Hc    = -Hc
         Hmx   = -Hmx
         hmrc  = -hmrc
         Tstop = T1 + 20._10*hmrc
         Bint  = -Bint
      endif
      Ntab = -1
      Mtab = -1
c
c printout title page
      call BODHED
c
c setup asteroid perturbations (individual and ring)
      call BODAST(Astflg)
 
      return
      end
