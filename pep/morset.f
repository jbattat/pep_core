      subroutine MORSET(neqns)
 
      implicit none

c
c r.king and r.cappallo   july 1977   subroutine morset
c setup computations for simultaneous integration of lunar orbit
c and rotation.  based on subroutines monset (m.ash, april 67) and
c mrtset (king/cappallo, aug 72/may 76).
c
c arguments
      integer neqns
c number of equations being integrated for other bodies
c if positive, this indicates n-body integration with no partials

c array dimensions
      include 'globdefs.inc'

c commons 
      include 'adams.inc'
      include 'bddtaint.inc'
      include 'ellips.inc'
      include 'emmips.inc'
      include 'empcnd.inc'
      real*10 alpha,beta,gamma
      equivalence (Mrcond(8),alpha),(Mrcond(9),beta),
     .            (Mrcond(10),gamma)
      real*10 meqinc,dels,deld
      equivalence (Mrcond(11),meqinc)
      equivalence (Mcond(25),dels),(Mcond(26),deld)
      include 'ethhar.inc'
      include 'funcon.inc'
      include 'incon.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'harmor.inc'
      real*10 mnz2,mnc22,mnc33,mnc31
      equivalence (Zhar,mnz2),(Char(2),mnc22),(Char(3),mnc31),
     .            (Char(5),mnc33)
      include 'metuna.inc'
      include 'monhar.inc'
      include 'morstf.inc'
      include 'namtim.inc'
      include 'nmstrt.inc'
      include 'orblun.inc'
      include 'param.inc'
      include 'plndta.inc'
      include 'prtcod.inc'
      include 'rstart.inc'
      include 'stint.inc'
      include 'xprcom.inc'
      integer*4 zxprcm/2982/   !r8=1494,r10=2982

c external functions
      real*10 LEGSCL,DOT

c local variables
      real*10 aukm,cmr2,emmr,eta,gme,gmm,quick
      real*10 mass3i,tsav,tsavc,tstrt,xac12,xh,xn,xsign
      integer*2 mplnti,corics(6)
      integer i,icoff,iroff,iparsv,j,j1,j2,jd0,jj,k,kim1,king,kintg,
     . kmsav,kong,kount,ktop,ktopa,l,lintg,llib,lnut,lrel,n1,n2
      character*72 errmsg
c
c some setup already done if this is part of n-body integration
      Knbd=neqns
      if(Knbd.gt.0) goto 100

c setup flags for simultaneous integration
      Orbint = .true.
      Rotint = .true.
      if(Jdmn0.eq.0 .or. Jdmn0.le.-3000000) Orbint = .false.
      if(Jdmr0.eq.0 .or. Jdmr0.le.-3000000) Rotint = .false.
      Corint = Rotint .and. Mrcond(23).ne.0._10
      Icor = Kkmr(11)
      jd0 = Jdmn0
      mplnti=Mplnt
      if(Rotint .and. Jdmn0.gt.0 .and. Jdmn0.ne.Jdmr0) call SUICID(
     . 'DIFFERENT JD0 FOR ORBIT AND ROTATION, STOP IN MORSET',13)
      if( .not. Orbint) then
         jd0  = Jdmr0
         Jdm1 = Jdmr1
         Jdm2 = Jdmr2
         mplnti=Mplntr
c if rotation alone integrated stuff kmr integration controls into
c km storage locations, but always need moon coordinates
         do i = 31, 39
            Km(i) = Kmr(i)
         end do
         Km(40)=1
         do i = 87, 92
            Km(i) = Kmr(i)
         end do
         Kkm(7)=Kkmr(7)
      endif
      Korb = 0
c
c see if this is checkpoint-restart
      Jdstrt = 0
      if(jd0.le.0) then
         Jdstrt = -jd0
         iroff  = 0
         icoff  = 0
         if(Rotint) then
            call RESTRT(Ilib, Mplntr, Mcentr, Intmn, Jdm1, jd0, Jdm2,
     .                  Mrcond, Mrcond(7), Mrcon1, Epsm, Kmr, Numkir,
     .                  Kir)
            if(.not.(Orbint.or.Corint)) goto 100
c need to clear out low order part of rstrt for orbit and/or core data
c iroff gives the starting location of the transferred rotation data
            iroff = i_mxeqn - Iparst
            do i = 1, Iparst
               do j = 1, 6
                  Rstrt(j,iroff+i) = Rstrt(j,i)
               end do
            end do
            tsav = Jdstrt+Frstrt
            if(Corint) then
c the core rotation must have been done with the same setup as the
c external rotation -- the restart uses the same locations for both
               iparsv=Iparst
               call RESTRT(Icor,Mplntr,Mcentr,Intmn,Jdm1,jd0,Jdm2,
     .          Mrcond,Mrcond(7),Mrcon1,Epsm,Kmr,Numkir,Kir)
c two minimal compatibility checks
               tsavc=Jdstrt+Frstrt
               if(Iparst.ne.iparsv .or. tsavc.ne.tsav) call SUICID(
     .'RESTART OF MISMATCHING LUNAR CORE AND MANTLE ROTATIONS, STOP MORS
     .ET ', 17)
               if(.not.Orbint) goto 100
c need to clear out low order part of rstrt again
               icoff = i_mxeqn - Iparst - iparsv
               do i = 1,Iparst
                  do j = 1,6
                     Rstrt(j,icoff+i) = Rstrt(j,i)
                  end do
               end do
            endif
         endif
 
c orbital restart ics
         call RESTRT(Imn,Mplnt,Mcentr,Intmn,Jdm1,jd0,Jdm2,Mcond,
     .               Mcond(7),Mcon1,Epsm,Km,Numkim,Kim)
c parallel integration restart -be sure the restart initial
c point falls at the beginning of a record for both o&r
         tstrt = Jdstrt + Frstrt
         if(Rotint .and. tstrt.ne.tsav) call SUICID(
     .'RESTART OF PARALLEL INTEGRATION MUST BEGIN ON A MULTIPLE OF 20 DA
     .YS FROM THE ORIGINAL EPOCH ', 23)
      endif
c
c setup for second harmonic of sun (already done if n-body integration)
      if(Km(63).ge.0) call SJ2SET
c
c setup of constants
  100 Gama   = Gauss**2
      Gamat  = Gama
      Gamem  = Gama*Mass(3)
      Gamtem = Gamem
      Mmoon  = Mass(3)*Mass(10)
      Masse  = 1._10 - Mass(10)
      Merth  = Masse*Mass(3)
      Gama3  = Gamem*Masse
      Gamat3 = Gama3
      i      = prmter(97) + 0.5001_10
      Tvary0 = i
      eta    = 4._10*Betapr - Gamapr - 3._10
      Delpls = Econ1(10) + Mcon1(10)
      Delmns = Econ1(10) - Mcon1(10)
      Delsum = dels + eta*Delpls
      Deldif = deld + eta*Delmns
 
c derived constants used only for printout
      emmr   = 1._10/Mass(10) - 1._10
      gme    = Gama*Masse*(Aultsc**3)*(Ltvel**3)/(86400._10**2)*Mass(3)
      gmm    = gme/emmr
      aukm   = Ltvel*Aultsc
      Cvel   = 86400._10/Aultsc
      mass3i = 1._10/Mass(3)
      Index1 = 3
      Jndex1 = 3
      Index2 = 3
      Ntab   = -1
      Ntabr  = -1
      Intxpr = Kkm(7)
 
c read first records of perturbing planet dataset
      lnut = 0
      if(Km(81).ge.0 .or. Km(84).ge.0 .or. Kmr(82).ge.0) lnut = 1
      llib = 0
      if(Km(82).ge.0) llib = 1
      lrel = 0
      if(Km(61).ge.0 .or. Kmr(83).ge.0) lrel = 1
c be sure solar position is calculated if used
      Dosun=Km(33)
      if(Rotint .and. (Kmr(81).ge.0 .or. Kmr(83).ge.0)) Dosun = 1
      if(Orbint .and. (Km(61).ge.0 .or. Km(85).ge.0)) Dosun = 1
      do i=1,9
         if(Km(i+30).ge.0) Dosun=1
      end do
      kmsav = Km(33)
      Km(33)=Dosun
      if(Knbd.eq.0) call PRTRD1(0,mplnti,Mcentr,Km(31),lrel,lnut,llib)
      Km(33) = kmsav
 
      do i = 1, i_mxplprt
         Icmtrl(i) = 0
         Icrtrl(i) = 0
      end do
 
c zero out harmonic vectors
         do i = 1, 19
            Imzone(i) = 0
            Izone(i) = 0
         end do
         do i = 1, 54
            Imcos(i) = 0
            Imsin(i) = 0
            Icos(i) = 0
            Isin(i) = 0
         end do
 
c set all initial conditions equal to zero
      do j=Knbd+1,6*i_mxeqn
         V0(j) = 0._10
      end do

      Iparm=0
      Iparmr=0

c setup for target planet quantities
      do i = 1, i_mxtrg
         Ntrg(i)= 0
      end do
      Numtar = 0
      i=1
      do while (i.le.i_mxtrg .and. Km(i).gt.0 .and. Knbd.eq.0)
         king = Km(i)
         if(king.gt.9) then
            write(errmsg,110) king
  110       format('TARGET BODY',i3,' NOT A PLANET, STOP IN MORSET ')
            call SUICID(errmsg,11)
         endif
         Ntrg(i) = king
         Numtar  = i
         if(Km(30+king).le.0) then
            write(errmsg,120) king,30+king,Km(30+king)
  120       format('TARGET BODY',i3,' NOT INCLUDED IN PARTIALS, KM(',i2,
     .       ')=',i3,', STOP IN MORSET ')
            call SUICID(errmsg,17)
         endif
         i=i+1
      end do

c set up for extra print on KOUT
      if(Kout.gt.0 .and. Kkm(7).gt.0) then
         call ZFILL(Relacc,zxprcm)
         Xprnam(1)=Aplnt(-2)
         Xprnam(3)=Aplnt(-3)
         write(Kout,125)
  125    format(' ACCELERATIONS AND OTHER INTERMEDIATE QUANTITIES')
      endif
c
c-----------------------------------------------------------------------
c
c           setup for orbit
c

      if(Orbint) then
c
c setup initial condition for equations of motion and
c equations for partial derivatives w.r.t. initial osculating
c eliptic orbital elements
         if(Icnd(-2).eq.0) then
            call INITL(Gauss, Mass(3), Mcond, 0)
            Km1 = Kim(1)
            Km1 = iabs(Km1)
            call ELIPT(Km1, 0._10)
            if(Kim(1).lt.0) Km1 = 0
            do j = 1, 6
               X0m(j, 1) = Ylpt(j)
               do i = 1, 6
                  Dx0m(i, j, 1) = Dylpt(i, j)
               end do
            end do
         else if(Icnd(-2).eq.-1) then
c initial conditions are cartesian coordinates, cannot do partials
            do j=1,6
               X0m(j,1)=Mcond(j)
               Ylpt(j)=Mcond(j)
               if(Kim(j+1).ge.0 .and. Kim(1).ne.0) call SUICID(
     .'CANNOT INTEGRATE PARTIALS FOR CARTESIAN ICS, STOP MORSET',14)
            end do
            if(Km(100).ge.0) call SUICID(
     .'CANNOT USE ENCKE INTEGRATION FOR CARTESIAN ICS, STOP MORSET ',15)
         else
            call SUICID('INVALID ICND, STOP MORSET  ',7)
         endif
         Nrec  = 0
         Iparm = 1
         N     = Knbd+6
         if(Km(100).lt.0) then
 
            do j = 1, 6
               V0(Knbd+j) = Ylpt(j)
            end do
            if(Kim(1).gt.0) then
               do j = 1, 6
                  if(Kim(j+1).ge.0) then
                     do i = 1, 6
                        N     = N + 1
                        V0(N) = Dylpt(i, j)
                     end do
                  endif
               end do
            endif
         else if(Km(100).eq.0) then
 
            if(Kim(1).gt.0) then
               do i = 2, 7
                  if(Kim(i).ge.0) N = N + 6
               end do
            endif
         else
 
            kim1 = Kim(1)
            call LUNSET(Gauss, Mass(3), Jdmn0, 0._10, kim1)
            do j = 1, 6
               V0(Knbd+j) = Ylpt(j) - Ylun(j)
            end do
            if(Kim(1).gt.0) then
               do j = 1, 6
                  if(Kim(j+1).ge.0) then
                     do i = 1, 6
                        N     = N + 1
                        V0(N) = Dylpt(i, j) - Dylun(i, j)
                     end do
                  endif
               end do
            endif
         endif
 
         if(Kim(1).ne.0) then
            do i = 2, 7
               if(Kim(i).ge.0) Iparm = Iparm + 1
            end do
         endif
c
c
c           setup for equations for partial derivatives with respect to
c           non-initial conditions
c
c           setup of moon  logic controls for partial derivatives
         kount = 0
         if(Kim(1).gt.0) then
            do i = 2, 7
               if(Kim(i).ge.0) then
                  kount = kount + 1
                  Icmtrl(kount) = -30 - (i - 1)
               endif
            end do
         endif
         i = 7
         do while(i.lt.Numkim)
            i = i + 1
            if(Kim(i).eq.0) goto 150
            kong = Kim(i)/100
            if(kong.ne.0) then
               king = iabs(Kim(i) - kong*100)

c initial conditions and parameters
               if(king.ge.1 .and. king.le.30) then
                  goto 130
                  
c zonal harmonics
               else if(king.eq.31) then
                  if(kong.ne.10) goto 180
                  n1 = Kim(i+1) - 1
                  n2 = Kim(i+2) - 1
                  do j = n1, n2
                     kount = kount + 1
                     Icmtrl(kount) = Kim(i)
                     Iparm = Iparm + 1
                     N = N + 6
                     Imzone(j) = kount
                  end do
                  i = i + 2
                  goto 140
 
c tesseral cosine harmonics
               else if(king.eq.41) then
                  if(kong.ne.10) goto 180
                  n1 = (Kim(i+1)*(Kim(i+1)-1))/2 + Kim(i+2) - 1
                  n2 = (Kim(i+3)*(Kim(i+3)-1))/2 + Kim(i+4) - 1
                  do j = n1, n2
                     kount = kount + 1
                     Icmtrl(kount) = Kim(i)
                     Iparm = Iparm + 1
                     N = N + 6
                     Imcos(j) = kount
                  end do
                  i = i + 4
                  goto 140
 
c tesseral sine harmonics
               else if(king.eq.51) then
                  if(kong.ne.10) goto 180
                  n1 = (Kim(i+1)*(Kim(i+1)-1))/2 + Kim(i+2) - 1
                  n2 = (Kim(i+3)*(Kim(i+3)-1))/2 + Kim(i+4) - 1
                  do j = n1, n2
                     kount = kount + 1
                     Icmtrl(kount) = Kim(i)
                     Iparm = Iparm + 1
                     N = N + 6
                     Imsin(j) = kount
                  end do
                  i = i + 4
                  goto 140
c anything else not implemented
               else
                  goto 180
               endif
            endif
  130       kount = kount + 1
            Icmtrl(kount) = Kim(i)
            Iparm = Iparm + 1
            N     = N + 6
 
  140    end do
  150    if(kount.gt.i_mxplprt) call SUICID(
     .       'TOO MANY ORBIT PARTIALS IN MORSET   ', 9)
         Korb = N
         if(Knbd.eq.0) call EMPRD1(0)
c
c setup non-zero initial conditions for equations for
c partials
         if(Jdstrt.le.0) then
            do j = 1, i_mxplprt
               if(Icmtrl(j).eq.-8) then
                  jj = 6*j
                  if(Km(61).lt.0) call SUICID(
     .       'GEODESIC PRECESSION PARTIAL WITHOUT REL, STOP IN MORSET '
     .       , 14)
                  call PRTCRD(jd0, 0._10)
                  call PLPCRD(jd0, 0._10)
                  Rpert2(3) = DOT(Xpert(1,3), Xpert(1,3))
                  Rpert(3)  = SQRT(Rpert2(3))
                  Rpert3(3) = Rpert2(3)*Rpert(3)
                  call CROSS(Xpert(4,3), Xpert(1,3), Xpert3)
                  call CROSS(Xpert3, X0m, V0(jj+4))
                  quick = -(Gamapm + 0.5_10)*Gamat/Cvel**2/Rpert3(3)
                  do l = 4, 6
                     V0(jj + l) = V0(jj + l)*quick
                  end do
               else if(Icmtrl(j).eq.3) then
                  jj    = 6*j
                  quick = 0.5_10/Mass(3)
                  do l = 4, 6
                     V0(jj + l) = quick*Ylpt(l)
                  end do
               endif
            end do
         endif
      endif
c-----------------------------------------------------------------------
c
c setup for rotation
c
      if(Rotint) then
 
         if( .not. Orbint) call EMPRD1(0)
c
c setup initial condition for equations of motion and
c equations for partial derivatives w.r.t. initial conditions
         j1 = 1
         if(Orbint) j1 = 2
         do j = 1, 6
            X0m(j, j1) = Mrcond(j)
            if(Corint) X0m(j,j1+1) = Mrcond(j+16)
            do i = 1, 6
               Dx0m(i, j, j1) = 0._10
               Dx0m(i, j, j1+1) = 0._10
            end do
            Dx0m(j, j, j1) = 1._10
            Dx0m(j, j, j1+1) = 1._10
         end do
         do j = 1, 6
            V0(j + Korb) = X0m(j, j1)
         end do
         N = Korb + 6
         Iparmr = 1
         if(Kir(1).gt.0) then
            do j = 1, 6
               if(Kir(j+1).ge.0) then
                  Iparmr = Iparmr + 1
                  do i = 1, 6
                     N     = N + 1
                     V0(N) = Dx0m(i, j, j1)
                  end do
               endif
            end do
         endif
         kount = 0
         if(Kir(1).gt.0) then
            do i = 2, 7
               if(Kir(i).ge.0) then
                  kount = kount + 1
                  Icrtrl(kount) = -30 - (i - 1)
c cross-reference to same partial derivative in orbit list, if present
                  do j=1,Iparm-1
                     if(Icmtrl(j).eq. -1000 - (i-1)) then
                        Icrref(kount)=j
                        Icmref(j)=kount
                     endif
                  end do
               endif
            end do
         endif
c
c setup for equations for partial derivatives with respect to
c non-initial conditions
 
c start of loop i=8,numkir
         i = 7
         do while( i.lt.Numkir )
            i = i + 1
            if(Kir(i).eq.0) goto 200
            kong = Kir(i)/100
            if(kong.gt.0) then
               if(kong.gt.10) goto 170
               king = Kir(i) - kong*100

c initial conditions and parameters
               if(king.ge.1 .and. king.le.30) then
                  kount = kount + 1
                  Icrtrl(kount) = Kir(i)
                  Iparmr = Iparmr + 1
                  N = N + 6
                  k=Kir(i)
                  if(kong.eq.10) then
                     k=-30-king
                     if(king.gt.6) k=6-king
                  endif
                  do j=1,Iparm-1
                     if(Icmtrl(j).eq.k) then
                        Icrref(kount)=j
                        Icmref(j)=kount
                     endif
                  end do
                  goto 160
                  
c zonal harmonics
               else if(king.eq.31) then
                  if(kong.ne.10) goto 170
                  n1 = Kir(i+1) - 1
                  n2 = Kir(i+2) - 1
                  do j = n1, n2
                     kount = kount + 1
                     Icrtrl(kount) = Kir(i)
                     Iparmr = Iparmr + 1
                     N = N + 6
                     Izone(j) = kount
                     if(Imzone(j).gt.0) then
                        Icrref(kount)=Imzone(j)
                        Icmref(Imzone(j))=kount
                     endif
                  end do
                  i = i + 2
                  goto 160
 
c tesseral cosine harmonics
               else if(king.eq.41) then
                  if(kong.ne.10) goto 170
                  n1 = (Kir(i+1)*(Kir(i+1)-1))/2 + Kir(i+2) - 1
                  n2 = (Kir(i+3)*(Kir(i+3)-1))/2 + Kir(i+4) - 1
                  do j = n1, n2
                     kount = kount + 1
                     Icrtrl(kount) = Kir(i)
                     Iparmr = Iparmr + 1
                     N = N + 6
                     Icos(j) = kount
                     if(Imcos(j).gt.0) then
                        Icrref(kount)=Imcos(j)
                        Icmref(Imcos(j))=kount
                     endif
                  end do
                  i = i + 4
                  goto 160
 
c tesseral sine harmonics
               else if(king.eq.51) then
                  if(kong.ne.10) goto 170
                  n1 = (Kir(i+1)*(Kir(i+1)-1))/2 + Kir(i+2) - 1
                  n2 = (Kir(i+3)*(Kir(i+3)-1))/2 + Kir(i+4) - 1
                  do j = n1, n2
                     kount = kount + 1
                     Icrtrl(kount) = Kir(i)
                     Iparmr = Iparmr + 1
                     N = N + 6
                     Isin(j) = kount
                     if(Imsin(j).gt.0) then
                        Icrref(kount)=Imsin(j)
                        Icmref(Imsin(j))=kount
                     endif
                  end do
                  i = i + 4
                  goto 160
c anything else not implemented
               else
                  goto 170
               endif
            endif
 
c usual parameter
            kount = kount + 1
 
c must be a valid parameter for which we have var'l. eqns.
            if(Kir(i).ne.-3 .and. Kir(i).ne.-4 .and. Kir(i).ne.-6
     .          .and. Kir(i).ne.-7 .and.
     .       (Kir(i).gt.-11 .or. Kir(i).lt.-18) .and.
     .       (Kir(i).lt.1 .or. Kir(i).gt.50)) goto 170
            Icrtrl(kount) = Kir(i)
            Iparmr = Iparmr + 1
            N = N + 6
            k=Kir(i)
            if(k.lt.0) k=k-1006
            do j=1,Iparm-1
               if(Icmtrl(j).eq.k) then
                  Icrref(kount)=j
                  Icmref(j)=kount
               endif
            end do
 
c end of loop i=8,numkir
  160    end do

         if(i.gt.Numkir .and. Numkir.gt.0) call SUICID(
     .    'KI EXCEEDS INPUT LIST IN MORSET ',8)
 
  200    if(kount.gt.i_mxplprt) call SUICID(
     .    'TOO MANY ROTATION PARTIALS IN MORSET', 9)
         Krot = N

         if(Corint) then
            do j = 1, 6
               V0(j + Krot) = X0m(j, j1+1)
            end do
            N = Krot + 6
            if(Kir(1).gt.0) then
               do j = 1, 6
                  if(Kir(j+1).ge.0) then
                     N = N + 6
                  endif
               end do
            endif

c get non-zero initial values for partials (only w.r.t. initial conditions)
            do j=1,6
               corics(j)=0
            end do
c start of loop i=8,numkir
            i = 7
            do while( i.lt.Numkir )
               i = i + 1
               if(Kir(i).ge.0 .or. Kir(i).lt.-16) goto 250
               if(Kir(i).le.-11 .and. Kir(i).ge.-16) then
                  j=-Kir(i)-10
                  corics(j)=N
                  do jj = 1, 6
                     V0(N+jj) = Dx0m(jj,j,j1+1)
                  end do
               endif
               N=N+6
            end do
c partials for core rotation are exactly the same in number and order
c as the partials for external rotation
  250       N= Krot + (Krot-Korb)
         endif

      endif

      if(Kkm(60).gt.0)
     . write(6,99222) (Icmtrl(j),Icrtrl(j),Icmref(j),Icrref(j),j=1,15)
99222 format('Icmtrl,Icrtrl,Icmref,Icrref'/(3(4i5,2x)))
c
c-----------------------------------------------------------------------
c
c numerical integration setup
      if(Knbd.eq.0) then
         L1    = 0
         M     = Km(90)
         T1    = Jdm1
         T2    = Jdm2
         Nsign = Jdm2 - Jdm1
         Nsign = ISIGN(1, Nsign)
         xsign = Nsign
         xh    = Km(91)
         kintg = Km(91)
         if(xh.le.0) then
            Hc = 2._10**kintg*xsign
         else
            Hc = xh*xsign
         endif
         xn    = Km(92)
         lintg = Km(92)
         if(xn.le.0) then
            Hmn = 2._10**lintg
         else
            Hmn = xn
         endif
         if(Intmn.lt.0) then
            Hmx = (2._10**Intmn)*xsign
         else if(Intmn.eq.0) then
            Hmx = 0.5_10*xsign
         else
            Hmx = Intmn*Nsign
         endif
         Epsi  = Epsm(3)
         T0    = jd0
         Tlpt0 = T0
         Tstop = T2 + 8._10*Hmx
 
c moon tape has 8 tabular points per record, rotation tape has 5
         Mint  = 7._10*Hmx
         Mint5 = 4._10*Hmx
c because only jd and not jd,fract is written at the start of each
c 8 tabular point moon record, intm must be greater or equal to -3
         Msign = ISIGN(1, Nsign*Ibdsgn)
c
c extra setup for adams-moulton integration
         Npredt = Km(89)
         Ncoret = Km(89)
         Neq    = N
         Nh     = Km(87)
         if(Orbint) Intmr = Intmn
         Inthmx = Intmn
         Itype  = Km(88)
c
c change initial conditions if this is checkpoint-restart
         Both  = 0._10
         Iboth = -1
         if(Jdstrt.gt.0) then
            T0 = Jdstrt + Frstrt
            T1 = T0
            if(Orbint) then
               if(Km(100).ge.0) then
                  T = Jdstrt - jd0
                  if(Km(100).gt.0) then
                     call LUNORB(Jdstrt, Frstrt, -1)
                  else
                     call ELIPT(1, T)
                  endif
               endif
               do i = 1, 6
                  X0m(i, 1) = Rstrt(i, 1)
                  if(Km(100).lt.0) then
                     V0(i) = Rstrt(i, 1)
                  else if(Km(100).eq.0) then
                     V0(i) = Rstrt(i, 1) - Ylpt(i)
                  else
                     V0(i) = Rstrt(i, 1) - Ylun(i)
                  endif
               end do
               j1 = 1
               l  = 6
               if(Kim(1).ne.0) then
                  do j = 1, 6
                     if(Kim(j+1).ge.0) then
                        j1 = j1 + 1
                        do i = 1, 6
                           Dx0m(i, j, 1) = Rstrt(i, j1)
                           l = l + 1
                           if(Kim(1).lt.0 .or. Km(100).lt.0) then
                              V0(l) = Rstrt(i, j1)
                           else if(Km(100).eq.0) then
                              V0(l) = Rstrt(i, j1) - Dylpt(i, j)
                           else
                              V0(l) = Rstrt(i, j1) - Dylun(i, j)
                           endif
                        end do
                     endif
                  end do
               endif
               do while( j1.lt.Iparst )
                  j1 = j1 + 1
                  do i = 1, 6
                     l     = l + 1
                     V0(l) = Rstrt(i, j1)
                  end do
               end do
            endif
            if(Rotint) then
 
c setup checkpoint restart i.c.s for rotation
               j1 = 1
               if(Orbint) j1 = 2
               do i = 1, 6
                  X0m(i,j1)   = Rstrt(i,iroff+1)
                  V0(Korb+i) = Rstrt(i,iroff+1)
               end do
               j2 = 1
               l  = 6
               if(Kir(1).ne.0) then
                  do j = 1, 6
                     if(Kir(j+1).ge.0) then
                        j2 = j2 + 1
                        do i = 1, 6
                           Dx0m(i, j, j1) = Rstrt(i, iroff + j2)
                           l = l + 1
                           V0(Korb + l) = Rstrt(i, iroff + j2)
                        end do
                     endif
                  end do
               endif

c transfer the rest of the variational equation ics
               if(j2.ne.Iparmr) then
                  j2 = j2 + 1
                  do j = j2, Iparmr
                     do i = 1, 6
                        l = l + 1
                        V0(Korb + l) = Rstrt(i, iroff + j)
                     end do
                  end do
               endif
               if(Corint) then
c corresponding setup for core rotation
                  do i = 1, 6
                     X0m(i,j1+1) = Rstrt(i,icoff+1)
                     V0(Krot+i) = Rstrt(i,icoff+1)
                  end do
                  do j = 1, 6
                     j2 = corics(j)/6+1
                     if(j2.gt.1) then
                        do i = 1, 6
                           Dx0m(i,j,j1+1) = Rstrt(i,icoff+j2)
                        end do
                     endif
                  end do

c transfer the rest of the variational equation ics
                  l=6
                  do j = 2,Iparmr
                     do i = 1, 6
                        l = l + 1
                        V0(Krot+l) = Rstrt(i,icoff+j)
                     end do
                  end do
               endif
            endif
c
c alter numerical integration setup for integrating forward
c and backward from initial epoch
         else if(jd0.ne.Jdm1) then
            Iboth = 0
            Nsign = -Nsign
            Hc    = -Hc
            Hmx   = -Hmx
            Mint  = -Mint
            Mint5 = -Mint5
            Msign = -Msign
            Tstop = T1 + 16._10*Hmx
         endif
      endif
c
c-----------------------------------------------------------------------
c           earth and moon gravitational potential harmonic setup
c
c           setup earth harmonic controls
      Nerzon = 0
      Keh    = Km(81)
      if(Kmr(82).ge.0 .and. Keh.lt.2) Keh = 2
      if(Keh.ge.0) then
         if(Keh.gt.20) call SUICID(
     .    'TOO MANY EARTH HARMONIC COEFFICIENTS, STOP IN MORSET',13)
         Nerzon = Keh
         Nezon1 = Nerzon - 1
         Nezonp = 1
         if(Nezon1.gt.0) then
            do i = 1, Nezon1
               Erzhar(i) = Ezhar(i)
            end do
         endif
      endif
      Erad = Econd(7)/(Aultsc*Ltvel)
c
c
c setup moon harmonic controls
      Kmh = Km(82)
      if(Kmh.ge.0 .or. Rotint) then
         if(Nmzone.gt.20 .or. Nmtess.gt.10) call SUICID(
     .    'TOO MANY MOON HARMONIC COEFFICIENTS, STOP IN MORSET ',13)
         Jzone = 1031
         Jcos  = 1041
         Jsin  = 1051
         Nzone = Nmzone
         Nzon1 = Nzone - 1
         Nzonp = 1
         do i = 1, Nzon1
            Zhar(i) = Mzhar(i)
         end do
         Ntess = Nmtess
         Ntes1 = Ntess - 1
         Ntesp = 1
         k = 0
         do j = 2, Ntess
            do i = 1, j
               k = k + 1
               Char(k) = Mchar(k)*LEGSCL(j, i)
               Shar(k) = Mshar(k)*LEGSCL(j, i)
            end do
         end do
         Mrad = Mcond(7)/(Aultsc*Ltvel)
c
c no. of harmonics included in partials integration
         Ntop = max0(Nzon1, Ntes1)
         if(Km(86).ge.0 .or. Kmr(86).ge.0) then
            ktop = Km(86)
            if(Kmr(86).gt.ktop) ktop = Kmr(86)
            ktopa = ktop/100
            Ntesp = min0(ktopa - 1, Ntes1)
            ktop  = ktop - ktopa*100 - 1
            Nzonp = min0(ktop, Nzon1)
         endif
c
c calculate derived constants for integration
         alpha = (beta - gamma)/(1._10 - beta*gamma)
         Mirat(1)=alpha
         Mirat(2)=beta
         Mirat(3)=gamma
 
c c/(m*r**2)
         cmr2  = 2._10*mnz2*(1._10+beta)/(2._10*beta-gamma+beta*gamma)
         mnc22 = gamma*cmr2/4._10
 
c ratios (lunar mass)/(moment of inertia) (1/au**2)
         Mmc = 1E30_10
         if(cmr2.ne.0._10) Mmc = 1._10/cmr2/Mrad**2
         Mma = Mmc*(1._10 + beta)/(1._10 - beta*gamma)
         Mmb = Mmc*(1._10 + beta)/(1._10 + gamma)
         Mmaw= Mma
         Mmbw= Mmb
         Mmcw= Mmc
 
c if meqinc value input dont derive
         if(meqinc.le.0._10) then
            meqinc = (5558.5_10 + 347.6_10*(beta-.00063_10)/.00003_10 -
     .               (gamma-.00022_10)/.00002_10)*Convds
            if(Ntess.gt.2) then
 
c xac12= .4/( .5*(a/m*r**2 + b/m*r**2)  )
               xac12  = (.8_10*Mrad**2)/(1._10/Mma + 1._10/Mmb)
               meqinc = meqinc +
     .                  (-20.303_10*mnc31/0.286E-4_10 + 5.134_10*mnc33/
     .                  0.27E-5_10)*xac12*Convds
            endif
         endif

         if(Rotint) then
c
c form rigid body inertia tensor and its inverse
            do i = 1, 3
               do j = 1, 3
                  I0(i,j)  = 0._10
                  I0i(i,j) = 0._10
               end do
            end do
            I0(1,1)  = Mmoon/Mma
            I0i(1,1) = Mma/Mmoon
            I0(2,2)  = Mmoon/Mmb
            I0i(2,2) = Mmb/Mmoon
            I0(3,3)  = Mmoon/Mmc
            I0i(3,3) = Mmc/Mmoon
            if(Kmr(83).ge.0) Index1 = 6
            Awhole = I0(1,1)
            Bwhole = I0(2,2)
            Cwhole = I0(3,3)
            Amantl = Awhole
            Bmantl = Bwhole
            Cmantl = Cwhole
            if(Corint) then
c substitute mantle inertia tensor for whole moon version
               Amantl = Awhole - Mrcond(24)
               if(Amantl.le.0._10 .or. Mrcond(24).lt.0._10) call SUICID(
     .   'UNPHYSICAL LUNAR CORE MOMENT OF INERTIA, STOP MORSET',13)
               Bmantl = Bwhole - Mrcond(24)
               Cmantl = Cwhole - Mrcond(24)
               Mma    = Mmoon/Amantl
               Mmb    = Mmoon/Bmantl
               Mmc    = Mmoon/Cmantl
               I0(1,1)= Amantl
               I0(2,2)= Bmantl
               I0(3,3)= Cmantl
               I0i(1,1)= 1._10/Amantl
               I0i(2,2)= 1._10/Bmantl
               I0i(3,3)= 1._10/Cmantl
c substitute mantle ratios for whole moon versions
               Mirat(1)=Mirat(1)*Awhole/Amantl
               Mirat(2)=Mirat(2)*Bwhole/Bmantl
               Mirat(3)=Mirat(3)*Cwhole/Cmantl
c core inertia tensor
               do i=1,3
                  do j=1,3
                     I0c(i,j)  = 0._10
                     I0ci(i,j) = 0._10
                  end do
               end do
               I0c(1,1) = Mrcond(24)
               I0c(2,2) = Mrcond(24)
               I0c(3,3) = Mrcond(24)
               I0ci(1,1) = 1e-10_10
               if(Mrcond(24).gt.0._10) I0ci(1,1) = 1._10/Mrcond(24)
               I0ci(2,2) = I0ci(1,1)
               I0ci(3,3) = I0ci(1,1)
            endif
         endif
c
c calculate partials matrix for libration parameters
c must set up J2 partial even if not integrating rotation
         if(Knbd.eq.0) call MORPAR
         if( .not. Orbint) goto 400
      endif
c
c-----------------------------------------------------------------------
c tidal friction setup
      if(Km(84).ge.0) Index1 = 6
      if(Km(84).gt.0) Index2 = 6
c
c general relativity setup
      if(Km(61).ge.0) then
         Index1 = 6
         Relmon = Relfct
         if(Km(97).gt.0) Relmon = Mcond(30)
         Jndex1 = 6
         if(Km(61).gt.0) Index2 = 6
         do i=1,9
            Kgr(i) = Km(i+30)
            Xmass(i) = Mass(i)
         end do
         Kgr(3) = 1
         Kgr(10) = 1
         Kgr(11) = 1
         Xmass(3) = Merth
         Xmass(10)= 1._10
         Xmass(11)= Mmoon
         Cvel2  = Cvel**2
         B2g2   = 2._10*(Betapm + Gamapm)
         B2g21  = B2g2 + 1._10
         B2m1   = 2._10*Betapm - 1._10
         G21    = 2._10*Gamapm + 1._10
         G22    = G21 + 1._10
         G11    = Gamapm + 1._10
         G7     = (4._10*Gamapm + 3._10)/2._10
         A44    = (Gauss*Aultsc/86400._10)**2
      endif
c
c print out data and control constants for moon integration
  400 if(Knbd.eq.0) call MORHED
      call PAGCHK(60, 12, 0)
      write(Iout, 500) mass3i, emmr, Aultsc, Mcond(7), beta, gamma,
     .                 mnz2
  500 format(/'-FUNDAMENTAL QUANTITIES FOR MOON INTEGRATION:'/
     .       ' MASS(SUN/E+M)= ', f17.10, 7x, 'MASS(E/M)= ', f17.14 /
     .       ' AU= ', f17.13, ' LT SEC', 7x, 'MRAD= ', f17.12,
     .       ' KM       BETA=', 1pd22.15/ ' GAMMA=', 1pd22.15,
     .       '      MNZ2=', 1pd22.15)
      write(Iout, 600) aukm, gme, gmm, cmr2, mnc22
  600 format('-DERIVED QUANTITIES:'/ ' AU=', 1pd22.15,
     .   ' KM   GM(EARTH)=', 0pf18.10, 4x, ' KM**3/SEC**2   GM(MOON)=',
     .   f18.12, 4x, ' KM**3/SEC**2  '/ ' MOON C/(M*R**2)=',
     .   0pd22.16, ' MNC22 (UNSCALED) =', 1pd22.15)
c
c print set up time
      if(Knbd.eq.0) call TIMRIT(
     . '  SETUP FOR MOON ORBIT AND ROTATION INTEGRATION ',12)
 
      return
  170 call SUICID('INVALID PARAMETER FOR MOONROT VAREQNS:MORSET',11)
  180 call SUICID('INVALID PARAMETER FOR MOON VAR.EQNS: MORSET ',11)
      end
