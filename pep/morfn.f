      real*10 function MORFN(k,j,s)

      implicit none

c r.king and r.cappallo   july 1977   real*10 function morfn
c evaluation of right side of differential equations for moon orbit
c and rotation.  this routine based on monfn, written by m.ash and
c m.slade oct 68, and mrtfn, written by king and cappallo may 76.
c
c arguments
      integer k,j
      real*10 s
c           k  =   equation number
c           j  =   iteration number
c           s  =   time (julian date+0.5_10)
c
c array dimensions
      include 'globdefs.inc'

c        common
      include 'ellips.inc'
      include 'emmips.inc'
      include 'empcnd.inc'
      real*10 alpha,beta,gamma,meqinc
      equivalence (Mrcond(8),alpha),(Mrcond(9),beta),
     .            (Mrcond(10),gamma),(Mrcond(11),meqinc)
      include 'harmor.inc'
      real*10 mmabc(3)
      equivalence (mmabc,Mma)
      include 'inodta.inc'
      include 'intstf.inc'
      include 'metuna.inc'
      include 'mnrtlb.inc'
      real*10 w(3)
      equivalence (w,W1)
      include 'morcrd.inc'
c
c     ecor  = earth relative to sun (re is distance) = rvec(.,3,10)
c     mcor  = moon relative to sun (rm is distance) = rvec(.,11,10)
c     mecor = moon relative to earth (rem is distance) = rvec(.,11,3)
c     mcor etc. are equatorial coordinates
c     smcor etc. are selenodetic coordinates
c
      include 'morstf.inc'
      real*10 dhdxes(3,3,2),h(3,2)
      equivalence (dhdxes(1,1,1),Dhedx),(h(1,1),He)
      include 'orblun.inc'
      include 'output.inc'
      include 'param.inc'
      include 'precmn.inc'
      include 'prtcod.inc'
      include 'xprcom.inc'
      include 'yvectplp.inc'
      include 'tapdtplp.inc'

c solar radiation pressure model
c include effect of absorption and reradiation
c assume negligible thermal inertia
c solar luminosity is given in equivalent solar mass per year
c possible refinements: thermal inertia for earth
c some day, the solar wind may be nonnegligible
      real*10 ertalb,monalb,sollum,twothirds
      parameter (ertalb=0.39_10,monalb=0.11_10,sollum=6.6E-14_10,
     . twothirds=0.666666666666666666_10)
      real*10 ertradp(3),monradp(3),radp(3),radp0,radp0e,radp0m
c
c quantities internal to this routine
      real*10 rpe2(9),rpe3(9),
     .       rpe5(9),rpm2(9),rpm3(9),rpm5(9),plf(3,9)
      real*10 hxe(3),hxs(3),pterms(9),mmterm
      equivalence (pterms,Ppsi2)
      real*10 dydp(6),dycdp(6),temp(3),temq(3),n0(3),fnr(3),dhdc(3,3),
     . dcordc(3,3),corpsi(3),corthe(3),corphi(3),ptrmq(3,3),
     . dfdy(6,3),dfdp(3),n0c(3),fnc(3),dwc(3),cmtrq(3),
     . mctrq(3),dfdyc(6,3),dfcdy(6,3),dfcdyc(6,3)
      equivalence (pterms,ptrmq(1,1)),(dcordc(1,1),corpsi),
     .            (dcordc(1,2),corthe),(dcordc(1,3),corphi)
      real*10 cor(3),dadx(3,3),dadxj2s(3,3),dadx3(3,3),dadxp(3,3,9),
     . daij,deldp,delsp,dfdx(6,3),dndx(3,3),
     . dphi,dphic,dpsy,dpsyc,dtheta,dthetac,dw(3),dw1,dw2,dw3,
     . gama1,gdpfc2,gdpfct,gdpomg(3),gdppar(3),gdpprv(3),
     . gfacte,gfacte2,gfacte5,gfactm,gfactm2,gfactm5,masfct,
     . pepe,pmpm,psum,psumh,st0,sum,sume(3),sumhe(3),sumhm(3),
     . summ(3),sump(3),sums(3),tde,tdm,
     . tdmp,tdep,termc,termre,termrm,vary,x(3),xcr(3,3)
      equivalence (dw1,dw),(dw2,dw(2)),(dw3,dw(3))
      equivalence (dadx3,dadxp(1,1,3))
      integer i,icm,icmkkk,ip1,ip2,is,itg,jj,kkk1,kt,l,l1,n1

c external functions
      real*10 DOT,DOTN

      real*10 DBL10
      real*4 x4
      DBL10(x4)=x4

      Fn(1) = 0._10
      Fn(2) = 0._10
      if(Kkm(60).gt.1)
     . write(6,99333) k,j,s-2440000.5_10
c
      if(k-Knbd.eq.10 .and. k.lt.Korb) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c  computations for motion partials done only once for a given
c  iteration of a given step
c
c determine if second harmonic of sun is included in partials
         if(Km(63).ge.1) then
            gfactm5 = 15._10*gfactm/Rm2
            gfacte5 = 15._10*gfacte/Re2
            do i=1,3
               termrm = Masse*(1._10+Mmoon)*(-5._10*Mcor(i)*gfactm2/Rm2
     .          + gfactm5*(C3(i)-gfactm*Mcor(i)/Rm))/Rm4
               termre = Mass(10)*(1._10+Merth)*(-5._10*Ecor(i)*gfacte2/
     .          Re2 + gfacte5*(C3(i)-gfacte*Ecor(i)/Re))/Re4
               termc = Masse*(1._10+Mmoon)*(gfactm5*Mcor(i)
     .          - 3._10*C3(i)/Rm)/Rm4
     .          + Mass(10)*(1._10+Merth)*(gfacte5*Ecor(i)
     .          - 3._10*C3(i)/Re)/Re4
               do jj=1,3
                  dadxj2s(jj,i)=termrm*Mcor(jj)+termre*Ecor(jj)+
     .             termc*C3(jj)
                  if(i.eq.jj) dadxj2s(jj,i)=dadxj2s(jj,i)
     .             + Masse*(1._10+Mmoon)*gfactm2/Rm4
     .             + Mass(10)*(1._10+Merth)*gfacte2/Re4
               end do
            end do
         endif

c determine partial derivatives of motion relative to central body
         n1 = Iparm - 1
         if(Orbint) then
            do l1 = 1,n1
               l = l1*6+Knbd
               do i = 1,Index2
                  l = l + 1
                  Dmcor(i,l1) = Y(l,j)
                  if(Icmtrl(l1).le.-31) then
                     if(Km(100).lt.0) then
                     else if(Km(100).eq.0) then
                        icm = -30 - Icmtrl(l1)
                        Dmcor(i,l1) = Dmcor(i,l1) + Dylpt(i,icm)
                     else
                        icm = -30 - Icmtrl(l1)
                        Dmcor(i,l1) = Dmcor(i,l1) + Dylun(i,icm)
                     endif
                  endif
               end do
            end do
            if(Km(61).gt.0) call MONREL(-2)
c
c compute gravity gradient wrt moon coordinates relative to earth (dadx)
c and wrt embary coordinates relative to sun (dadx3, part of dadxp)
c and wrt planet coordinates relative to sun (dadxp)
            do i=1,3
               do jj=1,3
                  dadx(i,jj)=0._10
                  dadx3(i,jj)=0._10
               end do
            end do
c contribution of planets (not earth) to gravity gradient
            do l=1,9
               if(Km(l+30).gt.0 .and. l.ne.3) then
                  tdmp=Gamat*Mass(l)/rpm3(l)
                  tdep=Gamat*Mass(l)/rpe3(l)
                  tdm=tdmp*Masse
                  tde=tdep*Mass(10)
                  dadxp(2,1,l)=0._10
                  dadxp(3,1,l)=0._10
                  dadxp(3,2,l)=0._10
                  do i=1,3
                     dadx(i,i)=dadx(i,i) - (tde+tdm)
                     dadx3(i,i)=dadx3(i,i) + tdep - tdmp
                     dadxp(i,i,l)= tdmp-tdep
                     do jj=1,i
                        pepe=pecor(jj,l)*pecor(i,l)/rpe2(l)
                        pmpm=pmcor(jj,l)*pmcor(i,l)/rpm2(l)
                        dadx(i,jj)=dadx(i,jj) + 3._10*(
     .                   tde*pepe+tdm*pmpm)
                        daij=3._10*(tdep*pepe-tdmp*pmpm)
                        dadx3(i,jj)=dadx3(i,jj) - daij
                        dadxp(i,jj,l)=dadxp(i,jj,l) + daij
                     end do
                  end do
c finish symmetric matrix
                  dadxp(1,2,l)=dadxp(2,1,l)
                  dadxp(1,3,l)=dadxp(3,1,l)
                  dadxp(2,3,l)=dadxp(3,2,l)
               endif
            end do

c contribution of sun to gravity gradient
c contrary to description of k(30+ncentr) in general, sun is included in
c equations for lunar partials whenever it's included in equations of
c motion
            if(Km(33).ge.0) then
               tdm=Gamat*Masse/Rm3
               tde=Gamat*Mass(10)/Re3
               tdmp=Gamat/Rm3
               tdep=Gamat/Re3
               do i=1,3
                  do jj=1,i
                     pepe=Ecor(jj)*Ecor(i)/Re2
                     pmpm=Mcor(jj)*Mcor(i)/Rm2
                     dadx(i,jj)=dadx(i,jj) + 3._10*(
     .                tde*pepe+tdm*pmpm)
                     dadx3(i,jj)=dadx3(i,jj) + 3._10*(
     .                -tdep*pepe+tdmp*pmpm)
                  end do
                  dadx(i,i)=dadx(i,i) - (tde+tdm)
                  dadx3(i,i)=dadx3(i,i) + (tdep-tdmp)
               end do
            endif

c contribution of earth to gravity gradient
            tde=Gamtem/Rem3
            do i=1,3
               do jj=1,i
                  dadx(i,jj)=dadx(i,jj)+3._10*tde*Mecor(jj)
     .             *Mecor(i)/Rem2
               end do
               dadx(i,i)=dadx(i,i) - tde
            end do

c finish symmetric matrices
            dadx(1,2)=dadx(2,1)
            dadx(1,3)=dadx(3,1)
            dadx(2,3)=dadx(3,2)
            dadx3(1,2)=dadx3(2,1)
            dadx3(1,3)=dadx3(3,1)
            dadx3(2,3)=dadx3(3,2)

c save for extra printout
            if(Kkm(7).gt.0) then
               do i=1,3
                  Pcorsav(i,1)=Mecor(i)
                  Pcorsav(i+3,1)=Mecor(i+3)
                  Pcorsav(i,3)=Xpert(i,3)
                  Pcorsav(i+3,3)=Xpert(i+3,3)
                  Parsav(i,1)=Dmcor(i,1)
                  Parsav(i+3,1)=Dmcor(i+3,1)
                  Indpsav(i,1,1)=DOT(Dmcor(1,1),dadx(1,i))
                  if(Parnum(3).gt.0) then
                     Parsav(i,3)=Dccor(i,Parnum(3))
                     Parsav(i+3,3)=Dccor(i+3,Parnum(3))
                     Indpsav(i,3,1)=DOT(Dccor(1,Parnum(3)),dadx3(1,i))
                  endif
                  do jj=1,3
                     Dadxsav(i,jj,1,1,1)=dadx(i,jj)
                     Dadxsav(i,jj,1,3,1)=dadx3(i,jj)
                  end do
                  if(Numtar.gt.0) then
                     itg=Nplpt(1)
                     Pcorsav(i,2)=Xpert(i,itg)
                     Pcorsav(i+3,2)=Xpert(i+3,itg)
                     if(Parnum(2).gt.0) then
                        Parsav(i,2)=Dtcor(i,Parnum(2),1)
                        Parsav(i+3,2)=Dtcor(i+3,Parnum(2),1)
                        Indpsav(i,2,1)=DOT(Dtcor(1,Parnum(2),1),
     .                   dadxp(1,i,itg))
                     endif
                     do jj=1,3
                        Dadxsav(i,jj,1,2,1)=dadxp(i,jj,itg)
                     end do
                  endif
               end do
            endif

         endif
         goto 300
      endif

      if(k-Knbd.ne.4) goto 300
      if(j.ne.2 .or. Knbd.gt.0) then
c
c-----------------------------------------------------------------------
c
c  computations done only once for a given step
c     note that computations depending on planets need to be redone
c     each iteration when the moon is included in the n-body list
c
c update time-dependent gravitation constant
c (truly only once per step)
         if(Km(62).ge.0 .and. j.ne.2) then
            Tvary  = s - Tvary0
            vary   = (1._10 + Gmvary*Tvary)
            Gamat  = Gama*vary
            Gamat3 = Gama3*vary
            Gamtem = Gamem*vary
         endif
c
c determine perturbing planet coordinates
         if(Knbd.eq.0) then
            Jd    = s
            Fract = Jd
            Fract = s - Fract
            call PRTCRD(Jd,Fract)
            do kt=0,Numtar
               call PLPCRDT(Jd,Fract,kt)
               if(Kkm(60).eq.-2) write(6,99112) kt,(Dtcor(i,7,kt),i=1,6)
99112          format('D',i1,'cor7=',1p3d20.12/3d20.12)
            end do
ccccc special insert to test moon code in nbody environment with
ccccc interpolated planets instead of simultaneously integrated
         else if(Ipert.gt.0 .and. Mcon1(9).eq.9._10) then
            call PRTCRD(Jd,Fract)
         endif
         if(Kkm(60).eq.-3) write(6,99113) s,((Xpert(i,jj),i=1,6),
     .    jj=1,9)
99113    format('s=',f13.4,' Xpert='/(1p3d22.15))

         if(Orbint) then

c determine general relativity quantities
            if(Km(61).ge.0) call MONREL(0)

c calculate expected geodesic precession rate vector
            gdpfct    = -(Gamapm + 0.5_10)*Gamat/Cvel**2
            Rpert2(3) = DOT(Xpert(1,3),Xpert(1,3))
            Rpert(3)  = SQRT(Rpert2(3))
            Rpert3(3) = Rpert2(3)*Rpert(3)
            do i = 1,3
               Xpert3(i,3) = Xpert(i,3)/Rpert3(3)
               gdppar(i)    = Xpert3(i,3)*gdpfct
            end do
            call CROSS(Xpert(4,3),gdppar,gdpomg)
            gdpfc2 = 3._10*DOT(Xpert(4,3),Xpert(1,3))/Rpert2(3)
c calculate partials of em barycenter coordinates wrt
c embarycent initial conditions and gmvary
c setup for moon cartesian coordinates (already done if part of n-body)
            if(Knbd.eq.0) then
               if(Km(100).lt.0) then
               else if(Km(100).eq.0) then
c
c determine elliptic orbit coordinates for encke's method
                  Tlpt = s - Tlpt0
                  call ELIPT(Km1,Tlpt)
               else
c
c determine brown mean lunar orbit coordinates
                  call LUNORB(Jd,Fract,Km1)
                  goto 20
               endif
            endif
c
c if rotation not integrated and brown mean lunar orbit not
c used, it is necessary to calculate aa matrix if effect of
c moon harmonics included
            if(.not. (Kmh.lt.0 .or. Rotint))
     .       call ECLPRC(Jd,Fract,0)
         endif
c
c determine precession nutation matrix if effect of
c earth harmonics or tidal friction included
c (truly only once per step)
   20    if((Keh.ge.0 .or. Km(84).ge.0) .and. j.ne.2) then
            call PRCES(s - 0.5_10)
            call PRENUT
         endif
c
c if orbit only, determine lunar rotation matrix from
c libration quantities on perturbing planet data set
c (truly only once per step)
         if((.not.(Rotint.or.Kmh.lt.0)) .and. j.ne.2) then
            if(Librat(1,2).ge.0.1) then

c euler angles on p-p tape
               call MONROT(0,DBL10(Librat(1,1)),DBL10(Librat(1,2)),
     .          DBL10(Librat(1,3)),0._10,0._10,0._10)
            else
               call MONLIB(0,0,0)
            endif
         endif
      endif

c ------------------------------------------------------------
c
c computations done only once for a given iteration of a given step
c
c  determine lunar rotation matrix from integrated euler angles
      if(Rotint) then
         call MONROT(0,Y(Korb+1,j),Y(Korb+2,j),Y(Korb+3,j),
     .    Y(Korb+4,j),Y(Korb+5,j),Y(Korb+6,j))
         dpsy   = Y(Korb+4,j)
         dtheta = Y(Korb+5,j)
         dphi   = Y(Korb+6,j)
c
c
c determine angular velocities in body-fixed frame
         W1 = dtheta*Cphi + dpsy*Sphi*Stheta
         W2 = -dtheta*Sphi + dpsy*Cphi*Stheta
         W3 = dpsy*Ctheta + dphi
c
c call monrot to print out classical libration elements
         if(Kout.ne.0) then
            if(Fract.eq.0._10) then
               if(j.ne.1) then
                  call ECLPRC(Jd,Fract,1)

c store jd temporarily in /precmn/ for monrot printout
                  Tpr = s - .5_10
                  call MONROT(1,Y(Korb+1,j),Y(Korb+2,j),Y(Korb+3,j),
     .             Y(Korb+4,j),Y(Korb+5,j),Y(Korb+6,j))
               endif
            endif
         endif
         if(Corint) then
            call MNCROT(Y(Krot+1,j),Y(Krot+2,j),Y(Krot+3,j),
     .       Y(Krot+4,j),Y(Krot+5,j),Y(Krot+6,j))
            dpsyc  = Y(Krot+4,j)
            dthetac= Y(Krot+5,j)
            dphic  = Y(Krot+6,j)

c get rotation matrix relating core-fixed to body-fixed
c multiply a body-fixed vector by mcmrot to get a core-fixed vector
            call PRODCT(Mcrtlb,Mrotlb,Mcmrot,3,-3,3)
c
c determine angular velocities in core-fixed frame
            Wc(1) = dthetac*Cphic + dpsyc*Sphic*Sthetac
            Wc(2) = -dthetac*Sphic + dpsyc*Cphic*Sthetac
            Wc(3) = dpsyc*Cthetac + dphic
            if(Kkm(61).gt.0.or.ABS(dpsyc).gt.1e-2_10)
     .       write(6,44401) s,j,Wc,dpsyc,dthetac,dphic,Mcmrot
44401       format(' T,j,Wc,dth,rot=',f16.6,i3,2x1p3e24.17/(3e20.13))
         endif
      endif
c
c determine coords of earth and sun in reference system
c
c 1-3 for newtonian, but need 1-6 for relativity
      if(.not. Orbint) then
         do i = 1,Index1
            Mecor(i) = Xpert(i,10)
            Rvec(i,3,11)= -Mecor(i)
         end do
      else
         do i = 1,Index1
            if(Km(100).lt.0) then
               Mecor(i) = Y(Knbd+i,j)
            else if(Km(100).eq.0) then
               Mecor(i) = Y(Knbd+i,j) + Ylpt(i)
            else
               Mecor(i) = Y(Knbd+i,j) + Ylun(i)
            endif
            Rvec(i,3,11)= -Mecor(i)
         end do
         call CROSS(gdpomg,Mecor,gdppar)
         call CROSS(gdpomg,Mecor(4),gdpprv)
      endif
      Rem2 = DOT(Mecor,Mecor)
      Rem  = SQRT(Rem2)
      Rab(3,11)= Rem
      Rem3 = Rem2*Rem
      Rem4 = Rem2*Rem2
      Rem5 = Rem2*Rem3
      if(Dosun.ge.0) then
         do i = 1,Index1
            if(i.le.Jndex1) then
               Ecor(i) = Xpert(i,3) - Mass(10)*Mecor(i)
               Mcor(i) = Xpert(i,3) + Masse*Mecor(i)
               Rvec(i,10,3) = -Ecor(i)
               Rvec(i,10,11)= -Mcor(i)
            endif
         end do
         Re2  = DOT(Ecor,Ecor)
         Re   = SQRT(Re2)
         Rab(10,3)= Re
         Re3  = Re2*Re
         Re4  = Re2*Re2
         Re5  = Re2*Re3
         Rm2  = DOT(Mcor,Mcor)
         Rm   = SQRT(Rm2)
         Rab(10,11) = Rm
         Rm3  = Rm2*Rm
         Rm4  = Rm2*Rm2
         Rm5  = Rm2*Rm3
      endif
c
c determine perturbing planets relative to earth and moon
      if(Orbint) then
         do l = 1,9
            if(l.ne.3) then
               if(Km(l+30).ge.0) then
                  do i = 1,3
                     pecor(i,l) = Xpert(i,l) - Ecor(i)
                     Rvec(i,3,l)= -pecor(i,l)
                     pmcor(i,l) = Xpert(i,l) - Mcor(i)
                     Rvec(i,11,l)= -pmcor(i,l)
                  end do
                  rpe2(l) = DOT(pecor(1,l),pecor(1,l))
                  rpe(l)  = SQRT(rpe2(l))
                  Rab(3,l)= rpe(l)
                  rpe3(l) = rpe2(l)*rpe(l)
                  rpe5(l) = rpe2(l)*rpe3(l)
                  rpm2(l) = DOT(pmcor(1,l),pmcor(1,l))
                  rpm(l)  = SQRT(rpm2(l))
                  Rab(11,l)= rpm(l)
                  rpm3(l) = rpm2(l)*rpm(l)
                  rpm5(l) = rpm2(l)*rpm3(l)
               endif
            endif
         end do
      endif
c
c determine selenodetic coordinates
      if(Rotint .or. Kmh.gt.0) then
         call PRODCT(Mrotlb,Mecor,Smecor,3,3,1)
         call PRODCT(Mrotlb,Mcor,Smcor,3,3,1)
         do i = 1,3
            Smecor(i) = -Smecor(i)
            Smcor(i)  = -Smcor(i)
         end do
         Srm2  = DOT(Smcor,Smcor)
         Srm   = SQRT(Srm2)
         Srem2 = DOT(Smecor,Smecor)
         Srem  = SQRT(Srem2)
      endif
c
c determine tidal friction quantities
      if(Orbint) then
         if(Km(84).ge.0) call MONTID(-1,s)
c
c determine general relativity quantities
         if(Km(61).ge.0) call MONREL(-1)
      endif
c
c determine quantities for harmonics of earth and moon
      call MORHAR(-1,j)
c
c determine quantities for solar j2
      if(Orbint .and. Km(63).ge.0) then
         gfacte  = DOT(C3,Ecor)/Re
         gfacte2 = (7.5_10*gfacte**2 - 1.5_10)/Re
         gfactm  = DOT(C3,Mcor)/Rm
         gfactm2 = (7.5_10*gfactm**2 - 1.5_10)/Rm
      endif

c set up for solar radiation pressure
      if(Orbint .and. Km(86).gt.0) then
         radp0=sollum*Cvel*0.25_10/365.25_10
         radp0m=radp0*Mrad**2/Rm3/Mmoon
         radp0e=radp0*Erad**2/Re3/(Masse*Mass(3))
         if(s.gt.2440005._10 .and. s.lt.2440007._10)
     .    write(6,99939) s,k,j,radp0,radp0m,radp0e
99939    format(f10.2,2i3,1p3d20.12)
         do i=1,3
            ertradp(i)=radp0e*(1._10 + twothirds*ertalb
     .       + (1._10-ertalb)*twothirds)*Ecor(i)
            monradp(i)=radp0m*(1._10 + (0.7_10*twothirds+0.3_10)*monalb
     .       + (1._10-monalb)*twothirds)*Mcor(i)
            radp(i)=monradp(i)-ertradp(i)
         end do
      endif
c
c determine additional quantities for right side of
c equations for rotation
      if(.not. Rotint) goto 300
      do i = 1,9
         pterms(i) = 0._10
      end do
      do i=1,3
         do l=1,6
            Ddwdx(i,l)=0._10
         end do
      end do

c partials of selenodetic transform matrix wrt euler angles
c Drtdy is transpose for convenience: Drtdy(i,j,k)=dMrotlb(j,i)/dy(k)
      if(Iparmr.gt.1) then
         call DMONRT(Drtdy,Dwdy,Spsi,Cpsi,Stheta,Ctheta,Sphi,Cphi,
     .    dpsy,dtheta,dphi)
         if(Corint)
     .    call DMONRT(Drtdyc,Dwcdyc,Spsic,Cpsic,Sthetac,Cthetac,
     .    Sphic,Cphic,dpsyc,dthetac,dphic)

      endif
c
c the following calculations done twice, once for the torque due
c to the earth (is=1), once for the torque due to the sun (is=2)
c
      is    = 1
      gama1 = Gamat3
      masfct= gama1
      do i = 1,3
         cor(i) = Smecor(i)
         x(i)   = -Mecor(i)
c        h(i)   = He(i)
      end do
c
c form harmonic cross product terms
      call CROSS(cor,h,hxe)
      do while( .true. )
c
c obtain partials of selenodetic coordinates wrt euler angles
         if(Iparmr.gt.1) then
            if(is.eq.2 .and. Kmr(81).lt.1) goto 200
            call PRODCT(Drtdy,x,dcordc,-9,3,1)
c save partials of earth coordinates for other calculations
            if(is.eq.1) then
               do l=1,3
                  do i=1,3
                     Dsxdy(i,l)=dcordc(i,l)
                  end do
               end do
            endif
c
c obtain partials of harmonic terms wrt euler angles
            call PRODCT(dhdxes(1,1,is),dcordc,dhdc,3,3,3)
c
c calculate harmonic terms for partials of rhs wrt euler angles
c ptrmq: (.,1) = ppsi2,4,5  (.,2)=pthe2,4,5  (.,3)=pphi2,4,5
            do l = 1,3
               call CROSS(dcordc(1,l),h(1,is),temp)
               call CROSS(cor,dhdc(1,l),temq)
               do i = 1,3
                  ptrmq(i,l) = ptrmq(i,l) + gama1*mmabc(i)
     .                          *(temp(i) + temq(i))
               end do
            end do

c obtain partials of rhs wrt orbit state
            if(Iparm.gt.1) then
c cross product is equivalent to multiplication by a matrix
               xcr(1,1)=0._10
               xcr(1,2)=+cor(3)
               xcr(1,3)=-cor(2)
               xcr(2,1)=-cor(3)
               xcr(2,2)=0._10
               xcr(2,3)=+cor(1)
               xcr(3,1)=+cor(2)
               xcr(3,2)=-cor(1)
               xcr(3,3)=0._10

               call PRODCT(xcr,dhdxes(1,1,is),dndx,3,3,3)
               dndx(1,2)=dndx(1,2)-h(3,is)
               dndx(1,3)=dndx(1,3)+h(2,is)
               dndx(2,1)=dndx(2,1)+h(3,is)
               dndx(2,3)=dndx(2,3)-h(1,is)
               dndx(3,1)=dndx(3,1)-h(2,is)
               dndx(3,2)=dndx(3,2)+h(1,is)
               do i=1,3
                  do l=1,3
                     Ddwdx(i,l)=Ddwdx(i,l)+dndx(i,l)*masfct*mmabc(i)
                  end do
               end do
            endif
         endif
         if(is.ne.2 .and. Kmr(81).ge.0) then
            is = 2
c
c enter loop for is=2 here
            gama1 = Gamat
            masfct= gama1*Masse
            do i = 1,3
               cor(i) = Smcor(i)
               x(i)   = -Mcor(i)
c               h(i)   = Hs(i)
            end do
            call CROSS(cor,h(1,is),hxs)
            goto 100
         endif
         if(Iparm.gt.1) then
c transform partials to inertial coordinates and zero out the velocity
c dependence (rigid-body approximation)

            call PRODCT(Ddwdx,Mrotlb,dndx,3,3,3)
            do i=1,3
               do l=1,3
                  Ddwdx(i,l)=dndx(i,l)
                  Ddwdx(i,l+3)=0._10
               end do
            end do
         endif
         goto 200
  100 end do
c
c determine rigid body torques
  200 do i     = 1,3
         n0(i) = Gamat3*Mmoon*hxe(i)

c include solar torque
         if(Kmr(81).ge.0) n0(i) = n0(i) + Gamat*Mmoon*hxs(i)
      end do
      if(Kkm(60).gt.2) write(6,99334)
99334 format('calling MORFFI')

c go get torque contribution from the figure-figure interaction
      if(Kmr(82).ge.0) call MORFFI(n0)

c get core-mantle torque
      if(Corint) then
         call PRODCT(Mcmrot,Wc,temp,-3,3,1)
         do i=1,3
            cmtrq(i)=temp(i)-w(i)
            n0(i)=n0(i)+Mrcond(23)*cmtrq(i)
         end do
         if(Kkm(61).gt.0.or.ABS(dpsyc).gt.1e-2_10)
     .    write(6,44402) temp,w,cmtrq,n0
44402    format('temp,w,cmtrq,n0='/(1p3e24.17))
      endif

c get kinematical cross terms
      call PRODCT(I0,W1,temp,3,3,1)
      call CROSS(W1,temp,temq)
      do i = 1,3
         N0mkin(i) = n0(i) - temq(i)
      end do

c rigid body angular acceleration vector put in dw
      call PRODCT(I0i,N0mkin,dw,3,3,1)

c include effects of lunar elasticity and dissipation if wanted
      if(Kmr(83).ge.0) then
         if(Kkm(60).gt.2) write(6,99335) 'calling MORED',dw
99335    format(a,' dw=',1p3e22.15)
         call MORED(s,dw,Kmr(83))
         if(Kkm(60).gt.2) write(6,99335) 'returned',dw
      endif
c
c obtain terms for rhs of equations of motion
      fnr(1) = (dw1*Sphi + dw2*Cphi + dtheta*dphi - dpsy*dtheta*Ctheta)
     .         /Stheta
      fnr(2) = dw1*Cphi - dw2*Sphi - dpsy*dphi*Stheta
      fnr(3) = dw3 - fnr(1)*Ctheta + dpsy*dtheta*Stheta
c
c obtain partials of rhs wrt euler angles and rates
      if(Iparmr.gt.1) then
c
c first get partials of angular accelerations
c ptrmq already has the partials of the harmonic terms of torque
c and is equivalenced with the first half of the array.
c Clear the second half of the array and add the partials of the
c kinematic terms. For historical reasons,
c the angle partials are denoted 2,4,5, and the angular velocity
c partials 1,3,6

         do i=1,3
            mmterm=Mirat(i)
            if(i.eq.2) mmterm=-mmterm
            ip1=i+1
            if(ip1.gt.3) ip1=1
            ip2=ip1+1
            if(ip2.gt.3) ip2=1
            do jj=4,6
               ddwdy(i,jj)=0._10
            end do
            do jj=1,6
               ddwdy(i,jj)=ddwdy(i,jj)-mmterm*(w(ip1)*Dwdy(ip2,jj)
     .          +w(ip2)*Dwdy(ip1,jj))
            end do
         end do

c include partials of core-mantle torque
         if(Corint) then
c first get dependence on angles, filling the first half of the array
c use the array for core partials for temporary storage

            call PRODCT(Mcrtlb,Wc,temp,-3,3,1)
            call PRODCT(Drtdy,temp,Ddwdyc,-9,3,1)

            do i=1,3
               do jj=4,6
                  Ddwdyc(i,jj)=0._10
               end do
               do jj=1,6
                  ddwdy(i,jj)=ddwdy(i,jj) + (Ddwdyc(i,jj)-Dwdy(i,jj))
     .             *Mrcond(23)*I0i(i,i)
               end do
            end do
         endif

c correct partials of angular acceleration w.r.t. state
      if(Kkm(60).gt.2) write(6,99336)
99336 format('calling MORFFI2')
         if(Kmr(82).ge.0) call MORFFI2(Iparm)
      if(Kkm(60).gt.2) write(6,99337)
99337 format('calling MORED2')
         if(Kmr(83).ge.0) call MORED2(s,Kmr(83),Iparm)
         
c dfdy(i,j) means df(j)/dy(i)
         dfdy(1,1) = (Sphi*Ppsi2 + Cphi*Ppsi4)/Stheta
         dfdy(1,2) = Cphi*Ppsi2 - Sphi*Ppsi4
         dfdy(1,3) = Ppsi5 - dfdy(1,1)*Ctheta
         dfdy(2,1) = (-Ctheta/Stheta*(dw1*Sphi+dw2*Cphi+dtheta*dphi)
     .               + Sphi*Pthe2 + Cphi*Pthe4
     .               + dpsy*dtheta/Stheta)/Stheta
         dfdy(2,2) = Cphi*Pthe2 - Sphi*Pthe4 - dpsy*dphi*Ctheta
         dfdy(2,3) = Pthe5 + fnr(1)*Stheta
     .               + (dpsy*dtheta - dfdy(2,1))*Ctheta
         dfdy(3,1) = (Sphi*(-dw2+Pphi2) + Cphi*(dw1+Pphi4))/Stheta
         dfdy(3,2) = Cphi*(-dw2+Pphi2) - Sphi*(dw1+Pphi4)
         dfdy(3,3) = Pphi5 - dfdy(3,1)*Ctheta
c
c obtain partials wrt time derivatives of euler angles
         dfdy(4,1) = (Sphi*Ppsi1 + Cphi*Ppsi3 - dtheta*Ctheta)/Stheta
         dfdy(4,2) = Cphi*Ppsi1 - Sphi*Ppsi3 - dphi*Stheta
         dfdy(4,3) = Ppsi6 - Ctheta*dfdy(4,1) + dtheta*Stheta
         dfdy(5,1) = (Sphi*Pthe1 + Cphi*Pthe3 + dphi - dpsy*Ctheta)
     .               /Stheta
         dfdy(5,2) = Cphi*Pthe1 - Sphi*Pthe3
         dfdy(5,3) = Pthe6 - dfdy(5,1)*Ctheta + dpsy*Stheta
         dfdy(6,1) = (Sphi*Pphi1 + Cphi*Pphi3 + dtheta)/Stheta
         dfdy(6,2) = Cphi*Pphi1 - Sphi*Pphi3 - dpsy*Stheta
         dfdy(6,3) = Pphi6 - Ctheta*dfdy(6,1)

         if(Corint) then
c also need partials of angular acceleration wrt core state

            call PRODCT(Mcmrot,Dwcdyc,Ddwdyc,-3,3,6)
            do jj=1,3
               call PRODCT(Drtdyc(1,1,jj),Wc,temp,3,3,1)
               call PRODCT(Mrotlb,temp,temq,3,3,1)
               do i=1,3
                  Ddwdyc(i,jj)=Ddwdyc(i,jj)+temq(i)
               end do
            end do

            do jj=1,6
               do i=1,3
                  Ddwdyc(i,jj)=Ddwdyc(i,jj)*Mrcond(23)*I0i(i,i)
               end do
               dfdyc(jj,1) = (Sphi*Ddwdyc(1,jj) + Cphi*Ddwdyc(2,jj))
     .          /Stheta
               dfdyc(jj,2) = Cphi*Ddwdyc(1,jj) - Sphi*Ddwdyc(2,jj)
               dfdyc(jj,3) = Ddwdyc(3,jj) - dfdyc(jj,1)*Ctheta
            end do
         endif

         if(Iparm.gt.1) then
c partials of rhs wrt orbit
            do i=1,6
               dfdx(i,1) = (Sphi*Ddwdx(1,i) + Cphi*Ddwdx(2,i))/Stheta
               dfdx(i,2) = Cphi*Ddwdx(1,i) - Sphi*Ddwdx(2,i)
               dfdx(i,3) = Ddwdx(3,i) - dfdx(i,1)*Ctheta
            end do
         endif
      endif

c quantities for right-hand side of equations of motion of lunar core
      if(.not.Corint) goto 300

c get mantle-core torque
      call PRODCT(Mcmrot,W,temp,3,3,1)
      do i=1,3
         mctrq(i)=temp(i)-Wc(i)
         n0c(i)=Mrcond(23)*mctrq(i)
      end do

c get kinematical cross terms
c (but these are zero in the case of a spherically symmetric core)
c      call PRODCT(I0c,Wc,temp,3,3,1)
c      call CROSS(Wc,temp,temq)
c      do i = 1,3
c         n0c(i) = n0c(i) - temq(i)
c      end do

c core angular acceleration vector put in dwc
      call PRODCT(I0ci,n0c,dwc,3,3,1)
      if(Kkm(61).gt.0.or.ABS(dpsyc).gt.1e-2_10)
     . write(6,44403) temp,Wc,mctrq
44403 format('temp,Wc,mctrq'/(1p3e24.17))

c obtain terms for rhs of equations of motion
      fnc(1) = (dwc(1)*Sphic + dwc(2)*Cphic + dthetac*dphic -
     . dpsyc*dthetac*Cthetac)/Sthetac
      fnc(2) = dwc(1)*Cphic - dwc(2)*Sphic - dpsyc*dphic*Sthetac
      fnc(3) = dwc(3) - fnc(1)*Cthetac + dpsyc*dthetac*Sthetac
      if(Kkm(61).gt.0.or.ABS(dpsyc).gt.1e-2_10) write(6,
     . '(''n0c,dwc,fnc,s/cpsi,s/cth,s/cphi=''/(1p3d22.15))')
     . n0c,dwc,fnc,Spsic,Cpsic,Sthetac,Cthetac,Sphic,Cphic

      if(Iparmr.le.1) goto 300

c first get partials of core angular acceleration wrt mantle state

      call PRODCT(Mcmrot,Dwdy,Ddwcdy,3,3,6)
      do jj=1,3
         call PRODCT(Drtdy(1,1,jj),W1,temp,3,3,1)
         call PRODCT(Mcrtlb,temp,temq,3,3,1)
         do i=1,3
            Ddwcdy(i,jj)=Ddwcdy(i,jj)+temq(i)
         end do
      end do

      do jj=1,6
         do i=1,3
            Ddwcdy(i,jj)=Ddwcdy(i,jj)*Mrcond(23)*I0ci(i,i)
         end do
         dfcdy(jj,1) = (Sphic*Ddwcdy(1,jj) + Cphic*Ddwcdy(2,jj))
     .    /Sthetac
         dfcdy(jj,2) = Cphic*Ddwcdy(1,jj) - Sphic*Ddwcdy(2,jj)
         dfcdy(jj,3) = Ddwcdy(3,jj) - dfcdy(jj,1)*Cthetac
      end do

c get partials of core angular acceleration wrt core state
c first get dependence on angles, filling the first half of the array
      call PRODCT(Mrotlb,w,temp,-3,3,1)
      call PRODCT(Drtdyc,temp,Ddwcdyc,-9,3,1)

      do i=1,3
         do jj=4,6
            Ddwcdyc(i,jj)=0._10
         end do
         do jj=1,6
            Ddwcdyc(i,jj)=(Ddwcdyc(i,jj)-Dwcdyc(i,jj))
     .       *Mrcond(23)*I0ci(i,i)
         end do
      end do

c dfcdyc(i,j) means dfc(j)/dyc(i)
      dfcdyc(1,1) = (Sphic*Ddwcdyc(1,1) + Cphic*Ddwcdyc(2,1))/Sthetac
      dfcdyc(1,2) = Cphic*Ddwcdyc(1,1) - Sphic*Ddwcdyc(2,1)
      dfcdyc(1,3) = Ddwcdyc(3,1) - dfcdyc(1,1)*Cthetac
      dfcdyc(2,1) = (-Cthetac/Sthetac*(dwc(1)*Sphic+dwc(2)*Cphic
     .    +dthetac*dphic) + Sphic*Ddwcdyc(1,2) + Cphic*Ddwcdyc(2,2)
     .    + dpsyc*dthetac/Sthetac)/Sthetac
      dfcdyc(2,2) = Cphic*Ddwcdyc(1,2) - Sphic*Ddwcdyc(2,2)
     .    - dpsyc*dphic*Cthetac
      dfcdyc(2,3) = Ddwcdyc(3,2) + fnc(1)*Sthetac
     .    + (dpsyc*dthetac - dfcdyc(2,1))*Cthetac
      dfcdyc(3,1) = (Sphic*(-dwc(2)+Ddwcdyc(1,3))
     .            + Cphic*(dwc(1)+Ddwcdyc(2,3)))/Sthetac
      dfcdyc(3,2) = Cphic*(-dwc(2)+Ddwcdyc(1,3))
     .           - Sphic*(dwc(1)+Ddwcdyc(2,3))
      dfcdyc(3,3) = Ddwcdyc(3,3) - dfcdyc(3,1)*Cthetac
c
c obtain partials wrt time derivatives of euler angles
      dfcdyc(4,1) = (Sphic*Ddwcdyc(1,4) + Cphic*Ddwcdyc(2,4)
     . - dthetac*Cthetac)/Sthetac
      dfcdyc(4,2) = Cphic*Ddwcdyc(1,4)-Sphic*Ddwcdyc(2,4)-dphic*Sthetac
      dfcdyc(4,3) = Ddwcdyc(3,4) - Cthetac*dfcdyc(4,1) + dthetac*Sthetac
      dfcdyc(5,1) = (Sphic*Ddwcdyc(1,5) + Cphic*Ddwcdyc(2,5) + dphic
     .     - dpsyc*Cthetac)/Sthetac
      dfcdyc(5,2) = Cphic*Ddwcdyc(1,5) - Sphic*Ddwcdyc(2,5)
      dfcdyc(5,3) = Ddwcdyc(3,5) - dfcdyc(5,1)*Cthetac + dpsyc*Sthetac
      dfcdyc(6,1) = (Sphic*Ddwcdyc(1,6) + Cphic*Ddwcdyc(2,6)
     .    + dthetac)/Sthetac
      dfcdyc(6,2) = Cphic*Ddwcdyc(1,6)-Sphic*Ddwcdyc(2,6)-dpsyc*Sthetac
      dfcdyc(6,3) = Ddwcdyc(3,6) - Cthetac*dfcdyc(6,1)
c
c-----------------------------------------------------------------------
c
c determine class of equation k
  300 Kkk = (k-Knbd-1)/6
      Kk  = k-Knbd - Kkk*6 - 3
      if(k.gt.Korb) Kkk = (k-Korb-1)/6
      if(k.gt.Krot) Kkk = (k-Krot-1)/6
      if(Kk.le.0) then
c time derivative of coordinate (or coordinate partial) is just the rate
c (or rate partial)
         Fn(1) = Y(k+3,j)
         if(k.le.Korb .and. Kkk.gt.0) then
            if(Icmtrl(Kkk).eq.-7) Fn(1) = Fn(1) - gdppar(Kk+3)
         endif
      else if(k.le.Korb) then
         if(Kkk.gt.0) then
c
c-----------------------------------------------------------------------
c
c           determine right side of equations for partial derivatives
c           of orbit with respect to parameters
c
c effect of sun and perturbing planets on partial derivatives
            icmkkk = Icmtrl(Kkk)
c check if parameter is mass of perturbing planet
            if(icmkkk.ge.1 .and. icmkkk.le.9 .and. icmkkk.ne.3)
     .       Fn(1)=Fn(1)+Gamat*plf(Kk,icmkkk)
c check if parameter is moon mass / embary mass
            if(icmkkk.eq.10) Fn(1)=Fn(1)-DOT(Mecor,dadx3(1,Kk))
c check if parameter's partial is available for embary/target orbits
            do kt=0,Numtar
               kkk1=Kpt(Kkk,kt)-1
               if(kkk1.gt.0) then
                  itg=Nplpt(kt)
                  Fn(1)=Fn(1)+DOT(Dtcor(1,kkk1,kt),dadxp(1,Kk,itg))
               endif
            end do
c check if parameter is gmvary
            if(icmkkk.eq.32) Fn(1) = Fn(1) + sump(Kk)*Gama*Tvary
c
c effects of harmonics of earth and moon on partial derivatives
            if(Kmh.gt.0 .or. Keh.gt.0) call MORHAR(k,j)
c
c effect of second harmonic of sun
            if(Km(63).ge.1) then
               psumh = DOT(dadxj2s(1,Kk),Dmcor(1,Kkk))
               Fn(1) = Fn(1) + Gamat*Shar2*psumh
c
c check if alpha is second harmonic of sun
               if(icmkkk.eq.33) then
                  Fn(1) = Fn(1) +Sunrd2*Gamat*((1._10+Mmoon)*sumhm(Kk)-
     .             (1._10+Merth)*sumhe(Kk))
c
c check if alpha is time variable gravitational constant
               else if(icmkkk.eq.32) then
                  Fn(1) = Fn(1) +Tvary*Shar2*Gama*(
     .             (1._10+Mmoon)*sumhm(Kk)-(1._10+Merth)*sumhe(Kk))
c
c check if alpha is mass of moon
               else if(icmkkk.eq.10) then
                  Fn(1) = Fn(1) +Shar2*Gamat*Mass(3)*(
     .             sumhm(Kk)+sumhe(Kk))
c
c check if alpha is mass of embary
               else if(icmkkk.eq.3) then
                  Fn(1) = Fn(1) +Shar2*Gamat*(Mass(10)*sumhm(Kk)-
     .             Masse*sumhe(Kk))
               endif
            endif
c
c tidal effects on partial derivatives
            if(Km(84).gt.0) then
               call MONTID(k,s)
               Fn(1)=Fn(1) + DOT(Dmcor(1,Kkk),Dadxt(1,Kk))
     .          + DOT(Dmcor(4,Kkk),Dadvt(1,Kk))
            endif
c
c effect of principle of equivalence on partial derivatives
            if(Km(85).gt.0) then
               if(icmkkk.ge.-20 .and. icmkkk.le.44) then
                  if(icmkkk.le.-19 .or. icmkkk.ge.43) then
                     delsp = Gamat/2._10*((Ecor(Kk)/Re3-Mcor(Kk)/Rm3)
     .                       - (Mass(10)+Masse)*Mass(3)*Mecor(Kk)
     .                       /Rem3 + (summ(Kk)-sume(Kk)))
                     deldp = Gamat/2._10*((Ecor(Kk)/Re3+Mcor(Kk)/Rm3)
     .                       - (Mass(10)-Masse)*Mass(3)*Mecor(Kk)
     .                       /Rem3 + (-summ(Kk)-sume(Kk)))
                     if(icmkkk.eq.-19) Fn(1) = Fn(1) + delsp
                     if(icmkkk.eq.-20) Fn(1) = Fn(1) + deldp
                     if(icmkkk.eq.43) Fn(1)  = Fn(1)
     .                    + 4._10*(Delpls*delsp + Delmns*deldp)
                     if(icmkkk.eq.44) Fn(1)  = Fn(1)
     .                    - Delpls*delsp - Delmns*deldp
                  endif
               endif
            endif
c
c effect of general relativity on partial derivatives
            if(Km(61).gt.0) call MONREL(k)
            if(icmkkk.eq.-8) Fn(1) = Fn(1) + 2._10*gdpprv(Kk)
     .          - gdpfc2*gdppar(Kk)
c
            if(Km(33).ge.0) then
c
c effect of sun on partial derivatives
               if(icmkkk.eq.32) Fn(1) = Fn(1) + sums(Kk)*Gama*Tvary
            endif
c
c effect of central force on partial derivatives
            psum = 0._10
c
c modify i.c. partials in Encke integrations
            if(icmkkk.le.-31) then
               kkk1=-icmkkk-30
               if(Km(100).lt.0) then
               else if(Km(100).eq.0) then

c subtract elliptic orbit acceleration partial
                  psum = psum -
     .             Gamem*(3._10*Ylpt(Kk)*DOT(Ylpt,Dylpt(1,kkk1))/Rylpt2
     .             - Dylpt(Kk,kkk1))/Rylpt3
               else

c subtract mean lunar orbit acceleration partial
                  psum = psum - Dylun(Kk+6,kkk1)
               endif
            endif

c add central force effect on partial to other effects
            Fn(1) = Fn(1) + psum

c add indirect term due to gravity gradient
            Fn(1)=Fn(1)+DOT(Dmcor(1,Kkk),dadx(1,Kk))

            if(icmkkk.eq.3) then
               Fn(1) = Fn(1) - Mecor(Kk)/Rem3*Gamat
            else if(icmkkk.eq.32) then
               Fn(1) = Fn(1) - Mecor(Kk)/Rem3*Gamem*Tvary
            endif
         else
c
c-----------------------------------------------------------------------
c
c           determine right side of equations of motion for orbit
c
c
c           effect of perturbing planets on motion of moon
            sump(Kk) = 0._10
            do l = 1,9
               if(l.ne.3 .and. Km(l+30).ge.0) then
                  plf(Kk,l)=pmcor(Kk,l)/rpm3(l) - pecor(Kk,l)/rpe3(l)
                  sump(Kk) = sump(Kk) + Mass(l)*plf(Kk,l)
               endif
            end do
            Fn(1) = Fn(1) + Gamat*sump(Kk)
c
c effect of harmonics of the earth and moon on orbital motion
            call MORHAR(k,j)
c
c effect of solar j2
            if(Km(63).ge.0) then
               sumhe(Kk)= (Ecor(Kk)*gfacte2-3._10*gfacte*C3(Kk))/Re4
               sumhm(Kk)= (Mcor(Kk)*gfactm2-3._10*gfactm*C3(Kk))/Rm4
               Fn(1)    = Fn(1) + Shar2*Gamat*((1._10+Mmoon)*sumhm(Kk)-
     .          (1._10+Merth)*sumhe(Kk))
            endif

c
c tidal effects on the motion of the moon
            if(Km(84).ge.0) call MONTID(k,s)

c
c effect of principle of equivalence on the motion of the moon
            if(Km(85).ge.0) then
               summ(Kk) = 0._10
               sume(Kk) = 0._10
               do l = 1,9
                  if(l.ne.3) then
                     if(Km(l+30).ge.0) then
                        summ(Kk) = summ(Kk)+Mass(l)*pmcor(Kk,l)/rpm3(l)
                        sume(Kk) = sume(Kk)+Mass(l)*pecor(Kk,l)/rpe3(l)
                     endif
                  endif
               end do
               Fn(1) = Fn(1)
     .                 - Gamat*((Delsum-Deldif)/2._10*Mcor(Kk)/Rm3 -
     .                 (Delsum+Deldif)/2._10*Ecor(Kk)/Re3 +
     .                 (Masse*Mass(3)*(Delsum-Deldif)/2._10
     .                 +Mass(10)*Mass(3)*(Delsum+Deldif)/2._10)
     .                 *Mecor(Kk)/Rem3 - summ(Kk)*(Delsum-Deldif)/2._10
     .                  + sume(Kk)*(Delsum+Deldif)/2._10)
            endif
c
c effect of general relativity on the motion of the moon
            if(Km(61).ge.0) call MONREL(k)

c solar radiation pressure
            if(Km(86).gt.0) Fn(1)=Fn(1)+radp(Kk)
c
c effect of the sun on the motion of the moon
            if(Km(33).ge.0) then

c if sun effect not included, ignorable underflows previous to this
               sums(Kk) = Ecor(Kk)/Re3 - Mcor(Kk)/Rm3
               Fn(1)    = Fn(1) + Gamat*sums(Kk)
            endif
c
c effect of central force on motion of the moon
            sum = -Gamtem*Mecor(Kk)/Rem3
            if(Km(100).lt.0) then
            else if(Km(100).eq.0) then

c subtract elliptic orbit acceleration
               sum = sum + Gamem*Ylpt(Kk)/Rylpt3
            else

c subtract mean lunar orbit acceleration
               sum = sum - Ylun(Kk + 6)
            endif

c add central force to other perturbing forces
            Fn(1) = Fn(1) + sum
         endif
c
c----------------------------------------------------------------------
c
c           determine right side of equations of motion for rotation
c
c
      else if(k.le.Krot) then
         if(Kkk.eq.0) then
            Fn(1) = fnr(Kk)
         else
c
c-----------------------------------------------------------------------
c
c           determine right side of equations for partial derivatives
c           of rotation with respect to parameters
            if(Kkm(60).gt.2) write(6,99340) Kk,Kkk
99340       format('kk,kkk=',2i10)
            if(Kk.ne.1) goto 600
            icmkkk = Icrtrl(Kkk)
            ddw1dp = 0._10
            ddw2dp = 0._10
            ddw3dp = 0._10
            do l = 1,6
               dydp(l) = Y(k+l-4,j)
               dycdp(l) = Y(k+l-4+Krot-Korb,j)
            end do
c
c check if paramter is core-mantle rotation coupling constant
            if(icmkkk.eq.-17) then
               ddw1dp=cmtrq(1)*I0i(1,1)
               ddw2dp=cmtrq(2)*I0i(2,2)
               ddw3dp=cmtrq(3)*I0i(3,3)
               goto 500
c
c check if paramter is mass of earth-moon barycenter
            else if(icmkkk.eq.3) then
               ddw1dp = Gamat*Masse*Mma*hxe(1)
               ddw2dp = Gamat*Masse*Mmb*hxe(2)
               ddw3dp = Gamat*Masse*Mmc*hxe(3)
               if(Corint) goto 400
               goto 500
c
c check if paramter is mass of moon
            else if(icmkkk.eq.10) then
               ddw1dp = -Gamtem*Mma*hxe(1)
               ddw2dp = -Gamtem*Mmb*hxe(2)
               ddw3dp = -Gamtem*Mmc*hxe(3)
               if(Kmr(81).ge.0) then
c Include the term due to the monthly modulation of the
c Sun-Moon distance.
                  call PRODCT(Dhsdx,Smecor,temp,3,3,1)
                  call CROSS(temp,Smcor,temq)
                  call CROSS(Hs,Smecor,temp)
                  do i=1,3
                     Ddwdp(i)=Ddwdp(i)+Gamat*mmabc(i)*(temp(i)+temq(i))
                  end do
               endif
               if(Corint) goto 400
               goto 500
c
c check if paramter is gmvary
            else if(icmkkk.eq.32) then
               ddw1dp = Gama3*Tvary*Mma*hxe(1)
               ddw2dp = Gama3*Tvary*Mmb*hxe(2)
               ddw3dp = Gama3*Tvary*Mmc*hxe(3)
               if(Kkm(60).eq.-5) write(6,11441) k,j,s-2440000.5_10,
     .          dw,(n0(i)*I0i(i,i),i=1,3),ddwdp
11441          format('k,j,s,dw ',2i3,f9.4,(1p3d22.15))
c include solar torque
               if(Kmr(81).ge.0) then
                  ddw1dp = ddw1dp + Gama*Tvary*Mma*hxs(1)
                  ddw2dp = ddw2dp + Gama*Tvary*Mmb*hxs(2)
                  ddw3dp = ddw3dp + Gama*Tvary*Mmc*hxs(3)
                  if(Kkm(60).eq.-5) write(6,11442) ddwdp
11442             format(1p3d22.15)
               endif
               if(Kkm(60).eq.-5) write(6,11442) Mrotlb,Mecor,Mcor,Ecor,
     .          (Xpert(i,4),i=1,3)
               goto 500

c see if parameter is k2(moon) or lunar response lag
            else if(icmkkk.eq.-6 .or. icmkkk.eq.-7) then
               if(Kkm(60).gt.2) write(6,99341)
99341          format('calling MOREDP')
               call MOREDP
               goto 550
            endif
c
c check if parameter is lunar harmonic coefficient
            if(icmkkk.ge.Jzone .and. icmkkk.le.Jsin)
     .       call MORHAR(k,j)
c
c check if parameter is beta, gamma, j2, or core moment
            if(icmkkk.eq.-3 .or. icmkkk.eq.-4 .or. icmkkk.eq.-18 .or.
     .       (icmkkk.eq.Jzone.and.Kkk.eq.Izone(1))) goto 400 
            goto 500
         endif
         goto 700
c----------------------------------------------------------------------
c
c        determine right side of equations of motion for core rotation
c
c
      else if(Kkk.eq.0) then
         Fn(1) = fnc(Kk)
      else
c
c-----------------------------------------------------------------------
c
c     determine right side of equations for partial derivatives
c     of core rotation with respect to parameters

         if(Kk.ne.1) goto 380
         icmkkk = Icrtrl(Kkk)
         ddw1dp = 0._10
         ddw2dp = 0._10
         ddw3dp = 0._10
         do l = 1,6
            dycdp(l) = Y(k+l-4,j)
            dydp(l) = Y(k+l-4+Korb-Krot,j)
         end do
c
c check if paramter is core-mantle rotation coupling
         if(icmkkk.eq.-17) then
            ddw1dp=mctrq(1)*I0ci(1,1)
            ddw2dp=mctrq(2)*I0ci(2,2)
            ddw3dp=mctrq(3)*I0ci(3,3)
c check if parameter is (spherical) core moment of inertia
         else if(icmkkk.eq.-18) then
            ddw1dp=-dwc(1)*I0ci(1,1)
            ddw2dp=-dwc(2)*I0ci(1,1)
            ddw3dp=-dwc(3)*I0ci(1,1)
         endif
c calculate the direct partial
         dfdp(1) = (Sphic*ddw1dp + Cphic*ddw2dp)/Sthetac
         dfdp(2) = Cphic*ddw1dp - Sphic*ddw2dp
         dfdp(3) = ddw3dp - Cthetac*dfdp(1)
c add up the indirect partials due to core and mantle states
  380    Fn(1)   = Fn(1) + DOTN(dfcdy(1,Kk),dydp,6)
     .    +DOTN(dfcdyc(1,Kk),dycdp,6)
         if(icmkkk.gt.-31) Fn(1) = Fn(1) + dfdp(Kk)
      endif
      goto 700

  400 ddw1dp = ddw1dp + Gamat3*Dmma(Kkk)*hxe(1)
      ddw2dp = ddw2dp + Gamat3*Dmmb(Kkk)*hxe(2)
      if(Kkm(60).eq.-6) write(6,99342) k,kk,kkk,s-2440000.5_10,ddw3dp
99342 format('k,kk,kkk',3i3,f9.4,' dw3=',1pe22.15)
      ddw3dp = ddw3dp + Gamat3*Dmmc(Kkk)*hxe(3)
      if(Kkm(60).eq.-6) write(6,99343) 'ehar',ddw3dp,Dmmc(kkk)
99343 format(a,1p2e22.15)
      if(Kmr(81).ge.1) then
         ddw1dp = ddw1dp + Gamat*Dmma(Kkk)*hxs(1)
         ddw2dp = ddw2dp + Gamat*Dmmb(Kkk)*hxs(2)
         ddw3dp = ddw3dp + Gamat*Dmmc(Kkk)*hxs(3)
         if(Kkm(60).eq.-6) write(6,99343) 'shar',ddw3dp
      endif
      if(Corint) then
         ddw1dp = ddw1dp + Mrcond(23)*cmtrq(1)*Dmmam(Kkk)
         ddw2dp = ddw2dp + Mrcond(23)*cmtrq(2)*Dmmbm(Kkk)
         ddw3dp = ddw3dp + Mrcond(23)*cmtrq(3)*Dmmcm(Kkk)
         if(Kkm(60).eq.-6) write(6,99343) 'core',ddw3dp,Dmmcm(kkk)
      endif

      ddw1dp = ddw1dp - W2*W3*Dalpha(Kkk)
      ddw2dp = ddw2dp + W1*W3*Dbeta(Kkk)
      ddw3dp = ddw3dp - W1*W2*Dgamma(Kkk)
      if(Kkm(60).eq.-6) write(6,99343) 'kine',ddw3dp,Dgamma(kkk)

c go pick up indirect c22 contributions
      if(Kkm(60).gt.2) write(6,99338) Kk,Kkk
99338 format('calling MORHR2,kk,kkk=2i10')
      call MORHR2
c morhr2 is entry point in morhar
c
c correct angular acceleration partials for earth-figure torque
  500 if(Kmr(82).ge.0) call MORFFIP
c correct angular acceleration partials for elasticity and dissipation
      if(Kmr(83).ge.0) call MOREDPC(Kmr(83))
      if(Kkm(60).gt.2) write(6,99339)
99339 format('called MOREDPC')

c obtain explicit partial terms
  550 continue
      if(Kkm(60).eq.-6) write(6,99343) 'final',ddw3dp
      dfdp(1) = (Sphi*ddw1dp + Cphi*ddw2dp)/Stheta
      dfdp(2) = Cphi*ddw1dp - Sphi*ddw2dp
      dfdp(3) = ddw3dp - Ctheta*dfdp(1)
      if(icmkkk.eq.32 .and. Kkm(60).eq.-5) write(6,11443) fnr,dfdp,
     . (DOTN(dfdx(1,i),Dmcor(1,Icrref(Kkk)),6),i=1,3),
     . (DOTN(dfdy(1,i),dydp,6),i=1,3),
     . (DOTN(dfdyc(1,i),dycdp,6),i=1,3)
11443 format('fnr ',1p3d22.15/'dfdp',3d22.15/'dfdx',3d22.15/
     . 'dfdy',3d22.15/'dfdc',3d22.15)
  600 Fn(1)   = Fn(1) + DOTN(dfdy(1,Kk),dydp,6)
      if(Corint) Fn(1)=Fn(1)+DOTN(dfdyc(1,Kk),dycdp,6)
      if(icmkkk.gt.-31) Fn(1) = Fn(1) + dfdp(Kk)
c indirect term due to position and velocity dependence
      if(Icrref(Kkk).gt.0)
     . Fn(1)=Fn(1)+DOTN(dfdx(1,Kk),Dmcor(1,Icrref(Kkk)),6)
c
c
c-----------------------------------------------------------------------
c
c           return right side of differential equation to calling prog.
  700 MORFN = Fn(1)
      if(Kkm(60).gt.0)
     . write(6,99333) k,j,s-2440000.5_10,Fn(1)
99333 format('k,j,s,fn',2i3,f11.4,1pd22.15)
      return
      end
