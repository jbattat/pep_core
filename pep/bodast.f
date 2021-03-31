      subroutine BODAST(astflg)
 
      implicit none

c
c m.e.ash   may 1974   subroutine bodast
c calculate asteroid perturbations (indidual and ring)
c for n-body integration
c
      logical   astflg
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'bdctrl.inc'
c kkbdy(j) =0 asteroid with planet number j+10 does not perturb n-body
c             integraton (j=1,20)
c kkbdy(j) =1 such individual asteroid does perturb n-body integration
c
c kkbdy(80)=0 do not include asteroid ring in n-body integration
c kkbdy(80)=1 effect of asteroid ring included in n-body integration
c
      include 'empcnd.inc'
      include 'fmmips.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'namtimq.inc'
      include 'param.inc'
      include 'petina.inc'
      include 'smlbdy.inc'
      include 'smlstf.inc'
c
c for distributed asteroidal perturbation force
      real*10 astnrm(3), zz(32), za(32)
      real*10 sinnod, cosnod, cosinc, sininc
      real*10 er(3), naxer(3), astd2, sith, coth, etheta(3),
     .          leg(32), leg1(32), asta, astb, astd, aste,
     .          astws1, astws2, astws3, astws4, astws5, astws6, astds
      real*10 rs, tbcor(3), tsb2, tsb, ratio0, ratio1, cos45,
     .          leg45(32), leg145(32), tstnrm(3)
      real*10 gauss9(12)
      real*10 sbcor(3), rsb, rsb2, rsb3
      logical   flg45
      real*10 astsm,crosan,DOT,den,dum,rysb2,rysb3,rysml,rysml2,rysml3,
     . s,ys3(3),ysc3,yscor(3)
      integer*4 i,iden,ij,iquit,j,jquit,k,kk,l,mm,n,nc,nexpan
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c           setup individual asteroids and asteroid ring
c           called from bodset
c
c           setup for perturbing asteroid or satellite orbits
c           (elliptic orbits assumed for bodies between 11 and 30)
      Numast = 0
      astflg = .false.
      do j = 11, 30
         if(Kkbdy(j-10).gt.0) then
            Numast = Numast + 1
            if(Numast.gt.12) call SUICID(
     .       ' MORE THAN 12 PERTURBING ASTEROIDS, STOP IN BODAST  ',13)
            Kpast(Numast) = Kkbdy(j - 10)
            do k = 1, Numpln
               if(Nqlnt(k).eq.j) then
                  Kast = k
                  go to 20
               endif
               end do
            call SUICID(
     .            ' PERTURBING ASTEROID NOT INPUT, STOP IN BODAST  ',12)
   20       if(Jcnd(Kast).ne.0) call SUICID(
     .   ' WRONG TYPE OF PERTURBING ASTEROID I.C., STOP IN BODAST ',14)
            Nplast(Numast) = Nqlnt(Kast)
            Ncnast(Numast) = Npcent(Kast)
            if(Ncnast(Numast).gt.0) call SUICID(
     .   'PERTURBING ASTEROID CENTRAL BODY NOT SUN, STOP IN BODAST',14)
            if(Jdpl0(Kast).gt.0) then
               Tlpast(Numast) = Jdpl0(Kast)
            else
 
c integration time variable is julian date + 0.5
               Tlpast(Numast) = Con11(1, Kast + 4) + 0.5_10
            endif
            Astmab(Numast) = Mass(j)
            Astmac(Numast) = 1._10
            if(Ncnast(Numast).gt.0) then
               nc = Ncnast(Numast)
               gauss9(Numast) = Gauss*SQRT(Mass(nc)*(1._10+Mass(j)))
               Astmas(Numast) = Mass(j)*Mass(nc)
               Astmac(Numast) = Mass(nc)
            else
               gauss9(Numast) = Gauss*SQRT(1._10 + Mass(j))
               Astmas(Numast) = Mass(j)
            endif
            call JNITL(gauss9(Numast), Pcond(1,Kast), Ast999(1,Numast),
     .                 0, dum)
         endif
         end do
      if(Numast.gt.0) then
         astflg = .true.
         write(Iout, 50) (Nplast(i), i = 1, Numast)
   50    format('0INDIVIDUAL PERTURBING ASTEROIDS ARE', 12I4)
      endif
c
c setup for distributed asteroidal perturbing force
      if(Kkbdy(80).gt.0) then
         astflg = .true.
         write(Iout, 100)
  100    format('0PERTURBATIONS DUE TO ASTEROID RING ARE INCLUDED')
         sinnod    = SIN(Asbasc*Convd)
         cosnod    = COS(Asbasc*Convd)
         cosinc    = COS(Asbinc*Convd)
         sininc    = SIN(Asbinc*Convd)
         astnrm(1) = sinnod*sininc
         astnrm(2) = -cosnod*sininc
         astnrm(3) = cosinc
         call LEGNDR(0._10, 1._10, 32, 0, za, zz, 0._10, 0._10)
         do mm = 1, 31
            zz(mm + 1) = za(mm)
            end do
         zz(1) = 0._10
         rs    = Asba*SQRT(2._10)
c
c     asteroid ring (inside) from l.friedman's mit thesis.
c     may 1974  m.ash   removed from subroutine solprb, added to n-body
c     no partials in n-body integration, just perturbing accelerations
c
         flg45 = .false.
      endif
c
c set up limited asteroids as perturbing bodies
c the masses are assumed to be negligible in determining the orbits,
c and so the scale may be squeezed to get the mean motion right
      do j = 1, 3
         Sumsml(j) = 0._10
      end do
      if(Kbdy(25).ge.0 .and. Numsml.gt.0) then
         write(Iout,105)
  105    format('0PERTURBATIONS DUE TO SMALL ASTEROIDS INCLUDED')
         astflg = .true.
         do i = 1, Numsml
            den = 2.0
            if(Denpts(i).gt.0) den = Prmter(Denpts(i))
            Smlvol(i)=Scond(7,i)**3*(4._10/3._10)*Pi*1E15_10/2E33_10
            Smlmas(i)=Smlvol(i)*den
            call JNITL(Gauss,Scond(1,i),Elptsm(1,i),0,dum)
         end do
      endif
      return
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c           calculations at start of integration step
c           called from bodfn
      entry BDAST0(s)
c
c determine perturbing asteroids relative to sun
      if(Numast.gt.0) then
         do l = 1, Numast
            call JLIPT(s-Tlpast(l),Ast999(1,l),0,Yast(1,l),
     .       Rastc(l),Rastc2(l),Rastc3(l),dum)
            do i = 1, 3
               Yastc3(i, l) = Yast(i, l)/Rastc3(l)
            end do
 
         end do
      endif
 
c compute limited asteroid coordinates
      if(Kbdy(25).ge.0 .and. Numsml.gt.0) then
         do i = 1,Numsml
            call JLIPT((Jd-Jdsml0(i))+Fract,Elptsm(1,i),0,
     .       Ysml(1,i),rysml,rysml2,rysml3,dum)
            do j = 1, 3
               Ysml3(j,i) = Ysml(j,i)/rysml3
            end do
         end do
      endif
      return
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
c           compute acceleration due to asteroids on coordinate kk
c           called from bodfn
      entry BDASTF(kk,sbcor,rsb,rsb3,astsm)
      astsm = 0._10
      if(Numast.gt.0 .or. Kkbdy(80).gt.0 .or. Kbdy(25).ge.0) then
c
c * * * * * * * * * * * * * * *
c
c           computations at start of iteration of a given step
c           for a given body
         if(kk.le.1) then
c
c determine perturbing asteroids relative to given body
            if(Numast.gt.0) then
               do l = 1,Numast
                  do i = 1,3
                     Yastb(i,l) = Yast(i,l) - sbcor(i)
                  end do
                  Rastb2(l) = Yastb(1,l)**2+Yastb(2,l)**2+Yastb(3,l)**2
                  Rastb(l)  = SQRT(Rastb2(l))
                  Rastb3(l) = Rastb2(l)*Rastb(l)
               end do
            endif
c
c determine if distributed asteroidal pertubation is included
            if(Kkbdy(80).gt.0) then
               rsb2 = rsb*rsb
               tsb2 = 0._10
               do i = 1,3
                  tbcor(i) = sbcor(i) + Asba*astnrm(i)
                  tsb2     = tsb2 + tbcor(i)**2
               end do
               tsb    = SQRT(tsb2)
               ratio0 = rsb/Asba
               if(ratio0.ge.1._10) ratio0 = 1._10/ratio0
               ratio1 = tsb/rs
               do i = 1,3
                  tstnrm(i) = astnrm(i)
               end do
               if(ratio1.lt..1_10) then
c expansion around new origin ill-conditioned at new origin
c move the new origin to the other side of the old origin
                  do i = 1,3
                     tstnrm(i) = -astnrm(i)
                     tbcor(i)  = sbcor(i) + Asba*tstnrm(i)
                  end do
                  tsb2   = DOT(tbcor,tbcor)
                  tsb    = SQRT(tsb2)
                  ratio1 = tsb/rs
               endif
               nexpan = 0
               if(ratio1.gt.1._10) nexpan   = 1
               if(ratio1.ge.1._10) ratio1   = 1._10/ratio1
               if(ratio1.ge.ratio0) nexpan = -1
               if(nexpan.eq.-1) then
                  astd2 = ratio0**2
                  do ij = 1,3
                     er(ij) = sbcor(ij)/rsb
                  end do
                  coth = DOT(astnrm,er)
                  call CROSS(astnrm,er,naxer)
                  call CROSS(naxer,er,etheta)
 
c partial derivative coding from solprb removed
                  sith = SQRT(1._10 - coth**2)
                  call LEGNDR(coth,sith,32,0,leg,leg1,0._10,0._10)
                  asta  = 0._10
                  astb  = 0._10
                  astds = 1._10
 
c partial derivative coding from solprb removed
                  do mm = 2,32,2
                     astds = astds*astd2
                     if(rsb.le.Asba) then
                        astws1 = astds*mm*zz(mm)
                     else
                        astws1 = astds*(mm + 1)*zz(mm)
                     endif
                     astws3 = astws1*leg(mm - 1)
                     astws2 = astds*rsb*zz(mm)
                     astws4 = astws2*leg1(mm - 1)
                     asta   = asta + astws3
                     astb   = astb + astws4
                     if(ABS(astws3).lt.1.E-16_10*ABS(asta) .and.
     .                  ABS(astws3).lt.1.E-16_10*ABS(astb) )
     .                  go to 110
 
c partial derivative coding from solprb removed
                  end do
  110             astws5 = Asbmas/(rsb2*Asba)
                  astws6 = Asbmas/(rsb2*rsb)
               else
                  if( .not. flg45) then
                     cos45 = 1._10/SQRT(2._10)
                     call LEGNDR(cos45,cos45,32,0,leg45,leg145,
     .                           0._10,0._10)
                     flg45 = .true.
                  endif
                  do i = 1,3
                     er(i) = tbcor(i)/tsb
                  end do
                  coth = DOT(tstnrm,er)
                  sith = SQRT(1._10 - coth**2)
                  call CROSS(tstnrm,er,naxer)
                  call CROSS(naxer,er,etheta)
                  call LEGNDR(coth,sith,32,0,leg,leg1,0._10,0._10)
                  asta = cos45*ratio1*coth
                  astb = cos45*ratio1*tsb
                  astd = cos45*ratio1*tsb*2._10
                  aste = ratio1*coth
                  if(nexpan.eq.1) then
                     aste = aste*2._10
                     asta = asta*2._10
                     astd = -.5_10*astd
                  endif
                  astds = ratio1
                  do n = 2,32
                     astds = astds*ratio1
                     if(tsb.le.rs) then
                        astws1 = astds*n*leg45(n - 1)
                     else
                        astws1 = astds*(n + 1)*leg45(n - 1)
                     endif
                     astws3 = astws1*leg(n - 1)
                     astws2 = astds*tsb*leg45(n - 1)
                     astws4 = astws2*leg1(n - 1)
                     iquit  = 1
                     jquit  = 1
                     if(ABS(astws3).gt.1.E-16_10*ABS(asta)) then
                        asta = asta + astws3
                     else
                        iquit = 0
                     endif
                     if(ABS(astws4).gt.1.E-16_10*ABS(astb)) then
                        astb = astb + astws4
                     else
                        jquit = 0
                     endif
                     if(iquit+jquit.eq.0) go to 120
 
c partial derivative coding from solprb removed
                  end do
  120             astws5 = Asbmas/(tsb2*rs)
                  astws6 = Asbmas/(tsb2*tsb)
               endif
            endif
c
c determine if limited asteroids are included
            if(Kbdy(25).ge.0) then
               do j = 1, 3
                  Sumsml(j) = 0._10
               end do
               do i = 1, Numsml
                  do j = 1, 3
                     yscor(j) = Ysml(j,i) - sbcor(j)
                  end do
                  rysb2 = DOT(yscor, yscor)
                  rysb3 = rysb2*SQRT(rysb2)
                  do j = 1, 3
                     ysc3     = yscor(j)/rysb3
                     ys3(j)   = ysc3 - Ysml3(j,i)
                     Sumsml(j)= Sumsml(j) + Smlmas(i)*ys3(j)
                  end do
               end do
            endif
         endif
c
c * * * * * * * * * * * * * * *
c
c           calculate acceleration on given body coordinate
c
c calculate acceleration due to individual asteroids
         if(Numast.gt.0) then
            do l = 1,Numast
               astsm = astsm + Astmas(l)
     .                 *(Yastb(kk,l)/Rastb3(l) - Yastc3(kk,l))
            end do
         endif
c
c effect of distributed asteroidal perturbation
         if(Kkbdy(80).gt.0) then
            if( ABS(rsb-Asba) .lt. .1_10*Asba ) then
               crosan = ACOS(DOT(astnrm,sbcor)/rsb)/Convd
               write(6,130) rsb,crosan
  130          format(' ',2E16.8)
            endif
            if(nexpan.lt.0) then
               if(rsb.le.Asba) then
                  astsm = astsm +
     .                    astws5*(asta*sbcor(kk) - astb*etheta(kk))
               else
                  astsm = astsm +
     .             astws6*((-1._10-asta)*sbcor(kk) - astb*etheta(kk))
               endif
            else if(nexpan.eq.0) then
               astsm = astsm + astws5*(asta*tbcor(kk) - astb*etheta(kk))
            else
               astsm = astsm +
     .          astws6*((-1._10-asta)*tbcor(kk) - astb*etheta(kk))
            endif
         endif
c
c effect of limited asteroids
         if(Kbdy(25).ge.0 .and. Numsml.gt.0) then
            astsm = astsm + Sumsml(kk)
         endif

      endif
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c
      return
      end
