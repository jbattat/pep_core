      subroutine PLNORB(ncall)
 
      implicit none
c
c n.beebe/r.king/r.reasenberg  sept 1970  subroutine plnorb
c computation of forces acting on artificial planetary orbiter
c

c arguments
      integer*4 ncall

c        ncall=-2 setup once per iteration of a given step for partials
c        ncall=-1 setup once per iteration of a given step for motion
c        ncall= 0 setup once per step
c        ncall= k evaluate acceleration for equation k.gt.0
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'astroi.inc'
      include 'cnthar.inc'
      include 'drgprt.inc'
      include 'engstf.inc'
      include 'gleaks.inc'
      real*10 tzero2,tzero,ax,ay,az,bz,tw
      equivalence (tzero2,Con1(6))
      equivalence (tzero,Con1(4)),(ax,Con),(ay,Con(2)),
     .            (az,Con(3)),(bz,Con(4)),(tw,Con(8))
      include 'intstf.inc'
      include 'lothrf.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'sbstuf.inc'
      include 'sbthng.inc'
c
c local
      real*10 at0,at3,at4,atmexp,cosprj,delt1,
     .          delt2,drhodr,fx,fxx,fy,fyy,fz,fz1,fz2
      real*10 fzz1,fzz2,haryc,par,solcn5,vatm
      integer   i,i1,i2,idsh,ipan,j,kdrggo,kradgo,ltfgo,ndish

c conversion factor for cm/sec**2 to a.u./day**2
      real*10 harycn/4.9900168876E-4_10/
 
c ltcon- solar intensity at 1 a.u. (dyne/cm**2)
      real*10 ltcon/4.65E-5_10/, ltcon2
      real*10 rradpr(3,3),altprs(3),radprs(3),pf(3)
      real*10 solcn(3)/3*0._10/,solcon,solcn2,solcn3
      equivalence (solcn,solcon),(solcn(2),solcn2),
     .            (solcn(3),solcn3)
      real*10 rfla(3),rfli(3)
      logical*4 opsflg
      real*10 ssavd/0._10/, sdel/0._10/
      real*10 srfmvm(3)
      real*10 rlim/1._10/, day/8.6400E4_10/, sbccm(6), vrotsb(3),
     .          velatm(3),vwind(3)/3*0._10/,vnorm(3),unorm(3),
     .          velatu(3),dpitch,vncrsv(3),vlift(3),drgmu,drgmag,
     .          dadro0(3),rsbkm,atdrag(3),rho,dadsh(3),at1,at2,aucm,
     .          s,rvec(3),bvec(3),bvmag,avec(3),rmtx(3,3),glacc(3)
      equivalence (rvec,rmtx(1,3))
      integer*2 idrag,iset/0/
      real*10 ltconz
      real*4    epsdsh
      real*10 eapan(3),eadsh(3)
      real*10 yvec(3),pb3dr,ry
c
c external functions
      real*10 DOT

c        sh = scale height (km)
c        rhoz= density at height h0 (gm/cm**3)
c        drgeps = upper limit on exponent in atmden
c        cdrg1 = drag coefficient
c        cdrg2 = drag coefficient
c        rhzkm = distance from center of planet to height h0
c
      if(ncall .lt. 0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c*  start=200
c        set-up once per iteration of a given step for motion
         if(ncall .ge. -1) then
 
c set up unit vectors for gas leak ,rad. pres., and drag
            if(Kp(83) .gt. 0 .or. Kdirct .gt. 0 .or. Kp(82) .ne. 0)
     .          then
 
c zero effect of drag on acceleration partials
               do i1 = 1, 3
                  do i2 = 1, 3
                     Dadxl(i1, i2) = 0._10
                     Dadxm(i1, i2) = 0._10
                     end do
                  end do
 
               rvec(1) = Bcor(1)/Rb
               rvec(2) = Bcor(2)/Rb
               rvec(3) = Bcor(3)/Rb
               if(Kp(83) .le. 0 .and. Kdirct .le. 0) go to 40
               call CROSS(Canvec, rvec, bvec)
               bvmag   = SQRT(bvec(1)**2 + bvec(2)**2 + bvec(3)**2)
               bvec(1) = bvec(1)/bvmag
               bvec(2) = bvec(2)/bvmag
               bvec(3) = bvec(3)/bvmag
               call CROSS(bvec, rvec, avec)
c in moyer's tr 32-1527  avec=t, bvec=n, beta=k
c determine if gasleaks are included
               if(Kp(83) .gt. 0) then
                  do i = 1, 3
                     rmtx(i, 1) = avec(i)*Cbeta + bvec(i)*Sbeta
                     rmtx(i, 2) = -avec(i)*Sbeta + bvec(i)*Cbeta
 
c rmtx(i,3)=rvec(i)      by equivalence
                     end do
c*  start=300
c gas leak acceleration
c
                  delt1 = s - tzero
                  fxx   = haryc
                  fyy   = haryc
                  fzz1  = haryc
                  fzz2  = haryc
                  if(Tx .ne. 0._10) fxx  = EXP(-Tx*delt1)*fxx
                  if(Ty .ne. 0._10) fyy  = EXP(-Ty*delt1)*fyy
                  if(Tz .ne. 0._10) fzz1 = EXP(-Tz*delt1)*fzz1
                  if(tw .ne. 0._10) fzz2 = EXP(-tw*delt1)*fzz2
                  fx  = ax*fxx
                  fy  = ay*fyy
                  fz1 = az*fzz1
                  fz2 = bz*fzz2
                  fz  = fz1 + fz2
 
c gas leak acceleration in standard equatorial coordinates
                  do i = 1, 3
                     glacc(i) = rmtx(i, 1)*fx + rmtx(i, 2)
     .                          *fy + rmtx(i, 3)*fz
                     end do
               endif
               if(Kdirct .le. 0) go to 40
c*  start=400
c        radiation pressure setup
c
c        kp(81)= kreflt*10+ kdirct
c          kdirct=0 no radiation pressure
c                =1 direct radiation,exact shadow
c                =2 direct radiation,approx. shadow.(not coded)
c                =3 direct radiation, no shadow
c
c          kreflt=1 albedo radiation
c                =2 infared radiation
c                =3 both inrared and albedo
c        kradgo    =iabs(kkp(78) )
c        0         mariner-9 general-purpose 3-axis
c                  seven-term model
c        1         rad georgevic model put in by
c                  l.e. conrod and r.mckinnis    (1975)
c        2         mvm model of c christensen
c        3         viking model of c christensen
c        4         pvo - put in by r.b.g.  dec. 1978
c
c
               Lambda = 1._10
               if(Ncentr .gt. 0) then
                  if(Kdirct .ne. 3) then
                     if(Kdirct .eq. 2) then
 
c simple on-off shadow avoid call to shadow routine
                        if(Rb .le. Rc) go to 10
                        if(Rsb2 - DOT(Sbcor,rvec)**2 .lt. Crad**2)
     .                      Lambda = 0._10
                     else
                        call SHADOW
                     endif
 
c end simple shadow check
                     if(Lambda .le. 0._10) then
 
c lambda not zero. direct radiation pressure forces included
                        do i = 1, 3
                           radprs(i) = 0._10
                           end do
                        go to 40
                     endif
                  endif
               endif
   10          solcn5 = ltconz*Lambda/Rb2
               solcon = ltcon2*Lambda/Rb2
c
c rotation matrix for radiation pressure
               do i = 1, 3
                  rradpr(i, 1) = avec(i)*Cbeta2 + bvec(i)*Sbeta2
                  rradpr(i, 2) = -avec(i)*Sbeta2 + bvec(i)*Cbeta2
                  rradpr(i, 3) = rvec(i)
                  end do
c
c acceleration in s/c coordinates
c basic 3 axis model
               altprs(1) = Con(9)*solcon
               altprs(2) = Con(10)*solcon
               altprs(3) = Con(11)*solcon
               kradgo    = Kkp(78)
               idsh = 0
 
c mission specific models
               if(kradgo .eq. 1) then
c for rg mvm acceleration must have both kdirct on and kkp(78)
c not zero then to bypass the classical direct calculation,
c set con(9) through con(15) to zero
                  call MVMSRF(s, srfmvm)
                  do i = 1, 3
c factor 1E-6_10 converts r.g.'s units of force to newtons
c factor 1E5_10 converts newtons to dynes
                     altprs(i) = altprs(i) + srfmvm(i)*1E-1_10*haryc
                     end do
               else if(kradgo .eq. 2) then
 
c mvm-specific model of c christensen
                  call PANEL(s, eapan, ipan)
                  if(ipan .gt. 0) then
                     do i = 1, 3
                        altprs(i) = altprs(i) + eapan(i)*solcn5*Con(12)
                        end do
                  endif
 
c viking uses dish part of mvm model
                  ndish = 1
                  go to 20
               else if(kradgo .eq. 3) then
                  ndish = 1
                  go to 20
               else if(kradgo .eq. 4) then
                  ndish = 2
                  go to 20
               else
 
c time dependant terms
                  delt2     = s - tzero2
                  solcn2    = solcon*delt2
                  solcn3    = solcn2*delt2
                  altprs(1) = altprs(1) + Con(12)*solcn2
                  altprs(2) = altprs(2) + Con(13)*solcn2
                  altprs(3) = altprs(3) + Con(14)*solcn2 + Con(15)
     .                        *solcn3
               endif
               goto 30
c find epsdsh, the sun-s/c-axis angle
c assume dish axis points to earth
   20          pb3dr  = DOT(Pbcor(1,3), rvec)
               epsdsh = ACOS(pb3dr/Rpb(3))
               call DISH(epsdsh, eadsh, idsh, ndish)
               if(idsh .eq. 2) then
                  do i = 1, 3
                     altprs(i) = altprs(i) + eadsh(i)*solcn5*Con(14)
                     end do
 
c error in dish
               endif
c
c acceleration in standard equatorial coordinates
   30          call PRODCT(rradpr, altprs, radprs, 3, 3, 1)
 
               if(idsh .eq. 1) then
                  do i = 1, 3
                     yvec(i) = Pbcor(i, 3) - rvec(i)*pb3dr
                     end do
                  ry = SQRT(DOT(yvec(1),yvec(1)))
                  do i = 1, 3
                     yvec(i) = yvec(i)/ry
                     end do
                  do i = 1, 3
                     radprs(i) = radprs(i) + solcn5*Con(14)
     .                           *(eadsh(2)*yvec(i) + eadsh(3)*rvec(i))
                     end do
               endif
            endif
 
c setup for reflected radiation pressure
   40       if(Kreflt .le. 0) then
            endif
c
c*  start=700
c determine if atmospheric drag is included
            if(Kp(82) .gt. 0) then
 
c do not compute if probe is too far away
               if(Rsb .gt. rlim) then
 
c signal later calls of this iteration not to compute
                  idrag = -1
               else
                  idrag = 0
c get spacecraft engineering data
c sbeng maintains ssavd - gives time increment to next call
c do units conversions
                  do i = 1, 3
                     sbccm(i)     = Sbcor(i)*aucm
                     sbccm(i + 3) = Sbcor(i + 3)*aucm/day
                     end do
                  rsbkm = Rsb*aucm/1E5_10
c compute rotational velocity of planet and cross with radial
c vector to probe
                  call CROSS(Omega, sbccm, vrotsb)
c compute wind velocity
c
c compute velocity of probe wrt atmosphere
                  do i = 1, 3
                     velatm(i) = sbccm(i + 3) - vrotsb(i) - vwind(i)
                     end do
c
c compute unit normal and velocity vectors
                  at0  = velatm(1)**2 + velatm(2)**2 + velatm(3)**2
                  vatm = SQRT(at0)
                  do i = 1, 3
                     velatu(i) = velatm(i)/vatm
                     end do
 
c find unit vector perpendicular to principal area
                  kdrggo = Kkp(79)
                  if(kdrggo .eq. 1) then
 
c s/c faces sun
                     do i = 1, 3
                        vnorm(i) = -rvec(i)
                        end do
                  else if(kdrggo .eq. 2) then
 
c s/c faces in direction of canvec
                     do i = 1, 3
                        vnorm(i) = Canvec(i)
                        end do
                  else
 
c spherical s/c
                     do i = 1, 3
                        vnorm(i) = velatu(i)
                        end do
                     dpitch = 1
                     go to 50
                  endif
 
c if probe moving away from sun, flip normal vector
                  dpitch = DOT(vnorm, velatu)
                  if(dpitch .lt. 0) then
                     dpitch = -dpitch
                     do i = 1, 3
                        vnorm(i) = -vnorm(i)
                        end do
                  else if(dpitch .eq. 0) then
                     idrag = -1
                     go to 60
                  endif
 
c get density of atmosphere
   50             call ATMDEN(rho, drhodr, atmexp)
                  at0 = at0*day**2/(2._10*Sbmass*aucm)
                  at1 = at0*dpitch*Sbarea(2)
                  at2 = 0._10
                  at3 = -rho*Cdrg1*day/(2._10*Sbmass)
                  at4 = at3*Sbarea(2)*dpitch
                  if(Cdrg2 .ne. 0._10) then
                     call CROSS(vnorm, velatu, vncrsv)
                     call CROSS(vncrsv, velatu, vlift)
                     at2 = Cdrg2/SQRT(vncrsv(1)**2 + vncrsv(2)
     .                     **2 + vncrsv(3)**2)
                  endif
 
c compute atmospheric drag acceleration
                  do i = 1, 3
                     dadro0(i) = at1*(-Cdrg1*velatu(i) + at2*vlift(i))
                     atdrag(i) = rho*dadro0(i)
                     end do
c
c partial w.r.t velocity term for principal area
                  do i2 = 1, 3
                     Dadxm(i2, i2) = vatm
                     do i1 = 1, 3
                        Dadxm(i2, i1) = (Dadxm(i2,i1) + velatm(i1)*
     .                                  velatu(i2))*at4
                        end do
                     end do
 
                  if(Sbarea(3) .gt. 0._10) then
 
c find unit vector perpendicular to secondary area vnorm
                     call CROSS(velatu, vnorm, vncrsv)
                     call CROSS(vnorm, vncrsv, unorm)
                     dpitch = DOT(unorm, velatu)
                     at1    = dpitch*Sbarea(3)*at0
                     at4    = at3*dpitch*Sbarea(3)
 
c lift on secondary area neglected
                     do i = 1, 3
                        dadro0(i) = dadro0(i) - at1*Cdrg1*velatu(i)
                        atdrag(i) = rho*dadro0(i)
                        end do
                  endif
c
c find effect of drag on acceleration partials
c effect of lift on partial w.r.t. velocity term ignored
                  if(ABS(rho) .gt. 1E-40_10) then
                     do i2 = 1, 3
                        Dadxm(i2, i2) = Dadxm(i2, i2) + vatm*at4
                        cosprj = Sbcor(i2)/Rsb
                        do i1 = 1, 3
 
c compute dadxl-transpose for convenience in sbfn1
                           Dadxl(i2, i1) = (atdrag(i1)*drhodr*cosprj)
     .                        /rho
                           Dadxm(i2, i1) = Dadxm(i2, i1) + velatm(i1)
     .                        *velatu(i2)*at4
                           end do
                        end do
                  endif
c
c*  start=1000
c        start of code that establishes rlim, the s/c height above
c        which drag is of no interest.
c        iset
c             0 initial
c            -1 drag was insignificant
c             1 drag was significant
c             2 drag has changed from insig. to sig. or vice-versa
c        idrag
c            -1 insignificant
c             0 significant
c
                  if(iset .ne. 2) then
c if acceleration insignificant, set distance limit for
c future calls
c drgmu=-gamat3/rsb2*1.0E-16_10
                     drgmu = -Gamat3/Rsb2*Atmeps
c loop to prevent underflows that result from exponential
c atm. density and high eccentricity
                     do i = 1, 3
                        if(ABS(atdrag(i)) .le. 1E-38_10) atdrag(i)
     .                      = 0._10
                        end do
                     drgmag = SQRT(atdrag(1)**2 + atdrag(2)
     .                        **2 + atdrag(3)**2)
                     if(drgmag .lt. drgmu) then
 
c acceleration insignificant
                        idrag = -1
                        rlim  = Rsb
                        if(iset .ne. 1) then
                           iset = -1
                           go to 60
                        endif
 
c acceleration significant
                     else if(iset .ne. -1) then
                        iset = 1
                        go to 60
                     endif
 
c time to freeze rlim, the distance limit
                     rlim = rlim*1.1_10
                     iset = 2
                  endif
               endif
            endif
c
c*  start=1200
c determine if distributed asteroidal pertubation is included
   60       if(Nbelt .gt. 0) call ASBFNC(ncall)
c once per step procedures that follow the once per iteration
c procedured during the first iteration
            if(opsflg) then
               opsflg = .false.
 
c effect of visible reflected light
               if(Kreflt .eq. 1 .or. Kreflt .eq. 3)
     .             call ALBEDO(s, Sbcor, rfla, Kreflt)
c
c effect of emitted infrared radiation
c kreflt also signals whether to recalculate orbit elements
               if(Kreflt .eq. 2 .or. Kreflt .eq. 3)
     .             call INFRED(s, Sbcor, rfli, Kreflt)
            endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=2000
c
c        set-up once per  iteration   for partials
c
c        determine if atmospheric drag included
         else if(Kp(82) .gt. 0) then
            if(idrag .ge. 0) then
               at1 = (rsbkm - Rhzkm)/Sh**2
 
c compute partials wrt scale height
               do i = 1, 3
                  dadsh(i) = atdrag(i)*at1
                  end do
c
c           determine if low-thrust force parameters are included in
c           partials
c        calculations done with motion expressions
c
c           determine if distributed asteroids is included in partials
c          calculations done above along with motion expressions
            endif
         endif
      else if(ncall .eq. 0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c*  start=100
c        set-up once per step
         s = Jd + Fract
         if(Ifleng .ne. 0) then
            if(ABS(s-ssavd) .ge. sdel) call SBENG(ssavd, sdel, s)
            ltconz = ltcon*harycn/Sbmass
            ltcon2 = ltconz*Sbarea(1)
            haryc  = harycn/Sbmass
            aucm   = Aultsc*Ltvel*1E5_10
         endif
c opsflg---once per step flag---see statement label=1350
c for once per step procedures that follow the once per
c iteration procedures.
         opsflg = .true.
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=2500
c
c        compute acceleration for motion
      else if(Kkk .gt. 0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c*  start=3000
c
c        compute accelerations for partials
c
c        atmospheric drag
         if(Kp(82) .gt. 0 .and. idrag .ge. 0) then
c
c check if alpha is atmospheric density
            if(Icntrl(Kkk) .eq. 100*Ncentr + 28) then
               Fn(1) = Fn(1) + dadro0(Kk)*atmexp
               return
c
c check if alpha is atmospheric scale height
            else if(Icntrl(Kkk) .eq. 100*Ncentr + 27) then
               Fn(1) = Fn(1) + dadsh(Kk)
               return
            endif
         endif
c
c*  start=3500
c           effect of low-thrust forces on partial derivatives
c          the effect on position and velocity is not included
c        see if partial is w/r/t a probe con
         if(Icntrl(Kkk) .gt. -1) then
c
c*  start=5000
c effect of distributed asteroidal perturbation on partials
c the effect on position and velocity is not included
            if(Nbelt .gt. 0) call ASBFNC(ncall)
         else
            pf(1) = 0._10
            pf(2) = 0._10
            pf(3) = 0._10
            ltfgo = -Icntrl(Kkk)
 
            if(ltfgo .le. 8) then
 
c gas leaks con(1)- con(8)
               if(ltfgo .eq. 2) then
                  pf(2) = fyy
               else if(ltfgo .eq. 3) then
                  pf(3) = fzz1
               else if(ltfgo .eq. 4) then
                  pf(3) = fzz2
               else if(ltfgo .eq. 5) then
                  pf(1) = -delt1*fx
                  pf(2) = delt1*fy
               else if(ltfgo .eq. 6) then
                  pf(1) = -delt1*fx
                  pf(2) = -delt1*fy
                  pf(3) = -delt1*fz1
               else if(ltfgo .eq. 7) then
                  pf(1) = delt1*fx
                  pf(2) = delt1*fy
                  pf(3) = -delt1*fz1
               else if(ltfgo .eq. 8) then
                  pf(3) = -delt1*fz2
               else
                  pf(1) = fxx
               endif
               par   = rmtx(Kk, 1)*pf(1) + rmtx(Kk, 2)*pf(2) + rmtx(Kk,
     .                 3)*pf(3)
               Fn(1) = Fn(1) + par
            else if(ltfgo .gt. 15) then
 
c reflected radiation partials
               if(ltfgo .eq. 16) Fn(1) = Fn(1) + rfla(Kk)
     .             *haryc*Sbarea(1)
               if(ltfgo .eq. 17) Fn(1) = Fn(1) + rfli(Kk)
     .             *haryc*Sbarea(1)
            else
c*  start=4000
c radiation pressure partials con(9)- con(15)
               if(ltfgo .eq. 15) ltfgo = 17
               i = 1 + mod(ltfgo, 3)
               j = ltfgo/3 - 2
 
c pf(i)=solcn(j)
               Fn(1) = Fn(1) + rradpr(Kk, i)*solcn(j)
            endif
         endif
      else
c
c atmospheric drag
         if(Kp(82) .gt. 0) then
            if(idrag .ge. 0) Fn(1) = Fn(1) + atdrag(Kk)
         endif
 
c direct sunlight pressure
         if(Kdirct .gt. 0) Fn(1) = Fn(1) + radprs(Kk)
 
c reflected radiation
         if(Kreflt .gt. 0) then
            if(Kreflt .eq. 1 .or. Kreflt .eq. 3) Fn(1) = Fn(1)
     .          + rfla(Kk)*haryc*Sbarea(1)*Con(16)
            if(Kreflt .eq. 2 .or. Kreflt .eq. 3) Fn(1) = Fn(1)
     .          + rfli(Kk)*haryc*Sbarea(1)*Con(17)
         endif
 
c gas leak
         if(Kp(83) .gt. 0) Fn(1) = Fn(1) + glacc(Kk)
c
c effect of distributed asteroidal perturbation
         if(Nbelt .gt. 0) call ASBFNC(ncall)
      endif
 
      return
      end
