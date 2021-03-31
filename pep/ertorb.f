      subroutine ERTORB(ncall)
 
      implicit none
c
c        ash / amuchastegui - october 1969 - subroutine ertorb
c     computations for forces on earth satellite motion due to
c     earth shadow, radiation pressure (direct and reflected),
c     air drag, charge drag, satellite thrustors and sensors, etc.
c
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
c        common
      include 'cnthar.inc'
      include 'engstf.inc'
      include 'gleaks.inc'
      include 'intstf.inc'
      include 'lothrf.inc'
      include 'output.inc'
      include 'petuna.inc'
      real*10 radcon(3)
      equivalence (radcon,Con(9))
      include 'sbstuf.inc'
      include 'sbthng.inc'
c
c quantities internal to this routine
      real*10 haryc,psum,s,solcon,v,ymag
      integer   i,idrag,j,k,kj,ll
      logical*4 opsflg
      real*10 rvec(3),rfla(3),rfli(3),ltcon2,radprs(3),
     .          vvec(3),ltconz,svec(3),yvec(3),qrad(3,3),q(3,3)
      real*10 ssavd/0._10/, sdel/0._10/
 
c conversion between cm/sec**2 to a.u./day**2
      real*10 harycn/4.9900168876E-4_10/
 
c ltcon- solar intensity at 1 a.u. (dyne/cm**2 )
      real*10 ltcon/4.65E-5_10/

c external functions
      real*10 DOT
c
c is this les-8 or les-9
      if(Kkp(91).ne.0) then
         call LESTHP(ncall)
 
      else if(ncall.lt.0) then
c
c
c-----------set up once per iteration of a given step for motion------
c
         if(ncall.lt.-1) then
c
c
c-----------set up once per iteration of a given step for partials------
c
            if(Kkp(78).gt.0) call NDSRAD(ncall,solcon,Sbcor,Rsb,
     .          Bcor, Rb, radcon, Lambda, radprs, qrad, q)
         else
c
c set up unit vectors for gas leak ,rad. pres., and drag
            if(Kp(83).gt.0 .or. Kdirct.gt.0 .or. Kp(82).ne.0)
     .          then
 
               if(Ncentr.eq.0)
     .              call SUICID('NCENTR=0 NOT ALLOWED IN ERTORB, STOP',
     .          9)
               svec(1) = Bcor(1)/Rb
               svec(2) = Bcor(2)/Rb
               svec(3) = Bcor(3)/Rb
               if(Kdirct.gt.0) then
c
c        radiation pressure setup
c
c     kp(81)= kreflt*10+ kdirct
c          kdirct=0 no radiation pressure
c                =1 direct radiation,exact shadow
c                =2 direct radiation,approx. shadow.(not coded)
c                =3 direct radiation, no shadow
c
c          kreflt=1 albedo radiation
c                =2 infared radiation
c                =3 both inrared and albedo
c
c
                  Lambda = 1._10
                  if(Ncentr.gt.0) then
                     if(Kdirct.ne.3) then
                        if(Kdirct.eq.2) then
 
c simple on-off shadow avoid call to shadow routine
                           if(Rb.le.Rc) goto 10
                           if(Rsb2 - DOT(Sbcor,svec)**2.lt.Crad**2)
     .                         Lambda = 0._10
                        else
                           call SHADOW
                        endif
 
c end simple shadow check
                        if(Lambda.le.0._10) then
                           do i = 1, 3
                              radprs(i) = 0._10
                           end do
                           goto 20
                        endif
                     endif
                  endif
   10             solcon = ltcon2*Lambda/Rb2
c
c spherical model
                  if(Kkp(78).gt.0) then
c
c nds model for gps satellites
                     call NDSRAD(ncall,solcon,Sbcor,Rsb,Bcor,Rb,
     .                           radcon, Lambda, radprs, qrad, q)
                  else
                     do i = 1, 3
                        radprs(i) = Con(9)*solcon*svec(i)
                     end do
 
c radiation pressure may have ad hoc y-axis component
                     do i = 1, 3
                        rvec(i) = Sbcor(i)/Rsb
                     end do
                     call CROSS(svec,rvec,yvec)
                     ymag = SQRT(DOT(yvec,yvec))
                     do i = 1, 3
                        yvec(i)   = yvec(i)/ymag
                        radprs(i) = radprs(i) + Con(10)*solcon*yvec(i)
                     end do
                  endif
               endif
c
c setup for atmospheric drag
   20          if(Kp(82).gt.0) idrag = 0
            endif
c
c setup for ad hoc along-track acceleration
            if(Kp(84).gt.0) then
               v = SQRT(DOT(Sbcor(4),Sbcor(4)))
               do i = 1, 3
                  vvec(i) = Sbcor(i + 3)/v
               end do
            endif
c
c calculate vectors for general relativity
            if(Kp(61).ge.0) call ERTREL(ncall)
c
c once per step procedures that follow the once per iteration
c procedures during the first iteration
c
            if(opsflg) then
               opsflg = .false.
 
c effect of visible reflected light
               if(Kreflt.eq.1 .or. Kreflt.eq.3)
     .             call ALBEDO(s,Sbcor,rfla,Kreflt)
c
c effect of emitted infrared radiation
c kreflt also signals whether to recalculate orbit elements
               if(Kreflt.eq.2 .or. Kreflt.eq.3)
     .             call INFRED(s,Sbcor,rfli,Kreflt)
            endif
         endif
      else if(ncall.eq.0) then
c
c
c-----------set up once per step----------------------------------------
c
         s = Jd + Fract
         if(Ifleng.ne.0) then
            if(ABS(s-ssavd).ge.sdel) call SBENG(ssavd,sdel,s)
            ltconz = ltcon*harycn/Sbmass
            ltcon2 = ltconz*Sbarea(1)
            haryc  = harycn/Sbmass
         endif
c
c
c        opsflg---once per step flag---see statement label=380
c        for once per step procedures that follow the once per
c        iteration procedures.
         opsflg = .true.
      else
c
c
c-----------compute acceleration for motion------------------------
c
         k = ncall
         if(Kkk.le.0) then
c
c atmospheric drag
            if(Kp(82).gt.0) then
               if(idrag.lt.0) then
               endif
            endif
c
c direct sunlight pressure
            if(Kdirct.gt.0) Fn(1) = Fn(1) + radprs(Kk)
c
c reflected radiation
            if(Kreflt.gt.0) then
               if(Kreflt.eq.1 .or. Kreflt.eq.3) Fn(1) = Fn(1)
     .             + rfla(Kk)*haryc*Sbarea(1)*Con(16)
               if(Kreflt.eq.2 .or. Kreflt.eq.3) Fn(1) = Fn(1)
     .             + rfli(Kk)*haryc*Sbarea(1)*Con(17)
            endif
c
c general relativity
            if(Kp(61).ge.0) call ERTREL(ncall)
c
c ad hoc along-track acceleration
            if(Kp(84).gt.0) Fn(1) = Fn(1) + vvec(Kk)*Con(18)
c
c
c-----------compute acceleration for partials------------------------
c
c
c         direct radiation pressure
         else if(Kdirct.gt.0) then
            if(Kkp(78).gt.0) then
c
c nds gps model
               if(Icntrl(Kkk).le.-9 .and. Icntrl(Kkk).ge.-11)
     .             then
                  ll   = -Icntrl(Kkk) - 8
                  psum = 0._10
                  do j = 1, 3
                     kj   = k - 4 + j
                     psum = psum + q(ll,j)*Y(kj,1)
                  end do
                  Fn(1) = Fn(1) + qrad(Kk,ll) + psum
                  return
               endif
c
c check if partial w.r.t. overall coefficient
            else if(Icntrl(Kkk).ne.-9) then
c
c check if partial w.r.t. y-axis acceleration coefficient
               if(Icntrl(Kkk).eq.-10) then
                  Fn(1) = Fn(1) + solcon*yvec(Kk)
                  return
               endif
            else
               Fn(1) = Fn(1) + solcon*svec(Kk)
               return
            endif
c
c
c ad hoc along-track accleration
            if(Kp(84).gt.0) then
               if(Icntrl(Kkk).eq.-18) Fn(1) = Fn(1) + vvec(Kk)
            endif
         endif
      endif
c
c
      return
      end
