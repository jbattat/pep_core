      subroutine SPOTPL(jdt,frt,nin,ict66,if,nvel)
 
      implicit none

c
c        R. Goldstein   July 1977
c        compar link implementation of planetary rotations
c        Computes the spot coordinates referred to equatorial inertial
c        frame from the body fixed coordinates.
c        Entry SPOTPR computes the partials of these inertial
c        spot coordinates wrt the free parameter of the rotation model
c        (CON(6)-CON(17)).
c        The array ROT is set up  to
c        allow SPOTCD to calculate partials of the inertial spot
c        coords wrt, spot lat, long, and radius.
c
c arguments
      real*10 frt
      integer jdt,nin,ict66,if,nvel

c        JDT       Julian day number of time when coord. desired
c        FRT       fraction of Julian day of same
c        YSPCD     input body fixed spot coords.
c        XSPCD     output equatorial inertial coords. of spot
c        NI        input index (1 or 2) of XSPCD and YSPCD
c        ROT       output rotation matrix
c        ICT66     input I*4 version of ICT(66) documented in PRMRED,
c                  mainly controls debug printout
c        IF        0: interpolation used to calculate INC, PSI,A,B
c                  1: full calculation of above for every observation
c
c        See "The rotation of Mars" by King & Reasenberg for
c        further documentation.
c
c        modified Oct 1978 R.B. Goldstein: pdot and seasonal terms
c             for variation of phi added
c

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'leon.inc'
      include 'mnsprt.inc'
      include 'number.inc'
      include 'obscrd.inc'
      include 'param.inc'
      include 'partcm.inc'
      include 'rotcom.inc'
      real*10 i0rd,psi0rd,sini0,phi0rd,omeg3,cosi0,
     . sn,sj,cn,cj,p,mu,pdotp0,anom
      equivalence (i0rd,Trig(13)),(psi0rd,Trig(16)),
     .            (sini0,Trig(14)),(phi0rd,Trig(19)),
     .            (omeg3,Trig(25)),(cosi0,Trig(15)),(sn,Trig(8)),
     .            (sj,Trig(11)),(cn,Trig(9)),(cj,Trig(12)),
     .            (p,Trig(26)),(mu,Trig(27)),(pdotp0,Trig(28)),
     .            (anom,Trig(29))

c local
      real*10 xsav(3,3),xpsi(3),xinc(3),xphi(3)
      equivalence (xsav(1,1),xpsi),(xsav(1,2),xinc),
     .            (xsav(1,3),xphi)
      real*10 rv(3,3),arm(4),tlim/1.E-5_10/,inc/0._10/,stold/0._10/
      equivalence (arm,sm),(arm(2),cm),(arm(3),s2m),
     .            (arm(4),c2m)
      integer*2 lgo,lpl(1),mpl(1)
      logical*1 lrotdp(12)/12*.false./
      logical   pcall,jump,skip,anok
      real*10 d,psav(3,10),prv(3,12),fudge(12),aui,f1,pbi0,d2
      integer*2 lsav
      integer*2 map(12)/1,2,6,7,3,4,5,8,9,10,11,12/
      character*6 cepch(2)/'1950.0','2000.0'/
      real*10 acf,bcf,c2m,cinc,cm,cphi,cpinc,cpsi,ctmt,dphi,eanom
      real*10 pdpsi,phi,phip,phipm,pinc,pnod,
     .          psi,psip,psipe,s2m,sinc,sm,sphi,spinc,spsi,st
      integer   i,iii,iptang,iptct,iptlog,iptpar,
     .          iptprv,iptrfx,iptrv,iq,ir,irv,jj,k,kick,ni
 
      ni = nin
      if(Nk1.lt.0) then
 
c initialize at beginning of series
         if(Jct(16).ne.0) ctmt = 0._10
         pcall    = .false.
         Rfx(1,1) = cn
         Rfx(1,2) = -sn*cj
         Rfx(1,3) = sn*sj
         Rfx(2,1) = sn
         Rfx(2,2) = cn*cj
         Rfx(2,3) = -cn*sj
         Rfx(3,1) = 0._10
         Rfx(3,2) = sj
         Rfx(3,3) = cj
         iptct     = MOD(ict66/2,2)
         iptang    = MOD(ict66/4,2)
         iptrv     = MOD(ict66/8,2)
         iptlog    = MOD(ict66/16,2)
         iptrfx    = MOD(ict66/32,2)
         iptprv    = MOD(ict66/64,2)
         iptpar    = MOD(ict66/128,2)
         if(iptang.eq.1) then
            cpinc = Coblq1*cj + Soblq1*sj*cn
            spinc = SQRT(1._10 - cpinc**2)
            pinc  = ACOS(cpinc)
            pnod  = ASIN(sj*sn/spinc)
            pdpsi = ASIN(Soblq1*sn/spinc)
            if(Line.gt.57) call OBSPAG
            write(Iout,20) Oblq1,pnod,pinc
   20       format(' SPOTPL: ECLIPTIC OBLIQ,NODE,INC =',3F12.9,
     .             ' (RAD)')
            Line = Line + 1
         endif
         skip = Jct(16).ne.0
         do i = 20,23
            if(Pcond(i,Klan).ne.0._10) then
               skip = .false.
               go to 100
            endif
         end do
      endif
c
c* start=100
c
c        see if the time has changed so that calculation
c        for inc, psi, a and b must be redone
c        tlim is set so that they are calculated once per
c        observation. when they are calculated, it is by
c        interpolation between values calculated at the
c        beginning and end of a "pass"
c
  100 st   = (jdt - Pcom(1)) + (frt - 0.5_10)
      jump = (ABS(st-stold).lt.tlim)
      if(jump .and. (if.eq.0)) go to 300
c
c calculate inc,psi,acf,bcf
c
      stold = st
      anom  = 2.957634951_10 + Mmot*((jdt-2433283) + frt)
      anok  = .false.
      if( .not. skip) then
         sm   = SIN(anom)
         cm   = COS(anom)
         s2m  = 2._10*sm*cm
         c2m  = 1._10 - 2._10*sm*sm
         anok = .true.
         dphi = Pcond(20,Klan)*sm + Pcond(21,Klan)*cm +
     .          Pcond(22,Klan)*s2m + Pcond(23,Klan)*c2m
         if(Jct(16).eq.0) then
            eanom = anom + Pcond(2,Klan)*sm
            ctmt  = 11.4E-3_10*SIN(eanom)/Secday
            if(iptct.ne.0) then
               if(Line.gt.57) call OBSPAG
               write(Iout,510) jdt,frt,anom,eanom,ctmt
  510          format(' SPOTPL: JD,FRACT,ANOM,EANOM,CT-MT(DAYS)=',i8,
     .                1p,4D15.6)
               Line = Line + 1
            endif
         endif
      endif
      call ROTIPS(st,acf,bcf,inc,psi,cinc,sinc,cpsi,spsi)
      if(iptang.eq.1) then
         if(Line.gt.57) call OBSPAG
         psip  = MOD(psi + Pi,Twopi)
         psipe = psip - pdpsi
         if(psipe.lt.0._10) psipe = psipe + Twopi
         write(Iout,250) acf,bcf,inc,psip,psipe
  250    format(' SPOTPL: ACF,BCF,INC,PSIP,PSIPE=',1p,5D18.9)
         Line = Line + 1
      endif
c
c* start=200
  300 st  = st - ctmt
      phi = phi0rd + omeg3*(1._10 + pdotp0*st)*st - psi*cosi0
      if( .not. skip) phi = phi + dphi
      cphi = COS(phi)
      sphi = SIN(phi)
      if(iptang.eq.1) then
         if(Line.gt.57) call OBSPAG
         phip  = phi + Pi
         phipm = MOD(phip,Twopi)
         write(Iout,350) phip,phipm
  350    format(' SPOTPL: PHI=',1pd20.12,' =',0pf14.10)
         Line = Line + 1
      endif
      call ROT3(cphi,sphi,Yspcd(1,ni),xphi)
      call ROT1(cinc,sinc,xphi,xinc)
      call ROT3(cpsi,spsi,xinc,xpsi)
      if(iptrv.eq.1) then
         if(Line.gt.57) call OBSPAG
         write(Iout,400) (Yspcd(k,ni),k = 1,3),xpsi
  400    format(' SPOTPL: X(FIXED),X(1978.0)=',1p,6D15.6)
         Line = Line + 1
      endif
      call PRODCT(Rfx,xpsi,Xspcd(1,ni),3,3,1)
      if(iptrfx.eq.1) then
         if(Line.gt.57) call OBSPAG
         write(Iout,450) cepch(Jct(13)+1),xpsi,(Xspcd(k,ni),k=1,3)
  450    format(' SPOTPL: X(1978.0),X(',a6,')=',1p,6D15.6)
         Line = Line + 1
      endif
c
c        set up rot to send to the rest of pep. rot is
c        needed to calculate partial wrt spot coords. the rot that
c        is needed is the transpose of the rot in the above
c        cited references due to fortran row and col. convention.
c
      rv(1,1) = cpsi*cphi - spsi*cinc*sphi
      rv(1,2) = -cpsi*sphi - spsi*cinc*cphi
      rv(1,3) = -spsi*sinc
      rv(2,1) = cphi*spsi + cpsi*cinc*sphi
      rv(2,2) = -sphi*spsi + cpsi*cinc*cphi
      rv(2,3) = cpsi*sinc
      rv(3,1) = -sinc*sphi
      rv(3,2) = -sinc*cphi
      rv(3,3) = cinc
 
      call PRODCT(Rfx,rv,Rot,3,3,-3)
 
      return
c
c
c
      entry SPOTPR(lgo,kick,lpl,mpl)
c
c* start=600
      if( .not. pcall) then
         if(Index.eq.6) call SUICID(
     .       'NEW ROT. CODE CANNOT CALC TIME DERIV OF PARTIALS',12)
         pcall = .true.
         lsav  = lgo
         call ROTLOG(lpl,lrotdp,mpl,Iabs1)
         if(iptlog.eq.1) then
            if(Line.gt.57) call OBSPAG
            write(Iout,520) lrotdp
  520       format(' SPOTPL: LROTDP=',5x,12L3)
            Line = Line + 1
         endif
c
c        set up fudge factors
c        the minus sign in these factors are due to the following:
c             obs = k*x(s/c-earth)
c             x(s/c-earth) =  terms  -x(s/c-mars)
c        the construction of the full partials are done in cpartc
c
         aui = 1._10/Aultsc
         f1  = aui*Convd
         fudge(1)  = -1._10
         fudge(2)  = -1._10
         fudge(3)  = +1._10
         fudge(4)  = -1._10
         fudge(5)  = -f1
         fudge(6)  = -f1
         fudge(7)  = -aui
         fudge(8)  = -aui
         fudge(9)  = -aui
         fudge(10) = -aui
         fudge(11) = -aui
         fudge(12) = -aui
      endif
c
c
      if(lgo.eq.lsav) then
c
c        calculations of partials of rv
c
c        at this point, the vector lrotdp has been set up which
c        indicates which partials to calculate in the order that
c        they must be calculated. it does not necessarily correspond to
c        lpl since partials wrt some parameters depend on calculations
c        performed when calculating partials wrt other parameters. psav
c        contains saved results that are used in later partials.
c             order of partials: phi0, p, i0, psi0, mu, delta0, alpha0
c
c        psav:results of operating on yspcd with the following operators
c             1: prot3(phi)
c             2: rot1(i)*prot3(phi)
c             3: rot3(psi)*rot1(i)*prot3(phi)
c             4: prot1(i)*rot3(phi)
c             5: rot3(psi)*prot1(i)*rot3(phi)
c             6: prot3(psi)*rot1(i)*rot3(phi)
c
c        from here to fortran statement 210 is a drop through sequence
c        controlled by lrotdp that results in the partial of rv wrt
c        the various rotation model parameters.
c
c* start=700
c        partial of rv(x) wrt phi0
         if( .not. lrotdp(1)) return
         call PROT3(cphi,sphi,Yspcd(1,ni),psav(1,1))
         call ROT1(cinc,sinc,psav(1,1),psav(1,2))
         call ROT3(cpsi,spsi,psav(1,2),psav(1,3))
         do jj = 1,3
            prv(jj,1) = psav(jj,3)
         end do
c
c partial of rv(x) wrt p
         if(lrotdp(2)) then
            d = -st*omeg3/p
            do jj = 1,3
               prv(jj,2) = psav(jj,3)*d
            end do
         endif
c
c partial of rv(x) wrt i0
c this partial modified sept 1977 to include term
c proportional to mu
         if(lrotdp(3)) then
            call PROT1(cinc,sinc,xphi,psav(1,4))
            call ROT3(cpsi,spsi,psav(1,4),psav(1,5))
            d    = psi*sini0
            pbi0 = bcf/sini0/cosi0
            d2   = 1._10 + mu*pbi0
            do jj = 1,3
               prv(jj,5) = d2*psav(jj,5) + d*psav(jj,3)
            end do
         endif
c
c partial of rv(x) wrt psi0
         if(lrotdp(4)) then
            call PROT3(cpsi,spsi,xinc,psav(1,6))
            d = -cosi0
            do jj = 1,3
               prv(jj,6) = psav(jj,6) + psav(jj,3)*d
            end do
         endif
c
c partial of rv(x) wrt mu
         if(lrotdp(5)) then
            d = -acf*cosi0
            do jj = 1,3
              prv(jj,7) = acf*psav(jj,6) + bcf*psav(jj,5) + d*psav(jj,3)
           end do
         endif
c
c partial of rv(x) wrt delta0
c use qmat and saved results
         if(lrotdp(6)) then
            do jj = 1,3
               prv(jj,3) = prv(jj,6)*Qmat(2,1) + prv(jj,5)*Qmat(2,2)
            end do
         endif
c
c partial of rv(x) wrt alpha0 - use qmat and saved results
         if(lrotdp(7)) then
            do jj = 1,3
               prv(jj,4) = prv(jj,6)*Qmat(1,1) + prv(jj,5)*Qmat(1,2)
            end do
         endif
c
c partial of rv(x) wrt pdot
         if(lrotdp(8)) then
            d = -omeg3*st*st/(2._10*p)
            do jj = 1,3
               prv(jj,8) = d*psav(jj,3)
            end do
         endif
c
c* start=1000
c partials of rv(x) wrt ut terms a,b,c,d
         do iq = 9,12
            if( .not. lrotdp(iq)) go to 550
            ir = iq - 8
            if(skip .and. ( .not. anok)) then
               sm   = SIN(anom)
               cm   = COS(anom)
               s2m  = 2._10*sm*cm
               c2m  = 1._10 - 2._10*sm*sm
               anok = .true.
            end if
            do jj = 1,3
               prv(jj,iq) = arm(ir)*psav(jj,3)
            end do
  550    end do
 
         if(iptprv.eq.1) then
            do i = 1,12
               iii = i + 5
               if(lrotdp(map(i))) then
                  if(Line.gt.57) call OBSPAG
                  write(Iout,560) iii,(prv(k,i),k = 1,3)
  560             format(' SPOTPR: PARTIAL RV WRT CON(',i2,')=',1p,
     .                   3D15.6)
                  Line = Line + 1
               endif
            end do
         endif
      endif
c
c* start=1100
c
c        apply rfix and calculate final partial derivatives. fill
c        derem...note...velocity part of derem not calculated as of
c        july, 1977. this implies that only phase delay method of
c        doppler calculation can be used  (ict(21)=-1).
c
      irv = lgo - 5
      call PRODCT(Rfx,prv(1,irv),Derem,3,3,1)
      do i = 1,3
         Derem(i,1) = fudge(irv)*Derem(i,1)
         Derem(i,2) = Derem(i,1)
      end do
      Ivze(1) = 1
 
      if(iptpar.eq.1) then
         if(Line.gt.57) call OBSPAG
         write(Iout,600) cepch(Jct(13)+1),lgo,(Derem(k,1),k = 1,3)
  600    format(' SPOTPR: PARTIAL X(',a4,') WRT CON(',i2,')=',1p,
     .          3D15.6)
         Line = Line + 1
      endif
 
      call CPARTC(kick)
c
c* start=9000
      return
      end
