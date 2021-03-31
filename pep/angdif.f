      subroutine ANGDIF(kathy)
 
      implicit none
 
c*** start of declarations inserted by spag
c*** end of declarations inserted by spag
 
 
c j.f.chandler - 1979 february - subroutine angdif
c calculate optical difference observables

c arguments
      integer*4 kathy
c          kathy= 1  photographic observation (topocentric, referred to
c                    equinox-equator of date, no aberration correction)
c          kathy= 2  photographic observation (topocentric with
c                    all aberration removed referred to the mean
c                    equinox and equator of ref. epoch) (astrographic)

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdateq.inc'
      include 'coord.inc'
      include 'coordxeq.inc'
      include 'coordoeq.inc'
      real*10 xpta(3),xptd(3),xz(3)
      equivalence (Raddum,xpta),(Raddum(4),xptd)
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'ltrapobs.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*4    acctim,accdst,accprc
      equivalence (acctim,Estf),(accdst,Estf(2)),(accprc,Estf(3))
      real*10 dist,distz
      equivalence (dist,Dstf(4)),(distz,Dstf(5))
      integer*2 itime,ntime
      equivalence (Kob(17),itime),(Kob(18),ntime)
c nice=-1 observation of right ascension only
c nice= 0 observation of right ascension and declination
c nice= 1 observation of declination only
      include 'jdfrequv.inc'
      include 'param.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'

c external functions
      real*10 DOT

c local
      real*10 adb,dum,frz,oa,od,off(3),opu2,pu,qxn2,
     .          refa,refr,refs,rfc,rfp,rltsc,rxa,rz,rzl,rzsq,xstp(3)
      integer   i,iter,iter1,j,jdz,ksu,lsw,lswp
      real*10 xdel(3),xalp(3),xpol(3)
c
c
c begin iteration to get light time correction
c for second body
      lswp = 1
      lsw  = 1
      iter = 0
  100 iter = iter + 1
      call TIMINC(Jd,fr(1),jdz,frz,-distz/Secday)
      if(Klans1.ne.0) then
         call SCTRP(1,jdz,frz,0,lsw,1)
         if(jdz.le.0) then
            Jd = 0
            return
         elseif( Ncs1.le.0 ) then
            goto 200
         endif
      endif
c
c second body requires "klan" position
      call PLTRP(1,jdz,frz,0,lswp)
      if(jdz.le.0) then
         Jd = 0
         return
      else
         lswp = -1
         if(Klans1.eq.0) then
 
c second body is central body of first
            do i = 1,3
               Xsc(i,1) = Xp(i)
            end do
         else
 
c second body is satellite of planet 'klan'
            do i = 1,3
               Xsc(i,1) = Xsc(i,1) + Xp(i)
            end do
         endif
      endif
 
c get vector to second body
  200 do i     = 1,3
         xz(i) = Xsc(i,1) + x(i,2) - Xsitau(i)
      end do
      rzsq = DOT(xz,xz)
      rz   = SQRT(rzsq)
      rzl  = rz*Aultsc
 
c determine if further iteration is needed
      if(ABS(distz - rzl).gt.accdst) then
 
c set up for next iteration
         distz = rzl
         lsw   = -1
         if(iter.lt.20) goto 100
         call SUICID(
     .'MORE THAN 20 ITERATIONS NEEDED FOR LIGHT TIME, STOP IN ANGDIF   '
     .,16)
      endif
c
c end of second body iteration
      Nit(19) = Nit(19) + iter
      distz   = rzl
c ought to save planet pqind
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c begin iteration for observed body
      lsw   = 1
      iter1 = 0
      do while( .true. )
         iter1 = iter1 + 1
         call TIMINC(Jd,fr(1),Jdy,fr(3),-dist/Secday)
         if(Klan.gt.0 .and. Klan.le.u_mxpl) then
c
c read planet tape, interpolate to get planet position
            call PLTRP(1,Jdy,fr(3),0,lswp)
            if(Jdy.le.0) then
               Jd = 0
               return
            else
 
c get planet relative to earth
               do j = 1,3
                  x(j,4) = Xp(j) + x(j,2)
               end do
            endif
         endif
c
c determine if probe or satellite is observed body
         if(Klanb.gt.0) then
c
c determine satellite position
            call SBTRP(1,Jdy,fr(3),0,lsw)
            if(Jdy.le.0) then
               Jd = 0
               return
c
c determine probe or satellite relative to earth
            elseif( Ncp0.ne.3 ) then
               ksu = 1
               if(Ncp0.ne.10) then
                  ksu = 4
                  if(Ncp0.eq.0) ksu = 2
               endif
               do j = 1,3
                  x(j,1) = Xsb(j)*Cmfct + x(j,ksu)
               end do
            else
               do j = 1,3
                  x(j,1) = Xsb(j)
               end do
            endif
         endif
c
c topocentric correction
         do j = 1,3
            x(j,Kst) = x(j,Kst) - Xsitau(j)
         end do
c
c see if there is to be further iteration
         r     = SQRT(DOT(x(1,Kst),x(1,Kst)))
         rltsc = r*Aultsc
         if(ABS(dist - rltsc).gt.accdst) then
c
c set up quantities for next iteration
            dist = rltsc
            lsw  = -1
            lswp = -1
            if(iter1.lt.20) goto 300
            call SUICID(
     .'MORE THAN 20 ITERATIONS NEEDED FOR LIGHT TIME, STOP IN ANGDIF   '
     .,16)
         endif
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c termination calculations
         Nit(20) = Nit(20) + iter1
c
c test for dummy observation below horizon
         if(Idumob.eq.1) then
            do i = 1,3
               xstp(i) = -x(i,Kst)
            end do
            call HORIZN(Xsitau,xstp,1,Jd)
            if(Jd.le.0) return
         endif
         goto 400
  300 end do
 
  400 dist = rltsc
 
c aberration correction
      do j = 1,3
         if(kathy.ne.2) then
            x(j,Kst) = x(j,Kst) + r*Ves(j)
            xz(j)     = xz(j) + rz*Ves(j)
         endif
         xo(j) = x(j,Kst)
         Xsitep(j,1) = -xo(j)
         Xsitep(j,2) = -xz(j)
      end do
 
      call UVECTR(3,Xsitep(1,1),Rsitp(1),Xsitp0(1,1),dum)
      call UVECTR(3,Xsitep(1,2),Rsitp(2),Xsitp0(1,2),dum)
 
c get apparent offset vector
      adb = DOT(Xsitp0(1,1),Xsitp0(1,2))
      do i = 1,3
         off(i) = Xsitp0(i,2) - Xsitp0(i,1)/adb
      end do
c
c correct for differential refraction
      if(kathy.ne.2) then
         call REFRCT(rfc,1)
         qxn1 = -DOT(Sitnrm,Xsitp0(1,2))
         do i = 1,3
            xnp(i) = Sitnrm(i,1) + Xsitp0(i,2)*qxn1
         end do
         qxn2 = qxn1**2
         refr = rfc/qxn2/(1._10 - qxn2)
         refs = refr*DOT(Sitnrm,off)
         do i = 1,3
            off(i) = off(i) - xnp(i)*refs
         end do
 
c decide if equator of date or of reference epoch
         if(kathy.ne.2) then
 
c using equator of date - get from nutpr
            do i = 1,3
               xpol(i) = Nutpr(3,i)
            end do
            call CROSS(Xsitp0(1,2),xpol,xalp)
            refa = refr*DOT(Sitnrm,xalp)
            do i = 1,3
               xalp(i) = xalp(i) - xnp(i)*refa
            end do
            rxa = SQRT(DOT(xalp,xalp))
            do i = 1,3
               xalp(i) = xalp(i)/rxa
            end do
            call CROSS(xalp,Xsitp0(1,2),xdel)
            goto 500
         endif
      endif
 
c using equator of reference epoch and no refraction
      xpol(1) = 0._10
      xpol(2) = 0._10
      xpol(3) = 1._10
      pu   = Xsitp0(3,2)
      opu2 = SQRT(1._10 - pu**2)
 
c get meridional unit vector
      do i = 1,3
         xdel(i) = (xpol(i) - pu*Xsitp0(i,2))/opu2
      end do
      call CROSS(Xsitp0(1,2),xdel,xalp)
  500 if(Nice.le.0) then
 
c calculate right ascension offset
         oa = DOT(off,xalp)
         Deriv(2,1) = oa/Convds
         Angdum(8)   = Deriv(2,1)
 
c apply plate scale bias, if any
         if(Neqnox.gt.0) Deriv(2,1) = Deriv(2,1)*Sclplt
      endif
      if(Nice.ge.0) then
 
c calculate declination offset
         od = DOT(off,xdel)
         Deriv(2,2) = od/Convds
         Angdum(9)   = Deriv(2,2)
 
c apply plate scale bias, if any
         if(Neqnox.gt.0) Deriv(2,2) = Deriv(2,2)*Sclplt
      endif
      if(Nphase.gt.0) call DPHCOR(kathy - 1)
c
c decide if partials are to be computed
      if(Ict(3).ge.0 .or. Idumob.ne.1) then
         if(Ict(1).gt.0) then
            rr1 = -Rsitp(1)*adb*Convds
            if(Neqnox.gt.0) rr1 = rr1/Sclplt
 
c set up for ra partial
            if(Nice.le.0) then
               do i = 1,3
                  xpta(i) = xalp(i) + oa*Xsitp0(i,2)
               end do
               if(kathy.ne.2) then
                  rfp = refr*DOT(Sitnrm,xalp)
                  do i = 1,3
                     xpta(i) = xpta(i) - rfp*xnp(i)
                  end do
               endif
            endif
 
c set up for dec partial
            if(Nice.ge.0) then
               do i = 1,3
                  xptd(i) = xdel(i) + od*Xsitp0(i,2)
               end do
               if(kathy.ne.2) then
                  rfp = refr*DOT(Sitnrm,xdel)
                  do i = 1,3
                     xptd(i) = xptd(i) - rfp*xnp(i)
                  end do
               endif
               cake = rr1
            endif
         endif
      endif
 
      return
      end
