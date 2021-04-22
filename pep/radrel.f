      subroutine RADREL(ngo)
 
      implicit none
 
c       oct. 1976  r.goldstein
c       effect of general relativity on tmdly and doppler shift
c       this routine replaces an older one of the same name written by
c       mike ash in sept. 1967. the old routine was only good for
c       planetary radar.

c arguments
      integer*4 ngo
c       values of ngo are:
c         1: time delay
c         2: doppler (instantaneous only)
c
c     note: this routine calculates a correction to a value of tmdly
c           that has been calculated and sent through common. for
c           doppler, dop is computed using a relativistic formulation
c           and previous "classical" calculations are overwritten.
c
c       redesigned 1978 june - j.f.chandler
c       now all radar routines should call radrel for both time delay
c       and doppler (skip doppler call if no doppler).  the logic to
c       determine whether apply relativistic corrections is performed
c       here.  some of the setup needed for phase delay doppler is done
c       in radrel(1).
c

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'coord.inc'
      include 'difnct.inc'
      real*10 posr(2)
      equivalence (posr,Rs1m)
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'ltrapobs.inc'
      real*10 tmdly, dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'namtim.inc'
      include 'number.inc'
      include 'obscrd.inc'
      include 'kobequiv.inc'
      include 'param.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'
 
c local
      real*10 psi(2),qr(2,2),dqr(2,2),a1,a2,b1,b2,da1,da2,db1,db2
      equivalence (qr(1,1),a1),(qr(2,1),b1),(qr(1,2),a2),(qr(2,2),b2),
     1 (dqr(1,1),da1),(dqr(2,1),db1),(dqr(1,2),da2),(dqr(2,2),db2)
      real*10 dop1,dop2,dop3,dop4,dop5,dop6,drsc3m,dxp,re3m(2),rs1m2,
     . rs2m2,rsc3m,rsc3m2,tmp,tmq,u,umv,vs1m2,vs2m2,x,ymx,qcr(2,2)
      integer i,jj

c external functions
      real*10 DOT
c
c
c
c
      if(ngo.eq.2) then
c
c
c---------------------------------doppler-------------------------
c        this code is for instaneous doppler only (ict(21)=1).
c        when doing instantaneous doppler, time delay is assumed to be
c        also calculated, and the upper half of this routine (for tmdly)
c        is assumed to have been executed.
c
c
c
         if(Jct(67).le.0) then
            x = DOT(Xemlsc(4,2),Xsitp0(1,2))
            u = DOT(Xsbsun(4),Xsitp0(1,1))
         endif
         umv = -DOT(Xsitep(4,1),Xsitp0(1,1))
         ymx = -DOT(Xsitep(4,2),Xsitp0(1,2))
c
c save umv and ymx for partials calculation
         Quan1(1) = -umv
         Quan1(2) = -ymx
c
c classical calculation of doppler shift
         dop1 = umv + ymx
         if(Jct(67).le.0) then
            dop2 = umv*(ymx + u) - ymx*x
            dop3 = umv*u*(ymx + u) + ymx*x*(x - umv)
            dop1 = dop1 + dop2 + dop3
         endif
         dop4 = 0._10
         if(ndop.ge.0) then
 
c special relativistic corrections
            vs1m2 = DOT(Xemlsc(4,1),Xemlsc(4,1))
            vs2m2 = DOT(Xemlsc(4,2),Xemlsc(4,2))
            dop4  = -.5_10*(vs2m2 - vs1m2)
            if(ndop.gt.0) then
 
c general relativistic corrections
               drsc3m = DOT(Xsbsun(4),Xsbsun(1))/rsc3m
               do i = 1, 2
                  psi(i)    = Gmc2/posr(i)
                  tmp       = drsc3m + DOT(Xemlsc(4,i),Xemlsc(1,
     .                        i))/posr(i)
                  dxp       = DOT(Xsitep(4,i),Xsitp0(1,i))
                  dqr(1,i) = tmp + dxp
                  dqr(2,i) = tmp - dxp
               end do
c
c  the next calculations  give delf/f (doppler).dop1,2,3 are the first,
c  second and third order parts of the 'classical' calc.except for a 3rd
c  order coef.which is included with a gen.rel. term in dop 4.   terms
c  from inclusion of gen.rel.(i.e. grav.effects) are in dop5 and dop6.
c  this expression is acc to order e-14. for e-12 accuracy dop4=dop6=0.
               dop5 = -2._10*(da2/a2 - db2/b2 + da1/a1 - db1/b1)
               dop6 = dop5*dop1
               Raddum(6) = dop5 + dop6
               dop4 = dop4 + .5_10*(psi(2)-psi(1)) + Raddum(6)*Reltrm(2)
            endif
         endif
         dop = dop1 + dop4
c
c
         return
      else
c
c---------------------------time delay--------------------------
c
         if(ntmdly.le.0 .and. ndop.le.0 .and. lopler.ne.-1)
     .       return
c
c unit vectors, e31 and e32 were set up by calling routine
c
         rs2m2  = DOT(Xemlsc(1,2),Xemlsc(1,2))
         Rs2m   = SQRT(rs2m2)
         rs1m2  = DOT(Xemlsc,Xemlsc)
         rsc3m2 = DOT(Xsbsun,Xsbsun)
         Rs1m   = SQRT(rs1m2)
         rsc3m  = SQRT(rsc3m2)
         if(ntmdly.le.0 .and. ndop.le.0) return
         do i = 1, 2
            tmp = posr(i) + rsc3m
            qr(1,i) = tmp + Rsitp(i)
            qr(2,i) = tmp - Rsitp(i)
         end do
         if(ntmdly.le.0) return
         Raddum(1) = 2._10*LOG((a2/b2)*(a1/b1))
         if(MOD(ntmdly/2,2).eq.1) then
c calculate delay due to earth
            do jj=1,2
               tmp = 0._10
               do i = 1, 3
                  tmp = tmp + (Xsitep(i,jj) - Xsite(i,jj))**2
               end do
               re3m(jj) = SQRT(tmp)
            end do
            tmp  = Rsite(2) + re3m(2)
            tmq  = Rsite(1) + re3m(1)
            Raddum(1) = Raddum(1) + 2._10*Mass(3)*(1._10 - Mass(10))
     .                  *LOG(((tmp+Rsitp(2))/(tmp-Rsitp(2)))
     .                  *((tmq+Rsitp(1))/(tmq-Rsitp(1))))
         endif
         if(MOD(ntmdly/4,2).eq.1) then
c calculate delay due to planet
            tmq = SQRT(DOT(Xsbpl,Xsbpl))
            do i=1,2
               tmp = SQRT(DOT(Xsitec(1,i),Xsitec(1,i))) + tmq
               qcr(1,i)=tmp+Rsitp(i)
               qcr(2,i)=tmp-Rsitp(i)
            end do
            Raddum(1)=Raddum(1)+2._10*Mass(Nplnt(Klan))*
     .       LOG((qcr(1,2)/qcr(2,2))*(qcr(1,1)/qcr(2,1)))
         endif

         tmdly = tmdly + Raddum(1)*Reltrm(1)

         if(mod(Jct(6)/2048,2).eq.1) then
            if(Line.gt.56) call OBSPAG
            write(Iout,100) Jd,Dstf(1),Raddum(1)*Reltrm(1)
  100       format(1x,'RADREL: JD.F=',i7,f13.12,' DT=',1pd15.7)
            Line = Line + 1
         endif
         return
      endif
      end
