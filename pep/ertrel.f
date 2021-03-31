      subroutine ERTREL(ncall)
 
      implicit none
c
c r.king   june 1980   subroutine ertrel
c evaluation of general relativity terms in right side of
c equations of motion for the orbit of an earth satrellite.
c called by ertorb and based on m.slade's code in monrel.
c
c arguments
      integer*4 ncall
c         ncall= -2 setup once per iteration of a given step for partial
c              = -1 setup once per iteration of a given step for motion
c              =  0 setup once per step
c              =  k evaluate acceleration for equation k.gt.0
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'intstf.inc'
      include 'sbstuf.inc'
      include 'sbthng.inc'
      real*10 xm(3),xe(3),xme(3)
      equivalence (Bcor,xm),(Ccor,xe),(Sbcor,xme)
c
c quantities internal to this routine
      real*10 c,dot1,dot12,dot2,dot22,dot3,dot4,dot5,
     . dot6,f1,f2,f3,f4,f5,ratme,ratmu,ve2,vevm,vm2,xi,xi2,xi3
      integer   i
      real*10 relfrc(3),vm(3),ve(3),msun,mex,mux
c external functions
      real*10 DOT
c
c setup of  quanities for general relativity accel's
      if(ncall.le.0) then
c
c     calculation of post-newtonian forces for the earth-moon system
c     from t.r. 425,pp. 28-29
c     lengths in a.u., time in days
c
         c = SQRT(Cvel2)
         do i = 1,3
            ve(i) = Ccor(i + 3)/c
            vm(i) = Bcor(i + 3)/c
         end do
         vevm  = DOT(ve,vm)
         ve2   = DOT(ve,ve)
         vm2   = DOT(vm,vm)
         xi    = Rsb
         xi2   = xi**2
         xi3   = xi2*xi
         dot1  = ve(1)*xme(1) + ve(2)*xme(2) + ve(3)*xme(3)
         dot2  = vm(1)*xme(1) + vm(2)*xme(2) + vm(3)*xme(3)
         dot3  = xme(1)*xe(1) + xme(2)*xe(2) + xme(3)*xe(3)
         dot4  = xme(1)*xm(1) + xme(2)*xm(2) + xme(3)*xm(3)
         dot5  = xm(1)*vm(1) + xm(2)*vm(2) + xm(3)*vm(3)
         dot6  = xe(1)*ve(1) + xe(2)*ve(2) + xe(3)*ve(3)
         dot12 = dot1**2
         dot22 = dot2**2
 
c masses in a.u.
         msun = Gamat/c**2
 
c the following two cards define the earth orbiter case:
         ratme = Mass1(3)
         ratmu = 0._10
         mex   = msun*ratme
         mux   = msun*ratmu
         f1    = ((4.*msun+mux)/Rb+(4.*mex+3.5*mux)/xi+mex/Rc - vm2)/Rb3
         f2    = ((4.*msun+mex)/Rc+(4.*mux+3.5*mex)/xi+mux/Rb - ve2)/Rc3
         f3    = (2.*(2.*ratme+ratmu)*(2.*mux+mex)/xi + (4.*mex+mux)
     .           /Rb + (4.*mux+mex)/Rc - (2.*ratme+ratmu)
     .           *ve2 - (2.*ratmu+ratme)*vm2 + 4.*(ratme+ratmu)
     .           *vevm + 1.5*ratme*dot12/xi2 + 1.5*ratmu*dot22/xi2 -
     .           0.5*mex*dot3/Rc3 + 0.5*mux*dot4/Rb3)/xi3
         f4    = 4.*dot5/Rb3 + (4.*ratme + 3.*ratmu)*dot2/xi3 -
     .           (3.*ratme + 4.*ratmu)*dot1/xi3
         f5    = 4.*dot6/Rc3 + (4.*ratme + 3.*ratmu)*dot2/xi3 -
     .           (3.*ratme + 4.*ratmu)*dot1/xi3
c-----------------------------------------------------------------------
c
c
c effect of general relativity on the motion of the satellite
      else if(Kkk.le.0) then
         relfrc(Kk) = (xm(Kk)*f1 - xe(Kk)*f2 + xme(Kk)*f3 + vm(Kk)
     .                *f4 - ve(Kk)*f5)
         Fn(1) = Fn(1) + Gamat*Rlfact*relfrc(Kk)
c
c
c effect of general relativity on partial derivatives
      else if(Icntrl(Kkk).eq.31) then
         Fn(1) = Fn(1) + Gamat*relfrc(Kk)
      else if(Icntrl(Kkk).eq.32) then
         Fn(1) = Fn(1) + Gama*Tvary*relfrc(Kk)
      endif
 
      return
      end
