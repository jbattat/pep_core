      subroutine PLCMC(jd,fr,xp,nv)
 
      implicit none

c        Subroutine PLCMC - J.F.Chandler - 1977 December 1
c        Compute corrections to coordinates in 'XP' for planet center of
c        mass instead of system center -- type of coordinate indicated
c        by iabs(NV)
c        NV=0 - position, NV=1 - velocity, NV=2 - acceleration
c        If NV=0 or 1, then XP is a 6-vector of position and velocity,
c        but if NV=2, then XP is 3-vector of acceleration
c        Controls set up in PLCMS (called from CMPAR3)

c parameters
      integer*4 jd,nv
      real*10 fr,xp(6)

c JCT(55) is a packed bits indicator for the planet center of mass
c         computations
c JCT(55)=0  no center of mass corrections
c    1 bit = 1: get satellite coordinates from interpolation if possible
c    2 bit = 1: get satellite coordinates from elliptic if necessary
c    4 bit = 1: use only satellites that are on s-body tape
c    8 bit = 1: for elliptic method, obtain elements from s-body tape
c                (if possible).  Warning: this option should normally
c                be set because the initial conditions in /EMPCND/ are
c                likely to have been overwritten from tape anyway and
c                may not match the epoch in CON1(1).
c   16 bit = 1: allow calculation of barycenter of planet+observed
c                satellite instead of planet
c   32 bit = 1: use precessing elliptic orbit instead of stationary

c array dimensions
      include 'globdefs.inc'

c common 
      include 'b2dta.inc'
      include 'b2ydta.inc'
      include 'bddta.inc'
      include 'bdydta.inc'
      include 'cmcke.inc'
      include 'comdat.inc'
      real*10 cmfct
      equivalence (Comcon(128),cmfct)
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'param.inc'
      include 'tabval.inc'
c local
      integer*2 kpl,np
      integer*4 i,inv,isgn,it,j,jf,jff,jy,k,n1,n2
      real*10 dimx,dum,pva(3),p(4),ry2,ry3,t,ti
c external functions
      real*10 D2TRPF,DTRPF,TERPF
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(Numpcm.le.0) return
      t   = jd + fr
      inv = -iabs(nv)
      do i = 1, Numpcm
         kpl = Kplcm(i)
         if(kpl.lt.0) then
            np = -kpl
         else
            np = Nplnt(kpl)
         endif

c skip correction for observed satellite
         if(np.ne.Nplnt0 .or. cmfct.eq.1._10) then
 
c if 'initial epoch' is non-zero, then do elliptic
            if(Tplc(i).le.0._10) then
 
c interpolate for coordinates from tape
               dimx  = Pcintx(i)
               Ntab1 = 1
               Ntab2 = 2
 
c test whether same tabular interval as last call
               ti = ((jd-2441000) + fr - Tlcm(i))/dimx
               if(ti.ge.-1._10 .and. ti.lt.2._10) then
 
c within one interval, can save old y vectors
                  it = (ti + 2._10)
                  it = it - 2
                  if(it.lt.0) then
c adjacent tabular interval, shift y vectors
c
c backwards in time
                     n2    = 2
                     n1    = 1
                     Ntab2 = 1
                  else if(it.eq.0) then
 
c same tabular interval as before
                     p(1) = ti
                     goto 10
                  else
 
c forwards in time
                     n2    = 1
                     n1    = 2
                     Ntab1 = 2
                  endif
 
c move y vector
                  do j = 1, 5
                     do k = 1, 3
                        Yplcm(j,n2,k,i) = Yplcm(j,n1,k,i)
                     end do
                  end do
               endif
               if(np.gt.10) then
                  Tgo(5) = 0._10
                  call B2REED(jd,fr,np,5)
                  isgn = Ib2sgn
               else
                  call BDREED(jd,np,6)
                  isgn = Ibdsgn
               endif
               if(jd.le.0) then
 
c error return
                  Tlcm(i) = 1E10_10
                  jd = 0
                  return
               endif
 
c figure index to desired tape values
               if(isgn.lt.0) then
                  jff     = Idxb2 + 10
                  Tlcm(i) = Ttb2(2) - 2441000._10 - dimx
               else
                  jff     = Idxb2 + 1
                  Tlcm(i) = Ttb2(2) - 2441000._10
               endif
               p(1) = ((jd-2441000) + fr - Tlcm(i))/dimx
               do jy = 1, 3
                  jf = jff
                  do j = 1, 10
                     if(np.gt.10) then
                        Tabvl(j) = Bod1(jy,jf,1)
                     else
                        Tabvl(j) = Merc(jy,jf)
                     endif
                     jf = jf + isgn
                  end do
                  call YCOFF(Yplcm(1,1,jy,i))
               end do
   10          p(2) = p(1)**2
               p(3) = 1._10 - p(1)
               p(4) = p(3)**2
               do j = 1, 3
                  if(inv.eq.0) then
                     pva(j) = TERPF(p,Yplcm(1,1,j,i))
                  else if(inv.eq.-1) then
                     pva(j) = DTRPF(p,Yplcm(1,1,j,i),dimx)
                  else if(inv.eq.-2) then
                     pva(j) = D2TRPF(p,Yplcm(1,1,j,i),dimx)
                  endif
               end do
            else if(Jflg(6)) then
               call PLIPT(t-Tplc(i),Elpcm(1,i),inv,pva,Rycm(i))
            else
               call JLIPT(t-Tplc(i),Elpcm(1,i),inv,pva,Rycm(i),
     .          ry2,ry3,dum)
            endif
 
c add on scaled satellite coordinates
            do j = 1, 3
               if(inv.ge.-1) then
                  k = j-3*inv
                  Xpcm(k,i) = pva(j)
               else
                  k = j
               endif
               xp(k) = xp(k) - pva(j)*Mass(np)
            end do
            if(Plcdbg(1-inv)) then
               if(Line.gt.56) call OBSPAG
               write(Iout,20) np,jd,fr,inv,pva
   20          format(' PLCMC(', i2, '): JD.F=', i7, f13.12, ' NV=',
     .                i2, 8x, 'X=', 1p3D23.15)
               Line = Line + 1
            endif
         endif
      end do
      return
      end
