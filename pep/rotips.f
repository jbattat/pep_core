      subroutine ROTIPS(st,a,b,inc,psi,cinc,sinc,cpsi,spsi)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   j
 
c*** end of declarations inserted by spag
 
 
      real*10 st, a, b, inc, psi, cinc, sinc, cpsi, spsi
c
c        r. b. goldstein  september 1977
c        input:  st in days from con1(1)
c        output: a,b,inc,psi and trig functions of inc,psi
c        algorithim:    calls rotab to calculate a,b at beginning and
c                       end of a "pass". all values of all returned
c                       parameters within the pass are obtained from
c                       linear interpolation.
c
c
      real*10 tsav/0._10/, as(2),bs(2),is(2),ps(2),cis(2),sis(2),
     .          cps(2),sps(2),ma,mb,mi,mp,mci,msi,mcp,msp,
     .          tsav0, tlim/0.5_10/, tl2/0.25_10/, tend, t
 
      include 'rotcom.inc'
      real*10 i0rd, psi0rd, mu
      equivalence (i0rd,Trig(13)),(psi0rd,Trig(16)),(mu,Trig(27))
 
      real*10 LINE, x, m, bd
      LINE(x,m,bd) = m*x + bd
c
c check if new endpoints are needed
c
      if(ABS(st-tsav).gt.tl2) then
         tsav0 = st
         tsav  = st + tl2
         tend  = st + tlim
         call ROTAB(as(1),bs(1),st)
         call ROTAB(as(2),bs(2),tend)
 
         do j = 1, 2
            is(j)  = i0rd + mu*bs(j)
            ps(j)  = psi0rd + mu*as(j)
            cis(j) = COS(is(j))
            sis(j) = SIN(is(j))
            cps(j) = COS(ps(j))
            sps(j) = SIN(ps(j))
         end do
 
         ma  = (as(2) - as(1))/tlim
         mb  = (bs(2) - bs(1))/tlim
         mi  = (is(2) - is(1))/tlim
         mp  = (ps(2) - ps(1))/tlim
         mci = (cis(2) - cis(1))/tlim
         msi = (sis(2) - sis(1))/tlim
         mcp = (cps(2) - cps(1))/tlim
         msp = (sps(2) - sps(1))/tlim
      endif
c
c
c
      t    = st - tsav0
      a    = LINE(t,ma,as(1))
      b    = LINE(t,mb,bs(1))
      inc  = LINE(t,mi,is(1))
      psi  = LINE(t,mp,ps(1))
      cinc = LINE(t,mci,cis(1))
      sinc = LINE(t,msi,sis(1))
      cpsi = LINE(t,mcp,cps(1))
      spsi = LINE(t,msp,sps(1))
 
      return
      end
