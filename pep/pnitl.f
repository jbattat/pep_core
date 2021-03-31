      subroutine PNITL(icnd1,cond1,cond2,setp)
 
      implicit none

c           pnitl - j.f.chandler - 1984 august
c           based on subroutine sate from duxbury and callahan 1982
c           compute position of satellite using precessing elliptic
c           elements in standard frame
c
c     setup is performed by calling pnitl:
c arguments 
c icnd1 type of initial conditions (should be 1, for longitudes)
c cond1 initial conditions (ignoring precession and inclination to the
c       laplacian plane)
c cond2 inclination to laplacian plane, angle rates, etc.
c setp  array to be filled with setup quantities
c
c     evaluation is performed by calling plipt:
c dt    time from initial epoch in days
c x     returned array of position or velocity
c ic    indicator: 0 -> position, -1 -> velocity
c
      integer*2 icnd1
      real*10 cond1(6),cond2(10),setp(33),dt,x(3)
      integer*4 ic
c
c commons
      include 'funcon.inc'
c
c local
      real*10 vbar(2),j,long,n,k,sm(3),cm(3)
      real*10 ang, anom, cj, ck, cn, cnck, cnsk, cv, e2, e3,
     .          oecv, omeg, per, r, rdotr, sj, sk, sn
      real*10 snck, snsk, sv, v, vdot, w
      integer   i

c external function
      real*10 DOT

c        use of array 'setp' --- storage for saved quantities from
c        setting up.  it is effectively equivalenced to the following:
c      b(3,2),tsv,a,e,anom0,motn,perrat,krat,nrat,jrat,mot2,fdga,fdgp,
c      1      7   8 9  10    11    12    13   14   15   16   17   18
c
c      k0,fct1,fct2,fct3,j0,n0,per0,quan1,nrfc,dvdm,xbar,ybar,rdfc,
c      19  20   21  22   23 24  25    26   27   28   29   30   31
c
c      sinc,cinc
c       32   33
c
c note: 'j' and 'n' are the node and inclination of the laplacian plane.
c       'i' is the inclination of the orbit to the laplacian plane.
c       'a' and 'e' are the usual semimajor axis and eccentricity.
c       'k' is the argument of node, not its longitude, while 'per' is
c       the longitude of periapse and 'long' is the mean longitude.
c       'fdga' and 'fdgp' represent a solar perturbation term for
c       deimos (wilken's theory)
c
c-----------------------------------------------------------------------
c           s e t t i n g   u p
c
c           get a,e,i from array
      setp(8)  = cond1(1)
      setp(9)  = cond1(2)
      ang      = cond2(1)*Convd
      setp(33) = COS(ang)
      setp(32) = SIN(ang)
 
c initialize 'previous time'
      setp(7) = 1E50_10
 
c get other angles from array
      do i = 2, 10
         setp(i+9) = cond2(i)*Convd
      end do
      do i = 3, 5
         setp(i+20) = cond1(i)*Convd
      end do
      setp(10) = cond1(6)*Convd
 
c convert arguments to longitudes, if necessary
      if(icnd1.ne.1) then
         setp(25) = setp(25) + setp(24)
         setp(10) = setp(10) + setp(25)
      endif
 
c set up true anomaly calculation
      e2 = setp(9)*setp(9)
      setp(26) = 1._10 - e2
      e3 = e2*setp(9)
      setp(20) = 2._10*setp(9) - e3/4._10
      setp(21) = 1.25_10*e2
      setp(22) = 13._10/12._10*e3
      return
c-----------------------------------------------------------------------
c e v a l u a t i o n
c
      entry PLIPT(dt,setp,ic,x,r)
 
      if(dt.ne.setp(7)) then
         setp(7) = dt
c
c get node, mean longitude and periapse (deg)
         k = MOD(setp(19) + dt*setp(13),Twopi)
         if(k.lt.0._10) k = k + Twopi
         ck   = COS(k)
         sk   = SIN(k)
         long = MOD(setp(10) + dt*(setp(11)+dt*setp(16)) + setp(17)
     .          *SIN(k-setp(18)),Twopi)
         per  = MOD(setp(25) + dt*setp(12),Twopi)
 
c compute laplacian pole (n=90+ra, j=90-dec)
         n  = MOD(setp(24) + dt*setp(14),Twopi)
         cn = COS(n)
         sn = SIN(n)
         j  = setp(23) + dt*setp(15)
         cj = COS(j)
         sj = SIN(j)
c
c compute mean and true anomalies in radians
         anom  = long - per
         sm(1) = SIN(anom)
         cm(1) = COS(anom)
         do i = 2, 3
            sm(i) = sm(1)*cm(i - 1) + cm(1)*sm(i - 1)
            cm(i) = cm(1)*cm(i - 1) - sm(1)*sm(i - 1)
         end do
         v = anom + DOT(setp(20),sm)
         setp(28) = 1._10 + setp(20)*cm(1) + 2._10*setp(21)*cm(2)
     .              + 3._10*setp(22)*cm(3)
         sv   = SIN(v)
         cv   = COS(v)
         oecv = 1._10 + setp(9)*cv
         setp(31) = setp(9)*sv/oecv
 
c angle from node on laplacian plane to satellite
         w = per - k - n + v
         r = setp(8)*setp(26)/oecv
 
c coordinates in orbital plane
         setp(29) = r*COS(w)
         setp(30) = r*SIN(w)
 
c setup rotation matrix
         cnck     = cn*ck
         cnsk     = cn*sk
         snck     = sn*ck
         snsk     = sn*sk
         setp(1)  = cnck - cj*snsk
         setp(2)  = snck + cj*cnsk
         setp(3)  = sj*sk
         setp(4)  = -setp(33)*(cnsk + cj*snck) + setp(32)*sn*sj
         setp(5)  = -setp(33)*(snsk - cj*cnck) - setp(32)*cn*sj
         setp(6)  = setp(33)*sj*ck + setp(32)*cj
         setp(27) = setp(33)*cj - setp(32)*sj*ck - 1._10
      endif
 
      if(ic.lt.0) then
c
c output velocities
         vdot    = (setp(11) + 2._10*dt*setp(16) - setp(12))*setp(28)
         rdotr   = vdot*setp(31)
         omeg    = vdot + setp(12) + setp(13)*(setp(33) - 1._10) + 
     .             setp(14)*setp(27)
         vbar(1) = setp(29)*rdotr - setp(30)*omeg
         vbar(2) = setp(30)*rdotr + setp(29)*omeg
         call PRODCT(setp,vbar,x,3,2,1)
 
         return
      else
 
c output coordinates
         call PRODCT(setp,setp(29),x,3,2,1)
         return
      endif
      end
