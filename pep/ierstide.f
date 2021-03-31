      subroutine IERSTIDE(nvel,gm,g,cphi,h2,h3,l2,l3,ri,a,menrad,
     . daidr,daidphi,daidlam,
     . xi,dxidh,dxidl,delr,delphi,dellam)
      implicit none
c
c IERS computation routine for vector solid-body tides
c
c Based on DEHANTTIDEINEL, ST1IDIU, ST1ISEM, and ST1L1 from the
c International Earth Rotation and Reference Systems Service (IERS)
c Conventions software collection.  The calculation is completely
c reorganized to segregate the offsets by perturber, rather than by
c type of term, and calculations are parametrized.  As such, this
c routine is a functional alternative to VETIDE with almost the same
c calling sequence, the differences being: instead of the local
c acceleration of gravity at the site, this routine expects G*M and
c mean radius, and it also expects both degree-2 and degree-3 Love
c numbers.
c
c NOTE: this routine contains some Earth-specific code that would have
c to be generalized to allow it to be used for the Moon.
c
c This subroutine calculates the tidal displacement due to the
c perturbation of the sun or moon, given the position vectors
c of the perturber (in inertial coordinates) and the earth site (in
c corotating body-fixed coordinates) and the transformation
c between body-fixed and inertial coordinates.
c
c  References:
c
c     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, ''Tidal station
c     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
c
c     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
c     IERS Technical Note No. 36, BKG (2010)

c parameters
      integer*4 nvel
      real*10 gm,g,cphi,h2,h3,l2,l3,ri(6),a(3),menrad,
     . daidr(3),daidphi(3),daidlam(3),
     . xi(6),dxidh(6),dxidl(6),delr,delphi,dellam
c
c  I N P U T   P A R A M E T E R S
c       nvel  =0 tide only
c             =1 tide and tide rate (not implemented, returns zero rate)
c         gm  big g * mass(perturber)
c          g  big g * mass(body)
c       cphi  cosine(site latitude)
c      h2,h3  love numbers - vertical
c      l2,l3  love numbers - horizontal
c         ri  vector from c.o.m. body to perturber's center (inertial)
c          a  vector from c.o.m. body to site (corotating)
c     menrad  mean radius of body
c      daidr  vectors of partial derivatives of the
c    daidphi    inertial site vector w.r.t. body-fixed
c    daidlam    spherical coordinates. (phi=latitude)
c
c  R E T U R N   P A R A M E T E R S
c         xi  tide vector (inertial)
c      dxidh  partials of tide vector (inertial)
c      dxidl   w.r.t. degree-2 love numbers
c       delr  tidal offsets divided by love number
c     delphi   in local body-fixed
c     dellam   frame
c
c all coordinates in light seconds

c local
      real*10 a1,r(6),dxtdh(6),dxtdl(6),r1,r2,sphi,c2lam,c2phi,
     . clam,dp2dh,dp2dl,dx2dl,fac2prt,fac3prt,mass_ratio,p2prt,p3prt,
     . s2lam,scprt,slam,x2prt,x3prt,top2eq(3,3),xt(3),rtopo(3),l1
      integer   i,nv

c external functions
      real*10 DOT

c
c     earth-specific coefficients
c coefficients, out-of-phase corrections induced by mantle anelasticity
c  in the diurnal band.
      real*10 dhid/-0.0025_10/,dlid/-0.0007_10/
c  and semidiurnal band
      real*10 dhis/-0.0022_10/,dlis/-0.0007_10/
c
c coefficients of latitude dependence of love numbers
      real*10 l1d/0.0012_10/,l1sd/0.0024_10/

      nv = 3
      if(nvel.ge.1) nv = 6
c velocities not implemented, returned as zero
      do i=4,nv
         xi(i)=0._10
         dxidh(i)=0._10
         dxidl(i)=0._10
      end do
c useful quantities
      r2 = DOT(ri,ri)
      r1 = SQRT(r2)
      a1 = SQRT(DOT(a,a))
c     sphi = a(3)/a1
c     slam = a(2)/cphi/a1
c     clam = a(1)/cphi/a1
c transformation from topocentric to equatorial in corotating system
      top2eq(1,1)=a(1)/a1
      top2eq(2,1)=a(2)/a1
      top2eq(3,1)=a(3)/a1
      clam=top2eq(1,1)/cphi
      slam=top2eq(2,1)/cphi
      sphi=top2eq(3,1)
      top2eq(1,2)=-sphi*clam
      top2eq(2,2)=-sphi*slam
      top2eq(3,2)=cphi
      top2eq(1,3)=-slam
      top2eq(2,3)=clam
      top2eq(3,3)=0._10
c other trigonometric functions
      c2phi = cphi*cphi - sphi*sphi
      c2lam=clam*clam-slam*slam
      s2lam=2._10*clam*slam

c transform to corotating body-fixed topocentric and equatorial
c cartesian coordinates (topocentric system is left-handed)
      rtopo(1)= DOT(ri,daidr)
      rtopo(2)= DOT(ri,daidphi)/a1
      rtopo(3)= -DOT(ri,daidlam)/a1/cphi
      call PRODCT(top2eq,rtopo,r,3,3,1)

c cosine of angle between site and perturbing body vectors
      scprt = rtopo(1)/r1

c coefficients for terms in direction of site vector
      dp2dh=1.5_10*scprt**2-0.5_10
      dp2dl=-3._10*scprt**2
      p2prt=h2*dp2dh+l2*dp2dl
      p3prt=2.5_10*(h3-3._10*l3)*scprt**3+1.5_10*(l3-h3)*scprt

c coefficients for terms in direction of perturber
      dx2dl=3._10*scprt
      x2prt=l2*dx2dl
      x3prt=1.5_10*l3*(5._10*scprt**2-1._10)

c calculate ordinary displacement and partials w.r.t. degree-2 love nos
      mass_ratio= gm/g
      fac2prt=mass_ratio*menrad*(menrad/r1)**3
      fac3prt=fac2prt*(menrad/r1)
      do i=1,3
         xt(i)=(fac2prt*x2prt+fac3prt*x3prt)*rtopo(i)/r1
         dxtdl(i)=fac2prt*dx2dl*rtopo(i)/r1         
         dxtdh(i)=0._10
      end do
      xt(1) = xt(1) + fac2prt*p2prt+fac3prt*p3prt
      dxtdl(1)=dxtdl(1) + fac2prt*dp2dl
      dxtdh(1)=fac2prt*dp2dh

c corrections for the out-of-phase part of love numbers (part h_2^(0)i  
c and l_2^(0)i )  

c diurnal band
      xt(1)=xt(1)+3d0*dhid*sphi*cphi*fac2prt*r(3)*rtopo(3)/r2
      xt(2)=xt(2)+3d0*dlid*c2phi*fac2prt*r(3)*rtopo(3)/r2
      xt(3)=xt(3)-3d0*dlid*sphi*fac2prt*r(3)*
     . (r(1)*clam+r(2)*slam)/r2

c semidiurnal band
      xt(1)=xt(1)-0.75_10*dhis*cphi**2*fac2prt*((r(1)**2-r(2)**2)*
     . s2lam-2._10*r(1)*r(2)*c2lam)/r2
      xt(2)=xt(2)+1.5_10*dlis*sphi*cphi*fac2prt*((r(1)**2-r(2)**2)*
     . s2lam-2._10*r(1)*r(2)*c2lam)/r2
      xt(3)=xt(3)-1.5_10*dlis*cphi*fac2prt*((r(1)**2-r(2)**2)*
     . c2lam+2._10*r(1)*r(2)*s2lam)/r2

c corrections for the latitude dependence of love numbers (part l^(1) )  

c diurnal band.
      l1=l1d*3._10
      xt(2)=xt(2)-l1*sphi**2*fac2prt*r(3)*(r(1)*clam+r(2)*slam)/r2
      xt(3)=xt(3)-l1*sphi*c2phi*fac2prt*r(3)*rtopo(3)/r2

c semidiurnal
      l1=l1sd*1.5_10
      xt(2)=xt(2)-l1*sphi*cphi*fac2prt*((r(1)**2-r(2)**2)*
     . c2lam+2d0*r(1)*r(2)*s2lam)/r2
      xt(3)=xt(3)-l1*sphi**2*cphi*fac2prt*((r(1)**2-r(2)**2)*
     . s2lam-2d0*r(1)*r(2)*c2lam)/r2

c convert to west longitude
      dellam=-xt(3)
      delphi= xt(2)
      delr  = xt(1)

c transform to inertial frame
      do i=1,3
         xi(i)=xt(1)*daidr(i)+xt(2)*daidphi(i)/a1
     .    -xt(3)*daidlam(i)/cphi/a1
         dxidh(i)=dxtdh(1)*daidr(i)+dxtdh(2)*daidphi(i)/a1
     .    -dxtdh(3)*daidlam(i)/cphi/a1
         dxidl(i)=dxtdl(1)*daidr(i)+dxtdl(2)*daidphi(i)/a1
     .    -dxtdl(3)*daidlam(i)/cphi/a1
      end do

c factor out love numbers
      if(h2.ne.0._10) delr=delr/h2
      if(l2.ne.0._10) then
         delphi=delphi/l2
         dellam=dellam/l2
      endif

      return
      end
      subroutine IERSTID2(a,cphi,fract,t,daidr,daidphi,daidlam,
     . xi,delr,delphi,dellam)
      implicit none

c     This routine is based on STEP2DIU and STEP2LON from the
c     International Earth Rotation and Reference Systems Service (IERS)
c     Conventions software collection.
c     
c     It gives the in-phase and out-of-phase tidal corrections induced
c     by mantle anelasticity.
c
c  References:
c
c     Mathews, P. M., Dehant, V., and Gipson, J. M., 1997, "Tidal station
c     displacements," J. Geophys. Res., 102(B9), pp. 20,469-20,477
c
c     Petit, G. and Luzum, B. (eds.), IERS Conventions (2010),
c     IERS Technical Note No. 36, BKG (2010)

c parameters
      real*10 a(3),cphi,fract,t,daidr(3),daidphi(3),daidlam(3),
     . xi(3),delr,delphi,dellam
c a     - earth-fixed cartesian coordinates of the site
c cphi  - cos(latitude) of site
c fract - UTC fraction of a day of the desired epoch
c t     - coordinate time from J2000.0 in Julian centuries
c daidr - vectors of partial derivatives of the
c daidphi    inertial site vector w.r.t. body-fixed
c daidlam    spherical coordinates. (phi=latitude)
c     on return
c xi    - correction to site position in inertial equatorial system
c delr  - correction to site position in
c delphi   local body-fixed
c dellam   frame
c
c all coordinates given in light seconds

      include 'globdefs.inc'
c common
      include 'funcon.inc'
      include 'param.inc'

c local
      integer*4 i,j
      real*10 s,tau,pr,h,p,zns,ps,a1,sphi,s2phi,c2phi,
     . clam,slam,f3phi,zla,thetaf,sth,cth,xt(3)

c external functions
      real*10 DOT

c tidal terms, amplitudes in millimeters
c fundamental arguments:
c 1 s   (1 month)
c 2 h   (1 year)
c 3 p   (8.8 years)
c 4 zns (18.6 years)
c 5 ps  (20936 years)
      real*8 datdi(9,31),datlo(9,5)
c diurnal terms: 6&7 radial, 8&9 horizontal (in and out of phase)
      data ((datdi(i,j),i=1,9),j=1,31)/
     . -3d0, 0d0, 2d0, 0d0, 0d0,  -0.01d0, 0d0,    0d0,    0d0,
     . -3d0, 2d0, 0d0, 0d0, 0d0,  -0.01d0, 0d0,    0d0,    0d0,
     . -2d0, 0d0, 1d0,-1d0, 0d0,  -0.02d0, 0d0,    0d0,    0d0,
     . -2d0, 0d0, 1d0, 0d0, 0d0,  -0.08d0, 0d0,   -0.01d0, 0.01d0,
     . -2d0, 2d0,-1d0, 0d0, 0d0,  -0.02d0, 0d0,    0d0,    0d0,
     . -1d0, 0d0, 0d0,-1d0, 0d0,  -0.10d0, 0d0,    0d0,    0d0,
     . -1d0, 0d0, 0d0, 0d0, 0d0,  -0.51d0, 0d0,   -0.02d0, 0.03d0,
     . -1d0, 2d0, 0d0, 0d0, 0d0,   0.01d0, 0d0,    0d0,    0d0,
     .  0d0,-2d0, 1d0, 0d0, 0d0,   0.01d0, 0d0,    0d0,    0d0,
     .  0d0, 0d0,-1d0, 0d0, 0d0,   0.02d0, 0d0,    0d0,    0d0,
     .  0d0, 0d0, 1d0, 0d0, 0d0,   0.06d0, 0d0,    0d0,    0d0,
     .  0d0, 0d0, 1d0, 1d0, 0d0,   0.01d0, 0d0,    0d0,    0d0,
     .  0d0, 2d0,-1d0, 0d0, 0d0,   0.01d0, 0d0,    0d0,    0d0,
     .  1d0,-3d0, 0d0, 0d0, 1d0,  -0.06d0, 0d0,    0d0,    0d0,
     .  1d0,-2d0, 0d0,-1d0, 0d0,   0.01d0, 0d0,    0d0,    0d0,
     .  1d0,-2d0, 0d0, 0d0, 0d0,  -1.23d0,-0.07d0, 0.06d0, 0.01d0,
     .  1d0,-1d0, 0d0, 0d0,-1d0,   0.02d0, 0d0,    0d0,    0d0,
     .  1d0,-1d0, 0d0, 0d0, 1d0,   0.04d0, 0d0,    0d0,    0d0,
     .  1d0, 0d0, 0d0,-1d0, 0d0,  -0.22d0, 0.01d0, 0.01d0, 0d0,
     .  1d0, 0d0, 0d0, 0d0, 0d0,  12.00d0,-0.80d0,-0.67d0,-0.03d0,
     .  1d0, 0d0, 0d0, 1d0, 0d0,   1.73d0,-0.12d0,-0.10d0, 0d0,
     .  1d0, 0d0, 0d0, 2d0, 0d0,  -0.04d0, 0d0,    0d0,    0d0,
     .  1d0, 1d0, 0d0, 0d0,-1d0,  -0.50d0,-0.01d0, 0.03d0, 0d0,
     .  1d0, 1d0, 0d0, 0d0, 1d0,   0.01d0, 0d0,    0d0,    0d0,
     .  0d0, 1d0, 0d0, 1d0,-1d0,  -0.01d0, 0d0,    0d0,    0d0,
     .  1d0, 2d0,-2d0, 0d0, 0d0,  -0.01d0, 0d0,    0d0,    0d0,
     .  1d0, 2d0, 0d0, 0d0, 0d0,  -0.11d0, 0.01d0, 0.01d0, 0d0,
     .  2d0,-2d0, 1d0, 0d0, 0d0,  -0.01d0, 0d0,    0d0,    0d0,
     .  2d0, 0d0,-1d0, 0d0, 0d0,  -0.02d0, 0d0,    0d0,    0d0,
     .  3d0, 0d0, 0d0, 0d0, 0d0,   0d0,    0d0,    0d0,    0d0,
     .  3d0, 0d0, 0d0, 1d0, 0d0,   0d0,    0d0,    0d0,    0d0/
c long-period terms: 6&8 radial, 7&9 horizontal (in and out of phase)
      data ((datlo(i,j),i=1,9),j=1,5)/
     .  0d0, 0d0, 0d0, 1d0, 0d0,   0.47d0, 0.23d0, 0.16d0, 0.07d0,
     .  0d0, 2d0, 0d0, 0d0, 0d0,  -0.20d0,-0.12d0,-0.11d0,-0.05d0,
     .  1d0, 0d0,-1d0, 0d0, 0d0,  -0.11d0,-0.08d0,-0.09d0,-0.04d0,
     .  2d0, 0d0, 0d0, 0d0, 0d0,  -0.13d0,-0.11d0,-0.15d0,-0.07d0,
     .  2d0, 0d0, 0d0, 1d0, 0d0,  -0.05d0,-0.05d0,-0.06d0,-0.03d0/

c  compute the phase angles in degrees.
      s = 218.31664563_10
     .  + (481267.88194_10
     .  + (-0.0014663889_10
     .  + (0.00000185139_10)*t)*t)*t

      tau = fract*360._10
     .    + 280.4606184_10
     .    + (36000.7700536_10
     .    + (0.00038793_10
     .    + (-0.0000000258_10)*t)*t)*t
     .    + (-s)

      pr = (1.396971278_10
     .   + (0.000308889_10
     .   + (0.000000021_10
     .   + (0.000000007_10)*t)*t)*t)*t

      s = s + pr

      h = 280.46645_10
     .  + (36000.7697489_10
     .  + (0.00030322222_10
     .  + (0.000000020_10
     .  + (-0.00000000654_10)*t)*t)*t)*t

      p = 83.35324312_10
     .  + (4069.01363525_10
     .  + (-0.01032172222_10
     .  + (-0.0000124991_10
     .  + (0.00000005263_10)*t)*t)*t)*t

      zns = 234.95544499_10
     .    + (1934.13626197_10
     .    + (-0.00207561111_10
     .    + (-0.00000213944_10
     .    + (0.00000001650_10)*t)*t)*t)*t

      ps = 282.93734098_10
     .   + (1.71945766667_10
     .   + (0.00045688889_10
     .   + (-0.00000001778_10
     .   + (-0.00000000334_10)*t)*t)*t)*t

c reduce angles to the range -360 to 360
      s =  MOD(s,360._10)
      tau = MOD(tau,360._10)
      h =  MOD(h,360._10)
      p =  MOD(p,360._10)
      zns = MOD(zns,360._10)
      ps = MOD(ps,360._10)

      a1 = SQRT(DOT(a,a))
      sphi = a(3)/a1
      s2phi= 2._10*sphi*cphi
      c2phi= cphi*cphi-sphi*sphi
      f3phi= (sphi*sphi-c2phi)*0.5_10

      clam = a(1)/cphi/a1
      slam = a(2)/cphi/a1
      zla = ATAN2(a(2),a(1))

c initialize topocentric correction components
      xt(1)=0._10
      xt(2)=0._10
      xt(3)=0._10

c diurnal band terms      
      do j=1,31
         thetaf=(tau+datdi(1,j)*s+datdi(2,j)*h+datdi(3,j)*p+
     .    datdi(4,j)*zns+datdi(5,j)*ps)*Convd + zla
         sth=SIN(thetaf)
         cth=COS(thetaf)

         xt(1)=xt(1)+datdi(6,j)*s2phi*sth + datdi(7,j)*s2phi*cth
         xt(2)=xt(2)+datdi(8,j)*c2phi*sth + datdi(9,j)*c2phi*cth
         xt(3)=xt(3)+datdi(8,j)*sphi*cth - datdi(9,j)*sphi*sth
      end do

c long-period terms
      do j=1,5
         thetaf=(datlo(1,j)*s+datlo(2,j)*h+datlo(3,j)*p+
     .    datlo(4,j)*zns+datlo(5,j)*ps)*Convd
         sth=SIN(thetaf)
         cth=COS(thetaf)

         xt(1)=xt(1)+datlo(6,j)*f3phi*cth + datlo(8,j)*f3phi*sth
         xt(2)=xt(2)+datlo(7,j)*s2phi*cth + datlo(9,j)*s2phi*sth
      end do

c convert to light seconds
      xt(1)=xt(1)/1e6_10/Ltvel
      xt(2)=xt(2)/1e6_10/Ltvel
      xt(3)=xt(3)/1e6_10/Ltvel

c convert to west longitude
      dellam=-xt(3)
      delphi= xt(2)
      delr  = xt(1)

c transform to inertial frame
      do i=1,3
         xi(i)=xt(1)*daidr(i)+xt(2)*daidphi(i)/a1
     .    -xt(3)*daidlam(i)/cphi/a1
      end do

      return
      end
