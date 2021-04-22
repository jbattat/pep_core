      subroutine MOREDQ(ss,cm,didk,didq,did2q,diddk,diddq,didd2q,
     . didkmx,didqmx,did2qmx)
 
      implicit none
 
c
c rj cappallo   november 1979   sr moredq
c
c     moredq calculates useful inertia tensor perturbation quantities
c     for an elastic moon with constant q.
c     ss  time (julian date)
c     cm  3rd principal moment of inertia of the Moon
c     didk,didq    partial of delta i wrt k2 or k2/q
c     diddk,diddq    corresponding partials of delta i dot (time deriv)
c     didkm,didqm    portion of above depending on moon mass
c                    note that the calculations are explicitly dependent
c                    only on the moon mass, and not at all on the earth,
c                    except insofar as Mmoon=Mass(3)*Mass(10)
c     analytic expressions for elastic perturbations supplied by dh eckh

c parameters
      real*10 ss,cm,didk(3,3),didq(3,3),did2q(3,3),
     . diddk(3,3),diddq(3,3),didd2q(3,3),
     . didkm(3,3),didqm(3,3),diddkm(3,3),diddqm(3,3),
     . didkmx(3,3),didqmx(3,3),did2qmx(3,3)
      real*10 k,q,didy(3,3,6),diddy(3,3,6),didx(3,3,6),diddx(3,3,6)
      integer*2 iparm
      real*10 ddidp(3,3),ddiddp(3,3)

c array dimensions
      include 'globdefs.inc'

c commons
      include 'funcon.inc'
      include 'harmor.inc'
      include 'morstf.inc'
      include 'param.inc'

c local arrays and variables
      real*4 ARG,ARGDOT,sum1(6),sum2(6),sum3(6),sum4(6),l,lp,f,d,
     .       ldot,lpdot,fdot,ddot,cosphi,phi,phidot,sinphi,t1,t2
      real*10 t,a0,a1,a2,mfact,diref(3,3),direfd(3,3),x
      integer i,j,n,nt
      integer*2 map(3,3)/1, 2, 3, 2, 4, 5, 3, 5, 6/
c
c elastic perturbations received from dh eckhardt on 11/30/79
c delaunay argument combinations
      integer*4 id(4, 22)/
     . 0, 0, 0, 0,
     . 0, 0, 0, 2,
     . 0, 0, 2, 0,
     . 1, 0, 0,-2,
     . 1, 0, 0, 0,
     . 2, 0, 0, 0,
     . 0, 0, 1, 0,
     . 1, 0, 1, 0,
     . 1, 0,-1, 0,
     . 0, 1, 0,-2,
     . 1,-1, 0, 0,
     . 1, 0, 2, 0,
     . 1, 1, 0,-2,
     . 1, 1, 0, 0,
     . 2, 0, 0,-2,
     . 3, 0, 0, 0,
     . 0, 1, 0, 0,
     . 1, 0, 0, 2,
     . 0, 0, 1,-2,
     . 0, 0, 1, 2,
     . 1, 0,-1,-2,
     . 2, 0, 1, 0/
c
c sine coeffs. in order of i11,i12,i13,i22,i23,i33
      real*4 s(22,6)/22*0.0,
     . 0.0, -.058, 0.014, 0.085, -.426, -.049, 3*0.0, 0.004, -.004, 0.0,
     .   0.004, 0.003, 0.005, -.006, 0.014, -.014, 4*0.0,
     . 6*0.0, 0.455, 0.062, -.009, 9*0.0, -.007, 0.010, -.012, 0.007,
     . 66*0.0/
 
c cosine matrix coefficients
      real*4 c(22, 6)/ -3.025, -.079, -.026, -.076, -.422, -.058, 3*0.0,
     .   -.005, -.004, -.004, -.003, 0.004, 0.006, -.006, 6*0.0,
     . 44*0.0,
     . 0.852, 0.043, 0.0, 0.036, 0.211, 0.041, 8*0.0, -.008, 0.005, 0.0,
     .   0.012, 4*0.0,
     . 7*0.0, -.025, 0.025, 10*0.0, -.004, -.005, -.004,
     . 2.174, 0.035, 0.026, 0.040, 0.211, 0.017, 5*0.0, 0.005, 5*0.0,
     .   0.006, 4*0.0/
      ARG(a0, a1, a2, x) = MOD(a0 + x*(a1+a2*x), 1._10)*Twopi
      ARGDOT(a1, a2, x)  = (a1 + 2._10*a2*x)*Twopi
 
      nt = 22
 
c evaluate delaunay arguments at time t
      t     = ss - 2.4150205E6_10
      l     = ARG(0.8225128_10, 3.6291645685E-2_10, 1.9E-14_10, t)
      ldot  = ARGDOT(3.6291645685E-2_10, 1.9E-14_10, t)
      lp    = ARG(0.9957662_10, 2.737778519E-3_10, 0._10, t)
      lpdot = ARGDOT(2.737778519E-3_10, 0._10, t)
      f     = ARG(3.12525E-2_10, 3.6748195692E-2_10, -7.E-15_10, t)
      fdot  = ARGDOT(3.6748195692E-2_10, -7.E-15_10, t)
      d     = ARG(0.9742708_10, 3.3863192198E-2_10, -3.E-15_10, t)
      ddot  = ARGDOT(3.3863192198E-2_10, -3.E-15_10, t)
 
      do i = 1, 6
         sum1(i) = 0.0
         sum2(i) = 0.0
         sum3(i) = 0.0
         sum4(i) = 0.0
      end do
 
c loop over different harmonic terms
      do j = 1, nt
         phi    = id(1,j)*l + id(2,j)*lp + id(3,j)*f + id(4,j)*d
         phidot = id(1,j)*ldot + id(2,j)*lpdot + id(3,j)
     .            *fdot + id(4,j)*ddot
         cosphi = cos(phi)
         sinphi = sin(phi)
         do i = 1, 6
            t1 = s(j, i)*sinphi + c(j, i)*cosphi
            t2 = s(j, i)*cosphi - c(j, i)*sinphi
            sum1(i) = sum1(i) + t1
            sum2(i) = sum2(i) - t2
            sum3(i) = sum3(i) + t2*phidot
            sum4(i) = sum4(i) + t1*phidot
         end do
      end do
 
c plug sums into proper matrices
      do i = 1, 3
         do j = 1, 3
            n = map(i,j)
            didk(i,j)  = sum1(n)*Convds*cm
            didq(i,j)  = sum2(n)*Convds*cm
            diddk(i,j) = sum3(n)*Convds*cm
            diddq(i,j) = sum4(n)*Convds*cm
            didkm(i,j) = didk(i,j)
            didqm(i,j) = didq(i,j)
            diddkm(i,j)= diddk(i,j)
            diddqm(i,j)= diddq(i,j)
            didkmx(i,j)= didk(i,j)
            didqmx(i,j)= didq(i,j)
c not used in this model, assumed linear
            did2q(i,j) = 0._10
            didd2q(i,j)= 0._10
            did2qmx(i,j)= 0._10
         end do
      end do
 
c subtract off rotational time average deformation
c (in principle, this operation should be done at a higher level)
      didk(1, 1) = didk(1, 1) + 4.2042891E-24_10
      didk(2, 2) = didk(2, 2) + 4.2042891E-24_10
      didk(3, 3) = didk(3, 3) - 8.4085782E-24_10
      return

      entry MOREDQ2(ss,cm,k,q,didy,diddy,iparm,didx,diddx)
c
c calculate sensitivity of delta I and delta I dot to change in the
c rotation state, to be used in calculating the partial derivatives
c of the angular acceleration
c if iparm.gt.1 then also calculate sensitivity to orbit state

c the formulation of this model is explicitly in terms of the time
c argument alone and is therefore independent of the integrated
c rotation state or orbit state

      do i=1,3
         do j=1,3
            do n=1,6
               didy(i,j,n)=0._10
               diddy(i,j,n)=0._10
               didx(i,j,n)=0._10
               diddx(i,j,n)=0._10
            end do
         end do
      end do
c set up for partials
      do i=1,3
         do j=1,3
            diref(i,j)=k*didkm(i,j)+q*didqm(i,j) 
            direfd(i,j)=k*diddkm(i,j)+q*diddqm(i,j)
         end do
      end do
      return

      entry MOREDQP(ddidp,ddiddp)
c Calculate partials of elasticity/dissipation corrections wrt parameters
c  Information passed in common:
c Mmc,Dmmc= Moon mass divided by C, and partials thereof
c Mass= array of mass ratios
c  Output:
c ddidp  = partial of inertia tensor correction
c ddiddp = partial of inertia tensor rate correction

      if(Icrtrl(Kkk).eq.3 .or. Icrtrl(Kkk).eq.10) then
c partial wrt earth or moon mass
         mfact=Mass(3)
         if(Icrtrl(Kkk).eq.10) mfact=Mass(10)
         do i=1,3
            do j=1,3
               ddidp(i,j)=diref(i,j)/mfact
               ddiddp(i,j)=direfd(i,j)/mfact
            end do
         end do
      else
c partial wrt other parameters for constant-Q model
         do i = 1,3
            do j = 1,3
               ddidp(i,j)  = -diref(i,j)*Dmmc(Kkk)/Mmc
               ddiddp(i,j) = -direfd(i,j)*Dmmc(Kkk)/Mmc
            end do
         end do
      endif

      return
      end
