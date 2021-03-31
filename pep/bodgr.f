      subroutine BODGR(rvec,rab,v,a,sumr)
 
      implicit none
 
      real*10 rab(15,15),rvec(3,15,15),v(3,15),sumr(3)
      integer*4 a
 
c     J.F.Chandler - 1991 Dec
c     Based on 'SOLGR' by R.W.Babcock - April 1984
c
c Compute general relativistic correction to force on sun-centered body
c for n-body integration.  Based on eq.  (6.34) of Will, Theory and
c Experiment in Gravitational Physics, Cambridge Univ. Press, 1981.
c (Actually this equation is in isotropic rather than harmonic
c coordinates, but these are the same to 1st post-newtonian order).
c
c    calling parameters
 
c  RVEC  - VECTORS TO BODY 1 FROM BODY 2
c  RAB   - DISTANCES
c  V     - VELOCITIES (IN SSBC FRAME)
c  A     - INDEX OF INTEGRATED BODY IN LIST NPLBDY
c           A < Nsun   EVALUATE ACCELERATION FOR BODY A
c           A = Nsun   SETUP ONCE PER ITERATION
c  SUMR  - RETURNED VECTOR GR ACCELERATION ON BODY A
c
c    local arrays
c
c  VAB2  - DOT PRODUCTS OF VELOCITIES OF BODY 1 AND BODY 2
c  RDOTV - DOT PRODUCTS OF RVEC1,2 AND V1

c array dimensions
      include 'globdefs.inc'
c
      include 'bodstf.inc'
      include 'intstf.inc'
      include 'param.inc'
c
c local
c
      real*10 vab2(15,15),rdotv(15,15),sums(3)
      real*10 DOT, f1, f2, f3, f4, f5, f6, f7, tsum
      integer*4 b,c,i
c
c     note on subscripts in this routine:
c         a,b,c (integers) match formula (6.34) or (6.78) in will
c         a is integrated body
c         last body in list is sun
c
      if(a.eq.Nsun) then
c----------------------------------------------------------------------
c     SETUP ONCE PER ITERATION OF A STEP
c
c           SET UP DOT PRODUCTS OF POSITIONS AND VELOCITIES
         do b=1,Nsun
            do c=b,Nsun
               rdotv(c,b)=DOT(rvec(1,c,b),v(1,c))
               if(c.ne.b) rdotv(b,c)=DOT(rvec(1,b,c),v(1,b))
               vab2(b,c)=DOT(v(1,b),v(1,c))
               vab2(c,b)=vab2(b,c)
            end do
         end do
c initialize for calculating sun's acceleration
         do i=1,3
            sumr(i)=0._10
         end do
      else

c for bodies other than the sun, subtract the sun's acceleration
         do i=1,3
            sumr(i)=-sums(i)
         end do
      endif
c
      do b = 1,Nsun
         if(b.ne.a) then
            f4   = Mass1(b)/rab(a,b)
            f5   = Mass1(a)/rab(a,b)
            tsum = B2g2*f4 + B2g21*f5
            f1   = 0._10
            f2   = 0._10
            f3   = 0._10
            do c = 1,Nsun
               if(c.ne.a .and. c.ne.b) then
                  f1 = f1 + Mass1(c)/rab(b,c)
                  f2 = f2 + Mass1(c)/rab(a,c)
                  f3 = f3 + Mass1(c)*DOT(rvec(1,a,b),rvec(1,b,c))/
     .                     rab(b,c)**3
               endif
            end do
 
            tsum = (tsum + B2m1*f1 + B2g2*f2 - f3/2._10)*A44
            tsum = tsum+(G22*vab2(a,b)-G11*vab2(b,b)-Gamapm*vab2(a,a)
     .                 + 1.5*(rdotv(b,a)/rab(a,b))**2)/Cvel2
            do i = 1,3
               sumr(i)  = sumr(i) + tsum*(f4/rab(a,b)**2)*rvec(i,a,b)
            end do

            f6 = f4*A44
            do c = 1,Nsun
               if(c.ne.a .and. c.ne.b) then
                  f7 = f6*Mass1(c)/rab(b,c)**3
                  do i = 1,3
                     sumr(i)  = sumr(i) - f7*rvec(i,b,c)*G7
                  end do
               endif
            end do
         endif
 
      end do
 
      do b = 1,Nsun
         if(b.ne.a) then
            f1 = rdotv(a,b)
            f2 = rdotv(b,a)
            f3 = (G22*f1 + G21*f2)*(Mass1(b)/(rab(a,b)**3*Cvel2))
            do i = 1,3
               sumr(i)  = sumr(i) + f3*(v(i,a) - v(i,b))
            end do
         endif
      end do
c save sun's acceleration for reuse
      if(a.eq.Nsun) then
         do i=1,3
            sums(i)=sumr(i)
         end do
      endif
c
      return
c
      end
