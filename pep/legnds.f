      subroutine LEGNDS(z,zz,nzone,ntess,leg,leg1,gleg,gleg1)

      implicit none

c
c j.davis - october 1992 - subroutine legnds
c modified 1993 January - J.F.Chandler
c derived from: f amuchastegui - september 1969 - subroutine legndr
c evaluation of scaled legendre polynomials and scaled legendre functions
c using recursion formulas
c
c array dimensions
      include 'globdefs.inc'
c
      real*10 glegh1,leg1h1,dnph1,d2np1
      integer*4 n,h,l,m,nsize,nsize1

c dimension in calling program (sbfn)
      real*10 leg(2),leg1(2),leg2(2),gleg(2),gleg1(2),gleg2(2)
      real*10 z,zz
      integer*4 nzone,ntess
c
c Legendre function scaling factors
      real*10 d2nm1,dnph,dnmh,dnmh1,sqr15
c
      real*10 anh((u_mxtes*(u_mxtes+1))/2-1),
     . bnh((u_mxtes*(u_mxtes+1))/2-1),dnh((u_mxtes*(u_mxtes+1))/2-1)
c
c anh(k) = a(n,h), bnh(k) = b(n,h), dnh(k) = d(n,h)
c      where k = (n*(n-1))/2 + h - 1) for n=2,50 & h=1,n
c
c    a(n,h) = sqrt{((2*n + 1)*(n-h))/((2*n - 1)*(n+h))}
c    b(n,h) = sqrt{((2*n + 1)*(n+h-1))/((2*n - 1)*(n+h))}     h>1
c    b(n,1) = sqrt{((2*n + 1)*(n)*2)/(n + 1)}
c    d(n,h) = sqrt{(n + h)*(n - h + 1)}                       h>1
c    d(n,1) = sqrt{(n + 1)*(n)*2*(2*n + 1)}
c
c        legnds  evaluates the legendre polynomials, its
c        derivatives and also the associated legendre functions
c        the legendre polynomials are defined as:
c        pn(z) = (1 / (n!*2**n) ) * (nth derivat of (z*z-1)**n)
c        and the associated legendre functions as
c        p(n,h) = ((sqrt(1-z*z))**h) * ph(n) ,
c        where ph(n) is the  hth derivative of  pn(z)
c        we can see that for h=0  p(n,0)=pn(z)
c        The Legendre functions are scaled by the factor
c        sqrt({2*(2*n+1)*(n-h)!}/{(n+h)!}.  The Legendre polynomials
c        are left unscaled in this routine.
c
c        the formulae are from  "tables of integrals, series
c        and products" - gradsteyn & ryzhik,pags 1004-27
c$$$$$$$ the factor (-1)**h  for the legendre function used in this
c        reference has been set = 1 to use the convention of
c        the smithsonian astr. observatory journal
c        the formulae for  p'(z)  and  p''(z)  are from "differential
c        equations with applications", ritger & rose, page 223
c
      nsize = nzone - 1
      if(ntess.gt.nzone) nsize = ntess - 1
c
c        evaluation of  p(z)
c        the first polynomial in the array  leg(n)  is the
c        second order polynomial, since the first order is not
c        used in PEP
c
      leg(1) = 1.5_10*z**2 - 0.5_10
      leg(2) = 2.5_10*z**3 - 1.5_10*z
      do n = 3, nsize
         d2np1 = 2*n + 1
         leg(n) = (d2np1*z*leg(n-1) - n*leg(n-2))/(n+1)
      end do
c
c evaluation of  p'(z)
      leg1(1) = 3._10*z
      leg1(2) = 7.5_10*z**2 - 1.5_10
      do n = 3, nsize
         d2np1 = 2*n + 1
         leg1(n) = leg1(n - 2) + d2np1*leg(n - 1)
      end do
c
c        evaluation of p(n,h)                 -- scaled from LEGNDR routine
c        the order of the legendre functions in the array gleg,
c        as well as the order of the partials in gleg1, gleg2 is
c        p(2,1), p(2,2), p(3,1), p(3,2), p(3,3), ....... , p(n,n)
c     The following recursive formulas were used for p(n,h) based on the
c        modification that p(n,h) is scaled for h>0 and unscaled for h=0
c        p(n,h) = a(n,h)*z*p(n-1,h) + zz*b(n,h)*p(n-1,h-1)
c        p(n,n) = zz*b(n,n)*p(n-1,n-1)
c
      if(ntess.le.1) return
      gleg(1) = sqr15*z*zz
      gleg(2) = 0.5_10*sqr15*zz**2
      l = 2
      do n = 3, ntess
         m    = ((n-1)*(n-2))/2
         glegh1 = leg(n-2)
         do h = 1, n
            l    = l + 1
            gleg(l) = bnh(l)*zz*glegh1
            if(h.ne.n) gleg(l) = z*gleg(m)*anh(l) + gleg(l)
            glegh1 = gleg(m)
            m = m + 1
         end do
      end do
c
c evaluation of  p'(n,h)                 -- scaled from LEGNDR routine
c     The following recursive formulas were used for p'(n,h) based on the
c        modification that p'(n,h) is scaled for h>0 and unscaled for h=0
c        p'(n,h) = h*z/zz^2 * p(n,h) - d(n,h)/zz * p(n, h-1)
c
      gleg1(1) = sqr15*(zz-(z**2)/zz)
      gleg1(2) = -sqr15*z
      l = 2
      do n = 3, ntess
         glegh1 = leg(n-1)
         do h = 1, n
            l    = l + 1
            gleg1(l) = (h*z*gleg(l) - dnh(l)*zz*glegh1)/(zz**2)
            glegh1 = gleg(l)
         end do
      end do
      return
c
c second derivatives of legendre polynomials and
c legendre function which are needed for partials
c
      entry LEGNS2(z, zz, nzone, ntess, leg, leg1, leg2, gleg, gleg1,
     .             gleg2)
c
c evaluation of p''(z)
      if(nzone.gt.1) leg2(1) = 3._10
      if(nzone.gt.2) leg2(2) = 15._10*z
      nsize1  = nzone - 1
      do n = 3, nsize1
         d2np1   = 2*n + 1
         leg2(n) = leg2(n - 2) + d2np1*leg1(n - 1)
      end do
c
c evaluation of  p''(n,h)                 -- scaled from LEGNDR routine
c     The following recursive formulas were used for p''(n,h) based on the
c        modification that p''(n,h) is scaled for h>0 and unscaled for h=0
c        p''(n,h) = h*(1+z^2)/zz^4 * p(n,h) + h*z/zz^2 * p'(n,h)
c		      - d(n,h)/zz^3 * {z * p(n,h-1) + zz^2 * p'(n, h-1)}
c
      if(ntess.le.1) return
      gleg2(1) = -(3._10 + (z**2)/(zz**2))*sqr15*z/zz
      gleg2(2) = -sqr15
      l = 2
      do n = 3, ntess
         leg1h1 = leg1(n-1)
         glegh1 = leg(n-1)
         do h = 1, n
            l    = l + 1
            gleg2(l) = h*(gleg(l)+z*gleg1(l)) + 2._10*z*gleg1(l)
     .                  - dnh(l)*(zz*leg1h1 - z*glegh1/zz)
            gleg2(l) = gleg2(l)/(zz**2)
            leg1h1 = gleg1(l)
            glegh1 = gleg(l)
         end do
      end do

      return
c
c Compute scaling factors intially for use in calls to LEGNDS/LEGNS2
c
      entry LEGSET
c
      sqr15 = SQRT(15._10)
c
      nsize = 50
      l = 0
      do n = 2,nsize
         d2np1 = 2*n+1
         d2nm1 = 2*n-1
         do h = 1,n
            l = l + 1
            dnph  = n+h
            dnmh  = n-h
            dnph1 = n+h-1
            dnmh1 = n-h+1
c modify scaling for d(n,1) and b(n,1) because of unscaled P(n,0)
            if(h.eq.1) then
               dnph1 = dnph1*2._10*d2nm1
               dnmh1 = dnmh1*2._10*d2np1
            endif
c
            anh(l) = SQRT((d2np1*dnmh)/(d2nm1*dnph))
            bnh(l) = SQRT((d2np1*dnph1)/(d2nm1*dnph))
            dnh(l) = SQRT(dnph*dnmh1)
         end do
      end do
      return   
      end
