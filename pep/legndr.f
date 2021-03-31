      subroutine LEGNDR(z, zz, nzone, ntess, leg, leg1, gleg, gleg1)

      implicit none

c
c f amuchastegui - september 1969 - subroutine legndr
c evaluation of legendre polynomials and legendre
c functions using recursion formulas
c
      real*10 glegh1, leg1h1, a2, dnph1, d2np1
      integer*4 n, h, l, m, nsize, nsize1

c dimension in calling program (sbfn)
      real*10 leg(2), leg1(2), leg2(2), gleg(2), gleg1(2), gleg2(2)
      real*10 z, zz
      integer*4 nzone, ntess
c
c        legndr  evaluates the legendre polynomials, its
c        derivatives and also the associated legendre functions
c        the legendre polynomials are defined as:
c        pn(z) = (1 / (n!*2**n) ) * (nth derivat of (z*z-1)**n)
c        and the associated legendre functions as
c        p(n,h) = ((sqrt(1-z*z))**h) * ph(n) ,
c        where ph(n) is the  hth derivative of  pn(z)
c        we can see that for h=0  p(n,0)=pn(z)
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
c        evaluation of p(n,h)
c        the order of the legendre functions in the array gleg,
c        as well as the order of the partials in gleg1, gleg2 is
c        p(2,1), p(2,2), p(3,1), p(3,2), p(3,3), ....... , p(n,n)
c
      if(ntess.le.1) return
      gleg(1) = 3._10*z*zz
      gleg(2) = 3._10*zz**2
      l = 2
      do n = 3, ntess
         m    = ((n-1)*(n-2))/2
         glegh1 = leg(n-2)
         do h = 1, n
            l    = l + 1
            dnph1  = n + h - 1
            gleg(l) = dnph1*zz*glegh1
            if(h.ne.n) gleg(l) = z*gleg(m) + gleg(l)
            glegh1 = gleg(m)
            m = m + 1
            end do
         end do
c
c evaluation of  p'(n,h)
      gleg1(1) = -3._10*(z**2)/zz + 3._10*zz
      gleg1(2) = -6._10*z
      l = 2
      do n = 3, ntess
         glegh1 = leg(n-1)
         do h = 1, n
            l    = l + 1
            a2   = (n-h+1)*(n+h)
            gleg1(l) = (h*z*gleg(l) - a2*zz*glegh1)/(zz**2)
            glegh1 = gleg(l)
            end do
         end do
      return
c
c second derivatives of legendre polynomials and
c legendre function which are needed for partials
c
      entry LEGND2(z, zz, nzone, ntess, leg, leg1, leg2, gleg, gleg1,
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
c evaluation of  p''(n,h)
      if(ntess.le.1) return
      gleg2(1) = -(3._10 + (z**2)/(zz**2))*3._10*z/zz
      gleg2(2) = -6._10
      l = 2
      do n = 3, ntess
         leg1h1 = leg1(n-1)
         glegh1 = leg(n-1)
         do h = 1, n
            l    = l + 1
            a2   = (n-h+1)*(n+h)
            gleg2(l) = h*(gleg(l)+z*gleg1(l)) + 2._10*z*gleg1(l)
     .                  - a2*(zz*leg1h1 - z*glegh1/zz)
            gleg2(l) = gleg2(l)/(zz**2)
            leg1h1 = gleg1(l)
            glegh1 = gleg(l)
            end do
         end do

      return
      end
