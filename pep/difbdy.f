      subroutine DIFBDY(prmi, prms, l, n, ndif, cdf, name1, name2)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 dift
      integer   ip, ityp, n, nn, np
 
c*** end of declarations inserted by spag
 
 
c subr. difbdy - j.f.chandler - 1980 may
c collect difference vector of input and saved nominals
c
c parameters
      character*8 name1
      character*6 name2
      real*10 prmi(30), prms(30)
      integer*2 l(30)
      integer*4 ndif(2)
      logical*4 cdf(5)
c  prmi - input nominal parameter values
c  prms - nominals from saved normal equations
c  l    - input l-vector (flags, ptrs for adjustable parameters)
c  n    - size of arrays (if neg. then only flags)
c  ndif - count of differences in both adjustable and non-adjustable
c         (free and fixed) parameters
c  cdf  - control indicators
c  name1- first name of parameter group (8-char)
c  name2- second name (6-char)
c
c        common
      include 'restor.inc'
 
      nn = 6
      if(n.lt.0) nn = -n
 
c run through parameter vectors
      do ip = 1, nn
         if(cdf(1)) prmi(ip) = prms(ip)
         dift = prms(ip) - prmi(ip)
 
c see if adjustable
         ityp = 2
         if(l(ip).gt.0) then
 
c adjustable (free)
            Nsav = Nsav + 1
            Sav(Nsav) = dift
            ityp = 1
         endif
 
c note differences
         if(dift.ne.0._10) then
            ndif(ityp) = ndif(ityp) + 1
            if(cdf(ityp+1)) call DIFPRT(name1, name2, ip, prmi(ip),
     .          prms(ip), dift, ityp)
         endif
         end do
      np = n - nn
      if(np.gt.0) call DIFPRM(prmi(7), prms(7), l(7), np, ndif,
     .                        cdf, name1, 'CON   ')
      return
      end
