      subroutine DIFHAR(prmi, prms, l, n, m, ns, ndif, cdf, name1,
     .                  name2)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 dift
      integer   ip, ityp, m, n, nn, ns
 
c*** end of declarations inserted by spag
 
 
c subr. difhar - j.f.chandler - 1980 may
c collect difference vector of input and saved nominals
c
c parameters
      character*8 name1
      character*6 name2
      real*10 prmi(ns, n), prms(ns, n)
      integer*2 l(ns, n)
      integer*4 ndif(2)
      logical*4 cdf(5)
c  prmi - input nominal parameter values
c  prms - nominals from saved normal equations
c  l    - input l-vector (flags for adjustable parameters)
c  n    - size of input arrays
c  m    - size of saved array of nominals
c  ns   - multiplicity of arrays
c  ndif - count of differences in both adjustable and non-adjustable
c         (free and fixed) parameters
c  cdf  - control indicators
c  name1- first name of parameter group (8-char)
c  name2- second name (6-char)
c
c        common
      include 'restor.inc'
 
      nn = n
      if(m.lt.n) nn = m
      if(nn.le.0) return
 
c run through parameter vectors
      do ip = 1, nn
         if(cdf(1)) prmi(1, ip) = prms(1, ip)
         dift = prms(1, ip) - prmi(1, ip)
 
c see if adjustable
         ityp = 2
         if(l(1,ip).gt.0) then
 
c adjustable (free)
            Nsav = Nsav + 1
            Sav(Nsav) = dift
            ityp = 1
         endif
 
c note differences
         if(dift.ne.0._10) then
            ndif(ityp) = ndif(ityp) + 1
            if(cdf(ityp+1)) call DIFPRT(name1, name2, ip, prmi(1,ip),
     .          prms(1,ip), dift, ityp)
         endif
         end do
      return
      end
