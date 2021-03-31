      subroutine DIFPRM(prmi, prms, l, n, ndif, cdf, name1, name2)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 dift
      integer   il, ip, ityp, lil, n
 
c*** end of declarations inserted by spag
 
 
c subr. difprm - j.f.chandler - 1980 may
c collect difference vector of input and saved nominals
c
c parameters
      character*8 name1
      character*6 name2
      real*10 prmi(n), prms(n)
      integer*2 l(n)
      integer*4 ndif(2)
      logical*4 cdf(5)
c  prmi - input nominal parameter values
c  prms - nominals from saved normal equations
c  l    - input l-vector (ptrs for adjustable parameters)
c  n    - size of arrays
c  ndif - count of differences in both adjustable and non-adjustable
c         (free and fixed) parameters
c  cdf  - control indicators
c  name1- first name of parameter group (8-char)
c  name2- second name (6-char)
c
c        common
      include 'restor.inc'
 
      il  = 0
      lil = 0
 
c run through parameter vectors
      do ip = 1, n
         if(cdf(1)) prmi(ip) = prms(ip)
         dift = prms(ip) - prmi(ip)
         do while(.true.)
 
c see if adjustable
            if(il.gt.n) then
 
c non-adjustable (fixed)
               ityp = 2
            else if(lil.lt.ip) then
 
c later in list
               il  = il + 1
               lil = l(il)
               if(lil.le.0) il = n + 1
               go to 50
            else if(lil.eq.ip) then
 
c adjustable (free)
               Nsav = Nsav + 1
               Sav(Nsav) = dift
               ityp = 1
            else
               ityp = 2
            endif
 
c note differences
            if(dift.ne.0._10) then
               ndif(ityp) = ndif(ityp) + 1
               if(cdf(ityp+1))
     .             call DIFPRT(name1, name2, ip, prmi(ip), prms(ip),
     .             dift, ityp)
            endif
            go to 100
   50    end do
  100 end do
      return
      end
