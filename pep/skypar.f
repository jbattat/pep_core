      subroutine SKYPAR(kick)
 
      implicit none

c           subr. skypar - j.f.chandler - 1983 jan
c           compute partials of ra and dec w.r.t. star catalog errors
c           (should not exist for other than optical observations)

c arguments
      integer*4 kick

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdat.inc'
      include 'funcon.inc'
      real*10 dsfac(2)
      equivalence (dsfac,Convds)
      include 'ltrapx.inc'
      integer*2 kind
      equivalence (Numpar,kind)
      include 'mtrapx.inc'
      include 'shpcom.inc'
      real*10 sncs(20,2),csky(20),ssky(20)
      equivalence (Leg,csky,sncs(1,1)),(ssky,sncs(1,2))

c local
      integer   i,iflag,isc,j1,j2,l,lim,m

      if(Nskyc.le.0) return
      l = 0
      call PCOPS(m,'SKY ',Iabs1)
      lim = Nskyc
      do while( .true. )
         iflag = 0
         call PCOPY(l,lim,iflag,1,Lsky,Msky)
         if(iflag.gt.0) goto 100
 
c must compute partial
         i   = (l + 3)/4
         isc = 1 + mod(l - 1,2)
         j1  = 1 + mod((l-1)/2,2)
         j2  = 3 - j1
         Deriv(kind,j2) = 0._10
         Deriv(kind,j1) = sncs(i,isc)/dsfac(j2)
         if(l.ge.lim) goto 100
         end do
 
  100 return
      end
