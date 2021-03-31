      subroutine SHPPAR(kick,dlyrad)
 
      implicit none

c     cappallo/ash   august 1970  subroutine shppar
c        modified june,1978 by r.b. goldstein and j. chandler
c     shppar calculates partial derivatives of height(lslope=1) or
c     slope(lslope=2) of a planetary surface wrt spherical harmonic
c     coefficients.
c     it is assumed that /legend/ contains leg. poly. for this co-lat.

c arguments
      integer*4 kick
      real*10 dlyrad

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdat.inc'
      include 'ltrapx.inc'
      integer*2 kind
      equivalence (kind,Numpar)
      include 'mtrapx.inc'
      include 'shpcom.inc'
      include 'shphar.inc'

c local
      integer*4 iflag,lim,m,mn,n,nn,nn1,noff,nsh

      nsh = Nshp + 1
      if(nsh.eq.2) then
c
c*  start=2000
c
c        fourier partials.
c
         call PCOPS(1,'SFOR',Iabs1)
         lim = 122
         n   = 0
         do while( .true. )
            iflag = 0
            call PCOPY(n,lim,iflag,1,Lszhar,Mszhar)
            if(iflag.gt.0) return
            Deriv(kind,1) = 0._10
            if(kick.eq.1) then
               if(Ncode.le.2) Deriv(kind,1) = dlyrad*Gleg(n)
            endif
            Deriv(kind,2) = 0._10
c
c
            if(n.ge.lim) return
            end do
      else if(nsh.eq.3) then
c
c*  start=3000
c grid partials
c
         call GRDPAR
         return
      else
c
c
c*  start=1000
c        spherical harmonics
c
c
         if(Lslop.eq.2) return
c
c partials wrt zonal coefficients
c
         if(Nszone.ge.1) then
 
            lim = Nszone
            n   = 0
            call PCOPS(1,'SZHR',Iabs1)
            do while( .true. )
               iflag = 0
               call PCOPY(n,lim,iflag,1,Lszhar,Mszhar)
               if(iflag.gt.0) go to 50
               Deriv(kind,1) = 0._10
               if(kick.eq.1) then
                  if(Ncode.le.2) Deriv(kind,1) = dlyrad*Leg(n + 1)
               endif
               Deriv(kind,2) = 0._10
               if(n.ge.lim) go to 50
               end do
         endif
      endif
c*  start=1200
c partials of height wrt tesseral and sectoral coeffs.
c first cosine harmonics
   50 noff = 15
      lim  = Nstess
      call PCOPS(1,'SCHR',Iabs1)
  100 n   = 2
      nn1 = 0
      nn  = 2
      mn  = 0
      do while( .true. )
         iflag = 0
         if(noff.eq.15) then
            call PCOPY(mn,lim,iflag,1,Lschar,Mschar)
         else
            call PCOPY(mn,lim,iflag,1,Lsshar,Msshar)
         endif
         if(iflag.gt.0) go to 200
 
c compute partials. first get correct value of m
         do while( nn.lt.mn )
            n   = n + 1
            nn1 = nn
            nn  = nn + n
            end do
         m = mn - nn1
 
         Deriv(kind,1) = 0._10
         if(kick.eq.1) then
            if(Ncode.le.2) Deriv(kind,1) = Gleg(mn + 1)
     .          *sincos(m + noff)*dlyrad
         endif
         Deriv(kind,2) = 0._10
         if(mn.ge.lim) go to 200
         end do
c*  start=1400
c
c sine harmonics set up
c
  200 if(noff.ne.0) then
         noff = 0
         call PCOPS(1,'SSHR',Iabs1)
c*  start=1600
c slope data not yet able to be processed
         go to 100
      endif
c
c*  start=9000
      return
      end
