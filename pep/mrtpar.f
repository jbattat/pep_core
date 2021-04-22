      subroutine MRTPAR(kick)
 
      implicit none
 
c
c king/ash september 1972   subroutine mrtpar
c calculate partial derivatives w.r.t. lunar libration parameters
c
c

c parameters
      integer kick

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdat.inc'
      include 'ltrapx.inc'
      integer*2 kind
      equivalence (Numpar, kind)
      include 'mnsprt.inc'
      include 'mtrapx.inc'
      include 'number.inc'
      include 'pemctl.inc'
      include 'partcm.inc'
      integer kmrt
      equivalence (Kprt,kmrt)
      include 'plndta.inc'
      include 'rotdtamr.inc'
      include 'tapdta.inc'

c local 
      integer i,idiss,iflag,j,k,l,lmon,lmrt,m,n,nstart
      integer*2 lprml,nqlnt
      integer*2 klakmr/0/
 
      kmrt = 1
      lmrt = 2
      Kmon = Lparm
      lmon = 8
      m    = 7
      l    = 0
 
c begin loop l=1,u_nmbod
      call PCOPS(m,'MR  ',Iabs1)
  100 if(l.lt.6) then
         iflag = 0
         call PCOPY(l,6,iflag,1,Lmrx,Mmrx)
         if(iflag.gt.0) goto 400
         if(Nspot.le.0 .and. Nspot2.le.0) goto 200
         if(Nplnt0.ne.10) goto 200
 
c search for partials w.r.t. initial conditions on ilib
         if(Ilib.le.0) goto 200
         call PBDIC(Kimr,lmrt,kmrt,l,'MR',klakmr,0)
c interpolate partials from tape and calculate partial of observation
         call CPARTL(6,3,kick)
c add contribution from moon orbit partials (these have neg. ki flags)
         lprml = 1000+l
         iflag = 1
         call PBDPRM(Nkimn,Kimn,lmon,Kmon,lprml,iflag)
         if(iflag.le.0) call CPARTL(2,3,kick)
         goto 380
      else if(l.eq.6) then
         kmrt = Lparmr
         lmrt = 8
      endif
c
c see if saved partials w.r.t. libration parameters can be used
      iflag = 1
      call PCOPY(l,u_nmbod,iflag,1,Lmrx,Mmrx)
      if(iflag.gt.0) goto 500
      if((Nspot.gt.0 .or. Nspot2.gt.0) .and. Nplnt0.eq.10) then
 
c search for partials w.r.t. libration parameters on ilib
         if(Ilib.le.0) goto 300
         iflag = 1
         call PBDPRM(Nkimr,Kimr,lmrt,kmrt,Lmrx(l),iflag)
         if(iflag.gt.0) goto 300
         call CPARTL(6,3,kick)
c add contribution from moon orbit partials (these have neg. ki flags)
         lprml = 1006+Lmrx(l)
         iflag = 1
         call PBDPRM(Nkimn,Kimn,lmon,Kmon,lprml,iflag)
         if(iflag.le.0) call CPARTL(2,3,kick)
         goto 380
      endif
  200 Deriv(kind,1) = 0._10
      Deriv(kind,2) = 0._10
      goto 100
 
c calculate partials w.r.t. libration parameter not on tape
  300 n = 1
      if(Nspot2.gt.0) n = 2
      if(Lmrx(l).eq.3) then
         do k = 1, n
            do i = 1, Mouse
               do j = 1, Index
                  Derpr(j,i,k) = -Dxdbet(j,k)
               end do
            end do
         end do
      else if(Lmrx(l).eq.4) then
         do k = 1, n
            do i = 1, Mouse
               do j = 1, Index
                  Derpr(j,i,k) = -Dxdgam(j,k)
               end do
            end do
         end do
      else if(Lmrx(l).eq.5) then
         do k = 1, n
            do i = 1, Mouse
               do j = 1, Index
                  Derpr(j,i,k) = -Dxdth(j,k)
               end do
            end do
         end do
c partials w.r.t. analytic dissipation coefficients
      else if(Lmrx(l).eq.8 .or. Lmrx(l).eq.9 .or. Lmrx(l).eq.10) then
         idiss = Lmrx(l)-7
         do k = 1, n
            do i = 1, Mouse
               do j = 1, Index
                  Derpr(j,i,k) = -Dxdphi(j,k)*Disstrm(idiss)
               end do
            end do
         end do
      else
         goto 200
      endif
      Ivze(1) = 1
c
c calculate partials of observation
  380 call CPARTC(kick)
  400 if(l.lt.u_nmbod) goto 100
c end loop l=1,u_nmbod
c calculate partials w.r.t. gravitational harmonics
  500 if(Klan.eq.u_mxpl+1 .and. Klanb.le.0) then
         kmrt  = Lparmr
         nqlnt = -10
 
c zonal harmonics
         call HPARTL(kick,Nszone,Lszhar,Mszhar,nqlnt,31,8,Iabs1)
 
c tesseral cosine harmonics
         call HPARTM(Nstess,Lschar,Mschar,41)
 
c tesseral sine harmonics
         call HPARTM(Nstess,Lsshar,Msshar,51)
      endif
 
      return
      end
