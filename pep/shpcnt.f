      subroutine SHPCNT
 
      implicit none

c     ash/cappallo    august 1970    subroutine shpcnt
c     calculate partial derivative control vectors for shape harmonics
c     modified oct 1974 by r. king to calculate control vectors for
c     gravitational harmonics of the moon as an observed body

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'fcntrl.inc'
      include 'ltrapx.inc'
      integer*2 lshape(1000)
      equivalence (lshape,Lszhar)
      include 'monhar.inc'
      include 'mtrapx.inc'
      include 'number.inc'
      include 'shphar.inc'

c local
      integer*4 i,j,n2
      integer*2 izr2/0/
c
c determine observed body shape or gravitational potential harmonic
c controls
      Nszone = 0
      Nstess = 0
      Lnfour = 0
      Lngd   = 0
      do j = 1, 1000
         lshape(j) = 0
      end do
      Shpnit = .false.
      n2     = 0
 
c zero out controls that are set only for grid model
      Lnshob = 0
      Lixshp = 0
      do i = 1, 4
         Lshobs(i) = 0
      end do
c
c shall partial derivative controls be determined
      if(Ict(1).gt.0) then
         if(Klan.gt.0) then
c
c gravitational harmonics for moon as observed body
            if(Klan.gt.u_mxpl) then
               n2 = Nmzone - 1
               call LVTHAR(Lszhar,Mszhar,Lmzhar,1,1,1,n2,Nszone,
     .                     Mszone, n2)
               n2 = Nmtess*(Nmtess + 1)/2 - 1
               call LVTHAR(Lschar,Mschar,Lmchar,1,1,1,n2,Nstess,
     .                     Mstess, n2)
               call LVTHAR(Lsshar,Msshar,Lmshar,1,1,1,n2,Nstess,
     .                     Mstess, n2)
c
c planet shape is observed
            else if(Npshp1.le.0) then
 
c triaxial ellipsoid
               call LVTHAR(Lszhar,Mszhar,izr2,1,1,1,n2,Nszone,
     .                     Mszone, n2)
               call LVTHAR(Lschar,Mschar,izr2,1,1,1,n2,Nstess,
     .                     Mstess, n2)
               call LVTHAR(Lsshar,Msshar,izr2,1,1,1,n2,Nstess,
     .                     Mstess, n2)
c
c        in the  setups for fourier and grid, set up nszone and nstess
c        so that in cmpar1,cmpar3,nrmict & other places where t3 record
c        is read or written, the l/m vectors get put in the correct
c        place with no overwriting of the "zonal" part by the
c        "tesseral" part of the read.
c
            else if(Nshp.eq.0) then
 
c spherical harmonic expansion
               n2 = Nz - 1
               call LVTHAR(Lszhar,Mszhar,Lz,1,1,1,n2,Nszone,
     .                     Mszone, n2)
               n2 = Nt*(Nt + 1)/2 - 1
               call LVTHAR(Lschar,Mschar,Lc,1,1,1,n2,Nstess,
     .                     Mstess, n2)
               call LVTHAR(Lsshar,Msshar,Ls,1,1,1,n2,Nstess,
     .                     Mstess, n2)
c
c fourier expansion
            else if(Nshp.eq.1) then
               n2 = 122
               call LVTHAR(Lszhar,Mszhar,lfour,1,1,1,n2,Lnfour,
     .                     Mnfour, n2)
               Nszone = 24
               Nstess = 122 - Nszone
c
c grid model
            else if(Nshp.eq.2) then
               n2 = ngdpts
               call LVTHAR(Lszhar,Mszhar,lsh,1,1,1,n2,Lngd,Mngd,
     .                     n2)
               Nszone = 24
               Nstess = Lngd - (24 + 209)
               if(Nstess.gt.0 .and. Nstess.lt.209) Nstess = 209
               if(Nstess.le.0) Nstess = Lngd - 24
               if(Nstess.lt.0) Nstess = 0
            endif
         endif
      endif
 
      return
      end
