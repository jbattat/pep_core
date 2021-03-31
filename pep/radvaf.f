      subroutine RADVAF(cond,x0,dx0,goose,kind)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 caz, cdecl, cfl, cond4, cra, goose, saz, sdecl, sfl, sra
      integer   i, j, kind
 
c*** end of declarations inserted by spag
 
 
 
c m.e. ash     sept 1971      subroutine radvaf
      real*10 cond(6),x0(6),dx0(6,6)
 
      include 'funcon.inc'
c
c  input quantities
c       cond(1) = distance from central body
c       cond(2) = right ascension of position (reference epoch)
c                 (between 0 and 360 deg)
c       cond(3) = declination of body (reference epoch)
c                 (between -90 and 90 deg)
c       cond(4) = magnitude of velocity minus velocity in
c                 circular orbit of radius cond(1)
c       cond(5) = azimuth of velocity vector on plane normal to
c                 radius vector (between 0 and 360 deg)
c       cond(6) = flight path angle measured from radius vector to
c                 velocity vector (between 0 and 180 deg)
c       goose   = sqrt(gm)
c       kind    = 0 if only coordinates needed, 1 if also partials
c  output quantities
c       x0(1-6) = cartesian components of position and velocity
c                 referred to the mean equinox and equator of ref.epoch
c       dx0(i,j)= partial derivative of x0(i) with respect to cond(j)
c
      cra   = cond(2)*Convd
      sra   = SIN(cra)
      cra   = COS(cra)
      cdecl = cond(3)*Convd
      sdecl = SIN(cdecl)
      cdecl = COS(cdecl)
      caz   = cond(5)*Convd
      saz   = SIN(caz)
      caz   = COS(caz)
      cfl   = cond(6)*Convd
      sfl   = SIN(cfl)
      cfl   = COS(cfl)
 
      x0(1) = cond(1)*cra*cdecl
      x0(2) = cond(1)*sra*cdecl
      x0(3) = cond(1)*sdecl
      cond4 = cond(4) + goose/SQRT(cond(1))
      x0(4) = cond4*(cfl*cra*cdecl - sfl*(caz*cra*sdecl+saz*sra))
      x0(5) = cond4*(cfl*sra*cdecl - sfl*(caz*sra*sdecl-saz*cra))
      x0(6) = cond4*(cfl*sdecl + sfl*caz*cdecl)
      if(kind.eq.0) return
 
      do i = 1, 3
         dx0(i,1)   = x0(i)/cond(1)
         dx0(i+3,1) = -0.5_10*goose*x0(i+3)
     .                   /(cond4*cond(1)*SQRT(cond(1)))
         dx0(i+3,4) = x0(i + 3)/cond4
         do j = 4, 6
            dx0(i,j) = 0.0_10
         end do
      end do
 
      dx0(1,2) = -x0(2)
      dx0(2,2) = x0(1)
      dx0(3,2) = 0.0_10
      dx0(4,2) = -x0(5)
      dx0(5,2) = x0(4)
      dx0(6,2) = 0.0_10
 
      dx0(1,3) = -cond(1)*cra*sdecl
      dx0(2,3) = -cond(1)*sra*sdecl
      dx0(3,3) = cond(1)*cdecl
      dx0(4,3) = -cond4*(cfl*cra*sdecl + sfl*caz*cra*cdecl)
      dx0(5,3) = -cond4*(cfl*sra*sdecl + sfl*caz*sra*cdecl)
      dx0(6,3) = cond4*(cfl*cdecl - sfl*caz*sdecl)
 
      dx0(4,5) = cond4*sfl*(saz*cra*sdecl - caz*sra)
      dx0(5,5) = cond4*sfl*(saz*sra*sdecl + caz*cra)
      dx0(6,5) = -cond4*sfl*saz*cdecl
 
      dx0(4,6) = -cond4*(sfl*cra*cdecl + cfl*(caz*cra*sdecl+saz*sra))
      dx0(5,6) = -cond4*(sfl*sra*cdecl + cfl*(caz*cra*sdecl-saz*cra))
      dx0(6,6) = cond4*(cfl*caz*cdecl - sfl*sdecl)
 
      return
      end
