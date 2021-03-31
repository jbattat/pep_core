      subroutine VETIDE(nvel, gm, g, cphi, h, l, r, a, dadr, dadphi,
     .                  dadlam, x, dxdh, dxdl, delr, delphi, dellam)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 a1, cphi, dellam, delldt, delpdt, delphi, delr, delrdt,
     .          DOT, dudlam, dudldt, dudpdt, dudphi, dudt, g, gm, h, r1,
     .          r3, r5
      real*10 rho1, rho3, rho5, u
      integer   i, nv, nvel
 
c*** end of declarations inserted by spag
 
 
c
c rj cappallo  may 1977  sr vetide
c mathematical computation routine for vector solid-body tides
c
      real*10 l, r(6), a(6), dadr(6), dadphi(6), dadlam(6), x(6), 
     .          dxdh(6), dxdl(6), rho(6)
c
c       i n p u t   p a r a m e t e r s
c          nvel  =0 tide only
c                =1 tide and tide rate
c            gm  big g * mass(perturber)
c             g  local gravity acceleration
c          cphi  cosine(site latitude)
c             h  love number - vertical
c             l  love number - horizontal
c             r  vector from c.o.m. body to perturber's center
c             a  vector from c.o.m. body to site
c          dadr  vectors of partial derivatives of the
c        dadphi    site vector w.r.t. body-fixed
c        dadlam    spherical coordinates. (phi=latitude)
c
c
c       r e t u r n   p a r a m e t e r s
c             x  tide vector (and rate if nvel=1)
c          dxdh  partials of tide vector
c          dxdl   w.r.t. love numbers (6-vector if nvel=1)
c
c     calculate useful quantities a priori
      nv = 3
      if( nvel .ge. 1 ) nv = 6
      do i = 1, nv
         rho(i) = r(i) - a(i)
      end do
      rho1 = SQRT(DOT(rho,rho))
      rho3 = rho1**3
      r1   = SQRT(DOT(r,r))
      r3   = r1**3
      a1   = SQRT(DOT(a,a))
 
      u = gm*(1._10/rho1 - 1._10/r1 - DOT(r,a)/r3)
      dudphi = gm*(DOT(rho,dadphi)/rho3 - DOT(r,dadphi)/r3)
      dudlam = gm*(DOT(rho,dadlam)/rho3 - DOT(r,dadlam)/r3)
 
c calculate body-fixed displacements / love numbers
      delr   = u/g
      delphi = dudphi/g
      dellam = dudlam/g/cphi
c form terms for cartesian tide-comps
c and partials w.r.t. love numbers
      do i = 1, 3
         dxdh(i) = delr*dadr(i)
         dxdl(i) = (delphi*dadphi(i) + dellam*dadlam(i)/cphi)/a1
      end do
      if( nvel .gt. 0 ) then
 
c t i d e   r a t e  -useful quantities
         r5   = r3*r1*r1
         rho5 = rho3*rho1*rho1
 
c temporal rates of potential gradients
         dudt   = gm*((DOT(r,rho(4))-DOT(r(4),a))/r3 - DOT(rho,rho(4))
     .            /rho3 + 3._10*DOT(r,a)*DOT(r,r(4))/r5)
         dudpdt = gm*((DOT(rho(4),dadphi)+DOT(rho,dadphi(4)))
     .            /rho3 - 3._10*DOT(rho,rho(4))*DOT(rho,dadphi)
     .            /rho5 - (DOT(r(4),dadphi)+DOT(r,dadphi(4)))
     .            /r3 + 3._10*DOT(r,r(4))*DOT(r,dadphi)/r5)
         dudldt = gm*((DOT(rho(4),dadlam)+DOT(rho,dadlam(4)))
     .            /rho3 - 3._10*DOT(rho,rho(4))*DOT(rho,dadlam)
     .            /rho5 - (DOT(r(4),dadlam)+DOT(r,dadlam(4)))
     .            /r3 + 3._10*DOT(r,r(4))*DOT(r,dadlam)/r5)
 
c calculate body-fixed tide rates / love numbers
         delrdt = dudt/g
         delpdt = dudpdt/g
         delldt = dudldt/g/cphi
 
c form terms for tide-rate cartesian components
         do i = 4, 6
            dxdh(i) = delrdt*dadr(i - 3) + delr*dadr(i)
            dxdl(i) = (delpdt*dadphi(i-3) + delphi*dadphi(i)
     .                + (delldt*dadlam(i-3)+dellam*dadlam(i))/cphi)/a1
         end do
      endif
 
c combine terms for tides, stuff partials
      do i = 1, nv
         x(i) = h*dxdh(i) + l*dxdl(i)
      end do
 
      return
      end
