      subroutine GFIELD(x0,rot,j,c,s,gm,a,nz,nt,irq,gru,
     .                  dgrudx, dgrudj, dgrudc, dgruds)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 a, clat, cterm, dn1, dn3, DOT, gm, hypot, q1, q2, r, r2,
     .          r3, r5, slat, sterm
      integer   i, ir, irq, k, m, n, n1, nm, nmax, nt, nz
 
c*** end of declarations inserted by spag
 
 
c
c rj cappallo     july 1978   sr gfield
c
      real*10 x(3),j(19),c(54),s(54),gru(3),dgrudx(3,3),
     .          dgrudj(3,19),dgrudc(3,54),dgruds(3,54),ar(20),
     .          cmlon(10),smlon(10),rot(3,3),pn(19),pnp(19),
     .          pnm(54),pnmp(54),dsldx(3),dthdx(3),d2sldx(3,3),
     .          d2thdx(3,3),bth(54),bthp(54),bthpp(54),pnpp(54),
     .          pnmpp(54),x0(3),vector(3),victor(3),tensor(3,3)
c
c     gfield evaluates the gradient of a spherical harmonic potential
c     field, as well as the partial derivatives of the resulting
c     force field wrt the harmonic coefficients and the three
c     coordinate directions.
c
c     i n p u t   p a r a m e t e r s
c       x0  position vector of the point of evaluation
c       rot rotation matrix to go from fundamental coord.system to
c             to body-fixed coords.(ignored unless irq's ten-digit set)
c       j   zonal coefficient array, starts with j2
c       c   tesseral cosine coefficient array, starts with c21
c       s   tesseral   sine coefficient array, starts with s21
c       gm  big g times mass of the gravitating body
c       a   radius of body corresponding to coefficient values
c       nz  maximum index of zonal coefficient array being used
c       nt  highest degree for tesseral coefficients
c       irq request code:
c           irq defines the requested options. it is a packed base-ten
c           number, with digits of 1 or 0 indicating whether or not a
c           particular option is desired. example:irq=101 means include
c           1/r**2 term & supply harmonic partials.
c
c           (power of ten)      (value of 1 implies)
c                   1        include central body force(1/r**2).
c                  10        x0 is supplied in a reference coord.
c                            system: use rot to rotate into body-fixed
c                            coordinates, and rotate results back
c                            into the reference system.
c                 100        return harmonic coefficient partials.
c                1000        return the tensor representing the partial
c                            of the potential gradient wrt the three
c                            coordinate directions.
c
c     r e t u r n   p a r a m e t e r s
c       gru           minus the gradient of the potential, u
c       dgrudx(i,k)   partial of gru(i) wrt x(k)
c       dgrudj(i,n)   partial of gru(i) wrt j(n)
c       dgrudc(i,nm)  partial of gru(i) wrt c(nm)
c       dgruds(i,nm)  partial of gru(i) wrt s(nm)
c
c     kronecker delta function
      real*10 d(3,3)/1._10,3*0._10,1._10,3*0._10,1._10/
      real*10 jterm
      logical*1 noparx, noparh, norot, nocent
      noparx = .true.
      noparh = .true.
      norot  = .true.
      nocent = .true.
 
c d e c o d e  i r q
      i = irq/1000
      if(i.eq.1) noparx = .false.
      ir = irq - 1000*i
      i  = ir/100
      if(i.eq.1) noparh = .false.
      ir = ir - 100*i
      i  = ir/10
      if(i.eq.1) norot = .false.
      ir = ir - 10*i
      if(ir.eq.1) nocent = .false.
 
c get new positon coords., if needed
      if(norot) then
 
         do i = 1, 3
            x(i) = x0(i)
         end do
      else
         call PRODCT(rot,x0,x,3,3,1)
      endif
 
c calculate useful scalar quantities
      r2    = DOT(x,x)
      r     = SQRT(r2)
      r3    = r*r2
      r5    = r3*r2
      ar(1) = a/r
      nmax  = max0(nt,nz)
      do n = 2, nmax
         ar(n) = ar(n - 1)*ar(1)
      end do
      slat = x(3)/r
      clat = SQRT(1._10 - slat**2)
 
c evaluate the necessary legendre polynomials
      call LEGNDR(slat,clat,nz,nt,pn,pnp,pnm,pnmp)
      if(.not. noparx) call LEGND2(slat,clat,nz,nt,pn,pnp,
     .                                pnpp, pnm, pnmp, pnmpp)
      hypot = r*clat
      if(nt.ge.2) then
         cmlon(1) = x(1)/hypot
         smlon(1) = x(2)/hypot
 
c figure multiple angle sines and cosines
         nm = 0
         do n = 2, nt
            n1 = n - 1
            cmlon(n) = cmlon(n1)*cmlon(1) - smlon(n1)*smlon(1)
            smlon(n) = smlon(n1)*cmlon(1) + cmlon(n1)*smlon(1)
            if(.not. (noparx)) then
               do m = 1, n
                  nm = nm + 1
                  bth(nm)   = c(nm)*cmlon(m) + s(nm)*smlon(m)
                  bthp(nm)  = m*(s(nm)*cmlon(m) - c(nm)*smlon(m))
                  bthpp(nm) = -m**2*bth(nm)
               end do
            endif
         end do
      endif
c
c loop over cartesian index k
c predetermine quantities dependent on i and k
      do k = 1, 3
         dsldx(k) = (r2*d(3,k) - x(3)*x(k))/r3
         dthdx(k) = (x(1)*d(2,k) - x(2)*d(1,k))/hypot**2
         gru(k)   = 0._10
         if(.not. (noparx)) then
            do i = 1, k
               d2sldx(i,k) = 3._10*x(3)*x(i)*x(k)
     .                        /r5 - (d(3,k)*x(i) + d(3,i)*x(k) + d(i,k)
     .                        *x(3))/r3
               d2thdx(i,k) = ((x(1)*d(1,k)+x(2)*d(2,k))*(x(2)*d(1,i)-x(
     .                        1)*d(2,i)) + (x(1)*d(1,i)+x(2)*d(2,i))
     .                        *(x(2)*d(1,k)-x(1)*d(2,k)))/hypot**4
               dgrudx(i,k) = 0._10
            end do
         endif
         nm = 0
 
c central body term
         if(.not. (nocent)) then
            gru(k) = -x(k)/r3
            if(.not. (noparx)) then
               do i = 1, k
                  dgrudx(i,k) = (3._10*x(i)*x(k))/r5 - d(i,k)/r3
               end do
            endif
         endif
 
c loop over degree n
         do n = 2, nmax
            dn1 = n + 1
            dn3 = n + 3
            n1  = n - 1
 
c zonal harmonic contributions
            if(n.le.nz) then
               jterm  = ar(n)*(dn1*x(k)*pn(n1)/r3 - pnp(n1)*dsldx(k)/r)
               gru(k) = gru(k) + j(n1)*jterm
               if(.not. noparh) dgrudj(k,n1) = gm*jterm
               if(.not. (noparx)) then
 
c also want gradient of the gradient matrix
                  do i = 1, k
                     dgrudx(i,k) = dgrudx(i,k) + j(n1)*ar(n)
     .                              *(dn1*((d(i,k)*pn(n1)+pnp(n1)
     .                              *(x(k)*dsldx(i)+x(i)*dsldx(k)))
     .                              /r3-dn3*x(i)*x(k)*pn(n1)/r5)
     .                              - (pnp(n1)*d2sldx(i,k)+pnpp(n1)
     .                              *dsldx(i)*dsldx(k))/r)
                  end do
               endif
            endif
 
            if(n.le.nt) then
 
c tesseral cosine and sine terms
               do m = 1, n
                  nm     = nm + 1
                  q1     = ar(n)*(pnmp(nm)*dsldx(k)/r - 
     .                     dn1*x(k)*pnm(nm)/r3)
                  q2     = ar(n)*m*dthdx(k)*pnm(nm)/r
                  cterm  = cmlon(m)*q1 - smlon(m)*q2
                  sterm  = smlon(m)*q1 + cmlon(m)*q2
                  gru(k) = gru(k) + c(nm)*cterm + s(nm)*sterm
 
                  if(.not. (noparh)) then
                     dgrudc(k,nm) = gm*cterm
                     dgruds(k,nm) = gm*sterm
                  endif
 
                  if(.not. (noparx)) then
                     do i = 1, k
                        dgrudx(i,k) = dgrudx(i,k) + ar(n)
     .                                 *(-dn1*((bthp(nm)*pnm(nm)
     .                                 *(x(i)*dthdx(k)+x(k)*dthdx(i))
     .                                 +bth(nm)*pnmp(nm)
     .                                 *(x(i)*dsldx(k)+x(k)*dsldx(i))
     .                                 +d(i,k)*bth(nm)*pnm(nm))
     .                                 /r3-dn3*x(i)*x(k)*bth(nm)*pnm(nm)
     .                                 /r5)
     .                                 + (pnm(nm)*(bthpp(nm)*dthdx(i)
     .                                 *dthdx(k)+bthp(nm)*d2thdx(i,k))
     .                                 +bthp(nm)*pnmp(nm)
     .                                 *(dthdx(k)*dsldx(i)+dthdx(i)
     .                                 *dsldx(k))+bth(nm)
     .                                 *(pnmpp(nm)*dsldx(i)*dsldx(k)
     .                                 +pnmp(nm)*d2sldx(i,k)))/r)
                     end do
                  endif
               end do
            endif
         end do
         gru(k) = gm*gru(k)
         if(.not. (noparx)) then
            do i = 1, k
               dgrudx(i,k) = gm*dgrudx(i,k)
               dgrudx(k,i) = dgrudx(i,k)
            end do
         endif
      end do
      if(.not. (norot)) then
 
c rotate results back into original coord. system
         call PRODCT(rot,gru,vector,-3,3,1)
         do i = 1, 3
            gru(i) = vector(i)
         end do
         if(.not. (noparh)) then
 
c rotate harmonic coeff. partials
            if(nz.ge.2) then
               n1 = nz - 1
               do n = 1, n1
                  call PRODCT(rot,dgrudj(1,n),vector,-3,3,1)
                  do i = 1, 3
                     dgrudj(i,n) = vector(i)
                  end do
               end do
            endif
 
            if(nt.ge.2) then
               do n = 1, nm
                  call PRODCT(rot,dgrudc(1,n),vector,-3,3,1)
                  call PRODCT(rot,dgruds(1,n),victor,-3,3,1)
                  do i = 1, 3
                     dgrudc(i,n) = vector(i)
                     dgruds(i,n) = victor(i)
                  end do
               end do
            endif
         endif
 
c transform gradient partial wrt position tensor
         if(.not. (noparx)) then
 
c use similarity transform to rotate
            call PRODCT(rot,dgrudx,tensor,-3,3,3)
            call PRODCT(tensor,rot,dgrudx,3,3,3)
         endif
      endif
 
      return
      end
