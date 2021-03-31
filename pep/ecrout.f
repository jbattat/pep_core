      subroutine ECROUT(a, r, d, eps, ni, mi, nn, mm, ind)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10   a, d, dp, dt, eps, p, r, s, t
      integer   i, ic, ii, ind, j, k, km, kp, m, mi, mm, n, ni, nn
 
c*** end of declarations inserted by spag
 
 
c     subroutine dcrout(a,r,d,eps,ni,mi,nn,mm,ind)
c
c     real*10 matrix inversion by crout elimination.
c        extended precision  paul macneil oct., 1978
c
c     refer to dcd memo number 688.
c
c----------------------------------------------------------------------
c
c     implicit real*8 (a,r,d,e,s,t,p,c)
c     real*10   f(7, 7), deltav(7, 7), ratio(7, 7), aorig(7, 7)
      dimension a(nn, nn), r(nn, mm)
      ind = 0
      n   = ni
      m   = mi
      if(m.lt.0) then
         m = n
         do i = 1, n
            do j = 1, n
               r(i, j) = 0.0
               end do
            r(i, i) = 1.0
            end do
      endif
      ic = 0
      ii = 0
      t  = abs(a(1,1))
      do i = 2, n
         if(t.lt.abs(a(i,1))) then
            ii = i
            t  = abs(a(i,1))
         endif
         end do
      if(ii.ne.0) then
         ic = ic + 1
         if(m.ne.0) then
            do j = 1, m
               s = r(1, j)
               r(1, j)  = r(ii, j)
               r(ii, j) = s
               end do
         endif
         do j = 1, n
            s = a(1, j)
            a(1, j)  = a(ii, j)
            a(ii, j) = s
            end do
      endif
      p = a(1, 1)
      if(abs(p).gt.eps) then
         do j = 2, n
            a(1, j) = a(1, j)/p
            end do
         if(m.ne.0) then
            do j = 1, m
               r(1, j) = r(1, j)/p
               end do
         endif
         do k = 2, n
            km = k - 1
            t  = -1.0
            do i = k, n
               dp = a(i, k)
               do j = 1, km
                  dp = dp - a(i, j)*a(j, k)
                  end do
               a(i, k) = dp
               if(t.lt.abs(a(i,k))) then
                  t  = abs(a(i,k))
                  ii = i
               endif
               end do
            if(ii.ne.k) then
               ic = ic + 1
               if(m.ne.0) then
                  do j = 1, m
                     s = r(k, j)
                     r(k, j)  = r(ii, j)
                     r(ii, j) = s
                     end do
               endif
               do j = 1, n
                  s = a(k, j)
                  a(k, j)  = a(ii, j)
                  a(ii, j) = s
                  end do
            endif
            dt = a(k, k)
            if(abs(dt).le.eps) go to 100
            p = p*dt
            if(k.ne.n) then
               kp = k + 1
               do j = kp, n
                  dp = a(k, j)
                  do i = 1, km
                     dp = dp - a(k, i)*a(i, j)
                     end do
                  a(k, j) = dp/dt
                  end do
            endif
            if(m.ne.0) then
               do j = 1, m
                  dp = r(k, j)
                  do i = 1, km
                     dp = dp - a(k, i)*r(i, j)
                     end do
                  r(k, j) = dp/dt
                  end do
            endif
            end do
         if(mod(ic,2).ne.0) p = -p
         d = p
         if(m.ne.0) then
            ii = n
            do k = 2, n
               kp = ii
               ii = ii - 1
               do j = 1, m
                  dp = r(ii, j)
                  do i = kp, n
                     dp = dp - a(ii, i)*r(i, j)
                     end do
                  r(ii, j) = dp
                  end do
               end do
         endif
         return
      endif
  100 ind = 1
      d   = 0.0
 
      return
      end
