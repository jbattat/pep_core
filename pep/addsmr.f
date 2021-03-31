      subroutine ADDSMR(b, s, h, sh, g, np, npnp, nptr, dptr, smat,
     .                  bthts, nd, a6, wr, wi)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ij, ik, il, ipt, iptrow, j, jpt, jptrow, k, kj, kl,
     .          kpt, l, lpt, lptrow, m, nd
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1974  subroutine addsmr
c
c         smear the covariance of noisy parameters
c         for algorithms used in this routine see memo of 21 august 1973
c         by r. d. reasenberg on "kalman filter equations".
c
c        this routine replaces p**-1 by m**-1 as outline on pages 3-6 of
c        memo actual procedure is steps 1-5 on page 6 of memo
c         parameters
      integer*4 np
      integer*4 npnp, smat
      real*10   b(1)
      real*10   bthts(nd, npnp)
      real*10 a6(npnp, npnp), wr(npnp), wi(npnp)
      logical*4 las/.true./
      logical*4 ls/.false./
      logical*4 ltm/.true./
      real*10 s(20, 20)
      real*10   h(npnp, npnp), sh(npnp, npnp), g(np, npnp)
      integer*2 dptr(np)
      integer*2 nptr(npnp)
c
c local
      real*10   d, dc
      real*10   has(20, 20)
 
      if( ls ) then
         do i = 1, npnp
            do j = 1, npnp
               a6(i, j) = s(i, j)
            end do
         end do
         call PAGCHK(60, 7, 0)
         call EIGDIS(a6, npnp, 'S  ', 1, wr, wi)
      endif
c
c initialize h
      do i = 1, npnp
         do j = 1, npnp
            h(j, i) = 0.0
         end do
      end do
c
c
c         1.  h = (e + a*s)**-1
c             a.  h = a*s
c                 (save h in sh)
      do i = 1, npnp
         ipt    = nptr(i)
         iptrow = ipt*(ipt - 1)/2
         do j = 1, npnp
            do k = 1, npnp
               kpt = nptr(k)
               ik  = iptrow + kpt
               if( kpt .gt. ipt ) ik = kpt*(kpt - 1)/2 + ipt
               h(i, j) = h(i, j) + b(ik)*s(k, j)
            end do
            sh(i, j) = h(i, j)
c
c want matrix s*a which is transpose of a*s
c
            a6(j, i) = h(i, j)
         end do
      end do
c
c calculate coherence time (use g for work)
c
      if( ltm ) call FCOHER(a6, npnp, wr, wi, g)
c
c b.  h = h + e
c (save h in sh)
      do i = 1, npnp
         h(i, i)  = h(i, i) + 1.0
         sh(i, i) = h(i, i)
         do j = 1, npnp
            a6(i, j) = h(i, j)
         end do
      end do
      if( las ) then
         call PAGCHK(60, 7, 0)
         call EIGDIS(a6, npnp, 'E+AS', 1, wr, wi)
c
c c.  h = h**-1  (use g for work)
         call ADDINV(sh, h, g, npnp, np)
      endif
c
c initialize sh
      do i = 1, npnp
         do j = 1, npnp
            sh(j, i) = 0.0
         end do
      end do
c
c
c 2.  a = h*a
c a.  sh = h*a
      do i = 1, npnp
         do j = 1, npnp
            jpt    = nptr(j)
            jptrow = jpt*(jpt - 1)/2
            do k = 1, npnp
               kpt = nptr(k)
               kj  = jptrow + kpt
               if( kpt .gt. jpt ) kj = kpt*(kpt - 1)/2 + jpt
               sh(i, j) = sh(i, j) + h(i, k)*b(kj)
            end do
         end do
      end do
c
c b.  a = sh
      do i = 1, npnp
         ipt    = nptr(i)
         iptrow = ipt*(ipt - 1)/2
         do j = 1, i
            jpt = nptr(j)
            ij  = iptrow + jpt
            if( jpt .gt. ipt ) ij = jpt*(jpt - 1)/2 + ipt
            b(ij) = sh(i, j)
         end do
      end do
c
c form has (ha now in sh)
      do i = 1, npnp
         do j = 1, npnp
            has(i, j) = 0.0
            do k = 1, npnp
               has(i, j) = has(i, j) + sh(i, k)*s(k, j)
            end do
         end do
      end do
c
c
c if all params are noise, finished
      m = np - npnp
      if( m .gt. 0 ) then
c
c reinitialize sh
         do i = 1, npnp
            do j = 1, npnp
               sh(j, i) = 0.0
            end do
         end do
         do i = 1, npnp
            do j = 1, np
               g(j, i) = 0.0
            end do
         end do
c
c
c 3.  sh = s*h
         do i = 1, npnp
            do j = 1, npnp
               do k = 1, npnp
                  sh(i, j) = sh(i, j) + s(i, k)*h(k, j)
               end do
            end do
         end do
c
c
c 4.  c = c - bt*sh*b
         do i = 1, m
            ipt    = dptr(i)
            iptrow = ipt*(ipt - 1)/2
            do l = 1, i
               lpt    = dptr(l)
               lptrow = lpt*(lpt - 1)/2
               dc     = 0.0
               do j = 1, npnp
                  d = 0.0
                  do k = 1, npnp
                     kpt = nptr(k)
                     kl  = lptrow + kpt
                     if( kpt .gt. lpt ) kl = kpt*(kpt - 1)/2 + lpt
                     d = d + sh(j, k)*b(kl)
                  end do
                  jpt = nptr(j)
                  ij  = iptrow + jpt
                  if( jpt .gt. ipt ) ij = jpt*(jpt - 1)/2 + ipt
                  dc = dc + b(ij)*d
               end do
               il = iptrow + lpt
               if( lpt .gt. ipt ) il = lptrow + ipt
               b(il) = b(il) - dc
            end do
         end do
c
c
c 5.  bt = bt*ht
c a.  g = bt*ht
         do i = 1, m
            ipt    = dptr(i)
            iptrow = ipt*(ipt - 1)/2
            do j = 1, npnp
               do k = 1, npnp
                  kpt = nptr(k)
                  ik  = iptrow + kpt
                  if( kpt .gt. ipt ) ik = kpt*(kpt - 1)/2 + ipt
                  g(i, j) = g(i, j) + b(ik)*h(j, k)
               end do
            end do
         end do
c
c
c form and save bthts
         do i = 1, m
            do j = 1, npnp
               bthts(i, j) = 0.0
               do k = 1, npnp
                  bthts(i, j) = bthts(i, j) + g(i, k)*s(k, j)
               end do
            end do
         end do
c
c b.  bt = g
         do i = 1, m
            ipt    = dptr(i)
            iptrow = ipt*(ipt - 1)/2
            do j = 1, npnp
               jpt = nptr(j)
               ij  = iptrow + jpt
               if( jpt .gt. ipt ) ij = jpt*(jpt - 1)/2 + ipt
               b(ij) = g(i, j)
            end do
         end do
      endif
c
c restore has to sh
      do i = 1, npnp
         do j = 1, npnp
            sh(i, j) = has(i, j)
         end do
      end do
 
      return
      end
