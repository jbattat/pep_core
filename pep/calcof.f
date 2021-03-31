      subroutine CALCOF

      implicit none
c
c subroutine ''calcof''    w.b.smith, 11-15-68
c calcof to calculate the coefficients for numerical integration

c array dimensions
      include 'globdefs.inc'

c common
      include 'adams.inc'
      include 'coeff.inc'

c local
      real*10 xk,zk
      integer   i,ii,j,k,kk,kv,l,loop,m,nc,nn,np,npd,nterm
      real*10 r(2,15),xa(2,15),s(2,15),z(2,15),
     .          xb(2,15),v(2,15),
     .          temp(2),binom(15,15),h(15),
     .          bf(15),wf(15),f(15),g(15),ss(15)

      np  = Npredt+1
      nc  = Ncoret+1
      npd = Itype
      do ii = 1,90
         A(ii) = 0._10
         end do
      binom(1,1) = 1._10
      binom(2,1) = 1._10
      binom(2,2) = 1._10
      do i = 3,15
         binom(i,1) = 1._10
         loop = i-1
         do k = 2,loop
            binom(i,k) = binom(i-1,k-1)+binom(i-1,k)
            end do
         binom(i,i) = 1._10
         end do
      do k = 1,15
         xk = k
         call XLOAD8(1._10)
         call XDIV8(xk)
         call XSTORE(r(1,k))
         end do
      call XLOAD8(1._10)
      call XSTORE(xa(1,1))
      do k = 2,15
         nterm = k-1
         call XLOAD8(1._10)
         call XSTORE(xa(1,k))
         do l = 1,nterm
            m = k-l+1
            call XLOAD(r(1,m))
            call XMUL(xa(1,l))
            call XSTORE(temp)
            call XLOAD(xa(1,k))
            call XSUB(temp)
            call XSTORE(xa(1,k))
            end do
         end do
      do k = 1,13,2
         ss(k)   = 1._10
         ss(k+1) = -1._10
         call XLOAD8(1._10)
         call XSTORE(s(1,k))
         call XLOAD8(-1._10)
         call XSTORE(s(1,k+1))
         end do
      ss(15) = 1._10
      call XLOAD8(1._10)
      call XSTORE(s(1,15))
      do k = 1,np
         call XLOAD8(0._10)
         call XSTORE(z(1,k))
         j = k
         do i = j,np
            call XLOAD(xa(1,i))
            call XMUL8(binom(i,k))
            call XMUL(s(1,k))
            call XADD(z(1,k))
            call XSTORE(z(1,k))
            end do
         end do
      do k = 1,np
         call XLOAD(z(1,k))
         call STORND(C(k))
         end do
      call XLOAD8(1._10)
      call XSTORE(xb(1,1))
      do k = 2,15
         nterm = k-1
         call XLOAD8(0._10)
         call XSTORE(xb(1,k))
         do l = 1,nterm
            m = k-l+1
            call XLOAD(r(1,m))
            call XMUL(xb(1,l))
            call XSTORE(temp)
            call XLOAD(xb(1,k))
            call XSUB(temp)
            call XSTORE(xb(1,k))
            end do
         end do
      do k = 1,nc
         call XLOAD8(0._10)
         call XSTORE(v(1,k))
         j = k
         do i = j,nc
            call XLOAD(xb(1,i))
            call XMUL8(binom(i,k))
            call XMUL(s(1,k))
            call XADD(v(1,k))
            call XSTORE(v(1,k))
            end do
         end do
      do k = 1,nc
         call XLOAD(v(1,k))
         call STORND(D(k))
         end do
      do k = 1,np
         call XLOAD(xa(1,k))
         call STORND(A(k))
         end do
      do k = 1,nc
         call XLOAD(xb(1,k))
         call STORND(h(k))
         call STORND(B(k))
         end do
      do k = 1,np
         W(k) = 0._10
         if(k.ne.1) then
            kv = k
            do kk = 1,k
               W(k) = W(k)+h(kk)*h(kv)
               kv   = kv-1
               end do
            E(k) = E(k-1)+W(k)
         else
            W(1) = h(1)**2
            E(1) = W(1)
         endif
         end do
      do k = 1,np
         call XLOAD(s(1,k))
         call STORND(h(k))
         end do
      do k = 1,np
         Cc(k) = 0._10
         j     = k
         do i = j,np
            Cc(k) = Cc(k)+h(k)*binom(i,k)*E(i)
            end do
         end do
      do k = 1,nc
         Dd(k) = 0._10
         j     = k
         do i = j,nc
            Dd(k) = Dd(k)+h(k)*binom(i,k)*W(i)
            end do
         end do

c compute coefficients for y(n)' corrector
      do k = 1,npd
         zk = k-1
         if(k.ne.1) then
            Pp(k) = 1._10/zk
         else
            Pp(k) = 0._10
         endif
         end do
      do k = 1,npd
         P(k) = 0._10
         j    = k
         do i = j,npd
            P(k) = P(k)+h(k)*binom(i,k)*Pp(i)
            end do
         end do

c compute coefficients for y(n+1)' predictor
      do k = 1,npd
         zk = k-1
         if(k.ne.1) then
            Pd(k) = Pd(k-1)+1._10/zk
         else
            Pd(k) = 0._10
         endif
         end do
      do k = 1,npd
         Pf(k) = 0._10
         j     = k
         do i = j,npd
            Pf(k) = Pf(k)+h(k)*binom(i,k)*Pd(i)
            end do
         end do
      do k = 1,15
         bf(k) = B(k)*ss(k)
         end do
      do k = 1,15
         wf(k) = 0._10
         if(k.ne.1) then
            kv = k
            do kk = 1,k
               wf(k) = wf(k)+bf(kk)*bf(kv)
               kv    = kv-1
               end do
         else
            wf(1) = bf(1)**2
         endif
         end do
      do k = 1,15
         zk   = k+1
         f(k) = ss(k)/zk
         end do
      do k = 1,15
         g(k) = 0._10
         kv   = k
         do kk = 1,k
            g(k) = g(k)+f(kk)*wf(kv)
            kv   = kv-1
            end do
         end do
      do k = 1,np
         A(k) = 0._10
         j    = k
         nn   = 0
         do i = j,np
            nn   = nn+1
            A(k) = A(k)+ss(nn)*binom(i,k)*g(i)
            end do
         end do
      return
      end
