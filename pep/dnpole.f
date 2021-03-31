      subroutine DNPOLE(n,narg,csfarg,nt,nfarg,bothcs,farg,sum)
 
      implicit none
c
c      r.w.king         mar 1980        subroutine dnpole
c      calculate diurnal polar motion
c      copied with minor modifications from vlbi3 subroutine written ori
c      by d.s.robertson
c
c arguments
      integer*4 n,narg,nt
      logical bothcs
      integer*2 nfarg(6,n)
      real*10 csfarg(nt,n),sum(2,nt),farg(6)
 
c commons
      include 'funcon.inc'

c local
      real*10 cosf, sinf, sumarg
      integer   i, j, k
 
      do j = 1, nt
         sum(1,j) = 0.0_10
         sum(2,j) = 0.0_10
      end do
      do i = 1, n
         sumarg = farg(1)*nfarg(1,i)
         do j = 2, narg
            sumarg = sumarg + farg(j)*nfarg(j,i)
         end do
         sinf = SIN(sumarg)
         if(bothcs) cosf = COS(sumarg)
         do k = 1, nt
            sum(1,k) = sum(1,k) + csfarg(k,i)*sinf*Convds
            if(bothcs) sum(2,k) = sum(2,k) + csfarg(k,i)*cosf*Convds
         end do
      end do
      return
      end
