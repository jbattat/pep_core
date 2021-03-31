      subroutine MATCH(coeff,names,numpar,xnom,apcoef,apnams,w,
     .                 nappar,xap,diagon,covar,apest,offset,nstop,
     .                 iskale,rhs,sqrtb,sumsq)
 
      implicit none
c
c d. white  april 1973  subroutine match
c
c inserts apriori values in appropriate position in coefficient
c array and accumulates statistics
c
c parameters
      character*8 names(2,1),apnams(2,1)
      real*10 coeff(1),apcoef(1),w,sumsq
      real*10 xnom(1),xap(1),iskale(1),rhs(1)
      integer*4 numpar,nappar,nstop
      logical*4 apest,covar,diagon,offset,sqrtb

c common
      include 'inodta.inc'
c
c local variables
c dimensions of pr, pc, iu, and pointr depend on maximum number
c of apriori parameters in one matrix
      integer*4 i,ipt,ipt1,iptdia,iptrow,irow,j,jpt,k
      real*10 pr(100),pc(100),b(1),pb(1)
      integer*2 iu(100),pointr(100)
      real*10 isk2,ww,rhsadd,apci

c external function
      integer*4 LEG
 
      if(offset .and. nappar.le.1) then
         write(Iout,300)
  300 format('0OFFSET=.TRUE. REQUIRES MORE THAN ONE APRIORI PARAMETER, E
     .RROR IN MATCH')
         if(Mout.gt.0) write(Mout,300)
         nstop = nstop + 1
         return
      endif
 
      ww = w
      if(covar) ww = 1._10/w
 
c find apnames among names of adjustable params
      k = nappar
      do i = 1, nappar
         do j = 1, numpar
            if(LEG(16,1,apnams(1,i),1,names(1,j)).eq.0) then

c form difference vector  xap - xnom
               if(apest.or.offset) xap(i) = xap(i)/iskale(j) - xnom(j)
c
c see if offset or diagonal
               if(diagon .and. .not.offset) then
c
c if diagonal, first get index to jth parameter coeff
c and square of scale factor
                  k    = j*(j+1)/2
                  isk2 = iskale(j)*iskale(j)
                  apci = apcoef(i)
                  if(sqrtb) apci = apci*apci
c if covariance, invert to get coeff
                  if(covar) apci = 1._10/apci

c add ith apriori value to jth coeff
c and add estimate to rhs, and update effective sum-squared residual
                  coeff(k) = coeff(k) + apci*isk2*ww
                  if(apest) then
                     rhsadd = apci*isk2*ww*xap(i)
                     rhs(j) = rhs(j) + rhsadd
                     sumsq  = sumsq + xap(i)*rhsadd
                  endif
               else
c
c if non diagonal, make pointer to jth param for ith apriori par
                  pointr(i) = j
               endif
               goto 250
            endif
         end do
c
c have not found ith apname
         write(Iout,240) apnams(1,i),apnams(2,i)
  240    format('0', 2(a8,1x),'NOT BEING ADJUSTED')
         if(Mout.gt.0) write(Mout,240) apnams(1,i),apnams(2,i)
         nstop = nstop + 1
         k     = k - 1
         if(k.le.0) return
  250 end do
      if(.not.offset) then
c
c if diagonal apriori matrix then all finished
         if(diagon) return
c
c if nondiagonal covariance matrix then invert to get coeff
         if(covar) then
            call SYMINV(apcoef,b,nappar,0,pr,pc,iu,pb,i)
            if(i.ge.1) nstop = nstop + 1
         endif
c
c place apcoeff values in coeff array using pointers
         do i = 1, nappar
            ipt    = pointr(i)
            iptrow = ipt*(ipt - 1)/2
            irow   = i*(i - 1)/2
            do j = 1, i
               jpt = pointr(j)
               k   = iptrow + jpt
               if(jpt.gt.ipt) k = jpt*(jpt - 1)/2 + ipt
               coeff(k) = coeff(k) + apcoef(irow + j)*iskale(ipt)
     .                    *iskale(jpt)*ww
            end do
         end do
c
         if(apest) then
c
c add in rhs contribution
            do i = 1, nappar
               ipt  = pointr(i)
               irow = i*(i - 1)/2
               do j = 1, nappar
                  jpt = pointr(j)
                  k   = irow + j
                  if(j.gt.i) k = j*(j - 1)/2 + i
                  rhsadd = apcoef(k)*iskale(ipt)*iskale(jpt)*xap(j)*ww
                  rhs(ipt) = rhs(ipt) + rhsadd
                  sumsq = sumsq + rhsadd*xap(i)
               end do
            end do
         endif

      else
c
c a priori constraint is on offset from first parameter in packet
c always diagonal
c
         ipt1   = pointr(1)
         iptdia = ipt1*(ipt1 + 1)/2
         do i = 2, nappar
            apci   = apcoef(i)
            if(sqrtb) apci = apci**2
            if(covar) apci = 1._10/apci
            ipt    = pointr(i)
            iptrow = ipt*(ipt - 1)/2
            coeff(iptdia) = coeff(iptdia) + apci*iskale(ipt1)**2*ww
            k = iptrow + ipt1
            if(ipt1.gt.ipt) k = iptdia - ipt1 + ipt
            coeff(k) = coeff(k) - apci*iskale(ipt1)*iskale(ipt)*ww
            k = iptrow + ipt
            coeff(k)  = coeff(k) + apci*iskale(ipt)**2*ww
            rhsadd    = (xap(1)*iskale(ipt1)-xap(i)*iskale(ipt))*apci*ww
            rhs(ipt1) = rhs(ipt1) + rhsadd*iskale(ipt1)
            rhs(ipt) = rhs(ipt) - rhsadd*iskale(ipt)
            sumsq = sumsq + rhsadd*(xap(1)*iskale(ipt1)-
     .       xap(i)*iskale(ipt))
         end do
      endif
      return
      end
