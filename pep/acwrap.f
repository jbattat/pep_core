      subroutine ACWRAP(coeff,rhs,xnom,numpar,aptitl,anyest,alldg,
     .                  nstop,anyoff,sumsq)
 
      implicit none
c
c d. white  april 1973  subroutine acwrap
c jfc/kcl january, 1980
c
c write compressed coefficient matrix, rhs contribution,
c and apriori estimate vector in format of saved normal eqn.
c
c parameters
      character*8 aptitl(10)
      real*10 coeff(1),rhs(1),xnom(1),sumsq
      integer*4 numpar,nstop
      logical*4 anyest,anyoff,alldg
c
c common
      include 'aprtbf.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
 
c shared work area
      common/WRKCOM/ Pr(1000),Pc(1000),Pointr(1000)
      real*10 Pr,Pc
      integer*2 Pointr
      real*10 pb,buff(1)
      equivalence (Pr,buff)
c
c local
      integer   i,ic,ipt,iptrow,irow,j,k,l,l0,nsave,nnom
      real*10 dmnus1/-1._10/
      integer*4 minus1/-1/
c
c build pointer, compressed rhs & xnom
      nsave = 0
      do i = 1, numpar
         k = i*(i + 1)/2
         if(coeff(k) .ne. 0._10) then
            nsave = nsave + 1
            Pointr(nsave) = i
            rhs(nsave)    = rhs(i)
            xnom(nsave)   = xnom(i)
         endif
      end do
      if(nsave .eq. 0) return
c
c build compressed coeff array
      do i = 1, nsave
         ipt    = Pointr(i)
         iptrow = ipt*(ipt - 1)/2
         irow   = i*(i - 1)/2
         do j = 1, i
            coeff(irow + j) = coeff(iptrow + Pointr(j))
         end do
      end do
c
c write title
      write(Ibuf2) aptitl
c
c write a priori sum-squared residual contribution
      write(Ibuf2) minus1,sumsq
c
c write nsave rows plus pointer vector
      write(Ibuf2) nsave,(Pointr(i),i=1,nsave)
c
c write right hand side contribution  (compressed)
      write(Ibuf2) nsave,(rhs(i),i=1,nsave)
c
c write rows of compressed coeff array
      l0 = 0
      do i = 1, nsave
         l  = l0
         ic = 1
         do j = 1, nsave
            l = l + ic
            buff(j) = coeff(l)
            if(j .ge. i) ic = j
         end do
         l0 = l0 + i
         write(Ibuf2) i, (buff(j), j = 1, nsave)
      end do
c
c calculate apriori estimate vector:  invert coeff,
c multiply by rhs, add nominal vector
c not possible with offset apriori, matrix is singular
      nnom=nsave
      if(.not. anyoff) then
         if(anyest) then
            if(.not. alldg) then
               call SYMINV(coeff, rhs, nsave, 1, Pr, Pc, Pointr, pb, i)
               if(i.gt.0) then
                  nnom=-1
                  if(Ict(1).gt.1) then
                     write(Iout, 10)
   10                format('0SYMINV ERROR IN ACWRAP. CANNOT ITERATE.')
                     nstop=nstop+1
                     if(Mout.gt.0) write(Mout, 10)
                  endif
               else
                  do i = 1, nsave
                     xnom(i) = xnom(i) + rhs(i)
                  end do
               endif
            else
               l = 0
               do i = 1, nsave
                  l = l + i
                  xnom(i) = xnom(i) + rhs(i)/coeff(l)
               end do
            endif
         endif
c
c write out estimate vector
         write(Ibuf2) nnom,(xnom(i),i=1,nsave)
      else if(Ict(1).le.1) then
c
c write nsave=-1,xnom=-1 to indicate no nominals
         write(Ibuf2) minus1, dmnus1
      else
         write(Iout, 50)
   50    format('0CAN''T ITERATE WITH OFFSET A PRIORI, ERROR IN ACWRAP')
         nstop = nstop + 1
         if(Mout .gt. 0) write(Mout, 50)
      endif
      endfile Ibuf2
      rewind Ibuf2
      return
      end
