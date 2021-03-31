      subroutine FWSIG(unit, nsize, b, sigma)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, l, l0, li, mattst, nsize
 
c*** end of declarations inserted by spag
 
 
      integer*4 unit, nparam
      real*10 b(1), sigma(1)
c
c          r. goldstein      may, 1975
c          routine for forming and writing out sigma from b.
c          formally, this was done in-line in nrmict and nrmsav
c           negative nsize indicates no suppression of zero rows
c
      nparam = iabs(nsize)
 
      mattst = 1
 
c initialize row pointer
      l0 = 0
      do i = 1, nparam
         if( nsize .gt. 0 ) mattst = 0
         l  = l0
         li = 1
         do j = 1, nparam
            l = l + li
            sigma(j) = b(l)
            if( j .ge. i ) li = j
            if( sigma(j) .ne. 0.0_10 ) mattst = 1
         end do
         if( mattst .gt. 0 .or. i .eq. nparam ) write(unit) i,
     .      (sigma(j), j = 1, nparam)
         l0 = l0 + i
      end do
 
      return
      end
