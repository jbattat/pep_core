      subroutine NRMSCL(b, scale, nsize, side, msize, print)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, i1, j, ldj, msize, ngzr, nsize
 
c*** end of declarations inserted by spag
 
 
c subroutine nrmscl - j.f.chandler - 1980 jan
c derived from m.e.ash   nov 1969    subroutine nrmsav
c scale the coefficient matrix and right side of the normal eqs
c option with nsize .ne. nparam used only for call from pprctl
      real*10 b(1), scale(1), side(1)
      logical*4 print
c       b - normal equation matrix (of size nsize)
c       scale(i)=1/sqrt(diagonal element i)   (output)
c       nsize - number of normal equations
c       side - vector of right hand sides
c       msize - if 0 then no rhs
c       print - if true then print scale factors
c
      include 'inodta.inc'
 
      integer*2 ierr(2)
      character*8 zerneg(2)/'ZERO','NEGATIVE'/
c
c calculate scale factors for the normal equations
c*  start=1000
      ldj     = 0
      ierr(1) = 0
      ierr(2) = 0
      do i = 1, nsize
         ldj = ldj + i
         if(b(ldj).gt.0._10) then
            scale(i) = 1._10/SQRT(b(ldj))
         else
            if(b(ldj).lt.0._10) then
               ngzr = 2
            else
               ngzr     = 1
               scale(i) = 1._10
            endif
            ierr(ngzr) = ierr(ngzr) + 1
            call PAGCHK(60, 1, 0)
            write(Iout, 20) i, zerneg(ngzr)
   20       format(' DIAGONAL ELEMENT', i4, ' IS ', a8)
         endif
      end do
c
c*  start=2000
c printout scale factors
      if(ierr(1).gt.0) then
         call PAGCHK(60, 2, 0)
         write(Iout, 50) ierr(1)
   50    format(/i4, ' DIAGONAL ELEMENTS OF COEFFICIENT MATRIX ARE ZERO,
     . CORRESPONDING SCALE FACTORS ARE ASSUMED 1.0_10')
      endif
      if(ierr(2).le.0) then
         if(print) then
            call PAGCHK(60, 4 + (nsize+7)/8, 0)
            write(Iout, 60) nsize, (scale(i), i = 1, nsize)
   60       format('0THE', i4, ' SCALE FACTORS APPLIED TO THE COEFFICIEN
     .T MATRIX AND RIGHT SIDES OF THE NORMAL EQUATIONS ARE'//
     . (4X,1P,8D16.8))
         endif
c
c*  start=3500
c scale the coefficient matrix and right sides of the
c normal equations
         b(1) = 1._10
         if(nsize.gt.1) then
            ldj = 1
            do i = 2, nsize
               b(ldj + i) = 1._10
               i1 = i - 1
               do j = 1, i1
                  b(ldj + j) = b(ldj + j)*scale(j)*scale(i)
               end do
               ldj = ldj + i
            end do
         endif
         if(msize.gt.0) then
            do i = 1, nsize
               side(i) = side(i)*scale(i)
            end do
         endif
      else
         call PAGCHK(60, 2, 0)
         write(Iout, 100) ierr(2)
  100    format(/i4,
     .' DIAGONAL ELEMENTS OF COEFFICIENT MATRIX ARE NEGATIVE, CANNOT SCA
     .LE NORMAL EQUATIONS')
         do i = 1, nsize
            scale(i) = 1._10
         end do
      endif
c
c*  start=9900
      return
      end
