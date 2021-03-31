      subroutine NRMRIT(title, b, nsize, measmt, switch)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, l1, l2, ldj, lininc, measmt, nsize, nt, nzero
 
c*** end of declarations inserted by spag
 
 
      character*8 title(4)
      real*10 b(1)
      logical   switch
c
c           print out lower diagonal matrix
c title  - 32-byte name for the matrix
c b      - matrix (in lower diagonal form)
c nsize  - dimension of matrix
c measmt - number of observations included in the matrix
c switch - logical flag: true if matrix to be printed
c
      include 'inodta.inc'
 
      character*8 stitl(13)/'THE LOWE', 'R DIAGON', 'AL HALF', 4*' ',
     .          ' MATRIX ', 'OF ORDER', '     FOR', 'MED FROM', ' ',
     .          ' OBS'/
      character*1 lchar, skip/'0'/
      character*4 row,blank/'    '/
      equivalence(lchar, row)
c
c* start=500
      if( switch ) then
         do i = 1, 4
            stitl(i + 3) = title(i)
         end do
         call EBCDI(nsize, stitl(10), 4)
         call EBCDI(measmt, stitl(12), 8)
         nt = -19
         if( measmt .gt. 0 ) nt = -25
         call PAGSET(stitl, nt)
         if( nsize .gt. 8 ) call PAGE(0, 1)
         if( nsize .le. 8 ) call PAGHED(0)
c
c*  start=1000
c write out coefficient matrix of normal equations
         do i = 1, nsize
            ldj = (i*(i-1))/2
            call EBCDI(i, row, 4)
            lchar  = skip
            lininc = 2
            l2     = 0
   20       nzero  = 0
            do while( .true. )
               l1 = l2 + 1
               l2 = l1 + 7
               if( l2 .ge. i ) then
                  l2 = i
                  go to 40
               else
                  do j = l1, l2
                     if( b(ldj+j) .ne. 0.0_10 ) go to 40
                  end do
                  nzero = nzero + 1
               endif
            end do
   40       if( nzero .lt. 1 ) then
            else if( nzero .eq. 1 ) then
               l2 = l1 - 1
               l1 = l2 - 7
            else
               call PAGCHK(60, lininc, 1)
               write(Iout, 50) row, nzero
   50          format(a4, i5, ' ZERO LINES')
               row    = blank
               lininc = 1
            endif
            call PAGCHK(60, lininc, 1)
            write(Iout, 60) row, (b(ldj+j), j = l1, l2)
   60       format(a4, 1p, 8D16.8)
            row    = blank
            lininc = 1
            if( l2 .lt. i ) go to 20
         end do
      else
         call PAGCHK(60, 2, 0)
         write(Iout, 100) title
  100    format('0', 4A8, ' MATRIX NOT PRINTED OUT')
      endif
c
c
c* start=9000
      return
      end
