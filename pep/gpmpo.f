      subroutine GPMPO(m, g, n, insert, label, nlabel)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, insert, j, k, lpw, m, n, nc, nlabel, np, nr
 
c*** end of declarations inserted by spag
 
 
c
c d. white  november 1974  subroutine gpmpo
c
c  General Purpose Matrix Printer Outer
c
c                Note: extended precision version is XPMPO
c
c  M - its dimension
c  G - a double dimension array  G(M,N) or GX(M,N)
c  N - the second dimension (M is first)
c      N=0 => G is vector G(M)
c      N=1 => G is lower diag half of M*M symmetric matrix
c      N>1 => G(M,N)
c  INSERT - page code; =1 start a new page
c  LABEL - optional label for matrix
c  NLABEL - number of bytes in label; >1 => print the label
c
c  parameters
      real*10 g(1)
      real*10   gx(1)
      character*1 label(nlabel)
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
 
      logical*4 long
 
      long = .false.
      goto 10
 
 
      entry XPMPO(m, gx, n, insert, label, nlabel)
      long = .true.
c
c page logic
   10 if(insert.eq.1 .or. Line.gt.55) call NEWPG
c
c label
      if(nlabel.gt.1) then
         write(Iout, 50) label
   50    format('-', 132A1)
         Line = Line + 3
      endif
 
      if(n.lt.1) then
 
c linear array
         nr = 1
         nc = m
      else if(n.eq.1) then
 
c lower diagonal half of symmetric matric
         nr = m
         nc = m
      else
 
c rectangular matrix
         nr = n
         nc = m
      endif
c
c print 'nr' rows of 'nc' items
      k = 0
      do i = 1, nr
         if(n.eq.1) nc = i
         np = nc
         if(nc.gt.6) then
            do j = 1, nc
               if(long) then
                  if(gx(k+j).ne.0._10) goto 100
               else
                  if(g(k+j).ne.0._10) goto 100
               endif
            end do
 
c row is all zeroes: abbreviate
            np = 1
         endif
 
c print appropriate-sized row
  100    lpw = (np - 1)/6 + 2
         call PAGCHK(60, lpw, 0)
         if(long) then
            write(Iout, 150) (gx(k+j), j = 1, np)
         else
            write(Iout, 150) (g(k+j), j = 1, np)
         endif
  150    format('0', (t2,1p,6D22.15))
         k = k + nc
      end do
      return
      end
