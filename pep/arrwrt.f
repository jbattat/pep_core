      subroutine ARRWRT(inum, jnum, ntype)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, i1, ij, inum, j, jinc, jnum, ksav, l
 
c*** end of declarations inserted by spag
 
 
c        subroutine to write subarrays of b on temporary storage
c        for subroutine pprctl
c        p. macneil  may, 1977
c
      character*4 ntype
c
c        commons
      include 'aprtbf.inc'
      include 'nrmmat.inc'
      include 'restor.inc'
      real*10   qb(1001), qsav(500)
      equivalence(qb, B), (qsav, Sav)
 
      character*4 nc/'C   '/,ND/'D   '/,NF/'F   '/
      character*4 dq/'DQ  '/
c
c determine array type
      if( ntype .eq. nc ) then
c
c c - matrix
         ksav = Ibuf6
      else if( ntype .eq. nd ) then
c
c d - inverse matrix
         ksav = Ibuf8
      else if( ntype .eq. nf ) then
c
c f - matrix
         ksav = Ibuf7
c
c write out rectangular sub - array
         do i = 1, inum
            do j = 1, jnum
               l = (i - 1)*jnum + j
               Sav(j) = B(l)
            end do
 
c write out row of array
            write(ksav) (Sav(j), j = 1, jnum)
 
         end do
         go to 100
      else if( ntype .eq. dq ) then
c
c extended-precision d-inverse
         ksav = Ibuf8
 
         i1 = 1
         do i = 1, inum
            ij   = i1
            jinc = 1
            do j = 1, inum
               qsav(j) = qb(ij)
               if( j .ge. i ) jinc = j
               ij = ij + jinc
            end do
 
c write out row of array
            write(Ibuf8) (qsav(j), j = 1, inum)
            i1 = i1 + i
         end do
         go to 100
      else
         call SUICID('BAD NTYPE, STOP IN ARRWRT   ', 7)
      endif
c
c write out symmetric array held in lower diagonal form
      do i = 1, inum
         do j = 1, jnum
            if( j .gt. i ) then
 
c element is above diagonal, reverse i and j
               l = i + (j*(j-1))/2
            else
               l = j + (i*(i-1))/2
            endif
 
c load sav
            Sav(j) = B(l)
 
         end do
 
c write out row of array
         write(ksav) (Sav(j), j = 1, jnum)
 
      end do
c
c
  100 endfile ksav
      rewind ksav
 
      return
      end
