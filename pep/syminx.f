      subroutine SYMINX(a, b, n, m, pvrow, pvcol, iuse, pvrwb, ierinv)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 a, b, diag, dmax, pivot, pvcol, pvrow, pvrwb
      integer   i, i21, ierinv, ii, ij, ik, il, iter, j, j21, k, l, l21,
     .          m, n
 
c*** end of declarations inserted by spag
 
 
c
c n.brenner    june 1968    subroutine syminx
c m ash/f amuchastegui  -  december 1968 - extended precision
c
      integer*2 iuse
c     inversion of a symmetric matrix in place.  algorithm 150, comm.
c     acm, h. rutishauser.
c     a is half of a symmetric matrix, n by n, arranged shortest rows
c     first.  that is, the first row is 1 long, the next 2, ... the
c     last n.  element a(i,j) (i.ge.j) is at a(j+i*(i-1)/2).
c     b is constant array of dimension n by m.  the solution of the
c     matrix equation ax=b is left in b.
c     pvrow, pvcol and iuse are temporary arrays of length n.  pvrwb
c     is a temporary array of length m.
      real*10 stuff(2)
      dimension a(1), b(1), pvrow(1), pvcol(1), iuse(1), pvrwb(1)
c dimension a(n*(n+1)/2),b(n,m),pvrow(n),pvcol(n),iuse(n),pvrwb(m)
c initially, all rows are unused.
      ierinv = 0
      do i = 1, n
         iuse(i) = 0
         end do
      do iter = 1, n
 
c search for the largest diagonal element, for use as a pivot.
         dmax = 0.
         ii   = 1
         do i = 1, n
            diag = ABS(a(ii))
            if(diag.gt.dmax .and. iuse(i).le.0) then
               dmax     = diag
               k        = i
               pivot = a(ii)
            endif
 
            ii = ii + i + 1
            end do
         if(dmax.le.0) then
            write(6, 20) iter, n, (iuse(i), i = 1, n)
   20       format(
     .'0FAILURE TO FIND NON-ZERO DIAGONAL ELEMENT IN SYMINV ON ITERATION
     .', i3, ' FOR N = ', i3, '.  IUSE = '/(1x,50I2))
            ierinv = 1
            return
         else
            iuse(k) = 1
c get the elements of the pivot row and the pivot column, with the
c proper signs.
            ik = 1 + (k*(k-1))/2
            do i = 1, n
               i21 = i*2 - 1
               if(i.lt.k) then
 
c 60 pvrow(i)=a(ik)
                  call XLOAD8(a(ik))
                  call XSTORE(pvrow(i21))
                  call XDIV8(pivot)
                  call XSTORE(pvcol(i21))
 
c pvcol(i)=a(ik)/pivot
                  a(ik) = 0.
                  ik    = ik + 1
                  if(iuse(i).gt.0) then
 
c 70 pvcol(i)=-pvcol(i)
                     call XLOAD8(0._10)
                     call XSUB(pvcol(i21))
                     call XSTORE(pvcol(i21))
c for j.gt.i, a(i,j) = (-1)**(iuse(i)+iuse(j)) * a(j,i).  the iuse
c array is used as it is before this pivot step, so that iuse(k)=0.
                  endif
                  go to 40
               else if(i.eq.k) then
 
c 80 pvrow(i)=1.
                  call XLOAD8(1._10)
                  call XSTORE(pvrow(i21))
 
c pvcol(i)=-1./pivot
                  call XLOAD8(-1._10)
                  call XDIV8(pivot)
                  call XSTORE(pvcol(i21))
                  a(ik) = 0.
                  il    = i
                  if(m.gt.0) then
                     do l = 1, m
                        l21 = l*2 - 1
 
c pvrwb(l)=b(il)
                        call XLOAD8(b(il))
                        call XSTORE(pvrwb(l21))
                        b(il) = 0.
                        il    = il + n
                        end do
                  endif
               else
 
c 90 pvrow(i)=a(ik)
c pvcol(i)=a(ik)/pivot
                  call XLOAD8(a(ik))
                  call XSTORE(pvrow(i21))
                  call XDIV8(pivot)
                  call XSTORE(pvcol(i21))
                  a(ik) = 0.
                  if(iuse(i).gt.0) then
 
c 100 pvrow(i)=-pvrow(i)
                     call XLOAD8(0._10)
                     call XSUB(pvrow(i21))
                     call XSTORE(pvrow(i21))
                  endif
               endif
               ik = ik + i
   40          end do
 
c perform the pivot step.
            ij = 1
            do i = 1, n
               i21 = i*2 - 1
               il  = i
               if(m.gt.0) then
                  do l = 1, m
                     l21     = l*2 - 1
                     call XLOAD(pvcol(i21))
                     call XMUL(pvrwb(l21))
                     call XSTORE(stuff)
                     call XLOAD8(b(il))
                     call XSUB(stuff)
                     call STORND(b(il))
 
c b(il)=b(il)-pvcol(i)*pvrwb(l)
                     il = il + n
                     end do
               endif
               do j = 1, i
                  j21     = j*2 - 1
                  call XLOAD(pvcol(i21))
                  call XMUL(pvrow(j21))
                  call XSTORE(stuff)
                  call XLOAD8(a(ij))
                  call XSUB(stuff)
                  call STORND(a(ij))
 
c a(ij)=a(ij)-pvcol(i)*pvrow(j)
                  ij = ij + 1
                  end do
               end do
         endif
         end do
      return
      end
