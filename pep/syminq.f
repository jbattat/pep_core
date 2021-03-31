      subroutine SYMINQ(a,b,n,m,pvrow,pvcol,iuse,pvrwb,ierinv)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10   a, b, diag, dmax, pivot, pvcol, pvrow, pvrwb
      integer   i, ierinv, ii, ij, ik, il, iter, j, k, l, m, n
 
c*** end of declarations inserted by spag
 
 
      integer*2 iuse
      dimension a(1),b(1),pvrow(1),pvcol(1),iuse(1),pvrwb(1)
c dimension a(n*(n+1)/2),b(n,m),pvrow(n),pvcol(n),iuse(n),pvrwb(m)
c subroutine syminq - j.f.chandler - 1989 feb
c based on syminv, but modified to avoid extended-precision
c divisions as much as possible (only one needed per row)
      ierinv = 0
      do i = 1, n
         iuse(i) = 0
      end do
      do iter = 1, n
 
c search for the largest diagonal element, for use as a pivot.
         dmax = 0.0
         ii   = 1
         do i = 1, n
            diag = abs(a(ii))
            if(diag.gt.dmax) then
               if(iuse(i).le.0) then
                  dmax  = diag
                  k     = i
                  pivot = a(ii)
               endif
            endif
            ii = ii + i + 1
         end do
         if(dmax.le.0) then
            write(6,20) iter,n,(iuse(i),i = 1,n)
   20       format(
     .'0FAILURE TO FIND NON-ZERO DIAGONAL ELEMENT IN SYMINQ ON ITERATION
     .', i3, ' FOR N = ', i3, '.  IUSE = '/(1x,50I2))
            ierinv = 1
            write(6,40) (a(i),i = 1,6)
   40       format('0B=', 1p, 6D16.9, ' ...')
            return
         else
            iuse(k) = 1
            pivot   = 1.0/pivot
c get the elements of the pivot row and the pivot column, with the
c proper signs.
            ik = 1 + (k*(k-1))/2
            do i = 1, n
               if(i.lt.k) then
                  pvrow(i) = a(ik)
                  pvcol(i) = a(ik)*pivot
                  a(ik)    = 0.
                  ik = ik + 1
c for j.gt.i, a(i,j) = (-1)**(iuse(i)+iuse(j)) * a(j,i).  the iuse
c array is used as it is before this pivot step, so that iuse(k)=0.
                  if(iuse(i).gt.0) pvcol(i) = -pvcol(i)
                  goto 60
               else if(i.eq.k) then
                  pvrow(i) = 1.
                  pvcol(i) = -pivot
                  a(ik)    = 0.
                  il = i
                  if(m.gt.0) then
                     do l = 1, m
                        pvrwb(l) = b(il)
                        b(il)    = 0.
                        il = il + n
                     end do
                  endif
               else
                  pvrow(i) = a(ik)
                  pvcol(i) = a(ik)*pivot
                  a(ik)    = 0.
                  if(iuse(i).gt.0) pvrow(i) = -pvrow(i)
               endif
               ik = ik + i
   60       end do
 
c perform the pivot step.
            ij = 1
            do i = 1, n
               il = i
               if(m.gt.0) then
                  do l = 1, m
                     b(il) = b(il) - pvcol(i)*pvrwb(l)
                     il    = il + n
                  end do
               endif
               do j = 1, i
                  a(ij) = a(ij) - pvcol(i)*pvrow(j)
                  ij    = ij + 1
               end do
            end do
         endif
      end do
      return
      end
