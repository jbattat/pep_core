      subroutine ZFILL(array, length)
c
c  j.f.chandler - 1983 feb
c
c  subroutine for zeroing core
c
c  array  - name of array or array element
c  length - size in bytes of region to zero
c
c
      implicit none
 
      integer*4 length, i, izero/0/
      character*1 array(length), zero
 
      equivalence (izero,zero)
 
      do i = 1, length
         array(i) = zero
         enddo
 
      end
