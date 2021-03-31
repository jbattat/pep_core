      integer function LEG(n, m1, list1, m2, list2)
c
c This program compares "n" characters starting in the m1th position
c of "list1" with the m2th position of "list2".
c It returns an integer value of -1, 0, or +1 according to whether
c the first string is less than, equal to, or greater than the second.
c The comparison is done assuming all "n" characters exist.
 
      implicit none
 
      character*(*) list1, list2
      integer   n, m1, m2
 
      if( list1(m1:m1+n-1) .eq. list2(m2:m2+n-1) ) then
         LEG = 0
      else if( list1(m1:m1+n-1) .lt. list2(m2:m2+n-1) ) then
         LEG = -1
      else
         LEG = 1
      endif
 
      return
      end
