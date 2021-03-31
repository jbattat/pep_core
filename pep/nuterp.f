      subroutine NUTERP(ynlb, nlb, dmoon, kiss)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   kiss
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   aug 1968    subroutine nuterp
c perform nutation,libration fourth difference interpolation
c
      real*10 ynlb(3, 2, 2), dmoon
      real*4    nlb(2)
 
      include 'prtpin.inc'
 
      nlb(1) = P(1)
     .         *(ynlb(1,1,Mtab1) + P(2)*(ynlb(2,1,Mtab1)+P(2)*ynlb(3,1,
     .         Mtab1))) + P(3)
     .         *(ynlb(1,1,Ntab1) + P(4)*(ynlb(2,1,Ntab1)+P(4)
     .         *ynlb(3,1,Ntab1)))
      if( kiss .gt. 1 ) nlb(2) = (ynlb(1,2,Mtab1) + P(2)*(ynlb(2,2,
     .                           Mtab1)+P(2)*ynlb(3,2,Mtab1))
     .                           - ynlb(1,2,Ntab1) - P(4)
     .                           *(ynlb(2,2,Ntab1)+P(4)*ynlb(3,2,Ntab1)
     .                           ))/dmoon
 
      return
      end
