      subroutine CHKHMS(jds, ihr, imin, sec)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   jds
      real      sec
 
c*** end of declarations inserted by spag
 
 
c
c ash/forni  march 1968  subroutine chkhms
c
c to check and change ihr,imin,sec less than 24,60,60
c respectively
      integer*2 ihr, imin
      if( sec .ge. 60.0 ) then
         sec  = sec - 60.0
         imin = imin + 1
      endif
      if( imin .ge. 60 ) then
         imin = imin - 60
         ihr  = ihr + 1
      endif
      if( ihr .ge. 24 ) then
         ihr = ihr - 24
         jds = jds + 1
      endif
      return
      end
