      subroutine KP2KI(kp, numki, ki)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i
 
c*** end of declarations inserted by spag
 
 
c           subroutine kp2ki - j.f.chandler - 1982 january
c           convert old-style body integration controls to new (if
c           necessary).   kp(1-30) -> ki(1-n), kp(93-96) -> kp(1-4)
c           halt processing if already new form but too many ki's
c
      integer*2 kp(100), numki, ki(99)
 
      include 'maxkidat.inc'
 
      if( numki .gt. maxki ) call SUICID(
     .    'TOO MANY INTEGRATION PARTIAL CONTROLS, STOP IN KP2KI', 13)
 
c already new form?
      if( numki .gt. 0 ) return
 
c copy old kp(1-30) into ki
      do i = 1, 30
         ki(i) = kp(i)
         kp(i) = 0
         if( ki(i) .ne. 0 ) numki = i
      end do
 
c make sure that 'numki' reflects at least ic's + 1
      if( numki .lt. 8 ) numki = 8
 
c copy target body numbers from kp(93-96) to kp(1-4)
      do i = 1, 4
         kp(i) = kp(i + 92)
         kp(92 + i) = 0
      end do
      return
      end
