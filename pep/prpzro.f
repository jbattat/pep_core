      subroutine PRPZRO
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i
 
c*** end of declarations inserted by spag
 
 
c t.m. eubanks feb 1977
c subroutine to zero prpgat variables
c ncal,cal,scal,ical,sumcor
c
      include 'prpgat.inc'
 
      Sumcor(1) = 0.
      Sumcor(2) = 0.
      Ncal = 0
      do i = 1, 100
         Ical(i) = 0
         Cal(i)  = 0.0
         Scal(i) = 0.0
 
      end do
      return
      end
