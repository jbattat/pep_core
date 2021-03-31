      integer function ENDTST(a)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, LEG
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash    oct 1970    integer function endtst
c test for &end on input namelist card searching columns 2-72
c
c endtst =0 no &end detected
c endtst =1 there was an &end detected
      character*80 a
      character*4 ampend/'&END'/
c integer function leg compares characters
c (leg stands for less than, equal to, or greater than)
c
      do i = 2, 69
         if( LEG(4,i,a,1,ampend) .eq. 0 ) then
            ENDTST = 1
            return
         endif
      end do
      ENDTST = 0
 
      return
      end
