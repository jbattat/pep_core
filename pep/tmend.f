      subroutine TMEND
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, itotrl
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   sept 1967    subroutine tmend
c print timer information for end of program run
c
 
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'timstf.inc'
 
      character*4 eate(2), jdayr, czero
      real*4    sec(5), zero4/0./
      integer*4 ihr(5), imin(5), itot(2)
      equivalence (czero,zero4)
      real*10 f1, f2
c
c make sure all tapes are rewound
      do i = 1, Ntptop
         if(Itrwnd(i).gt.0) rewind i
         end do
c
c write out times for this just completed run
      jdayr   = Idayr
      eate(1) = czero
      call TIMSET(ihr, imin, sec, itot, eate)
      itotrl = Ireal0 - Itimst
      if(Idayr.ne.jdayr) itotrl = 8640000 + itotrl
      call CONVRT(ihr(3), imin(3), sec(3), itotrl)
      call CONVRT(ihr(4), imin(4), sec(4), Itotsk)
      call CONVRT(ihr(5), imin(5), sec(5), Ireal0)
      write(Iout, 100) eate, Idayr, ihr(5), imin(5), sec(5)
  100 format('0 TERMINATION DATE       = ',2a4/
     .       '  TERMINATION DAY OF YEAR= ',a3/
     .       '  TERMINATION CLOCK TIME =',   i3,'H',i3,'M',f6.2,'S')
      write(Iout, 200) (ihr(i), imin(i), sec(i), i = 1, 4)
  200 format(/28x, 'REAL TIMER', 6x, 'TASK TIMER'/
     .            9x, 'CLOCK INCREMENT =', 2(i3,'H',i3,'M',f6.2,'S ')/
     .       4x, 'TOTAL EXECUTION TIME =', 2(i3,'H',i3,'M',f6.2,'S '))
      f1 = itotrl
      f2 = Itotsk
      f2 = f2/f1
      write(Iout, 300) f2
  300 format(' (TASK TIME)/(REAL TIME) =', f8.5)
      write(Iout, 400)
  400 format('-'/'-'/'-'/'-'/'-'/
     .36x,'***************   ***               ***    **************  '/
     .36x,'***************   ****              ***    *************** '/
     .36x,'***************   *****             ***    ****************'/
     .36x,'***               ******            ***    ***          ***'/
     .36x,'***               *** ***           ***    ***          ***'/
     .36x,'***               ***  ***          ***    ***          ***'/
     .36x,'***               ***   ***         ***    ***          ***'/
     .36x,'***               ***    ***        ***    ***          ***'/
     .36x,'*************     ***     ***       ***    ***          ***'/
     .36x,'*************     ***      ***      ***    ***          ***'/
     .36x,'*************     ***       ***     ***    ***          ***'/
     .36x,'***               ***        ***    ***    ***          ***'/
     .36x,'***               ***         ***   ***    ***          ***'/
     .36x,'***               ***          ***  ***    ***          ***'/
     .36x,'***               ***           *** ***    ***          ***'/
     .36x,'***               ***            ******    ***          ***'/
     .36x,'***************   ***             *****    ****************'/
     .36x,'***************   ***              ****    *************** '/
     .36x,'***************   ***               ***    **************  ')
      return
      end
