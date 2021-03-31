      subroutine DIFPRT(name1, name2, ip, prmi8, prms8, dift, ityp)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 prmi, prms
      integer   ip, ityp, ndtot
 
c*** end of declarations inserted by spag
 
 
c           subr. difprt - j.f.chandler - 1980 june
c           print disparities between input and saved nominals
c           (alt. entry difpts for setup)
c           (alt. entry difpt4 for real*4 values)
c
c           parameters
      character*8 name1
      character*6 name2
      real*10 dift, prmi8, prms8
      real*4    prmi4, prms4
 
      include 'fcntrl.inc'
      include 'inodta.inc'
 
      character*8 type(2)/'FREE','FIXED'/
 
      prmi = prmi8
      prms = prms8
  100 if(Line.gt.58 .or. ndtot.le.1) then
         if(Line.gt.40) call NEWPG
         write(Iout, 150)
  150    format('0', 10x, 'NOMINAL PARAMETER DISPARITIES'/
     .          5x, 'PARAMETER NAME', 7x, 'INPUT', 13x, 'SAVED', 11x,
     .          'DIFFERENCE', 7x, 'TYPE')
         Line = Line + 3
      endif
      write(Iout, 200) name1, name2, ip, prmi, prms, dift, type(ityp)
  200 format(2x, a8, 1x, a6, i3, 1p, 3D18.8, 3x, a8)
      Line  = Line + 1
      ndtot = 2
      return
 
      entry DIFPTS
c
c set up for print
      ndtot = 1
      return
 
      entry DIFPT4(name1, name2, ip, prmi4, prms4, dift, ityp)
c
c copy values into r*8 locations
      prmi = prmi4
      prms = prms4
      go to 100
      end
