      subroutine LINCHK
 
      implicit none
 
c
c ash/forni  october 1967  subroutine linchk
c to write or not to write a blank line is the question
c
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'adjstf.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
 
      if(Il.le.1) then
         if(Line.lt.60) then
            write(Iout,20)
            if(Jout.gt.0) write(Jout,20)
   20       format(' ')
            Line = Line + 1
         else
            call PAGE(0,1)
         endif
         Il = 2
      endif
      return
      end
