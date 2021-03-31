      subroutine NMLFIX(a)

      implicit none


c*** start of declarations inserted by spag
      integer   i, ISCAN

c*** end of declarations inserted by spag

c subroutine NMLFIX - J.F.Chandler - 1992 July

c Replace any occurrences of '&' with the proper namelist character
c on the current system (defined by the NMLCHRDT include file).

      character*80 a

      include 'nmlchrdt.inc'

      if(Nmlchr.eq.'&') return

      do while (.true.)
         i = ISCAN(a,80,'&') + 1
         if(i.le.0) goto 100
         a(i:i) = Nmlchr
         end do

  100 return
      end
