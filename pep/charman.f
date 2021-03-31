c-----------------------------------------------------------------------
c
c function iscan(string,nlength,'char')  scans for char
c result is integer offset in string (0 up) if 'char' is there,
c or -999 if not there
c
c-----------------------------------------------------------------------
      integer function ISCAN(string, length, char)
 
      implicit none
 
      integer   i, length
      character*1 string(length)
      character*(*) char
 
      do i = 1, length
         if(string(i).eq.char(1:1)) then
            ISCAN = i - 1
            return
         endif
      enddo
 
      ISCAN = -999
 
      return
      end
 
 
c-----------------------------------------------------------------------
c
c function nscan is just like iscan, but scans for non-match
c
c-----------------------------------------------------------------------
      integer function NSCAN(string, length, char)
 
      implicit none
 
      integer   i, length
      character*1 string(length)
      character*(*) char
 
      do i = 1, length
         if(string(i).ne.char(1:1)) then
            NSCAN = i - 1
            return
         endif
      enddo
 
      NSCAN = -999
 
      return
      end
 
 
c-----------------------------------------------------------------------
c
c subroutine movebl(from, lfrom, to, lto)  move 'lfrom' characters
c from 'from' to 'to', then blank out the remaining part, if
c any, up to total length 'lto'
c
c warning: if lfrom > lto, some data may be clobbered
c
c-----------------------------------------------------------------------
      subroutine MOVEBL(from, lfrom, to, lto)
 
      implicit none
 
      character*(*) from, to
      integer   i, lfrom, lto
 
      to(1:lfrom) = from(1:lfrom)
 
      do i = lfrom + 1, lto
         to(i:i) = ' '
         enddo
 
      return
      end
