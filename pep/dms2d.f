      real*10 function DMS2D(deg, imin, sec)
 
      implicit none
 
c       J.F. Chandler   1977 Feb
c
c        Expanded 1977 Oct.
c        Seconds option added 1977 Nov.
c        Converted to Fortran 77 1992 Jun.
c
c  Purpose:
c        Convert astronomical notation (sign of degrees applies to all)
c
c    DEG:  CHARACTER*3 array of length 3 containing maybe + or -
c           (e.g., read from card as 3A1)
c           Note: invalid numeric characters are treated as zeroes.
c
c    IMIN: INTEGER*4 minutes
c    SEC:  REAL*8 seconds
c
c        sign (if any) in DEG is applied to entire result
c        e.g.    DEC = DMS2D('-01',5,13.0_10)
c                DEC = -1.08694444444444_10
c
c        In addition, a minus sign on any field  applies to all
c        e.g.    DMS2D('  0',-7,3._10) = DMS2D(' -0',7,-3._10)
c              = DMS2D('   ',-7,-3._10)
c
c        FUNCTION DMS2S is exactly the same except that the result is
c        in seconds rather than degrees.
c
c        Also, both functions work for hours instead of degrees.
c
c   output:      dms2d  = decimal degrees
c                         (double precision)
c
c
c-----------------------------------------------------------------------
 
      real*10 sec, DMS2S
      integer*4 imin, ideg, i, ibyte
      character*3 deg
      character*1 byte
      logical negate, secout
 
      secout = .false.
      goto 100
 
      entry DMS2S(deg, imin, sec)
      secout = .true.
 
  100 negate = imin.lt.0 .or. sec.lt.0._10
 
c convert character*3 degrees to an integer
 
      ideg=0
 
      do i=1,3
        byte = deg(i:i)
        if(byte.eq.'-') then
          negate= .true.
        else
          ibyte= ICHAR(byte)-ICHAR('0')
          if(ibyte.lt.0 .or. ibyte.gt.9) ibyte=0
          ideg=ideg*10 + ibyte
          end if
        end do
 
c convert degrees and minutes to seconds
 
      DMS2D = ideg*3600 + IABS(imin)*60
 
c add in seconds and choose a sign
 
      DMS2D = DMS2D + ABS(sec)
      if(negate) DMS2D = -DMS2D
 
c convert seconds to fraction if necessary
 
      if(.not.secout) DMS2D = DMS2D/3600._10
      return
      end
