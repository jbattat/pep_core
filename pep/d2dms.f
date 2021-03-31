      subroutine D2DMS(deg,imin,sec,ddeg,idigit)
 
      implicit none
 
c       J.F. Chandler   1977 Feb (originally part of DMS2D)
c
c        Expanded 1977 Oct.
c        Seconds option added 1977 Nov.
c        Converted to Fortran 77 1992 Jun.
c
c  Purpose:
c     Convert to astronomical notation (sign of degrees applies to all)
c
c Arguments
c  DEG:  CHARACTER*3 array of length 3 containing maybe formatted
c        number.  The minus sign (if any) is placed in the 1st col.
c  IMIN: INTEGER*4 minutes output
c  SEC:  REAL*8 seconds output
c  DDEG: REAL*8 degrees, to be converted to dms format as used
c        in DMS2D
c  IDIGIT: number of significant digits to the right of the decimal
c         point in seconds.  D2DMS will adjust the output so that
c         the 'seconds' are always less than 59.0 - 0.5*10**-IDIGIT.
c         Must be restricted to (-2 to +4).  Thus, if the seconds are
c         printed using format Fm.IDIGIT, they will come out strictly
c         less than 60.
c
c    Alternate entry point S2DMS is the same, except that the 4th
c    argument is
c  DSEC is REAL*8 seconds to be converted to dms format
c
c        Note: both entry points work for hours instead of degrees.
c
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 
      real*10 sec, dsec, ddeg, degin, valabs, round, roundd
      integer*4 imin, ideg, idigit
      character*3 deg
 
      degin = ddeg
      goto 100
 
      entry S2DMS(deg,imin,sec,dsec,idigit)
      degin = dsec/3600._10
 
  100 round = 10._10**(-idigit) / 2._10
      if(round.gt.50._10) round=30._10
      roundd = round/3600._10
 
c convert to character*3 degrees
 
      valabs= ABS(degin) + roundd
      ideg=valabs
      write(deg,110) ideg
  110 format(I3)
      if(degin.lt.0._10) deg(1:1)='-'
 
c get minutes
 
      valabs = (valabs - ideg) * 60._10
      imin = valabs
 
c get seconds and deduct rounding criterion
 
      sec = (valabs - imin) * 60._10 - round
      if(sec.lt.0._10) sec=0._10
 
      return
      end
