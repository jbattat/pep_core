      SUBROUTINE DAYTIM
      IMPLICIT NONE
C       PRINT OUT DATE AND TIME USING ROUTINE 'TODAY'
C       FORMAT IS YY/MM/DD AND HH:MM:SS
C
C         COMMON
      include 'inodtabc.inc'
      include 'plttyp.inc'
C
C         LOCAL VARIABLES
      CHARACTER*8 VERSN/'92/09/22'/
      CHARACTER*4 DBF(5),D(2),T(2)
      EQUIVALENCE (DBF(1),D(1)),(DBF(4),T(1))
C
C         START
      CALL TODAY(DBF)
      WRITE(PRINT,10) VERSN,LOCALE,PLTRNM
   10 FORMAT('1***** ABC  VERSION ',A8,' - DEVICE',I2,1X,A8,' *****'//)
      WRITE (PRINT, 1) D, T
    1 FORMAT ('-*** DATE ',2A4,' *** TIME ',2A4,' ***'/)
      RETURN
      END
