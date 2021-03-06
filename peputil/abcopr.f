      SUBROUTINE ABCOPR(MAXPT,DATA)
      IMPLICIT NONE
C       SUBR. ABCOPR - J.F.CHANDLER - 1981 FEB
C        CODE TAKEN FROM ABC MAIN ROUTINE
C
C        PERFORM REQUESTED ABC FUNCTIONS BY CALLING THE WORKING
C        SUBROUTINES
C
C        PARAMETERS
      INTEGER*4 MAXPT
      REAL*4 DATA(MAXPT)
C  MAXPT - SIZE OF DATA ARRAY (IN FULLWORDS)
C  DATA  - WORK AREA FOR PASSING OBSLIB DATA AMONG SUBROUTINES
C
C        COMMON
      include 'graphs.inc'
      include 'misc.inc'
C
C        LOCAL
      INTEGER*4 TTMODE
C
      TTMODE=TMODE
C
C         TAKE LOOK AT TAPE
      TMODE=TTMODE
      IF(LOOK.GE.1) CALL LOOKER(MAXPT,DATA)
      IF(TASK.LE.0) GOTO 4000
C
C         CHECK PARAMETERS
      CALL PARMO(MAXPT,DATA)
      TMODE=TTMODE
      IF(TMODE.EQ.3) TMODE=1

      IF(FFT.OR.TASK.LT.3) THEN
C
C SELECT DATA
         CALL SELECT(DATA,DATA)
C DISPLAY SELECTED DATA
         CALL DISPLY(DATA,DATA)
C TRANSFORM
         IF(TASK.GE.3) CALL TRANS(DATA)
      ELSE
C SELECT, DISPLAY, AND TRANSFORM IN ONE STEP
         CALL FGEN(DATA,DATA)
      ENDIF
C
      IF(TASK.LT.3) GOTO 4000
C
C         LIST AND GRAPH POWER SPECTRA
      CALL SHOW(DATA)
C
C         AUTOCORRELATION
      IF(TASK.LT.4) GOTO 4000
      CALL AUTOC(DATA)
 4000 RETURN
      END
