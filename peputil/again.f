      SUBROUTINE AGAIN
      IMPLICIT NONE
C         REINITIALIZE PLOTER STUFF
C
C         COMMON
      include 'inodtabc.inc'
      include 'graphs.inc'
      include 'pltfrm.inc'
C
C        LOCAL
      REAL*4 XREM
      INTEGER*4 I,IOB
C
C           JUST TO BE SAFE, RESET ALL TAPES
      TAPENO = 0
      DO 10 I = 1, 10
      IOB = OBSLIB(I)
      IF(IOB.LE.0) GO TO 20
      REWIND IOB
 10   CONTINUE
 20   CONTINUE
C
C         RESET PLOT PARAMS
      IF(IPLOT.EQ.0) RETURN
      XLNOM=0.
      IF(YEND+YLEN*FACT.GT.YFRAME) TFIRST=.TRUE.
      XREM=XFRAME-XEND
      IF(XREM.LT.XEND.AND.XREM.LT.XLEN*FACT) XEND=XFRAME
      RETURN
      END
