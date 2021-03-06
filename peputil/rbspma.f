      SUBROUTINE RBSPMA(NPLNT, K, NAMES, LVEC, SITE, SERIES)
C         GENERATES NAMES FOR RADAR BIASES
C
C         PARAMETERS
      CHARACTER*8 NAMES(2, 1)
      CHARACTER*4 SITE(2, 2), SERIES
      INTEGER*4  K, NPLNT
      INTEGER*2  LVEC(2)
C
C         LOCALS
      CHARACTER*3 RB /'TDS'/
      CHARACTER*8 TEMP
C
C         LOOP THRU BIASES
      DO 20 I = 1, 2
      IF (LVEC(I) .EQ. 0) GO TO 20
      K = K + 1
      CALL MVC (SITE(1, 1), 1, 4, NAMES(1, K), 1)
      CALL MVC (SITE(1, 2), 1, 4, NAMES(1, K), 5)
      CALL MVC (SERIES, 1, 4, TEMP, 1)
      TEMP(5:5)='0'
      N=6
      IF(NPLNT.GT.9.OR.NPLNT.LT.0) N=5
      CALL EBCDIX(NPLNT,TEMP,N,7-N)
      CALL MVC(RB,I,2,TEMP,7)
      NAMES(2, K) = TEMP
 20   CONTINUE
      RETURN
      END
