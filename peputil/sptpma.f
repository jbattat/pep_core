      SUBROUTINE SPTPMA(NPLNT, K, NAMES, LVEC, SNAME)
C         GENERATES NAMES FOR SPOT COORDINATES
C
C         PARAMETERS
      CHARACTER*8 NAMES(2, 1)
      INTEGER*4 K, NPLNT
      CHARACTER*4 SNAME
      INTEGER*2 LVEC(3)
C
C         LOCALS
      CHARACTER*8 COORD(3) /'RAD', 'LONG', 'LAT'/
      CHARACTER*4 BLANKS /'    '/
C
C         LOOP THRU COORDINATES
      DO 20 I = 1, 3
      IF (LVEC(I) .EQ. 0) GO TO 20
      K = K + 1
      CALL MVC (SNAME, 1, 4, NAMES(1, K), 1)
      CALL MVC (BLANKS, 1, 4, NAMES(1, K), 5)
      NAMES(2, K) = COORD(I)
 20   CONTINUE
      RETURN
      END
