      CHARACTER*8 FUNCTION BNAME(NPLNT)
      IMPLICIT NONE
C         GENERATE NAME FOR A BODY
C
C         PARAMETER
      INTEGER*4 NPLNT
C
C         LOCAL
      CHARACTER*8  BODIES(10) /'MERCURY', 'VENUS', 'EARTH', 'MARS',
     1  'JUPITER', 'SATURN', 'URANUS', 'NEPTUNE', 'PLUTO', 'MOON'/
      CHARACTER*8 TEMP,BODY/'BODY    '/
      CHARACTER*1 R/'R'/
      INTEGER*4 N
C
C         IF NPLNT > 10, USE TWO DIGITS
      N = IABS(NPLNT)
      IF (N .LE. 10) GO TO 10
      TEMP=BODY
      CALL EBCDIX(N,TEMP,5,2)
      IF(NPLNT.LT.0) TEMP(7:7)=R
      GOTO 20
C
C         STANDARD NAME
 10   CONTINUE
      TEMP = BODIES(N)
      IF(NPLNT.LT.0) TEMP(8:8) = R
   20 BNAME = TEMP
      RETURN
      END
