      SUBROUTINE TRANS(DATA)
      IMPLICIT NONE
C         CALLS OUT APPROPRIATE TRANSFORM ROUTINE FOR EACH SPAN
C         USING TIME AS WORK ARRAY FOR FOURG
C
C        PARAMETER
      REAL*4 DATA(1)
C           DATA  - STORAGE ARRAY CONTAINING DATA
C         COMMON
      include 'misc.inc'
      include 'span.inc'

C LOCAL
      INTEGER*4 I,ID,IT
C
C         STATEMENT FUNCTIONS FOR INDICES TO MASTER DATA ARRAY
      include 'indices.inc'
C
C         GET A POINTER FOR WORK
      IT=KD(NSPAN+1,1)
C
C         LOOP THRU SPANS
      DO I = 1, NSPAN
         ID = KD(I, 1)
         IF(USEF1) THEN
            CALL FOUR1(DATA(ID),LOG2N,-1)
         ELSE
            CALL FOURG(DATA(ID),NPOINT,-1,DATA(IT))
         ENDIF
      END DO
      RETURN
      END
