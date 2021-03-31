      SUBROUTINE WRKARR(KBYTES)
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C*          SUBR. WRKARR - J.F.CHANDLER - 1981 FEB
C*          CALLED FROM ABC MAIN PROGRAM WITH REQUEST IN KILOBYTES
C*          OBTAIN FREE STORAGE IN THAT AMOUNT (FAILURE MEANS ABEND)
C*          THEN CALL ABCOPR(SIZE,DATA)
C* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      COMMON/WRKARX/ STORAG(1024000)
      REAL*4 STORAG
      INTEGER*4 NSTOR/1024000/
C
      KBMAX=(NSTOR*4)/1024
      IF(KBYTES.GT.KBMAX) THEN
         WRITE(6,100) KBYTES,KBMAX
  100    FORMAT(' *** REQUEST FOR',I12,' KBYTES: MAX=',I5)
         STOP 10
      ENDIF   
      CALL ABCOPR(NSTOR,STORAG)
      RETURN
      END
