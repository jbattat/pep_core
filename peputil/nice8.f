      SUBROUTINE NICE8(TOP,BOT,RNG,NDIV,DIV)
 
      IMPLICIT NONE
      REAL*8 TOP,BOT,RNG,DIV
      INTEGER*4 NDIV
 
C           R. REASENBERG, 72/8
C           CALLING SEQUENCE CHANGED, J.F.CHANDLER, 80/10
C
C  TOP  IN: TOP OF NOMINAL RANGE, OUT: TOP OF NICE RANGE
C  BOT  DITTO FOR BOTTOM OF RANGE
C  RNG  OUT: NICE TOP-BOT
C  NDIV IN: EXPECTED NUMBER OF STEPS IN RANGE IS IABS(NDIV)
C  DIV OUT: NICE NUMBER OF STEPS TO USE, NORMALLY AN INTEGER.  IF NDIV
C           IS NEGATIVE ON INPUT, DON'T CHANGE TOP OR BOT AND GIVE A
C           POSSIBLY NON-INTEGRAL VALUE FOR DIV
C
      REAL*4 LOG2/0.301030/, LOG2P5/0.397940/, LOG5/0.698970/
C
      REAL*8 AMS,DS,DSL,DSLQ,DSLR,DSM,DSS,DS0,RNG0,RNG1,RNG2
      INTEGER*4 I,LDS,MS,NBOT,NBOTS,NTOP,NTOPS
 
C           MS IS THE 'MAXIMUM' NUMBER OF SCALE POINTS
      MS=10
      IF(IABS(NDIV).GT.1) MS=IABS(NDIV)
      AMS=MS
C
      RNG=1D58
      DSS=0.
      RNG0=TOP-BOT
      DS0=RNG0/AMS
      RNG1=RNG0
      IF(RNG1.LE.0.) RNG1=0.01
C        LOOK FOR SPECIAL CASES
      IF(BOT.GE.0. .AND. BOT.LT.6. .AND. TOP.GT.174. .AND. TOP.LE.180.)
     . THEN
         RNG=180.
         NBOTS=0
         NTOPS=NDIV
         IF(NDIV.LE.1 .OR. NDIV.EQ.5 .OR. NDIV.EQ.7 .OR. NDIV.EQ.8) 
     .    NTOPS=6
         IF(NDIV.GT.9) NTOPS=9
         DSS=RNG/NTOPS
         GOTO 70
      ELSE IF(BOT.GE.0. .AND. BOT.LT.6. .AND. TOP.GT.354. .AND.
     .    TOP.LE.360.) THEN
         RNG=360.
         NBOTS=0
         NTOPS=NDIV
         IF(NDIV.LE.1 .OR. NDIV.EQ.5) NTOPS=4
         IF(NDIV.EQ.7) NTOPS=8
         IF(NDIV.EQ.10 .OR. NDIV.EQ.11) NTOPS=9
         IF(NDIV.GT.12) NTOPS=12
         DSS=RNG/NTOPS
         GOTO 70
      ELSE IF(BOT.GE.0. .AND. BOT.LT.0.4 .AND. TOP.GT.23.6 .AND.
     .    TOP.LE.24.) THEN
         RNG=24.
         NBOTS=0
         NTOPS=NDIV
         IF(NDIV.LE.1 .OR. NDIV.EQ.5) NTOPS=4
         IF(NDIV.EQ.7 .OR. NDIV.EQ.9) NTOPS=8
         IF(NDIV.GT.9) NTOPS=12
         DSS=RNG/NTOPS
         GOTO 70
      ENDIF
C
C           REPEAT CALCULATION TO LOOK FOR BEST FIT
      DO 60 I=1,3
      DSM=RNG1/AMS
C
C           DSM IS THE SMALLEST POSSIBLE DS
      DSL=DLOG10(DSM)+100.0
      LDS=DSL-0.00001
      DSLQ=DSL-LDS
C
      IF(NDIV.LT.0) THEN
        DSLR=DINT(1D1**DSLQ + 0.99999)
      ELSE
        DSLR=2.0
        IF(DSLQ.GT.LOG2) DSLR=2.5
        IF(DSLQ.GT.LOG2.AND.LDS.EQ.100) DSLR=5.0
        IF(DSLQ.GT.LOG2P5) DSLR=5.0
        IF(DSLQ.GT.LOG5) DSLR=10.0
      ENDIF
      DS=DSLR*1D1**(LDS-100)
      IF(DS.EQ.DSS) GOTO 60
C
      NTOP=(TOP*0.99999+BOT*0.00001)/DS
      NBOT=(BOT*0.99999+TOP*0.00001)/DS
      IF(TOP.GT.0.) NTOP=NTOP+1
      IF(BOT.LT.0.) NBOT=NBOT-1
      RNG2=(NTOP-NBOT)*DS
C           LOOK FOR CLOSEST FIT OF RANGE
      IF(RNG2.GE.RNG) GOTO 60
C           SAVE BEST SO FAR
      NTOPS=NTOP
      NBOTS=NBOT
      DSS=DS
      RNG=RNG2
      IF(DS0.LE.0.) GOTO 70
   60 RNG1=RNG1+DS0
C
   70 IF(NTOPS.EQ.NBOTS) NTOPS=NBOTS+1
      DIV=NTOPS-NBOTS
      IF(NDIV.GE.0) THEN
        TOP=NTOPS*DSS
        BOT=NBOTS*DSS
      ELSE
         RNG=RNG0
         DIV=RNG0/DSS
      ENDIF
      RETURN
      END
