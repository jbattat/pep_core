      SUBROUTINE SQALE(PLMODE,A,N,K,LEN,EXTRM,SXTRM,INICE)
 
      IMPLICIT NONE
 
C         DETERMINE RANGE AND SCALING FOR ARRAY OF DATA
C             PLMODE - TYPE OF PLOTTING TO DO:
C                  -1  NON-LOOKER PLOTS
C                   0  IMMEDIATE LOOKER PLOT
C                   1  FIRST SERIES OF ONEAXE PLOT
C                   2  SUBSEQUENT SERIES FOR ONEAXE PLOT
C                   3  FINISH OFF ONEAXE PLOT
C             A - ARRAY (A8 - SAME, BUT REAL*8)
C             N - NUMBER OF POINTS
C             K - REPEAT FACTOR
C             LEN - LENGTH OF DESIRED AXIS IN INCHES
C             EXTRM - INPUT LIMITS FOR COORDINATES + RANGE MODIFIER
C                     (TO BE ADDED AT EACH END OF ACTUAL RANGE IF
C                      THE INPUT LIMITS ARE OMITTED)
C                     IF THE MODIFIER IS -MAX, THEN NEITHER CHANGE THE
C                     RANGE NOR ALLOW FOR CHARACTER PLOT WIDTHS.
C             SXTRM - LIMITS ADOPTED FOR PLOT
C             INICE - IF POSITIVE, CALL NICE WITH DS=INICE
C                     IF NEGATIVE, JUST EXPAND RANGE TO NEAREST INTEGERS
C                     IF THE DATA RANGE IS LESS THAN 1.
C
C         PARAMETERS
      INTEGER*4 PLMODE,N,K,INICE
      REAL*4 A(1),ERR(1),LEN
      REAL*8 A8(1),ERR8(1),EXTRM(3),SXTRM(4)
C
C         LOCALS
      REAL*8 MIN4/2D-38/
      REAL*8 MAX/1D37/, RANGE,DIV,PTOP,PBOT,VMIN,VMAX,SYMWID
      INTEGER*4 I,IANICE,L,NDIV
      LOGICAL DUBLP,ERRORP
C
      DUBLP=.FALSE.
      ERRORP=.FALSE.
      GOTO 1
C
      ENTRY SQALE8(PLMODE,A8,N,K,LEN,EXTRM,SXTRM,INICE)
      DUBLP=.TRUE.
      ERRORP=.FALSE.
      GOTO 1
C
      ENTRY SQALR(PLMODE,A,ERR,N,K,LEN,EXTRM,SXTRM,INICE)
      DUBLP=.FALSE.
      ERRORP=.TRUE.
      GOTO 1
C
      ENTRY SQALR8(PLMODE,A8,ERR8,N,K,LEN,EXTRM,SXTRM,INICE)
      DUBLP=.TRUE.
      ERRORP=.TRUE.
C
C           KEEP SAME LIMITS FOR NON-SORT LOOKER ONEAXE
    1 IF(PLMODE.GT.1) GOTO 170
C           SET UP SCALING
      SXTRM(4)=10.
      SXTRM(1)=0.
      SXTRM(2)=0.
      NDIV=LEN
      IANICE=IABS(INICE)
      IF(IANICE.GT.1) NDIV=IANICE
      DIV=NDIV
C           TEST FOR USER-SUPPLIED LIMITS FOR LOOKER PLOTS
      IF(PLMODE.GE.0 .AND.
     .  (EXTRM(1).NE.0.0 .OR. EXTRM(2).GT.EXTRM(1))) GOTO 100
C           NO USER-SUPPLIED LIMITS, SCAN ARRAY (IF ANY)
      IF(N.LE.0) GOTO 30
C           INITIALIZE MIN,MAX
      SXTRM(1)=MAX
      SXTRM(2)=-MAX
C         FIND MAX AND MIN OF A
      L=N*K
      DO 10 I=1,L,K
         IF(DUBLP) THEN
            IF(ERRORP) THEN
               PTOP=A8(I)+ERR8(I)*.5D0
               PBOT=A8(I)-ERR8(I)*.5D0
            ELSE
               PTOP=A8(I)
               PBOT=A8(I)
            ENDIF
         ELSE
            IF(ERRORP) THEN
               PTOP=A(I)+ERR(I)*.5
               PBOT=A(I)-ERR(I)*.5
            ELSE
               PTOP=A(I)
               PBOT=A(I)
            ENDIF
         ENDIF
         IF(PTOP.GT.SXTRM(2)) SXTRM(2)=PTOP
         IF(PBOT.LT.SXTRM(1)) SXTRM(1)=PBOT
   10 CONTINUE
      IF(SXTRM(2)+EXTRM(3).GT.SXTRM(1)-EXTRM(3)) THEN
         SXTRM(1)=SXTRM(1)-EXTRM(3)
         SXTRM(2)=SXTRM(2)+EXTRM(3)
      ENDIF
C
   30 IF(INICE.GT.0) THEN
C ROUND OFF RANGE
         GOTO 140
      ELSE IF(INICE.LT.0) THEN
C EXPAND RANGE TO INTEGERS FOR LOG PLOTS
         IF(SXTRM(2)-SXTRM(1).GE.1.) GOTO 140
         PBOT=DINT(SXTRM(1))
         IF(PBOT.GT.SXTRM(1)) PBOT=PBOT-1.
         SXTRM(1)=PBOT
         PTOP=DINT(SXTRM(2))
         IF(PTOP.LT.SXTRM(2) .OR. PTOP.EQ.PBOT) PTOP=PTOP+1.
         SXTRM(2)=PTOP
         GOTO 150
      ELSE
C KEEP RANGE AS IS
         GOTO 150
      ENDIF
C
C           TAKE USER-SUPPLIED LIMITS (WITHOUT CHANGE)
  100 SXTRM(1)=EXTRM(1)
      SXTRM(2)=EXTRM(2)
      NDIV=-NDIV
      IF(NDIV.GE.0) GOTO 150
C
C           GET ROUNDED LIMITS AND DIVISIONS
  140 VMIN=SXTRM(1)
      VMAX=SXTRM(2)
      CALL NICE8(SXTRM(2),SXTRM(1),RANGE,NDIV,DIV)
      IF(LEN.LE.0.) RETURN
C--MAKE ALLOWANCE FOR WIDTH OF A PLOTTING SYMBOL
C--(SHOULD MATCH SYMBOL SIZE IN SUBROUTINE LINE)
      IF(EXTRM(3).EQ.0.) THEN
         SYMWID=RANGE*0.05/LEN
         VMIN=VMIN-SYMWID
         VMAX=VMAX+SYMWID
         IF(VMIN.LT.SXTRM(1) .OR. VMAX.GT.SXTRM(2)) THEN
            SXTRM(1)=DMIN1(SXTRM(1),VMIN)
            SXTRM(2)=DMAX1(SXTRM(2),VMAX)
            DIV=DIV*(SXTRM(2)-SXTRM(1))/RANGE
         ENDIF
      ENDIF
C
C           SET PLOTTING SCALE
  150 IF(LEN.LE.0.) RETURN
      SXTRM(3)=(SXTRM(2)-SXTRM(1))/LEN
      SXTRM(4)=10.*DIV/LEN
C
C           SET UP FOR CALL TO 'LINE'
  170 IF(DUBLP) THEN
         A8(N*K+1)  =SXTRM(1)
         A8(N*K+K+1)=SXTRM(3)
      ELSE
         A(N*K+1)  =SXTRM(1)
         A(N*K+K+1)=SXTRM(3)
      ENDIF
      RETURN
      END
