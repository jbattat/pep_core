      SUBROUTINE SORTER(PLMODE,DATA,N,K)
C           FILL PLOTTING SYMBOL FOR SERIES INTO ARRAY
C
C             PARAMETERS
C
C        PLMODE - TYPE OF PLOTTING TO DO:
C                  -1  NON-LOOKER PLOTS
C                   0  IMMEDIATE LOOKER PLOT
C                   1  FIRST SERIES OF ONEAXE PLOT
C                   2  SUBSEQUENT SERIES FOR ONEAXE PLOT
C                   3  FINISH OFF ONEAXE PLOT
C        DATA 2-D ARRAY WITH (Y,X,SYMBOL,...) FOR EACH OF 'N' POINTS
C        N    NUMBER OF POINTS = LENGTH OF BOTH X AND Y ARRAYS
C        K    REPEAT FACTOR (USE EVERY K-TH MEMBER OF THE ARRAYS)
C
C **********************************************************************
C
      INTEGER*4 PLMODE,N,K
      REAL*4 DATA(1)
C
C        COMMONS
      include 'graphs.inc'
      REAL*4 SYM
      EQUIVALENCE (NSYM,SYM)
C
C             LOCALS
  110 L=N*K
      DO 120 J=1,L,K
  120 DATA(J+2)=SYM
      RETURN
      END
