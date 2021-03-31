      SUBROUTINE DISPLY(DATA,DAT8)
C         PRINTS AND GRAPHS DATA SELECTED FOR TRANSFORMATION.
C
C        PARAMETERS
      REAL*4 DATA(9)
      REAL*8 DAT8(9)
C
C         COMMON
      include 'inodtabc.inc'
      include 'graphs.inc'
      include 'misc.inc'
      include 'pltlab.inc'
      include 'span.inc'
C
C         LOCAL VARIABLES
      CHARACTER*8 G/'JD.TIME '/
      REAL*8 FIRST
C
C         STATEMENT FUNCTIONS FOR INDICES TO MASTER DATA ARRAY
C         DATA - INDEXING INTO 2 INTERLEAVED TIERS OF REAL*4
      include 'indices.inc'
C
C         TIME - INDEXING INTO ONE TIER OF REAL*8
      KT(ISPAN,IPOINT)= (NSPAN+ISPAN-1)*(NPOINT+2)+IPOINT
C
C         WEIGHT - INDEXING INTO FIFTH TIER OF REAL*4
      KW(ISPAN,IPOINT)= (4*NSPAN+ISPAN-1)*(NPOINT+2)+IPOINT
C
C         HEADING
      WRITE(PRINT,10)
   10 FORMAT('1*** SELECTION SECTION ***'/)
C
C         FIRST DO STATISTICS -
C         CYCLE THRU DATA TO COUNT UP ZEROES
      N = NPOINT
      IF (PADZ) N = NPOINT/2
      NPTTOT=0
      DO I=1,NSPAN
         IDAT= KD(I, 1)
         NZERO = 0
         DO J=1,N
            IF(DATA(IDAT) .EQ. 0.0) NZERO = NZERO + 1
            IDAT= IDAT+ 2
         END DO
         NPTTOT=NPTTOT+NPTSPN(I)
         M = 100*(N - NZERO)/N
         WRITE(PRINT,15) I,NPTSPN(I),NZERO,M
   15    FORMAT (/' STATISTICS FOR SPAN', I4/
     .    ' NUMBER OF DATA POINTS INCLUDED', I6/
     1    ' NUMBER OF ZERO ENTRIES IN DATA', I6/
     2    ' PERCENT OF DATA ARRAY FILLED  ', I6)
      END DO
      WRITE(PRINT,25) NPTTOT
   25 FORMAT(/' TOTAL DATA POINTS INCLUDED IN ALL SPANS',I6)
C           SET UP PLOT LABELS
      LLX=21
      CALL MVC('Time relative to span',1,LLX,XLAB,1)
      CALL SETLAB(YLAB,LLY,NDV,OBSNM,OBSUN,OBSLU)
      LLT=26
      CALL MVC('Data selected for span ',1,24,PTITLE,1)
C
C         SEE IF FULL DISPLAY WANTED
      IF(MOD(IPRINT/2,2).EQ.1) THEN
C
C         PRINT DATA FOR EACH SPAN
         DO I = 1, NSPAN
            WRITE(PRINT,30) I
   30       FORMAT(//' DATA ARRAY FOR SPAN', I4)
            WRITE(PRINT,1001)
 1001       FORMAT(/4X, 'IND. VAR.', 15X, 'DEP. VAR.', 6X, 'WEIGHT')
            WRITE(PRINT,35) G,YLAB(1:LLY)
   35       FORMAT(4X,A8,16X,A)
            IDAT= KD(I, 1)
            ITIM= KT(I, 1)
            IWGT= KW(I, 1)
            DO J=1,N
               WRITE(PRINT, 1002) DAT8(ITIM), DATA(IDAT),DATA(IWGT)
 1002          FORMAT (F25.13, 1PE14.5, 0PF9.3)
               IDAT= IDAT+ 2
               ITIM= ITIM+ 1
               IWGT= IWGT+ 1
            END DO
         END DO
      ENDIF
C
C         SEE WHETHER TO PLOT
      IF(MOD(IPLOT/2,2).NE.1) RETURN
C
C         REPLACE ABSOLUTE TIMES BY TIMES RELATIVE TO FIRST IN SPAN
C         AND PLOT TIME VS DATA FOR EACH SPAN
      DO I = 1, NSPAN
         IDAT= KD(I, 1)
         ITIM= KT(I, 1)
         ITIM4=ITIM*2-1
         FIRST = DAT8(ITIM)
         DATA(ITIM4)=0.
         K=ITIM4
         DO J = 2, N
            K=K+2
            DATA(K)=DAT8(ITIM+J-1)-FIRST
         END DO
C
C         TITLE AND PLOT DATA
         CALL EBCDI(I,PTITLE(25:26),2)
         CALL PLOTER(-1,DATA(ITIM4),DATA(IDAT),N,2)
      END DO
      RETURN
      END
