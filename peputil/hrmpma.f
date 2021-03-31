      SUBROUTINE HRMPMA(LVEC,PREFIX, K, NAMES, TYPE, ORDER)
      IMPLICIT NONE
C         GENERATES NAMES FOR HARMONICS
C
C         PARAMETERS
      INTEGER*2 LVEC(1), ORDER
      INTEGER*4 TYPE, K
      CHARACTER*8 PREFIX, NAMES(2, 1)
C
C         LOCALS
      CHARACTER*8 ZHARM1 /'J( )'/, ZHARM2 /'J(  )'/
      CHARACTER*8 THARM1/'C( , )'/, THARM2/'C(  , )'/
      CHARACTER*8 FOURS(8)/'PA( )', 'PB( )', 'PC( )', 'PD( )', 'PAP( )',
     1 'PCP( )', 'PAPP', 'PBPP'/, GRIDS/'GRID'/
      CHARACTER*8 TEMP
      INTEGER*4 IP,IT,JT,L,N,NL,NP
C
      IF(TYPE.LT.0) GOTO 100
C
C         ZONALS
      IF(TYPE .GT. 1) GO TO 30
      TEMP=ZHARM1
      N=1
      DO 20 L=2,ORDER
      IF(LVEC(L-1) .EQ. 0) GO TO 20
      K = K + 1
      NAMES(1, K) = PREFIX
      IF(L.LT.10) GOTO 10
      TEMP=ZHARM2
      N=2
   10 CALL EBCDIX(L,TEMP,3,N)
 15   NAMES(2, K) = TEMP
 20   CONTINUE
      GO TO 900
C
C         TESSERALS
 30   CONTINUE
      IP=0
      IT=1
      TEMP=THARM1
      IF(TYPE.EQ.3) TEMP(1:1)='S'
      DO 80 L = 2, ORDER
      IF(L.NE.10) GOTO 40
      IT=2
      TEMP=THARM2
      IF(TYPE.EQ.3) TEMP(1:1)='S'
   40 CALL EBCDIX(L,TEMP,3,IT)
      JT=1
      DO 80 N=1,L
      IP=IP+1
      IF(N.NE.10) GOTO 50
      JT=2
      TEMP(8:8)=')'
   50 IF(LVEC(IP).LE.0) GOTO 80
      K = K + 1
      NAMES(1, K) = PREFIX
      CALL EBCDIX(N,TEMP,IT+4,JT)
      NAMES(2, K) = TEMP
 80   CONTINUE
      GOTO 900
C
C           FOURIER SHAPE MODEL
  100 IF(TYPE.LT.-1) GOTO 200
      IP=0
      NP=4
      DO 120 N=1,6
      TEMP=FOURS(N)
      IF(N.EQ.5) NP=5
      NL=1
      DO 120 L=1,20
      IP=IP+1
      IF(L.NE.10) GOTO 110
      NL=2
      TEMP(NP+NL:NP+NL)=')'
  110 IF(LVEC(IP).LE.0) GOTO 120
      K=K+1
      NAMES(1,K)=PREFIX
      CALL EBCDIX(L,TEMP,NP,NL)
      NAMES(2,K)=TEMP
  120 CONTINUE
      DO 130 N=7,8
      IP=IP+1
      IF(LVEC(IP).LE.0) GOTO 130
      K=K+1
      NAMES(1,K)=PREFIX
      NAMES(2,K)=FOURS(N)
  130 CONTINUE
      GOTO 900
C
C           GRID SHAPE MODEL
  200 CONTINUE
      TEMP=GRIDS
      DO 210 L=1,4
      K=K+1
      NAMES(1,K)=PREFIX
      CALL EBCDIX(L,TEMP,5,1)
      NAMES(2,K)=TEMP
  210 CONTINUE
C
  900 RETURN
      END
