      SUBROUTINE FNAXE(PLMODE,N)
      IMPLICIT NONE
C           FIND LENGTH OF X AXIS:
C           IF ONEAXE=F, LESSER OF XLEN AND N/XDEN (AS LONG AS BOTH
C                        XLEN AND XDEN NON-ZERO)
C           IF ONEAXE=T, XLEN
C     ALSO, RE-INITIALIZE YL FROM YLEN
C
C        PARAMETERS
      INTEGER*4 PLMODE,N
C             PLMODE - TYPE OF PLOTTING TO DO:
C                  -1  NON-LOOKER PLOTS
C                   0  IMMEDIATE LOOKER PLOT
C                   1  FIRST SERIES OF ONEAXE PLOT
C                   2  SUBSEQUENT SERIES FOR ONEAXE PLOT
C                   3  FINISH OFF ONEAXE PLOT
C             N - NUMBER OF DATA POINTS TO PLOT
C
C        COMMONS
      include 'graphs.inc'
      include 'pltfrm.inc'
      include 'pltlab.inc'
      include 'plttyp.inc'
C
C        LOCAL
      REAL*4 TNOM,XLT
C
C           LENGTH OF Y AXIS
      YL=YLEN
C           FIND LENGTH OF X AXIS
      XL=XLEN
      IF(PLMODE.GT.0) GOTO 100
      XLT=XLEN
      IF(XDEN.GT.0.0) XLT=N/XDEN
      IF(XL.LE.0.0.OR.XLT.LT.XL) XL=XLT
  100 IF(XL.GT.XFRAME-XMRG) XL=XFRAME-XMRG
      IF(XL.LT.1.0) XL=1.0
      XL=AINT(XL)
      XLNOM=XL
      TNOM=(LLX*0.2*PROPOR+XL)/2.0
      IF(TNOM.GT.XLNOM) XLNOM=TNOM
      RETURN
      END
