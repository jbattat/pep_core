      SUBROUTINE SYMBOL(X,Y,HEIGHT,STRING,THETA,NCHAR)
C
C This is the standard FORTRAN call to subroutine SYMBOL:
C
C (X,Y)   are the starting position (lower left corner) of the first
C         character of the output string.  X, Y are measured in inches
C         from the origin.  If X or Y is 999, then the corresponding
C         coordinate is taken from the position at the end of the last
C         call to SYMBOL or NUMBER.
C HEIGHT  is the height of the characters to be output in inches.  If
C         HEIGHT is not positive, the height and angle of the previous
C         call are used instead.
C STRING  is a numerical array containing the characters to be drawn
C         or an integer index pointing to a single character within the
C         font (depending on NCHAR)
C         If the integer index is -1 or 16, the symbol is a dot.
C THETA   is the angle of the baseline along which the characters
C         are to be plotted, measured in degrees counter-clockwise
C         from the horizontal (x-axis) axis.
C NCHAR   is the number of characters to be plotted.  If non-positive,
C         STRING is a right-justified integer character code.  If NCHAR
C         is less than -1, the movement to (X,Y) is made with the pen
C         down.
C
C Separate file to define the character font (could be ASCII or EBCDIC)
      COMMON/SYMBVT/ NFONT,IPTR(2,257),IVECT(10)
      INTEGER*4 NFONT
      INTEGER*2 IPTR,IVECT
C
      INCLUDE '../peputil/plttyp.inc'
C
      REAL*4 THHW/0./, HTHW/0.25/, DX/1./, DY/0./
      INTEGER*4 ICHR
      CHARACTER*1 STRING(4),CHR(4)
      EQUIVALENCE (ICHR,CHR)
      LOGICAL HW
C DOT COORDINATES
c     REAL*4 VECT(2,5)/2.,2., 2.,2.04, 1.97,1.98, 2.03,1.98, 2., 2.04/
      REAL*4 VECT(2,5)/2.,2., 2.,2.5, 1.5,2., 2.,1.5, 2.5,2./
C
      IPEN=3
C Divide by 7 for convenient sub-units
      HT=HEIGHT/7.
      NCHT=NCHAR
      HW=.TRUE.
C
C Marking symbols: integer, maybe connect with lines, maybe bigger
      IF(NCHAR.LE.0) THEN
         HW=.FALSE.
         IF(NCHAR.LT.-1) IPEN=2
         NCHT=1
         DO I=1,4
            CHR(I)=STRING(I)
         END DO
         IF(ICHR.EQ.-1) ICHR=16
         ICHR=MOD(ICHR,256)
         IF(ICHR.LT.0) ICHR=ICHR+256
         IF(ICHR.LE.13) THEN
            HT=HEIGHT/4.
         ENDIF
      ENDIF
C Normal call with positive height: set relative X and Y units
C Start at level 1 (neither super- nor sub-script)
      IF(HT.GT.0.) THEN
         HTHW=HEIGHT
         THHW=THETA
         THETAR=THHW*1.7453292E-2
         DY=SIN(THETAR)
         DX=COS(THETAR)
         LEVEL=1
         HTSV=HT
         DXU=HT*DX
         DYU=HT*DY
      ENDIF
      IF(LOCALE.NE.8) HW=.FALSE.
C Set reference point 2 units down and left
      IF(HW) THEN
         XORG=X
         YORG=Y
      ELSE
         XORG=X-2.*DXU+2.*DYU
         YORG=Y-2.*DXU-2.*DYU
      ENDIF
      IF(X.EQ.999.) XORG=XLAST
      IF(Y.EQ.999.) YORG=YLAST
      XLAST=XORG
      YLAST=YORG
      XW=XLAST
      YW=YLAST
      IF(HW) THEN
         CALL PLOT(XORG,YORG,IPEN)
         IF(NCHAR.GT.0) THEN
            CALL SYMBHW(STRING,NCHAR,HTHW,THHW)
         ELSE
            CALL SYMBHW(CHAR(ICHR),1,HTHW,THHW)
         ENDIF
         XLAST=XORG+NCHT*HTHW*PROPOR*DX
         YLAST=YORG+NCHT*HTHW*PROPOR*DY
         RETURN
      ENDIF
      DO 1000 JCHR=1,NCHT
         IF(NCHAR.GT.0) THEN
            ICHR=ICHAR(STRING(JCHR))
         ENDIF
         IF(ICHR.EQ.16) GOTO 420
         N=MOD(ICHR,NFONT)+1
         IVPT1=IPTR(1,N)
         IVPT2=IVPT1-1+IPTR(2,N)
         DO 400 IVPTR=IVPT1,IVPT2,2
            IVX=IVECT(IVPTR)
            IVY=IVECT(IVPTR+1)
            IF(IVX.LT.15) THEN
               XW=XW+DXU*IVX-DYU*IVY
               YW=YW+DYU*IVX+DXU*IVY
               CALL PLOT(XW,YW,IPEN)
               IPEN=2
               XW=XLAST
               YW=YLAST
            ELSE
               IF(IVY.EQ.0) IPEN=3
               IF(IVY.GT.0) GOTO 500
            ENDIF
  400    CONTINUE
         GOTO 450
  420    DO IVPTR=1,5
            CALL PLOT(XW+DXU*VECT(1,IVPTR)-DYU*VECT(2,IVPTR),
     .       YW+DYU*VECT(1,IVPTR)+DXU*VECT(2,IVPTR),IPEN)
            IPEN=2
         END DO
  450    XW=XW+7.*PROPOR*DXU
         YW=YW+7.*PROPOR*DYU
         GOTO 600
C
C Special operations flagged by IVX=15
  500    IF(IVY.EQ.1) THEN
C IVY=1 - enter superscript mode or leave subscript mode
            IF(LEVEL.EQ.0) THEN
               LEVEL=1
               DXU=HTSV*DX
               DYU=HTSV*DY
               XW=XW-2.*DYU
               YW=YW+2.*DXU
            ELSE IF(LEVEL.EQ.1) THEN
               XW=XW-4.*DYU
               YW=YW+4.*DXU
               LEVEL=2
               DXU=DXU*0.7
               DYU=DYU*0.7
            ENDIF
C IVY=2 - leave superscript mode or enter subscript mode
         ELSE IF(IVY.EQ.2) THEN
            IF(LEVEL.EQ.2) THEN
               LEVEL=1
               DXU=HTSV*DX
               DYU=HTSV*DY
               XW=XW+4.*DYU
               YW=YW-4.*DXU
            ELSE IF(LEVEL.EQ.1) THEN
               LEVEL=0
               XW=XW+2.*DYU
               YW=YW-2.*DXU
               DXU=DXU*0.7
               DYU=DYU*0.7
            ENDIF
C IVY=3 - linefeed
         ELSE IF(IVY.EQ.3) THEN
            XORG=XORG+12.*DYU
            YORG=YORG-12.*DXU
            XW=XORG
            YW=YORG
C IVY=4 - backspace
         ELSE
            XW=XW-6.*DXU
            YW=YW-6.*DYU
         ENDIF
  600    CONTINUE
         IPEN=3
         XLAST=XW
         YLAST=YW
 1000 CONTINUE
      RETURN
      END
