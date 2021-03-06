      REAL*10 FUNCTION LAGRNG(T,TTAB,FTAB,IND)
C        F.AMUCHASTEGUI - 1970 APRIL - FUNCTION LAGRNG
C           J.F.CHANDLER - 1986 SEP - VARIABLE ARRAY SIZE
C        PERFORMS LAGRANGE INTERPOLATION FOR UN-EQUALLY SPACED
C        VALUES OF THE FUNCTION
      REAL*10 T,TTAB(99),FTAB(99),TT(99),A
C T    - INDEPENDENT VARIABLE FOR EVALUATION (IGNORED IF IND>0)
C TTAB - TABLE OF INPUT INDEPENDENT VARIABLES (IGNORED IF IND>0)
C FTAB - TABLE OF CORRESPONDING FUNCTION VALUES
C IND  - POS => COEFFICIENTS ALREADY SET UP, ZERO OR NEG => MUST SET UP
C        ARRAY SIZE IS IABS(IND), EXCEPT 0 OR 1 IMPLIES SIZE=4
C
      NTAB=4
      IF(IND.NE.0.AND.IND.NE.1) NTAB=IABS(IND)
      IF(IND.GT.0) GOTO 30
C
C           EVALUATE LAGRANGE INTERPOLATION COEFFICIENTS
      DO 20 I=1,NTAB
      A=1._10
      DO 10 J=1,NTAB
      IF(J.EQ.I) GOTO 10
      A=A*(T-TTAB(J))/(TTAB(I)-TTAB(J))
   10 CONTINUE
   20 TT(I)=A
C
C        LAGRANGE INTERPOLATION (ARRAY ALREADY SET UP)
   30 LAGRNG=0._10
      DO 40 J=1,NTAB
   40 LAGRNG=LAGRNG+TT(J)*FTAB(J)
C
      RETURN
      END
