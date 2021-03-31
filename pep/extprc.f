      SUBROUTINE EXTPRC
C
C  THE CALL TO EXTPRC SETS A FLAG: -1 MEANS USE DOUBLE PRECISION
C                                  0  MEANS USE HARDWARE EXTENDED
C                                  +1  SOFTWARE EXTENDED PRECISION
C                                  +2   USE SOFTWARE EXTENDED PRECISION
C                                       FOR DIVIDE, HARDWARE EXTENDED
C                                    PRECISION FOR THE OTHER OPERATIONS
C
C  IF FLAG=-1, IT IS ASSUMED THAT ONLY THE 1ST 64 BITS OF THE CALLING
C  ARGUMENTS ARE SIGNIFICANT (TREATED AS REAL*8 VARIABLES), BUT ALL
C  128 BITS ARE STILL STORED/FETCHED AS INDICATED.
C
C        CALL XLOAD(A)       THE 128 BITS AT "A" ARE LOADED INTO LOCAL
C                            FLOATING POINT STORAGE OF THE RIGHT TYPE.
C
C        CALL XSTORE(A)      THE 128 BITS AT "A" ARE SET FROM THE LOCAL
C                            FLOATING POINT QUANTITY (OR 64 IF FLAG=-1).
C
C        CALL XADD(A)        THE 128 BITS AT "A" ARE ADDED TO LOCAL
C                            FLOATING POINT STORAGE.
C
C        CALL XSUB(A)        THE 128 BITS AT "A" ARE SUBTRACTED FROM
C                            LOCAL FLOATING POINT STORAGE.
C
C        CALL XMUL(A)        THE 128 BIT FLOATING POINT NUMBER IN LOCAL
C                            FLOATING POINT STORAGE IS MULTIPLIED BY
C                            THE 128 BIT FLOATING POINT NUMBER "A".
C
C        CALL XDIV(A)        THE 128 BIT FLOATING POINT NUMBER IN LOCAL
C                            FLOATING POINT STORAGE IS DIVIDED BY
C                            THE 128 BIT FLOATING POINT NUMBER "A".
C
C        CALL STORND(B)      THE 128-BIT EXTENDED PRECISION NUMBER
C                            IN LOCAL FLOATING POINT STORAGE IS COPIED
C                            TO A 64-BIT VALUE AT LOCATION "B".  THE
C                            LOCAL STORAGE IS UNCHANGED.
C
C        E = STORNE(B)       THE DOUBLE PRECISION NUMBER B IS ROUNDED
C                            TO SINGLE PRECISION AND RETURNED AS THE
C                            FUNCTION VALUE.
C                            WITH CALLING PROGRAM CODED AS ILLUSTRATED,
C                            ROUNDED RESULT IS STORED IN THE SINGLE
C                            PRECISION NUMBER E. MUST SPECIFY
C                            REAL*4 STORNE IN CALLING PROGRAM.
C
C ENTRY POINTS XLOAD8, XADD8, XSUB8, XMUL8, AND XDIV8 ARE SIMILAR, BUT
C EXPECT REAL*8 ARGUMENTS AND EXTEND THE ARGUMENT TO 128 BITS.
C
C EXAMPLE 1:
C        G = (A*B-C)/D+E*F
C        real*10 A(2),B(2),C(2),D(2),E(2),F(2),G(2)
C        CALL XLOAD(A)
C        CALL XMUL(B)
C        CALL XSUB(C)
C        CALL XDIV(D)
C        CALL XSTORE(G)
C        CALL XLOAD(E)
C        CALL XMUL(F)
C        CALL XADD(G)
C        CALL XSTORE(G)
C
C
C EXAMPLE 2:
C        Y = A(1)+A(2)*X+A(3)*X**2+. . .+A(N+1)*X**N
C WE EMPLOY THE USUAL TRICK FOR EVALUATING POLYNOMIALS AND WRITE THIS
C IN THE FORM
C        Y = A(1)+X*(A(2)+X*(A(3)+. . .+X*(A(N)+X*A(N+1))). . . )
C        real*10 A(2,N+1),X(2),Y(2)
C        CALL XLOAD(A(1,N+1))
C        DO 11 I=1,N
C        J = N+1-I
C        CALL XMUL(X)
C        CALL XADD(A(1,J))
C     11 CONTINUE
C        CALL XSTORE(Y)
C
C
C EXAMPLE 3:
C        MATRIX MULTIPLICATION C=A*B
C        real*10 A(2,N,N),B(2,N,N),C(2,N,N),SUM(2)
C        DO 15 I=1,N
C        DO 13 J=1,N
C        SUM(1) = 0.0_10
C        SUM(2) = 0.0_10
C        DO 11 K=1,N
C        CALL XLOAD(A(1,I,K))
C        CALL XMUL(B(1,K,J))
C        CALL XADD(SUM)
C        CALL XSTORE(SUM)
C     11 CONTINUE
C        C(1,I,J) = SUM(1)
C        C(2,I,J) = SUM(2)
C     13 CONTINUE
C     15 CONTINUE
C
      real*10 VAL,ARG,XARG
      real*10 VAL8,ARG8,XARG8
      EQUIVALENCE (VAL,VAL8),(XARG,XARG8)
      INTEGER*2 INFLG
      INTEGER*4 FLAG
C
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C LOAD EXTENDED PRECISION FLOATING POINT REGISTER WITH CALLING VAR
      ENTRY XLOAD(ARG)
      VAL=ARG
      RETURN
C
      ENTRY XLOAD8(ARG8)
      IF(FLAG.LT.0) THEN
        VAL8=ARG8
      ELSE
        VAL=ARG8
      ENDIF
      RETURN
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C STORE EXTENDED PRECISION FLOATING POINT REGISTER IN CALLING VAR
      ENTRY XSTORE(ARG)
      ARG=VAL
      RETURN
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C EXTENDED PRECISION ADD CALLING VARIABLE TO FLOATING POINT REGISTER
      ENTRY XADD(ARG)
      IF(FLAG.LT.0) THEN
        XARG=ARG
        VAL8=VAL8+XARG8
      ELSE
        VAL=VAL+ARG
      ENDIF
      RETURN
C
      ENTRY XADD8(ARG8)
      IF(FLAG.LT.0) THEN
        VAL8=VAL8+ARG8
      ELSE
        VAL=VAL+ARG8
      ENDIF
      RETURN
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C EXTENDED PRECISION SUBTRACT CALLING VAR FROM FLOATING POINT REGISTER
      ENTRY XSUB(ARG)
      IF(FLAG.LT.0) THEN
        XARG=ARG
        VAL8=VAL8-XARG8
      ELSE
        VAL=VAL-ARG
      ENDIF
      RETURN
C
      ENTRY XSUB8(ARG8)
      IF(FLAG.LT.0) THEN
        VAL8=VAL8-ARG8
      ELSE
        VAL=VAL-ARG8
      ENDIF
      RETURN
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C EXTENDED PRECISION MULTIPLY FLOATING POINT REGISTER BY CALLING VAR
      ENTRY XMUL(ARG)
      IF(FLAG.LT.0) THEN
        XARG=ARG
        VAL8=VAL8*XARG8
      ELSE
        VAL=VAL*ARG
      ENDIF
      RETURN
C
      ENTRY XMUL8(ARG8)
      IF(FLAG.LT.0) THEN
        VAL8=VAL8*ARG8
      ELSE
        VAL=VAL*ARG8
      ENDIF
      RETURN
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C EXTENDED PRECISION DIVIDE FLOATING POINT REGISTER BY CALLING VAR
      ENTRY XDIV(ARG)
      IF(FLAG.LT.0) THEN
        XARG=ARG
        VAL8=VAL8/XARG8
      ELSE
        VAL=VAL/ARG
      ENDIF
      RETURN
C
      ENTRY XDIV8(ARG8)
      IF(FLAG.LT.0) THEN
        VAL8=VAL8/ARG8
      ELSE
        VAL=VAL/ARG8
      ENDIF
      RETURN
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C STORE ROUNDED EXT.PREC. INTO DOUBLE PREC. CALLING VARIABLE
      ENTRY STORND(ARG8)
      IF(FLAG.LT.0) THEN
        ARG8=VAL8
      ELSE
        ARG8=VAL
      ENDIF
      RETURN
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C
C      SET FLAG FOR EXTENDED PRECISION ARITHMETIC TYPE
      ENTRY EXTFLG(INFLG)
      FLAG=INFLG
      RETURN
      END
C * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      REAL FUNCTION STORNE(ARG8)
      real*10 ARG8
C
C ROUND DOUBLE PRECISION CALLING VAR. INTO SINGLE PRECISION VAR.
C
      STORNE=ARG8
      RETURN
      END
