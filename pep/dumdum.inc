C           INTEGRATION COEFFICIENTS
      COMMON/DUMDUM/ A(6*i_mxeqn),B(6*i_mxeqn),C(6*i_mxeqn),
     1 D(6*i_mxeqn),L4,N1,N2,L5,EE(6*i_mxeqn)
      REAL*10 A,B,C,D,EE
      INTEGER*4 L4,N1,N2,L5
      REAL*10 AA(6,i_mxeqn,5)
      EQUIVALENCE (AA(1,1,1),A(1))
