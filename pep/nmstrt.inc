C           STARTING VALUES FOR INTEGRATION
      COMMON/NMSTRT/ A0(6*i_mxeqn),B0(6*i_mxeqn),C0(6*i_mxeqn),
     1 D0(6*i_mxeqn),EE0(6*i_mxeqn),F0(6*i_mxeqn),
     2 HSAVE,BOTH,IBOTH,NSIGN
      REAL*10 A0,B0,C0,D0,EE0,F0,HSAVE,BOTH
      INTEGER*4 IBOTH,NSIGN
      REAL*10 AA0(6,i_mxeqn,5),FF0(6,i_mxeqn)
      EQUIVALENCE (AA0(1,1,1),A0(1)),(F0(1),FF0(1,1))
