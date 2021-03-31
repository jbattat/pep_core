C           NUMERICAL INTEGRATION RESULTS
      COMMON/OUTPUT/ Y(6*i_mxeqn,4),F(6*i_mxeqn,3),K1,K2
      REAL*10 Y,F
      INTEGER*4 K1,K2
      REAL*10 FF(6,i_mxeqn,3)
      EQUIVALENCE (FF(1,1,1),F(1,1))
