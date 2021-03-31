      SUBROUTINE SUICID(MESSAG,NN)                                      00000010
      CHARACTER*4 MESSAG(32), NORMAL(2)/' NOR','MAL '/
C                                                                       00000030
      N =IABS(NN)                                                       00000040
C                                                                       00000050
      WRITE(6,31) (MESSAG(I),I=1,N)                                     00000060
   31 FORMAT(//(1X,32A4))                                               00000070
      IF(NN.LT.0) RETURN                                                00000080
      IF(MESSAG(1).EQ.NORMAL(1).AND.MESSAG(2).EQ.NORMAL(2)) STOP        00000090
      STOP 20                                                           00000100
      END                                                               00000110
