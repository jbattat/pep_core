      REAL*10 JD(3),W(14)
      read(2,100) JD,W
  100 FORMAT(1P,3D25.17)
      write(1) JD
      write(1) W
      STOP
      END
