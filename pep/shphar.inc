C           COEFFICIENTS, L-VECTORS, AND CONTROLS FOR CURRENT BODY
      COMMON/SHPHAR/ NSHP,SCNTRL(9),NZ,NT,I4FILL,
     1 ZONE(24),CTESS(209),STESS(209),R8FILL(58),
     2 LZ(24),LC(209),LS(209),I2FILL(558),SHPZER,SHPNIT
      REAL*10 ZONE,CTESS,STESS,R8FILL
      REAL*4 SCNTRL
      INTEGER*4 NSHP,I4FILL
      INTEGER*2 NZ,NT,LZ,LC,LS,I2FILL
      LOGICAL*1 SHPZER,SHPNIT
      INTEGER*4 NGDPTS
      EQUIVALENCE (SCNTRL(9),NGDPTS)
      INTEGER*2 LSH(1000)
      EQUIVALENCE (LSH(1),LZ(1))
      INTEGER*2 LFOUR(122)
      EQUIVALENCE(LFOUR(1),LZ(1))