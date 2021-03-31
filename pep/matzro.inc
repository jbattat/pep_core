c A priori B matrix (dimension declared in BLKDATA8)
c
c **** MATZRO MUST BE INCLUDED BEFORE NRMMAT IF BOTH APPEAR ****
c
c Logically equivalence with NRMMAT for un-overlayed PEP.
c In this case, one extra doubleword is allocated.
      COMMON/NRMMAT/ Coeff(1)
c Maintain separate common for overlayed PEP.
C      COMMON/MATZRO/ Coeff(1)
C
      real*10    Coeff
