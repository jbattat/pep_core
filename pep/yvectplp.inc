C           CENTRAL/TARGET-BODY Y-VECTORS, COORDINATES, PARTIALS
      COMMON/SBDTA/ YPLP(5,6,2,i_mxplp,0:i_mxtrg),
     . YTP(6,i_mxplp,0:i_mxtrg)
      REAL*10 YPLP,YTP
      REAL*10 DCCOR(6,i_mxplp),DTCOR(6,i_mxplp,1)
      EQUIVALENCE (YTP(1,2,0),DCCOR),(YTP(1,2,1),DTCOR)
