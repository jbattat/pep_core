C           DATA FROM CENTRAL-BODY AND TARGET-BODY DATA SETS
C           NOTE: LAST DIMENSION OF CENTBD IS FROM GLOBDEFS, AS ARE LAST
C           TWO DIMENSIONS OF TRGBD AND SO ON
      COMMON/TAPDTP/ CENTBD(6,24,i_mxplprt+1),
     . TRGBD(6,24,i_mxplprt+1,i_mxtrg),TINT(0:i_mxtrg),
     1 TFRACT(3,0:i_mxtrg),JDT1(0:i_mxtrg),JDT2(0:i_mxtrg),
     2 INTTX(0:i_mxtrg),INTT5(0:i_mxtrg),IDIRT(0:i_mxtrg),
     3 IPART(0:i_mxtrg),JDT(3,0:i_mxtrg),IPVELT(3,0:i_mxtrg),
     4 IBADT(3,0:i_mxtrg),JPLNTG(0:i_mxtrg),LMVLT(0:i_mxtrg),
     5 NQT(0:i_mxtrg),KQT(i_mxplp,0:i_mxtrg),KRT(i_mxplp,0:i_mxtrg),
     6 NPLPT(0:i_mxtrg),KPT(i_mxeqn-1,0:i_mxtrg),MPPT(0:i_mxtrg),
     7 NPTSPR(0:i_mxtrg)
      REAL*10 CENTBD,TRGBD,TINT,TFRACT
      INTEGER*4 JDC1,JDT1,JDC2,JDT2,INTTX,INTT5,
     1 IDIRT,IPART,JDT,IPVELT,IBADT,
     2 JPLNTG,LMVLT,NQT
      INTEGER*2 KQT,KRT,NPLPT,KPT,MPPT,NPTSPR
      INTEGER*4 JPLNC
      INTEGER*2 NPLPC
      EQUIVALENCE (JPLNC,JPLNTG),(NPLPC,NPLPT)
