C           DATA FROM PLANET,MOON,EARTH ROTATION TAPES
C           NOTE: FOR MINIMUM SIZE OF KI.. ARRAYS SEE 'MAXKIDAT'
      COMMON/ROTDTA/ PLNMON(6,i_mxplprt+1,15),EARTH(6,i_mxplprt+1,15),
     1 FMR(3),FER(3),FRMR1,FRMR2,FRER1,FRER2,JMR,JER,JDMR(3),JDER(3),
     2 JDMR1,JDMR2,JDER1,JDER2,IPARMR,IPARER,LPARMR,LPARER,KMR(u_nmprm),
     3 KER(u_nmprm),IMRVEL(3),IERVEL(3),NKIMR,KIMR(99),NKIER,KIER(99)
      REAL*10 PLNMON,EARTH,FMR,FER,FRMR1,FRMR2,FRER1,FRER2
      INTEGER*4 JMR,JER,JDMR,JDER,JDMR1,JDMR2,JDER1,JDER2,
     1 IPARMR,IPARER,LPARMR,LPARER,IMRVEL,IERVEL
      INTEGER*2 KMR,KER,NKIMR,KIMR,NKIER,KIER