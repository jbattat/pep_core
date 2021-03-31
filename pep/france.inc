C           QUANTITIES THAT GO INTO IPLCON PERIPHERAL DATA SET
      COMMON/FRANCE/ DUMCON(12,u_mxpl+4),SECP(u_mxpl+4),
     1 DTCON(30,u_mxpl+4),DUMEPS(6,u_mxpl+4),JD1(u_mxpl+4),
     2 JD2(u_mxpl+4),INTP1(u_mxpl+4),INTP2(u_mxpl+4),
     3 KKK(u_nmprm,u_mxpl+4),INT(u_mxpl+4),IHRP(u_mxpl+4),
     4 IMINP(u_mxpl+4),KKP(100,u_mxpl+4),
     5 NDUMKI(u_mxpl+4),KDUMI(99,u_mxpl+4)
      REAL*10 DUMCON,SECP,DTCON
      REAL*4 DUMEPS
      INTEGER*4 JD1,JD2,INTP1,INTP2
      INTEGER*2 KKK,INT,IHRP,IMINP,KKP,NDUMKI,KDUMI
