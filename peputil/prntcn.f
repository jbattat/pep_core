      SUBROUTINE PRNTCN
C        J.F.CHANDLER - 1978 MAY
C  PRINT ALL INPUT CONSTANTS ON UNIT 'PRINT'
C
C         COMMON
      include 'burst.inc'
      include 'bursu.inc'
      include 'checks.inc'
      include 'graphs.inc'
      include 'hrmnsc.inc'
      REAL*4 HMS(3,4)
      EQUIVALENCE (HMS(1,1),STHR)
      include 'inodtabc.inc'
      include 'misc.inc'
      include 'namen.inc'
      include 'pltfrm.inc'
      include 'skip.inc'
      include 'sortv.inc'
      include 'span.inc'
      REAL*8 DYS(4)
      EQUIVALENCE (DYS(1),START)
C
      CHARACTER*8 CRITER(10)/'RSITE   ','SSITE   ','SERIES  ',
     1 'NSEQ    ','NCODF   ','NPLNT   ','NCENTB  ','SPOT    ',
     . 'SPOT2   ','FREQ    '/
      CHARACTER*8 STSTP(4)/'START   ','STOP    ','SIZSPN  ','DELTA   '/
      CHARACTER*2 ISTSTP(4)/'ST','SP','SZ','DE'/
      CHARACTER*8 BLANK/'        '/
      CHARACTER*4 BLNK4
      EQUIVALENCE (BLANK,BLNK4)
      CHARACTER*4 EX/'  X '/,WHY/'  Y '/
      CHARACTER*2 EX2(2),WHY2(2)
      EQUIVALENCE (EX2(1),EX),(WHY2(1),WHY)
C
      WRITE(PRINT,100)
  100 FORMAT('1*** ALL INPUT PARAMETERS ***'/)
C
      WRITE(PRINT,110) CONTRL,PRINT,SPOOL,OBSOUT,OBSLIB,IPRINT,REPEAT,
     1 TASK,LOOK,FFT,WTFT,KBYTES,ISEED,RNDMDT
  110 FORMAT('0   DATA SETS AND CONTROL SWITCHES'/ '  CONTRL=',I2,
     1 '   PRINT=',I3,'   SPOOL=',I3,'   OBSOUT=',I2,'   OBSLIB=',10I3/
     2 '  IPRINT=',I2,'   REPEAT=',I2,'   TASK=',I4,'   LOOK=',I4,
     3 '     FFT=',L2,'     WTFT=',l2/
     . ' KBYTES=',I4,'   ISEED=',I6,'   RNDMDT=',L2)
      WRITE(PRINT,120)NOPART,NTH,NTHBEG,DTMIN,DTMAX
  120 FORMAT('0   DATA SKIPPING CONTROLS'/'  NOPART=',L2,'   NTH=',I5,
     1 ' NTHBEG=',I4,'   DTMIN,DTMAX=',1P,2E15.6)
      IF(NCHKS.GT.0) GOTO 140
      WRITE(PRINT,130)
  130 FORMAT('0   NO INPUT SERIES SELECTION CRITERIA')
      GOTO 1000
  140 WRITE(PRINT,150) NCHKS
  150 FORMAT('0   SERIES SELECTION CRITERIA  (NCHKS=',I3,')' )
      ICRIT=1
      IF(RSITE(1).NE.BLANK) GOTO 180
C        NONE OF THIS CATEGORY
  160 WRITE(PRINT,170) CRITER(ICRIT)
  170 FORMAT('  NO VALUES INPUT FOR ',A8)
      GOTO (220,260,300,340,380,420,460,500,540,1000), ICRIT
  180 I1=NCHKS+1
      DO 190 I=1,NCHKS
      I1=I1-1
      IF(RSITE(I1).NE.BLANK) GOTO 200
  190 CONTINUE
  200 WRITE(PRINT,210) CRITER(ICRIT),(RSITE(I),I=1,I1)
  210 FORMAT(2X,A6,'=',(T10,10(2X,A8,2X)))
  220 ICRIT=2
      IF(SSITE(1).EQ.BLANK) GOTO 160
      I1=NCHKS+1
      DO 230 I=1,NCHKS
      I1=I1-1
      IF(SSITE(I1).NE.BLANK) GOTO 240
  230 CONTINUE
  240 WRITE(PRINT,210) CRITER(ICRIT),(SSITE(I),I=1,I1)
  260 ICRIT=3
      IF(SERIES(1).EQ.BLNK4) GOTO 160
      I1=NCHKS+1
      DO 270 I=1,NCHKS
      I1=I1-1
      IF(SERIES(I1).NE.BLNK4) GOTO 280
  270 CONTINUE
  280 WRITE(PRINT,290) CRITER(ICRIT),(SERIES(I),I=1,I1)
  290 FORMAT(2X,A6,'=',(T10,20(2X,A4)))
  300 ICRIT=4
      IF(NSEQ(1).EQ.-999) GOTO 160
      I1=NCHKS+1
      DO 310 I=1,NCHKS
      I1=I1-1
      IF(NSEQ(I1).NE.-999) GOTO 320
  310 CONTINUE
  320 WRITE(PRINT,330) CRITER(ICRIT),(NSEQ(I),I=1,I1)
  330 FORMAT(2X,A6,'=',(T10,20I6))
  340 ICRIT=5
      IF(NCODF(1).EQ.-999) GOTO 160
      I1=NCHKS+1
      DO 350 I=1,NCHKS
      I1=I1-1
      IF(NCODF(I1).NE.-999) GOTO 360
  350 CONTINUE
  360 WRITE(PRINT,330) CRITER(ICRIT),(NCODF(I),I=1,I1)
  380 ICRIT=6
      IF(NPLNT(1).EQ.-999) GOTO 160
      I1=NCHKS+1
      DO 390 I=1,NCHKS
      I1=I1-1
      IF(NPLNT(I1).NE.-999) GOTO 400
  390 CONTINUE
  400 WRITE(PRINT,330) CRITER(ICRIT),(NPLNT(I),I=1,I1)
  420 ICRIT=7
      IF(NCENTB(1).EQ.-999) GOTO 160
      I1=NCHKS+1
      DO 430 I=1,NCHKS
      I1=I1-1
      IF(NCENTB(I1).NE.-999) GOTO 440
  430 CONTINUE
  440 WRITE(PRINT,330) CRITER(ICRIT),(NCENTB(I),I=1,I1)
  460 ICRIT=8
      IF(SPOT(1).EQ.BLNK4) GOTO 160
      I1=NCHKS+1
      DO 470 I=1,NCHKS
      I1=I1-1
      IF(SPOT(I1).NE.BLNK4) GOTO 480
  470 CONTINUE
  480 WRITE(PRINT,290) CRITER(ICRIT),(SPOT(I),I=1,I1)
  500 ICRIT=9
      IF(SPOT2(1).EQ.BLNK4) GOTO 160
      I1=NCHKS+1
      DO 510 I=1,NCHKS
      I1=I1-1
      IF(SPOT2(I1).NE.BLNK4) GOTO 520
  510 CONTINUE
  520 WRITE(PRINT,290) CRITER(ICRIT),(SPOT2(I),I=1,I1)
  540 ICRIT=10
      IF(FREQ(1).EQ.-999._10) GOTO 160
      I1=NCHKS+1
      DO I=1,NCHKS
         I1=I1-1
         IF(FREQ(I1).NE.-999._10) GOTO 560
      END DO
  560 WRITE(PRINT,570) CRITER(ICRIT),(FREQ(I),I=1,I1)
  570 FORMAT(2X,A6,'=',(T10,1P10E12.5))
 1000 CONTINUE
      WRITE(PRINT,1010) (STSTP(J),DYS(J),ISTSTP(J),(HMS(I,J),I=1,3),
     1  J=1,4)
 1010 FORMAT('0   INDEPENDENT VARIABLE CONTROLS'/ (2X,A6,'=',1PD23.15,
     1 3X,A2,'(HR,MIN,SEC)=',0P,3F15.7) )
      WRITE(PRINT,1020) TSPAN,NSPAN,NPOINT,NSERIE,AVGSPN,PADZ,PERIV,
     1 (STRTVL(I),I=1,NSPAN)
 1020 FORMAT('  TSPAN=',F12.4,'   NSPAN=',I5,'   NPOINT=',I6,
     1 '   NSERIE=',I6,'   AVGSPN=',L2,'   PADZ=',L2,'   PERIV=',1P,
     2 2D18.10/ '  STRTVL=',5D23.14/(9X,5D23.14))
      WRITE(PRINT,1030) NIV,MINIV,MAXIV,TLAT,NORMIV
 1030 FORMAT('  NIV=',I5,'  MIN,MAXIV=',1P,2E12.4,7X,'TLAT=',2E12.4,
     1 '   NORMIV=',L2)
      WRITE(PRINT,1040) KTH,FILTER,TWOOBS
 1040 FORMAT('  KTH=',I5,'   FILTER=',I2,'   TWOOBS=',L2)
      WRITE(PRINT,1100) NCODE,NDV,NAMES,FNAME,LNAME,CUT,NORMDV,PERDV
 1100 FORMAT('0   DEPENDENT VARIABLE CONTROLS'/
     1 '  NCODE=',I3,'   NDV=',I5,'   NAMES=',L3,'   FNAME,LNAME=',
     2 2(1X,A8),'   CUT=',1P,2E15.6,'   NORMDV=',L2/'  PERDV=',2D18.10)
      WRITE(PRINT,1120) NBOX,SZBOX
 1120 FORMAT('0   HISTOGRAM CONTROLS'/ '  NBOX=',I4,
     1 '   SZBOX=',1PE15.7)
      WRITE(PRINT,1150) SORT,BGAP,BSPAN,BSHR,BSMIN,BSSEC,BURSTL
 1150 FORMAT('0   SORT CONTROLS'/ '  SORT=',L4,'  BGAP=',F10.4,
     1 '  BSPAN=',F6.2,'  BSHR=',F5.2,'  BSMIN=',
     1 F5.2,'  BSSEC=',F5.2,'  BURSTL=',F10.5)
      WRITE(PRINT,1200)
 1200 FORMAT('0   PLOTTING CONTROLS')
      WRITE(PRINT,1210) IPLOT,ONEAXE,OUTIN,ZLINE,SYMBEL,DOTPLT,NDEC,
     1 TMODE,FACT
 1210 FORMAT('  IPLOT=',I3,'   ONEAXE=',L2,'   OUTIN=',L3,'   ZLINE=',
     1 L3,'   SYMBEL=',L2,'   DOTPLT=',L2,'   NDEC=',I4,'   TMODE=',I3,
     2 '   FACT=',1PE15.7)
      WRITE(PRINT,1220) SYMOPT,PLTSYM,HIPREC,XDEN,JDREF,INDEX
 1220 FORMAT('  SYMOPT=',L2,'  PLTSYM=',I3,'   HIPREC=',L2,
     1 '   XDEN=',F10.4,'   JDREF=',I9 / '  INDEX=',20I4)
      WRITE(PRINT,1225) EX,XLEN, EX,XFRAME, EX,XMRG, EX2(2),EXTRMX,
     1                 WHY,YLEN,WHY,YFRAME,WHY,YMRG,WHY2(2),EXTRMY
 1225 FORMAT(A3,'LEN=',0PF10.4,A3,'FRAME=',F10.4,A3,'MRG=',F10.4,
     1 '  EXTRM',A1,'=',1P,3D15.7)
      WRITE(PRINT,1230) PROBNO,PROGNO,ASPECT
 1230 FORMAT('  PROBNO=',A8,'   PROGNO=',A8,'   ASPECT=',L2)
      RETURN
      END
