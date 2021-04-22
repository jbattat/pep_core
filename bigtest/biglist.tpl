1INPUT STREAM TO THE PLANETARY EPHEMERIS PROGRAM (PEP)                                                                   PAGE    1
0PEP VERSION= 20210302 PEP.PEPLOAD.PEP790                          
0TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS                 TITLE CARD
  &NMLST1                                                                         SOLAR SYSTEM PARAMETERS &NMLST1 
  EXTPRC=0,                 $ USE HARDWARE EXTENDED PRECISION                    
  ICT(1)=1,                 $ ONE ITERATION                                      
  ICT(3)=1,IOBS1=31,        $ WRITE DUMMIES ON OBSLIB                            
  ICT(15)=-1,               $ NO SOLUTION                                        
  ICT(20)=-4,               $ PRINT PARTIALS                                     
  ICT(28)=1                 $ CALC. NUT-PREC AT BOTH SEND AND RECV TIME          
  ICT(27)=-1,IBODY=90,NBODY=1,PRMTER(81)=1.D0 $USE S.S.B.C                       
  ICT(34)=1,                $ ANNUAL BUT NOT SEMI-ANNUAL TERMS IN A1-UT1         
  ICT(37)=-1,               $ DO NOT SKIP DUMOBS BELOW HORIZON                   
  ICT(39)=1,                $ USE ANY AVAILABLE INTEGRATIONS                     
  ICT(40)=1,                $ USE SBFN FOR EMBARY INTEGRATION                    
  ICT(43)=1,                $ USE NEW CODE                                       
  JCT(2)=11,1, JCAL(9)=22   $ PROPCO WITH EXTRA PRINT, DO SOLAR PLASMA           
  JCT(27)=1,                $ * COMMANDS                                         
  JCT(33)=93,               $ UT1 AND WOBBLE FROM FT93 AND FT94                  
  JCT(54)=1,                $ USE MOR ROUTINES                                   
  MASS(11)=6.7D-10,         $ CERES MASS                                         
  GMVARY=1D-35, SUNHAR=1D-8,                                                     
  PRMTER(47)= 0, 23.4433,   $ A-BELT ASC. & INCL.                                
  PRMTER(49)= 2.9, 1D-9,    $ A-BELT DIST & MASS                                 
  PRMTER(60)=10,            $ ELECTRON DENSITY AT 1 AU                           
  MDSTSC=0,                 $ GET STANDARD MOON UNITS                            
  LPRM=3,4,5,10,31,32,33,40,41,42,49,50,51,60                                    
  LOUT=6,                   $ DISPLAY PARAMETER NAMES                            
0*OBJECT EMBARY                                                                  *****
  COND= .9958267665339834, .09780299873377071, .04240314161324497,                INDIVIDUAL BODY PARAMETERS &NMLST2  
   -.002111811083818306, .01563879101388465, .006781326667821355,                
  ICND=-5,                                                                       
  CON1(1)=2442320.5,                                                             
  JD1=2442304, JD0=2442321, JD2=2442347,                                         
  K=1, 3*1,3*-1, 3,4,10,31,32,33,40,41,42, $ PARTIALS                            
  K(31)= -1,-1,3*1,4*-1,1,  $ INCLUDE MARS, JUP, MOON                            
  K(61)=3*1, K(70)=3*1,     $ INCLUDE RELFCT,GMVARY,SUNHAR,POE,BETA,GAMMA        
  K(79)=2*1,                $ INCLUDE A-BELT                                     
  K(82)=0,                  $ REACTION ON J2 BY SUN                              
  K(87)=2, INT=2,           $ INTERVALS                                          
  K(88)=2,6,                $ ADAMS-MOULTON METHOD, 7 TERMS                      
  K(91)=-3,-6, EPS(3)=1E-9  $ STARTING INTERVALS                                 
  K(98)=4,0,-1,             $ PRINT + TAPE; ORDINARY EQNS OF MOTION              
  L=3*1,3*0, 24,                                                                 
0*OBJECT EARTH-ROTATION                                                          *****
  NAME=' EROTAT ', ITAPE=0,                                                       INDIVIDUAL BODY PARAMETERS &NMLST2  
  CON(1)=6*0,               $ ZERO ALL AD HOC ROTATION ANGLES AND RATES          
  CON(10)=2442048.0         $ EPOCH FOR DRIFT IN A1-UT1 IS 0 JAN 1974            
  CON(22)=5026.75           $ NOMINAL PRECESSION CONSTANT                        
 $CON(23)=84404.84          $ NOMINAL OBLIQUITY AT 1950.0                        
  CON(24)=9.21              $ NOMINAL NUTATION CONSTANT                          
0*OBJECT MOON                                                                    *****
  CON(1)=1738.09,           $ MEAN RADIUS FOR MOON HARMONICS (JPL LLB-5)          INDIVIDUAL BODY PARAMETERS &NMLST2  
  CON(16)=8.93645681309210D-2,                                                   
  K=1, 1,3*-1,2*1, -16,-20,3,10,32,303, 1031,2,3  $ INTEGRATED PARTIALS          
  K(40)=1,K(61)=1,1,                                                             
  K(81)=4,3,3,1,1,                                                               
  K(87)=-2, INT=-1,         $ INTEGRATION STEP SIZE                              
  K(88)=2,6,                $ A-M INTEGRATOR, 7 TERMS                            
  K(91)=-4,-6,              $ STARTING INTERVALS                                 
  K(98)= 10,0,-1,           $ PRINT+TAPE, USUAL EQNS OF MOTION                   
  COND= .0026232835202694, -.0005304792874056, .0000239484232369,                
        .0000854824499752,  .0005195905088138, .0002166804074577,                
  ICND=-5,                  $ CONVERT TO ELLIPTIC                                
  JD0=2442321,JD1=2442321,JD2=2442332,                                           
  NZONE=3, NTESS=3,                                                              
  J2=2.03822000000000D-4, J3=8.78597916047804D-6,                                
  C2=0.,  3.57863661189565D-5,                                                   
  C3=3.02761364595177D-5, 1.42242064346298D-5, 1.35361385680119D-5,              
  S3=5.97661136879657D-6, 4.72438898591151D-6, -2.97262223597402D-6,             
  LJ2=1, LJ3=1,                                                                  
  LC3=1,1,                                                                       
  L=1,3*0,2*1, 16,20    $ I.C.'S + TIDAL FRICTION + ETA DELTA                    
  EPS(3)=1.0D-10,                                                                
0*OBJECT MOON-ROTATION                                                           *****
  NAME=' MROTAT ', ITAPE=21,                                                      INDIVIDUAL BODY PARAMETERS &NMLST2  
  K=1, 2*1,4*-1, -3,-6,-7,1031,2,3,1041,3,1,3,2,  $ PARTIALS                     
  K(31)=2*-1,1,6*-1,1,  $ INCLUDE ONLY EARTH AND MOON                            
  K(61)=-1,             $ NO RELATIVITY                                          
  K(81)=1,              $ SUN INCLUDED IN EQNS OF MOTION AND PARTIALS            
  K(82)=0,              $ F-F INTERACTION IN EOM                                 
  K(83)=0,              $ ELAS. & DISS. IN EQS. OF MOTION                        
  K(86)=303,            $ USE 3RD ORDER ZONAL AND TESS. HARMS IN PARTIALS        
  K(87)=-2,INT=-1,      $ INTEGRATION STEP SIZE                                  
  K(88)=2,6,            $ ADAMS-MOULTON INTEGRATION, 7 TERMS                     
  K(91)=-4,-6,          $ STARTING INTERVAL 1/128 DAY, MIN 1/256                 
  K(98)=10,0,           $ PRINTOUT OF PARTIALS & TAPE                            
  K(100)=-1,            $ INDICATE EULER ANGLE INTEGRATION                       
  CON(1)=1738.09,       $ MEAN RADIUS                                            
  CON(2)=0,             $ USE DERIVED VALUE OF ALPHA                             
  CON(3)=6.31716849185511D-4, CON(4)=2.27809893710516D-4,                        
  CON(5)=2.6920298D-2                                                            
  CON(6)=2.41511993653923D-2, 4.68318905574876D-3, $ K2,K2*DELTA                 
  L=2*1,4*0,  3,6,7,                                                             
  COND=.0640689056244628,  .4174732965826930, 442.772215182822606,               
      -.0001116475592635, -.0000391413652830,    .2300859910558332,              
  JD0=2442321,JD1=2442321,JD2=2442332                                            
0*OBJECT MARS                                                                    *****
  NCENTR=-1,                                                                      INDIVIDUAL BODY PARAMETERS &NMLST2  
  COND=-1.582279061283973, -.3680747075174319, -.1265083724103855,               
      .003872124750611298, -.01123881652406281,-.005262720398365393              
  ICND=-5,                                                                       
  CON(1)=3392.459,      $ RADIUS                                                 
  CON(6)=151.446, 1.025956, 52.6951, 317.3116 $ PHASE, PERIOD, DEC, RA           
  CON1(1)=2443509.5     $ J1978.0                                                
  JD1=2442321, JD0=2442321, JD2=2442336,                                         
  K=1, 3*-1,3*1, 3,5,11,31,40,41,49,50,506,                                      
  K(31)= -1,-1,3*1,5*-1, 1,                                                      
  K(61)=1, K(70)=2*1, K(79)=2*1,                                                 
  K(87)=1, INT=1,       $ INTERVALS                                              
  K(91)=-10,-14,        $ STARTING INTERVALS                                     
  K(93)=5,              $ TARGET BODY (FOR MEAN ANOMALY)                         
  K(98)=4,0,0,          $ PRINT+TAPE, ENCKE METHOD                               
  L=3*0,3*1, 6,7,8,9, 24                                                         
0*OBJECT JUPITER                                                                 *****
  ITAPE=0                                                                         INDIVIDUAL BODY PARAMETERS &NMLST2  
  L(6)=1,                                                                        
  CON1(1)=2440000.5,                                                             
0*OBJECT 11                                                                      *****
  NAME=' CERES 1', ITAPE=0,                                                       INDIVIDUAL BODY PARAMETERS &NMLST2  
 $ I.C. FROM ASTRONOMICAL PAPERS FOR USE OF AMERICAN EPHEMERIS VOL 16,PART 3     
  CON1(1)=2440000.5D0,                                                           
  COND(1)=2.766333786805837, COND(2)=7.862823976140470D-2,                       
  COND(3)=27.17884370592106, COND(4)=23.40674566424140,                          
  COND(5)=128.9601772696749, COND(6)=58.45329102597486,                          
0*SPOTS                                                                          *****
 OLYM  4    3396.015377   -133.71234       17.81234               1 1 1           SPOT OR STAR CARDS  
0*DLTREAD                                                                        *****
  1 10 HAYS LASR HAYS AP15 1.0E0 1.0E0 1.E-9 1       4.317800000000D+14 0 -5 1000 OVERRIDE ERRORS AND DUMMY OBS. CARDS
 2442326 10 01 30.5      2442327  2 03 00.0       1.0E-08            7200        
                                                                                 
  6 10 NICE ASTG           1.0E0 1.0E0 1.E-9 1       4.317800000000D+14 0 -5 1010
 2442326 10 01 30.5      2442327  2 03 00.0        .01     .1        7200        
                                                                                 
  1  4 AREC RADR AREC OLYM 1.0E0 1.0E0 1.E-9 1       2.380000000000D+09 0 -5 1020
 2442326 20 01 30.5      2442327 12 03 00.0       1.0E-08            7200        
                                                                                 
  3  4 AREC RADR AREC OLYM 1.0E0 1.0E0 1.E-9 1       2.380000000000D+09 0 -5 1030
 2442326 20 01 30.5      2442327  4 03 00.0       1.0E-08            7200        
                                                                                 
  5  4 PARI ASTM           1.0E0 1.0E0 1.E-9 1       4.317800000000D+14 0 -5 1040
 2442327 00 01 30.5      2442327 12 03 00.0        .01     .1        7200        
 **END                                                                           
0**END                                                                           *****
 -------------------------------------------------------------------------------- END OF INPUT STREAM
1PLANETARY EPHEMERIS PROGRAM INPUT DATA   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE    2
0DAY OF YEAR=  61    CLOCK TIME AT EXECUTION= 12H 30M 42.00S
-
-
-

                                  **   **           **   **********    **       **   ************
                                  **   ***          **   ***********   **       **   ************
                                  **   ****         **   **       **   **       **        ** 
                                  **   ** **        **   **       **   **       **        ** 
                                  **   **  **       **   **       **   **       **        ** 
                                  **   **   **      **   **       **   **       **        ** 
                                  **   **    **     **   ***********   **       **        ** 
                                  **   **     **    **   **********    **       **        ** 
                                  **   **      **   **   **            **       **        ** 
                                  **   **       **  **   **            **       **        ** 
                                  **   **        ** **   **            **       **        ** 
                                  **   **         ****   **            **       **        ** 
                                  **   **          ***   **            ***********        ** 
                                  **   **           **   **             *********         ** 
-
-
                                          **           **    **           **   **       **
                                          **           **    ***          **   **      ** 
                                          **           **    ****         **   **     **  
                                          **           **    ** **        **   **    **   
                                          **           **    **  **       **   **   **    
                                          **           **    **   **      **   *** **     
                                          **           **    **    **     **   *****      
                                          **           **    **     **    **   ** **      
                                          **           **    **      **   **   **  **     
                                          **           **    **       **  **   **   **    
                                          **           **    **        ** **   **    **   
                                          **           **    **         ****   **     **  
                                          ***********  **    **          ***   **      ** 
                                          ***********  **    **           **   **       **
- PEP VERSION= 20210302 PEP.PEPLOAD.PEP790                          
-
0SPECIFICATION OF THE UNITS OF MASS, LENGTH AND TIME
0MASS OF SUN                                   = 1.0000000000000
 SQRT(GRAVITATIONAL CONSTANT TIMES MASS OF SUN)= 0.0172020989500 ((AU)**3/2)/(DAY)
 INDEPENDENT TIME VARIABLE OF EPHEMERIDES      = A.1 + 32.15S
1PLANETARY EPHEMERIS PROGRAM INPUT DATA   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE    3
0PERIPHERAL DATA SET NUMBERS NOT BELONGING TO A SPECIFIC BODY OR OBSERVATION LIBRARY
      IN =  5    IOUT =  6    JOUT =  0    KOUT =  0   IPUNCH=  7   IGRAPH= 97   INTERN= 99    IBUF =  4    IMAT =  0    IPERT= 90
    JPERT=  0    KPERT=  0    IVCNT=  0   IPLCON=  1   IOBCON=  2    IOBS =  5   ITHING=  3    IMAT1=  0    IMAT2=  0    JMAT =  0
    IMAT3=  0    IMAT4=  0    ICTAT=  0   LIBHOC=  0
 NUMMT0=  0   IMAT0=
0CONTROL CONSTANTS FOR PROGRAM FLOW AND LEAST-SQUARES ANALYSIS     NPARAM=  44
  ICT( 1)=  1  ICT( 2)=  0  ICT( 3)=  1  ICT( 4)=  0  ICT( 5)= -1  ICT( 6)=  0  ICT( 7)=  0  ICT( 8)=  0  ICT( 9)=  0  ICT(10)= -2
  ICT(11)=  0  ICT(12)=  0  ICT(13)=  0  ICT(14)=  0  ICT(15)= -1  ICT(16)=  0  ICT(17)=  0  ICT(18)=  0  ICT(19)=  0  ICT(20)= -4
  ICT(21)=  0  ICT(22)=  0  ICT(23)=  0  ICT(24)=  0  ICT(25)=  0  ICT(26)=  1  ICT(27)= -1  ICT(28)=  1  ICT(29)=  0  ICT(30)=  0
  ICT(31)=  0  ICT(32)=  0  ICT(33)=  0  ICT(34)=  1  ICT(35)=  0  ICT(36)=  0  ICT(37)= -1  ICT(38)=  0  ICT(39)=  1  ICT(40)=  1
  ICT(41)=  0  ICT(42)=  0  ICT(43)=  1  ICT(44)=  0  ICT(45)=  0  ICT(46)= -1  ICT(47)=  0  ICT(48)=  0  ICT(49)=  0  ICT(50)=  0
  ICT(51)=  0  ICT(52)=  0  ICT(53)=  0  ICT(54)=  0  ICT(55)=  0  ICT(56)=  0  ICT(57)=  0  ICT(58)=  0  ICT(59)=  0  ICT(60)=  0
  ICT(61)=  0  ICT(62)=  0  ICT(63)=  0  ICT(64)=  0  ICT(65)=  0  ICT(66)=  0  ICT(67)=  0  ICT(68)=  0  ICT(69)=  0  ICT(70)=  0
  ICT(71)=  0  ICT(72)=  0  ICT(73)=  0  ICT(74)=  0  ICT(75)=  0  ICT(76)=  0  ICT(77)=  0  ICT(78)=  0  ICT(79)=  0  ICT(80)=  0
  EPS( 1)= 1.00000E+03  EPS( 2)= 1.00000E+02  EPS( 3)= 1.00000E+01  EPS( 4)= 1.00000E+01  EPS( 5)= 1.00000E+01  EPS( 6)= 1.00000E+01
  EPS( 7)= 1.00000E+02  EPS( 8)= 1.00000E+01  EPS( 9)= 1.00000E-01  EPS(10)= 1.00000E-02  EPS(11)= 1.00000E-03  EPS(12)= 1.00000E-06
  EPS(13)= 1.00000E-32  EPS(14)= 1.00000E+00  EPS(15)= 0.00000E+00  EPS(16)= 0.00000E+00  EPS(17)= 0.00000E+00  EPS(18)= 0.00000E+00
  EPS(19)= 0.00000E+00  EPS(20)= 0.00000E+00  EPS(21)= 0.00000E+00  EPS(22)= 0.00000E+00  EPS(23)= 0.00000E+00  EPS(24)= 0.00000E+00
  EPS(25)= 0.00000E+00  EPS(26)= 0.00000E+00  EPS(27)= 0.00000E+00  EPS(28)= 0.00000E+00  EPS(29)= 0.00000E+00  EPS(30)= 0.00000E+00
0CONTROL CONSTANTS FOR LEAST SQUARES PARAMETER ADJUSTMENT AND THE CORRESPONDING PARAMETER VALUES
 LPRM( 1)=  3 LPRM( 2)=  4 LPRM( 3)=  5 LPRM( 4)= 10 LPRM( 5)= 31 LPRM( 6)= 32 LPRM( 7)= 33 LPRM( 8)= 40 LPRM( 9)= 41 LPRM(10)= 42
 LPRM(11)= 49 LPRM(12)= 50 LPRM(13)= 51 LPRM(14)= 60 LPRM(15)=  0 LPRM(16)=  0 LPRM(17)=  0 LPRM(18)=  0 LPRM(19)=  0 LPRM(20)=  0
 LPRM(21)=  0 LPRM(22)=  0 LPRM(23)=  0 LPRM(24)=  0 LPRM(25)=  0 LPRM(26)=  0 LPRM(27)=  0 LPRM(28)=  0 LPRM(29)=  0 LPRM(30)=  0
 LPRM(31)=  0 LPRM(32)=  0 LPRM(33)=  0 LPRM(34)=  0 LPRM(35)=  0 LPRM(36)=  0 LPRM(37)=  0 LPRM(38)=  0 LPRM(39)=  0 LPRM(40)=  0
 LPRM(41)=  0 LPRM(42)=  0 LPRM(43)=  0 LPRM(44)=  0 LPRM(45)=  0 LPRM(46)=  0 LPRM(47)=  0 LPRM(48)=  0 LPRM(49)=  0 LPRM(50)=  0
 LPRM(51)=  0 LPRM(52)=  0 LPRM(53)=  0 LPRM(54)=  0 LPRM(55)=  0 LPRM(56)=  0 LPRM(57)=  0 LPRM(58)=  0 LPRM(59)=  0 LPRM(60)=  0
 LPRM(61)=  0 LPRM(62)=  0 LPRM(63)=  0 LPRM(64)=  0 LPRM(65)=  0 LPRM(66)=  0 LPRM(67)=  0 LPRM(68)=  0 LPRM(69)=  0 LPRM(70)=  0
 LPRM(71)=  0 LPRM(72)=  0 LPRM(73)=  0 LPRM(74)=  0 LPRM(75)=  0 LPRM(76)=  0 LPRM(77)=  0 LPRM(78)=  0 LPRM(79)=  0 LPRM(80)=  0
 LPRM(81)=  0 LPRM(82)=  0 LPRM(83)=  0 LPRM(84)=  0 LPRM(85)=  0 LPRM(86)=  0 LPRM(87)=  0 LPRM(88)=  0 LPRM(89)=  0 LPRM(90)=  0
 LPRM(91)=  0 LPRM(92)=  0 LPRM(93)=  0 LPRM(94)=  0 LPRM(95)=  0 LPRM(96)=  0 LPRM(97)=  0 LPRM(98)=  0 LPRM(99)=  0 LPRM(**)=  0
  MASS( 1)= 1.657848020429993D-07  MASS( 2)= 2.447848585877872D-06  MASS( 3)= 3.040436898620584D-06  MASS( 4)= 3.227159776680543D-07
  MASS( 5)= 9.547534692876814D-04  MASS( 6)= 2.857796067672611D-04  MASS( 7)= 4.361098996947231D-05  MASS( 8)= 5.192107995846314D-05
  MASS( 9)= 2.500000000000000D-07  MASS(10)= 1.215052064980984D-02  MASS(11)= 6.700000000000000D-10  MASS(12)= 0.000000000000000D+00
  MASS(13)= 0.000000000000000D+00  MASS(14)= 0.000000000000000D+00  MASS(15)= 0.000000000000000D+00  MASS(16)= 0.000000000000000D+00
  MASS(17)= 0.000000000000000D+00  MASS(18)= 0.000000000000000D+00  MASS(19)= 0.000000000000000D+00  MASS(20)= 0.000000000000000D+00
  MASS(21)= 0.000000000000000D+00  MASS(22)= 0.000000000000000D+00  MASS(23)= 0.000000000000000D+00  MASS(24)= 0.000000000000000D+00
  MASS(25)= 0.000000000000000D+00  MASS(26)= 0.000000000000000D+00  MASS(27)= 0.000000000000000D+00  MASS(28)= 0.000000000000000D+00
  MASS(29)= 0.000000000000000D+00  MASS(30)= 0.000000000000000D+00 RELFCT31 = 1.000000000000000D+00 GMVARY32 = 1.000000000000000D-35
 SUNHAR33 = 1.000000000000000D-08 PRMTR(34)= 0.000000000000000D+00 PRMTR(35)= 0.000000000000000D+00 PRMTR(36)= 0.000000000000000D+00
 PRMTR(37)= 0.000000000000000D+00 PRMTR(38)= 0.000000000000000D+00 SUNPOE39 =-7.800000000000000D-04 PRMTR(40)= 0.000000000000000D+00
 BETA  41 = 1.000000000000000D+00 GAMMA 42 = 1.000000000000000D+00 BETA' 43 = 1.000000000000000D+00 GAMMA'44 = 1.000000000000000D+00
 PRMTR(45)= 0.000000000000000D+00 PRMTR(46)= 0.000000000000000D+00 PRMTR(47)= 0.000000000000000D+00 PRMTR(48)= 2.344330000000000D+01
  ASBA 49 = 2.900000000000000D+00 ASBMAS50 = 1.000000000000000D-09 AULTSC51 = 4.990047800000000D+02 LTVARY52 = 0.000000000000000D+00
 RELDEL53 = 1.000000000000000D+00 RELDOP54 = 0.000000000000000D+00 PRMTR(55)= 0.000000000000000D+00 PRMTR(56)= 0.000000000000000D+00
 PRMTR(57)= 0.000000000000000D+00 PRMTR(58)= 0.000000000000000D+00 PRMTR(59)= 0.000000000000000D+00 PRMTR(60)= 1.000000000000000D+01
 PRMTR(61)= 0.000000000000000D+00 PRMTR(62)= 1.000000000000000D+00 PRMTR(63)= 1.000000000000000D+00 PRMTR(64)= 0.000000000000000D+00
 PRMTR(65)= 0.000000000000000D+00 PRMTR(66)= 0.000000000000000D+00 PRMTR(67)= 0.000000000000000D+00 PRMTR(68)= 0.000000000000000D+00
 PRMTR(69)= 0.000000000000000D+00 PRMTR(70)= 0.000000000000000D+00 PRMTR(71)= 0.000000000000000D+00 CTVARY72 = 0.000000000000000D+00
 PRMTR(73)= 0.000000000000000D+00 PRMTR(74)= 0.000000000000000D+00 PRMTR(75)= 0.000000000000000D+00 PRMTR(76)= 0.000000000000000D+00
 PRMTR(77)= 0.000000000000000D+00 PRMTR(78)= 0.000000000000000D+00 PRMTR(79)= 0.000000000000000D+00 PRMTR(80)= 0.000000000000000D+00
 PRMTR(81)= 1.000000000000000D+00 PRMTR(82)= 0.000000000000000D+00 PRMTR(83)= 0.000000000000000D+00 PRMTR(84)= 0.000000000000000D+00
 PRMTR(85)= 0.000000000000000D+00 PRMTR(86)= 0.000000000000000D+00 PRMTR(87)= 0.000000000000000D+00 PRMTR(88)= 0.000000000000000D+00
 PRTMR(89)= 0.000000000000000D+00 PRMTR(90)= 0.000000000000000D+00 ECINC 91 = 2.344578706750000D+01 SEQINC92 = 7.250000000000000D+00
 SEQASC93 = 7.506250000000000D+01 SUNRAD94 = 6.960000000000000D+05 PRMTR(95)= 6.960000000000000D+05 PRMTR(96)= 0.000000000000000D+00
 PRMTR(97)= 2.440000500000000D+06 MDSTAU98 = 4.263529034000000D-05 MDSTSC99 = 0.000000000000000D+00 LTVEL100 = 2.997924580000000D+05
1PLANETARY EPHEMERIS PROGRAM INPUT DATA   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE    4
0INVERSE MASSES CALCULATED FROM PREVIOUS PAGE
 IMASS( 1)= 6.031916000000000D+06 IMASS( 2)= 4.085220000000000D+05 IMASS( 3)= 3.289001000000000D+05 IMASS( 4)= 3.098700000000000D+06
 IMASS( 5)= 1.047390800000000D+03 IMASS( 6)= 3.499200000000000D+03 IMASS( 7)= 2.293000000000000D+04 IMASS( 8)= 1.926000000000000D+04
 IMASS( 9)= 4.000000000000000D+06 IMASS(10)= 8.230100000000000D+01 IMASS(11)= 1.492537313432836D+09 IMASS(12)= 0.000000000000000D+00
 IMASS(13)= 0.000000000000000D+00 IMASS(14)= 0.000000000000000D+00 IMASS(15)= 0.000000000000000D+00 IMASS(16)= 0.000000000000000D+00
 IMASS(17)= 0.000000000000000D+00 IMASS(18)= 0.000000000000000D+00 IMASS(19)= 0.000000000000000D+00 IMASS(20)= 0.000000000000000D+00
 IMASS(21)= 0.000000000000000D+00 IMASS(22)= 0.000000000000000D+00 IMASS(23)= 0.000000000000000D+00 IMASS(24)= 0.000000000000000D+00
 IMASS(25)= 0.000000000000000D+00 IMASS(26)= 0.000000000000000D+00 IMASS(27)= 0.000000000000000D+00 IMASS(28)= 0.000000000000000D+00
 IMASS(29)= 0.000000000000000D+00 IMASS(30)= 0.000000000000000D+00
0ALTERNATE PERTURBING PLANET DATA SETS, BUFFERS, ACCURACY CONSTANTS
   IPERT0= 90   IPERT1=  0   IPERT2=  0   JPERT0=  0   JPERT1=  0   JPERT2=  0   KPERT0=  0   KPERT1=  0   KPERT2=  0
    IBUF1=  0    IBUF2=  0    IBUF3=  0    IBUF4=  0    IBUF5=  0    IBUF6=  0    IBUF7=  0    IBUF8=  0    IBUF9=  0   IBUF10=  0
 EPSA( 1)= 0.00000E+00 EPSA( 2)= 0.00000E+00 EPSA( 3)= 0.00000E+00 EPSA( 4)= 0.00000E+00 EPSA( 5)= 0.00000E+00 EPSA( 6)= 0.00000E+00
 EPSA( 7)= 0.00000E+00 EPSA( 8)= 0.00000E+00 EPSA( 9)= 0.00000E+00 EPSA(10)= 0.00000E+00 EPSA(11)= 0.00000E+00 EPSA(12)= 0.00000E+00
 EPSA(13)= 0.00000E+00 EPSA(14)= 0.00000E+00 EPSA(15)= 0.00000E+00 EPSA(16)= 0.00000E+00 EPSA(17)= 0.00000E+00 EPSA(18)= 0.00000E+00
 EPSA(19)= 0.00000E+00 EPSA(20)= 0.00000E+00 EPSA(21)= 0.00000E+00 EPSA(22)= 0.00000E+00 EPSA(23)= 0.00000E+00 EPSA(24)= 0.00000E+00
 EPSA(25)= 0.00000E+00 EPSA(26)= 0.00000E+00 EPSA(27)= 0.00000E+00 EPSA(28)= 0.00000E+00 EPSA(29)= 0.00000E+00 EPSA(30)= 5.00000E-05
0ADDITIONAL PERIPHERAL DATA SET NUMBERS NOT BELONGING TO A SPECIFIC BODY OR OBSERVATION LIBRARY + SWITCHES
   JPUNCH=  0   KPUNCH=  0    LOUT =  6    MOUT =  0    NOUT =  0    IENG =  0    JENG =  0    KENG =  0
   TYPOUT=  0   EXTPRC=  0   NOPRNT=  0    ISEQ =  0
0ADDITIONAL CONTROL CONSTANTS FOR PROGRAM FLOW AND LEAST-SQUARES ANALYSIS
  JCT( 1)=  0  JCT( 2)= 11  JCT( 3)=  1  JCT( 4)=  0  JCT( 5)=  0  JCT( 6)=  0  JCT( 7)=  0  JCT( 8)=  0  JCT( 9)=  0  JCT(10)=  0
  JCT(11)=  0  JCT(12)=  0  JCT(13)=  0  JCT(14)=  0  JCT(15)=  0  JCT(16)=  0  JCT(17)=  0  JCT(18)=  0  JCT(19)=  0  JCT(20)=  0
  JCT(21)= -1  JCT(22)=  0  JCT(23)=  0  JCT(24)=  0  JCT(25)=  0  JCT(26)=  0  JCT(27)=  1  JCT(28)=  0  JCT(29)=  0  JCT(30)=  0
  JCT(31)=  0  JCT(32)=  0  JCT(33)= 93  JCT(34)=  0  JCT(35)=  0  JCT(36)=  0  JCT(37)=  0  JCT(38)=  0  JCT(39)=  0  JCT(40)=  0
  JCT(41)=  0  JCT(42)=  0  JCT(43)=  0  JCT(44)=  0  JCT(45)=  0  JCT(46)=  0  JCT(47)=  0  JCT(48)=  0  JCT(49)=  0  JCT(50)=  0
  JCT(51)=  0  JCT(52)=  0  JCT(53)=  0  JCT(54)=  1  JCT(55)=  0  JCT(56)=  0  JCT(57)=  0  JCT(58)=  0  JCT(59)=  0  JCT(60)=  0
  JCT(61)= -1  JCT(62)= -1  JCT(63)=  0  JCT(64)=  0  JCT(65)=  0  JCT(66)=  0  JCT(67)=  0  JCT(68)=  0  JCT(69)=  0  JCT(70)=  0
  JCT(71)=  0  JCT(72)=  0  JCT(73)=  0  JCT(74)=  0  JCT(75)=  0  JCT(76)=  0  JCT(77)=  0  JCT(78)=  0  JCT(79)=  0  JCT(80)=  0
  JCT(81)=  0  JCT(82)=  0  JCT(83)=  0  JCT(84)=  0  JCT(85)=  0  JCT(86)=  0  JCT(87)=  0  JCT(88)=  0  JCT(89)=  0  JCT(90)=  0
  JCT(91)=  0  JCT(92)=  0  JCT(93)=  0  JCT(94)=  0  JCT(95)=  0  JCT(96)=  0  JCT(97)=  0  JCT(98)=  0  JCT(99)=  0  JCT(**)=  0
0OBSERVATION LIBRARY PERIPHERAL DATA SET NUMBERS (NUMOBT= 1)
  IOBS0( 1)= 0  IOBS1( 1)=31  IOBS2( 1)= 0
0WEIGHTS APPLIED TO NORMAL EQUATIONS FORMED FROM OBSERVATION DATA SETS OR RESTORED FROM SAVED NORMAL EQUATIONS
 WGTOBS(1)= 1.000000000000000D+00
0INPUT CONTROL CONSTANTS FOR N-BODY INTEGRATION
   NBODY =  1   IBODY = 90   JDBDY1=       0   JDBDY0=       0   JDBDY2=       0   JVLBDY=     0   EPSBDY= 1.00000E-16   INTBDY=  2
 KBDY( 1)= -1 KBDY( 2)= -1 KBDY( 3)=  0 KBDY( 4)= -1 KBDY( 5)= -1 KBDY( 6)= -1 KBDY( 7)= -1 KBDY( 8)= -1 KBDY( 9)= -1 KBDY(10)= -1
 KBDY(11)= -1 KBDY(12)= -1 KBDY(13)= -1 KBDY(14)= -1 KBDY(15)= -1 KBDY(16)= -1 KBDY(17)= -1 KBDY(18)= -1 KBDY(19)= -1 KBDY(20)= -1
 KBDY(21)=  0 KBDY(22)= -1 KBDY(23)= -1 KBDY(24)= -1 KBDY(25)= -1 KBDY(26)= -1 KBDY(27)= -1 KBDY(28)=  3 KBDY(29)= 11 KBDY(30)=  6
 KBDY(31)=-12 KBDY(32)=-14 KBDY(33)= -1 KBDY(34)= -1 KBDY(35)= -1 KBDY(36)=  0 KBDY(37)= -1 KBDY(38)= -1 KBDY(39)=  0 KBDY(40)= -1
 NPLBDY( 1)= 1
1PLANETARY EPHEMERIS PROGRAM INPUT DATA   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE    5
0 EMBARY     NPLNT(-3)=  3    IPLNT= 13         JDEM1= 2442304       JDEM0= 2442321       JDEM2= 2442347     INT =  2   NCENTR=  0
  LEM( 1)=  1  LEM( 2)=  1  LEM( 3)=  1  LEM( 4)=  0  LEM( 5)=  0  LEM( 6)=  0  LEM( 7)= 24  LEM( 8)=  0  LEM( 9)=  0  LEM(10)=  0
  LEM(11)=  0  LEM(12)=  0  LEM(13)=  0  LEM(14)=  0  LEM(15)=  0  LEM(16)=  0  LEM(17)=  0  LEM(18)=  0  LEM(19)=  0  LEM(20)=  0
  LEM(21)=  0  LEM(22)=  0  LEM(23)=  0  LEM(24)=  0  LEM(25)=  0  LEM(26)=  0  LEM(27)=  0  LEM(28)=  0  LEM(29)=  0  LEM(30)=  0
 INITIAL EPOCH (COORD.TIME) 2442320.50000000 IHR= 0 IMIN= 0 SEC= 0.0000 INT1=         0 INT2=  0 FRACT= 0.0000000000000  ICND=  0
         A= 1.000002347957079D+00         E= 1.669041814744669D-02       INC= 2.344258335876651D+01       ASC= 8.391300825564327D-04
       PER= 1.022643964965586D+02      ANOM= 2.657488867161210D+02    RADIUS= 6.378166000000000D+03      FLAT= 3.352329869259135D-03
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 0.000000000000000D+00
   CON( 7)= 0.000000000000000D+00   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00  CON1( 1)= 2.442320500000000D+06  CON1( 2)= 0.000000000000000D+00
  CON1( 3)= 0.000000000000000D+00  CON1( 4)= 0.000000000000000D+00  CON1( 5)= 0.000000000000000D+00  CON1( 6)= 0.000000000000000D+00
  CON1( 7)= 0.000000000000000D+00  CON1( 8)= 0.000000000000000D+00  CON1( 9)= 0.000000000000000D+00  CON1(10)=-4.600000000000000D-10
  CON1(11)= 0.000000000000000D+00  CON1(12)= 0.000000000000000D+00
   K( 1)=   0   K( 2)=   0   K( 3)=   0   K( 4)=   0   K( 5)=   0   K( 6)=   0   K( 7)=   0   K( 8)=   0   K( 9)=   0   K(10)=   0
   K(11)=   0   K(12)=   0   K(13)=   0   K(14)=   0   K(15)=   0   K(16)=   0   K(17)=   0   K(18)=   0   K(19)=   0   K(20)=   0
   K(21)=   0   K(22)=   0   K(23)=   0   K(24)=   0   K(25)=   0   K(26)=   0   K(27)=   0   K(28)=   0   K(29)=   0   K(30)=   0
   K(31)=  -1   K(32)=  -1   K(33)=   1   K(34)=   1   K(35)=   1   K(36)=  -1   K(37)=  -1   K(38)=  -1   K(39)=  -1   K(40)=   1
   K(41)=  -1   K(42)=  -1   K(43)=  -1   K(44)=  -1   K(45)=  -1   K(46)=  -1   K(47)=  -1   K(48)=  -1   K(49)=  -1   K(50)=  -1
   K(51)=  -1   K(52)=  -1   K(53)=  -1   K(54)=  -1   K(55)=  -1   K(56)=  -1   K(57)=  -1   K(58)=  -1   K(59)=  -1   K(60)=  -1
   K(61)=   1   K(62)=   1   K(63)=   1   K(64)=  -1   K(65)=  -1   K(66)=  -1   K(67)=  -1   K(68)=  -1   K(69)=  -1   K(70)=   1
   K(71)=   1   K(72)=   1   K(73)=  -1   K(74)=  -1   K(75)=  -1   K(76)=  -1   K(77)=  -1   K(78)=  -1   K(79)=   1   K(80)=   1
   K(81)=  -1   K(82)=   0   K(83)=  -1   K(84)=  -1   K(85)=  -1   K(86)=  -1   K(87)=   2   K(88)=   2   K(89)=   6   K(90)=   6
   K(91)=  -3   K(92)=  -6   K(93)=   0   K(94)=   0   K(95)=   0   K(96)=   0   K(97)=   0   K(98)=   4   K(99)=   0   K(**)=  -1
 NUMKI= 16  KI= 1 1 1 1-1-1-1   3   4  10  31  32  33  40  41  42                                                                    
  EPS( 1)= 0.00000E+00  EPS( 2)= 0.00000E+00  EPS( 3)= 1.00000E-09  EPS( 4)= 0.00000E+00  EPS( 5)= 0.00000E+00  EPS( 6)= 0.00000E+00
0          ZONAL GRAVITATIONAL POTENTIAL HARMONIC COEFFICIENTS -  EARTH   (NPLNT=  3, NZONE= 11)                                     
       J2= 1.082637000000000D-03       J3=-2.541000000000000D-06       J4=-1.618000000000000D-06       J5=-2.280000000000000D-07
       J6= 5.520000000000000D-07       J7=-3.520000000000000D-07       J8=-2.050000000000000D-07       J9=-1.540000000000000D-07
      J10=-2.370000000000000D-07      J11= 3.120000000000000D-07
      LJ2=  0      LJ3=  0      LJ4=  0      LJ5=  0      LJ6=  0      LJ7=  0      LJ8=  0      LJ9=  0     LJ10=  0     LJ11=  0

0TESSERAL COSINE GRAVITATIONAL POTENTIAL HARMONIC COEFFICIENTS -  EARTH   (NPLNT=  3, NTESS=  5)                                     
    C2(1)= 0.000000000000000D+00    C2(2)= 2.379900000000000D-06    C3(1)= 1.997700000000000D-06    C3(2)= 7.783000000000000D-07
    C3(3)= 4.901100000000000D-07    C4(1)=-5.174800000000000D-07    C4(2)= 3.429600000000000D-07    C4(3)= 1.039000000000000D-06
    C4(4)=-1.051200000000000D-07    C5(1)=-5.366700000000000D-08    C5(2)= 5.986900000000000D-07    C5(3)=-5.842900000000000D-07
    C5(4)=-1.158300000000000D-07    C5(5)= 1.395600000000000D-07
    LC2(1)= 0    LC2(2)= 0    LC3(1)= 0    LC3(2)= 0    LC3(3)= 0    LC4(1)= 0    LC4(2)= 0    LC4(3)= 0    LC4(4)= 0    LC5(1)= 0
    LC5(2)= 0    LC5(3)= 0    LC5(4)= 0    LC5(5)= 0
0TESSERAL  SINE  GRAVITATIONAL POTENTIAL HARMONIC COEFFICIENTS -  EARTH   (NPLNT=  3, NTESS=  5)                                     
    S2(1)= 0.000000000000000D+00    S2(2)=-1.365600000000000D-06    S3(1)= 2.233700000000000D-07    S3(2)=-7.551900000000000D-07
    S3(3)= 1.528300000000000D-06    S4(1)=-4.814000000000000D-07    S4(2)= 6.717400000000000D-07    S4(3)=-1.192300000000000D-07
    S4(4)= 3.566100000000000D-07    S5(1)=-7.997300000000000D-08    S5(2)=-3.991000000000000D-07    S5(3)=-1.633800000000000D-07
    S5(4)=-4.539300000000000D-08    S5(5)=-8.684100000000000D-07
    LS2(1)= 0    LS2(2)= 0    LS3(1)= 0    LS3(2)= 0    LS3(3)= 0    LS4(1)= 0    LS4(2)= 0    LS4(3)= 0    LS4(4)= 0    LS5(1)= 0
    LS5(2)= 0    LS5(3)= 0    LS5(4)= 0    LS5(5)= 0
1PLANETARY EPHEMERIS PROGRAM INPUT DATA   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE    6
0  MOON      NPLNT(-2)= 10    IPLNT= 20         JDMN1= 2442321       JDMN0= 2442321       JDMN2= 2442332     INT = -1   NCENTR=  3
  LMN( 1)=  1  LMN( 2)=  0  LMN( 3)=  0  LMN( 4)=  0  LMN( 5)=  1  LMN( 6)=  1  LMN( 7)= 16  LMN( 8)= 20  LMN( 9)=  0  LMN(10)=  0
  LMN(11)=  0  LMN(12)=  0  LMN(13)=  0  LMN(14)=  0  LMN(15)=  0  LMN(16)=  0  LMN(17)=  0  LMN(18)=  0  LMN(19)=  0  LMN(20)=  0
  LMN(21)=  0  LMN(22)=  0  LMN(23)=  0  LMN(24)=  0  LMN(25)=  0  LMN(26)=  0  LMN(27)=  0  LMN(28)=  0  LMN(29)=  0  LMN(30)=  0
 INITIAL EPOCH (COORD.TIME) 2442320.50000000 IHR= 0 IMIN= 0 SEC= 0.0000 INT1=         0 INT2=  0 FRACT= 0.0000000000000  ICND=  0
         A= 2.584844478828888D-03         E= 4.663456366998220D-02       INC= 2.240083560206691D+01       ASC= 3.473239233192494D+02
       PER= 1.425397834496920D+02      ANOM= 2.222476213123667D+02    RADIUS= 1.738090000000000D+03      FLAT= 0.000000000000000D+00
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 0.000000000000000D+00
   CON( 7)= 0.000000000000000D+00   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 8.936456813092100D-02   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00  CON1( 1)= 0.000000000000000D+00  CON1( 2)= 0.000000000000000D+00
  CON1( 3)= 0.000000000000000D+00  CON1( 4)= 0.000000000000000D+00  CON1( 5)= 0.000000000000000D+00  CON1( 6)= 0.000000000000000D+00
  CON1( 7)= 0.000000000000000D+00  CON1( 8)= 0.000000000000000D+00  CON1( 9)= 0.000000000000000D+00  CON1(10)=-1.900000000000000D-11
  CON1(11)= 0.000000000000000D+00  CON1(12)= 0.000000000000000D+00
   K( 1)=   0   K( 2)=   0   K( 3)=   0   K( 4)=   0   K( 5)=   0   K( 6)=   0   K( 7)=   0   K( 8)=   0   K( 9)=   0   K(10)=   0
   K(11)=   0   K(12)=   0   K(13)=   0   K(14)=   0   K(15)=   0   K(16)=   0   K(17)=   0   K(18)=   0   K(19)=   0   K(20)=   0
   K(21)=   0   K(22)=   0   K(23)=   0   K(24)=   0   K(25)=   0   K(26)=   0   K(27)=   0   K(28)=   0   K(29)=   0   K(30)=   0
   K(31)=   1   K(32)=   1   K(33)=   1   K(34)=   1   K(35)=   1   K(36)=   1   K(37)=   1   K(38)=   1   K(39)=   1   K(40)=   1
   K(41)=  -1   K(42)=  -1   K(43)=  -1   K(44)=  -1   K(45)=  -1   K(46)=  -1   K(47)=  -1   K(48)=  -1   K(49)=  -1   K(50)=  -1
   K(51)=  -1   K(52)=  -1   K(53)=  -1   K(54)=  -1   K(55)=  -1   K(56)=  -1   K(57)=  -1   K(58)=  -1   K(59)=  -1   K(60)=  -1
   K(61)=   1   K(62)=   1   K(63)=  -1   K(64)=  -1   K(65)=  -1   K(66)=  -1   K(67)=  -1   K(68)=  -1   K(69)=  -1   K(70)=  -1
   K(71)=  -1   K(72)=  -1   K(73)=  -1   K(74)=  -1   K(75)=  -1   K(76)=  -1   K(77)=  -1   K(78)=  -1   K(79)=  -1   K(80)=  -1
   K(81)=   4   K(82)=   3   K(83)=   3   K(84)=   1   K(85)=   1   K(86)=  -1   K(87)=  -2   K(88)=   2   K(89)=   6   K(90)=   6
   K(91)=  -4   K(92)=  -6   K(93)=   0   K(94)=   0   K(95)=   0   K(96)=   0   K(97)=   0   K(98)=  10   K(99)=   0   K(**)=  -1
 NUMKI= 16  KI= 1 1-1-1-1 1 1 -16 -20   3  10  32 3031031   2   3                                                                    
  EPS( 1)= 0.00000E+00  EPS( 2)= 0.00000E+00  EPS( 3)= 1.00000E-10  EPS( 4)= 0.00000E+00  EPS( 5)= 0.00000E+00  EPS( 6)= 0.00000E+00
0          ZONAL GRAVITATIONAL POTENTIAL HARMONIC COEFFICIENTS -   MOON   (NPLNT= 10, NZONE=  3)                                     
       J2= 2.038220000000000D-04       J3= 8.785979160478040D-06
      LJ2=  1      LJ3=  1
0TESSERAL COSINE GRAVITATIONAL POTENTIAL HARMONIC COEFFICIENTS -   MOON   (NPLNT= 10, NTESS=  3)                                     
    C2(1)= 0.000000000000000D+00    C2(2)= 3.578636611895650D-05    C3(1)= 3.027613645951770D-05    C3(2)= 1.422420643462980D-05
    C3(3)= 1.353613856801190D-05
    LC2(1)= 0    LC2(2)= 0    LC3(1)= 1    LC3(2)= 1    LC3(3)= 0
0TESSERAL  SINE  GRAVITATIONAL POTENTIAL HARMONIC COEFFICIENTS -   MOON   (NPLNT= 10, NTESS=  3)                                     
    S2(1)= 0.000000000000000D+00    S2(2)= 0.000000000000000D+00    S3(1)= 5.976611368796570D-06    S3(2)= 4.724388985911510D-06
    S3(3)=-2.972622235974020D-06
    LS2(1)= 0    LS2(2)= 0    LS3(1)= 0    LS3(2)= 0    LS3(3)= 0
1PLANETARY EPHEMERIS PROGRAM INPUT DATA   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE    7
0 EROTAT     NPLNT(-1)= -3    IPLNT=  0         JDER1=       0       JDER0=       0       JDER2=       0     INT =  0   NCENTR=  0
  LER( 1)=  0  LER( 2)=  0  LER( 3)=  0  LER( 4)=  0  LER( 5)=  0  LER( 6)=  0  LER( 7)=  0  LER( 8)=  0  LER( 9)=  0  LER(10)=  0
  LER(11)=  0  LER(12)=  0  LER(13)=  0  LER(14)=  0  LER(15)=  0  LER(16)=  0  LER(17)=  0  LER(18)=  0  LER(19)=  0  LER(20)=  0
  LER(21)=  0  LER(22)=  0  LER(23)=  0  LER(24)=  0  LER(25)=  0  LER(26)=  0  LER(27)=  0  LER(28)=  0  LER(29)=  0  LER(30)=  0
 INITIAL EPOCH (COORD.TIME)      -0.50000000 IHR= 0 IMIN= 0 SEC= 0.0000 INT1=         0 INT2=  0 FRACT= 0.0000000000000  ICND=  0
         A= 0.000000000000000D+00         E= 0.000000000000000D+00       INC= 0.000000000000000D+00       ASC= 0.000000000000000D+00
       PER= 0.000000000000000D+00      ANOM= 0.000000000000000D+00    RADIUS= 0.000000000000000D+00      FLAT= 0.000000000000000D+00
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 0.000000000000000D+00
   CON( 7)= 0.000000000000000D+00   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 2.442048000000000D+06
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 5.026750000000000D+03
   CON(23)= 0.000000000000000D+00   CON(24)= 9.210000000000000D+00  CON1( 1)= 0.000000000000000D+00  CON1( 2)= 0.000000000000000D+00
  CON1( 3)= 0.000000000000000D+00  CON1( 4)= 0.000000000000000D+00  CON1( 5)= 0.000000000000000D+00  CON1( 6)= 0.000000000000000D+00
  CON1( 7)= 0.000000000000000D+00  CON1( 8)= 0.000000000000000D+00  CON1( 9)= 0.000000000000000D+00  CON1(10)= 0.000000000000000D+00
  CON1(11)= 0.000000000000000D+00  CON1(12)= 0.000000000000000D+00
   K( 1)=   0   K( 2)=   0   K( 3)=   0   K( 4)=   0   K( 5)=   0   K( 6)=   0   K( 7)=   0   K( 8)=   0   K( 9)=   0   K(10)=   0
   K(11)=   0   K(12)=   0   K(13)=   0   K(14)=   0   K(15)=   0   K(16)=   0   K(17)=   0   K(18)=   0   K(19)=   0   K(20)=   0
   K(21)=   0   K(22)=   0   K(23)=   0   K(24)=   0   K(25)=   0   K(26)=   0   K(27)=   0   K(28)=   0   K(29)=   0   K(30)=   0
   K(31)=   1   K(32)=   1   K(33)=   1   K(34)=   1   K(35)=   1   K(36)=   1   K(37)=   1   K(38)=   1   K(39)=   1   K(40)=  -1
   K(41)=  -1   K(42)=  -1   K(43)=  -1   K(44)=  -1   K(45)=  -1   K(46)=  -1   K(47)=  -1   K(48)=  -1   K(49)=  -1   K(50)=  -1
   K(51)=  -1   K(52)=  -1   K(53)=  -1   K(54)=  -1   K(55)=  -1   K(56)=  -1   K(57)=  -1   K(58)=  -1   K(59)=  -1   K(60)=  -1
   K(61)=   1   K(62)=  -1   K(63)=  -1   K(64)=  -1   K(65)=  -1   K(66)=  -1   K(67)=  -1   K(68)=  -1   K(69)=  -1   K(70)=  -1
   K(71)=  -1   K(72)=  -1   K(73)=  -1   K(74)=  -1   K(75)=  -1   K(76)=  -1   K(77)=  -1   K(78)=  -1   K(79)=  -1   K(80)=  -1
   K(81)=  -1   K(82)=  -1   K(83)=  -1   K(84)=  -1   K(85)=  -1   K(86)=  -1   K(87)=   0   K(88)=   3   K(89)=  11   K(90)=   6
   K(91)= -12   K(92)= -14   K(93)=   0   K(94)=   0   K(95)=   0   K(96)=   0   K(97)=   0   K(98)=   0   K(99)=   0   K(**)=   0
 NUMKI=  8  KI= 0 0 0 0 0 0 0   0                                                                                                    
  EPS( 1)= 0.00000E+00  EPS( 2)= 0.00000E+00  EPS( 3)= 1.00000E-16  EPS( 4)= 0.00000E+00  EPS( 5)= 0.00000E+00  EPS( 6)= 0.00000E+00
0 MROTAT     NPLNT( 0)=-10    IPLNT= 21         JDMR1= 2442321       JDMR0= 2442321       JDMR2= 2442332     INT = -1   NCENTR=  0
  LMR( 1)=  1  LMR( 2)=  1  LMR( 3)=  0  LMR( 4)=  0  LMR( 5)=  0  LMR( 6)=  0  LMR( 7)=  3  LMR( 8)=  6  LMR( 9)=  7  LMR(10)=  0
  LMR(11)=  0  LMR(12)=  0  LMR(13)=  0  LMR(14)=  0  LMR(15)=  0  LMR(16)=  0  LMR(17)=  0  LMR(18)=  0  LMR(19)=  0  LMR(20)=  0
  LMR(21)=  0  LMR(22)=  0  LMR(23)=  0  LMR(24)=  0  LMR(25)=  0  LMR(26)=  0  LMR(27)=  0  LMR(28)=  0  LMR(29)=  0  LMR(30)=  0
 INITIAL EPOCH (COORD.TIME) 2442320.50000000 IHR= 0 IMIN= 0 SEC= 0.0000 INT1=         0 INT2=  0 FRACT= 0.0000000000000  ICND=  0
         A= 6.406890562446280D-02         E= 4.174732965826930D-01       INC= 4.427722151828226D+02       ASC=-1.116475592635000D-04
       PER=-3.914136528300000D-05      ANOM= 2.300859910558332D-01    RADIUS= 1.738090000000000D+03      FLAT= 0.000000000000000D+00
   CON( 3)= 6.317168491855110D-04   CON( 4)= 2.278098937105160D-04   CON( 5)= 2.692029800000000D-02   CON( 6)= 2.415119936539230D-02
   CON( 7)= 4.683189055748760D-03   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00  CON1( 1)= 0.000000000000000D+00  CON1( 2)= 0.000000000000000D+00
  CON1( 3)= 0.000000000000000D+00  CON1( 4)= 0.000000000000000D+00  CON1( 5)= 0.000000000000000D+00  CON1( 6)= 0.000000000000000D+00
  CON1( 7)= 0.000000000000000D+00  CON1( 8)= 0.000000000000000D+00  CON1( 9)= 0.000000000000000D+00  CON1(10)= 0.000000000000000D+00
  CON1(11)= 0.000000000000000D+00  CON1(12)= 0.000000000000000D+00
   K( 1)=   0   K( 2)=   0   K( 3)=   0   K( 4)=   0   K( 5)=   0   K( 6)=   0   K( 7)=   0   K( 8)=   0   K( 9)=   0   K(10)=   0
   K(11)=   0   K(12)=   0   K(13)=   0   K(14)=   0   K(15)=   0   K(16)=   0   K(17)=   0   K(18)=   0   K(19)=   0   K(20)=   0
   K(21)=   0   K(22)=   0   K(23)=   0   K(24)=   0   K(25)=   0   K(26)=   0   K(27)=   0   K(28)=   0   K(29)=   0   K(30)=   0
   K(31)=  -1   K(32)=  -1   K(33)=   1   K(34)=  -1   K(35)=  -1   K(36)=  -1   K(37)=  -1   K(38)=  -1   K(39)=  -1   K(40)=   1
   K(41)=  -1   K(42)=  -1   K(43)=  -1   K(44)=  -1   K(45)=  -1   K(46)=  -1   K(47)=  -1   K(48)=  -1   K(49)=  -1   K(50)=  -1
   K(51)=  -1   K(52)=  -1   K(53)=  -1   K(54)=  -1   K(55)=  -1   K(56)=  -1   K(57)=  -1   K(58)=  -1   K(59)=  -1   K(60)=  -1
   K(61)=  -1   K(62)=  -1   K(63)=  -1   K(64)=  -1   K(65)=  -1   K(66)=  -1   K(67)=  -1   K(68)=  -1   K(69)=  -1   K(70)=  -1
   K(71)=  -1   K(72)=  -1   K(73)=  -1   K(74)=  -1   K(75)=  -1   K(76)=  -1   K(77)=  -1   K(78)=  -1   K(79)=  -1   K(80)=  -1
   K(81)=   1   K(82)=   0   K(83)=   0   K(84)=  -1   K(85)=  -1   K(86)= 303   K(87)=  -2   K(88)=   2   K(89)=   6   K(90)=   6
   K(91)=  -4   K(92)=  -6   K(93)=   0   K(94)=   0   K(95)=   0   K(96)=   0   K(97)=   0   K(98)=  10   K(99)=   0   K(**)=  -1
 NUMKI= 18  KI= 1 1 1-1-1-1-1  -3  -6  -71031   2   31041   3   1   3   2                                                            
  EPS( 1)= 0.00000E+00  EPS( 2)= 0.00000E+00  EPS( 3)= 1.00000E-16  EPS( 4)= 0.00000E+00  EPS( 5)= 0.00000E+00  EPS( 6)= 0.00000E+00
1PLANETARY EPHEMERIS PROGRAM INPUT DATA   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE    8
0  MARS      NPLNT( 1)=  4    IPLNT= 14         JDPL1= 2442321       JDPL0= 2442321       JDPL2= 2442336     INT =  1   NCENTR= -1
  LPL( 1)=  0  LPL( 2)=  0  LPL( 3)=  0  LPL( 4)=  1  LPL( 5)=  1  LPL( 6)=  1  LPL( 7)=  6  LPL( 8)=  7  LPL( 9)=  8  LPL(10)=  9
  LPL(11)= 24  LPL(12)=  0  LPL(13)=  0  LPL(14)=  0  LPL(15)=  0  LPL(16)=  0  LPL(17)=  0  LPL(18)=  0  LPL(19)=  0  LPL(20)=  0
  LPL(21)=  0  LPL(22)=  0  LPL(23)=  0  LPL(24)=  0  LPL(25)=  0  LPL(26)=  0  LPL(27)=  0  LPL(28)=  0  LPL(29)=  0  LPL(30)=  0
 INITIAL EPOCH (COORD.TIME) 2442320.50000000 IHR= 0 IMIN= 0 SEC= 0.0000 INT1=         0 INT2=  0 FRACT= 0.0000000000000  ICND=  0
         A= 1.523705402995026D+00         E= 9.330246794792663D-02       INC= 2.469315159246902D+01       ASC= 3.344609353024999D+00
       PER= 3.322243666631956D+02      ANOM= 2.255191494289617D+02    RADIUS= 3.392459000000000D+03      FLAT= 0.000000000000000D+00
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 1.514460000000000D+02
   CON( 7)= 1.025956000000000D+00   CON( 8)= 5.269510000000000D+01   CON( 9)= 3.173116000000000D+02   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00  CON1( 1)= 2.443509500000000D+06  CON1( 2)= 0.000000000000000D+00
  CON1( 3)= 0.000000000000000D+00  CON1( 4)= 0.000000000000000D+00  CON1( 5)= 0.000000000000000D+00  CON1( 6)= 0.000000000000000D+00
  CON1( 7)= 0.000000000000000D+00  CON1( 8)= 0.000000000000000D+00  CON1( 9)= 0.000000000000000D+00  CON1(10)=-8.400000000000000D-11
  CON1(11)= 0.000000000000000D+00  CON1(12)= 0.000000000000000D+00
   K( 1)=   5   K( 2)=   0   K( 3)=   0   K( 4)=   0   K( 5)=   0   K( 6)=   0   K( 7)=   0   K( 8)=   0   K( 9)=   0   K(10)=   0
   K(11)=   0   K(12)=   0   K(13)=   0   K(14)=   0   K(15)=   0   K(16)=   0   K(17)=   0   K(18)=   0   K(19)=   0   K(20)=   0
   K(21)=   0   K(22)=   0   K(23)=   0   K(24)=   0   K(25)=   0   K(26)=   0   K(27)=   0   K(28)=   0   K(29)=   0   K(30)=   0
   K(31)=  -1   K(32)=  -1   K(33)=   1   K(34)=   1   K(35)=   1   K(36)=  -1   K(37)=  -1   K(38)=  -1   K(39)=  -1   K(40)=  -1
   K(41)=   1   K(42)=  -1   K(43)=  -1   K(44)=  -1   K(45)=  -1   K(46)=  -1   K(47)=  -1   K(48)=  -1   K(49)=  -1   K(50)=  -1
   K(51)=  -1   K(52)=  -1   K(53)=  -1   K(54)=  -1   K(55)=  -1   K(56)=  -1   K(57)=  -1   K(58)=  -1   K(59)=  -1   K(60)=  -1
   K(61)=   1   K(62)=  -1   K(63)=  -1   K(64)=  -1   K(65)=  -1   K(66)=  -1   K(67)=  -1   K(68)=  -1   K(69)=  -1   K(70)=   1
   K(71)=   1   K(72)=  -1   K(73)=  -1   K(74)=  -1   K(75)=  -1   K(76)=  -1   K(77)=  -1   K(78)=  -1   K(79)=   1   K(80)=   1
   K(81)=  -1   K(82)=  -1   K(83)=  -1   K(84)=  -1   K(85)=  -1   K(86)=  -1   K(87)=   1   K(88)=   3   K(89)=  11   K(90)=   6
   K(91)= -10   K(92)= -14   K(93)=   0   K(94)=   0   K(95)=   0   K(96)=   0   K(97)=   0   K(98)=   4   K(99)=   0   K(**)=   0
 NUMKI= 16  KI= 1-1-1-1 1 1 1   3   5  11  31  40  41  49  50 506                                                                    
  EPS( 1)= 0.00000E+00  EPS( 2)= 0.00000E+00  EPS( 3)= 1.00000E-16  EPS( 4)= 0.00000E+00  EPS( 5)= 0.00000E+00  EPS( 6)= 0.00000E+00
0JUPITER     NPLNT( 2)=  5    IPLNT=  0         JDPL1=       0       JDPL0=       0       JDPL2=       0     INT = 20   NCENTR=  0
  LPL( 1)=  0  LPL( 2)=  0  LPL( 3)=  0  LPL( 4)=  0  LPL( 5)=  0  LPL( 6)=  1  LPL( 7)=  0  LPL( 8)=  0  LPL( 9)=  0  LPL(10)=  0
  LPL(11)=  0  LPL(12)=  0  LPL(13)=  0  LPL(14)=  0  LPL(15)=  0  LPL(16)=  0  LPL(17)=  0  LPL(18)=  0  LPL(19)=  0  LPL(20)=  0
  LPL(21)=  0  LPL(22)=  0  LPL(23)=  0  LPL(24)=  0  LPL(25)=  0  LPL(26)=  0  LPL(27)=  0  LPL(28)=  0  LPL(29)=  0  LPL(30)=  0
 INITIAL EPOCH (COORD.TIME)      -0.50000000 IHR= 0 IMIN= 0 SEC= 0.0000 INT1=         0 INT2=  0 FRACT= 0.0000000000000  ICND=  0
         A= 0.000000000000000D+00         E= 0.000000000000000D+00       INC= 0.000000000000000D+00       ASC= 0.000000000000000D+00
       PER= 0.000000000000000D+00      ANOM= 0.000000000000000D+00    RADIUS= 7.135000000000000D+04      FLAT= 0.000000000000000D+00
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 0.000000000000000D+00
   CON( 7)= 0.000000000000000D+00   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00  CON1( 1)= 2.440000500000000D+06  CON1( 2)= 0.000000000000000D+00
  CON1( 3)= 0.000000000000000D+00  CON1( 4)= 0.000000000000000D+00  CON1( 5)= 0.000000000000000D+00  CON1( 6)= 0.000000000000000D+00
  CON1( 7)= 0.000000000000000D+00  CON1( 8)= 0.000000000000000D+00  CON1( 9)= 0.000000000000000D+00  CON1(10)=-5.600000000000000D-07
  CON1(11)= 0.000000000000000D+00  CON1(12)= 0.000000000000000D+00
   K( 1)=   0   K( 2)=   0   K( 3)=   0   K( 4)=   0   K( 5)=   0   K( 6)=   0   K( 7)=   0   K( 8)=   0   K( 9)=   0   K(10)=   0
   K(11)=   0   K(12)=   0   K(13)=   0   K(14)=   0   K(15)=   0   K(16)=   0   K(17)=   0   K(18)=   0   K(19)=   0   K(20)=   0
   K(21)=   0   K(22)=   0   K(23)=   0   K(24)=   0   K(25)=   0   K(26)=   0   K(27)=   0   K(28)=   0   K(29)=   0   K(30)=   0
   K(31)=   1   K(32)=   1   K(33)=   1   K(34)=   1   K(35)=   1   K(36)=   1   K(37)=   1   K(38)=   1   K(39)=   1   K(40)=  -1
   K(41)=  -1   K(42)=  -1   K(43)=  -1   K(44)=  -1   K(45)=  -1   K(46)=  -1   K(47)=  -1   K(48)=  -1   K(49)=  -1   K(50)=  -1
   K(51)=  -1   K(52)=  -1   K(53)=  -1   K(54)=  -1   K(55)=  -1   K(56)=  -1   K(57)=  -1   K(58)=  -1   K(59)=  -1   K(60)=  -1
   K(61)=   1   K(62)=  -1   K(63)=  -1   K(64)=  -1   K(65)=  -1   K(66)=  -1   K(67)=  -1   K(68)=  -1   K(69)=  -1   K(70)=  -1
   K(71)=  -1   K(72)=  -1   K(73)=  -1   K(74)=  -1   K(75)=  -1   K(76)=  -1   K(77)=  -1   K(78)=  -1   K(79)=  -1   K(80)=  -1
   K(81)=  -1   K(82)=  -1   K(83)=  -1   K(84)=  -1   K(85)=  -1   K(86)=  -1   K(87)=   4   K(88)=   3   K(89)=  11   K(90)=   6
   K(91)= -12   K(92)= -14   K(93)=   0   K(94)=   0   K(95)=   0   K(96)=   0   K(97)=   0   K(98)=   0   K(99)=   0   K(**)=   0
 NUMKI=  8  KI= 0 0 0 0 0 0 0   0                                                                                                    
  EPS( 1)= 0.00000E+00  EPS( 2)= 0.00000E+00  EPS( 3)= 1.00000E-16  EPS( 4)= 0.00000E+00  EPS( 5)= 0.00000E+00  EPS( 6)= 0.00000E+00
1PLANETARY EPHEMERIS PROGRAM INPUT DATA   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE    9
0 CERES 1    NPLNT( 3)= 11    IPLNT=  0         JDPL1=       0       JDPL0=       0       JDPL2=       0     INT =  0   NCENTR=  0
  LPL( 1)=  0  LPL( 2)=  0  LPL( 3)=  0  LPL( 4)=  0  LPL( 5)=  0  LPL( 6)=  0  LPL( 7)=  0  LPL( 8)=  0  LPL( 9)=  0  LPL(10)=  0
  LPL(11)=  0  LPL(12)=  0  LPL(13)=  0  LPL(14)=  0  LPL(15)=  0  LPL(16)=  0  LPL(17)=  0  LPL(18)=  0  LPL(19)=  0  LPL(20)=  0
  LPL(21)=  0  LPL(22)=  0  LPL(23)=  0  LPL(24)=  0  LPL(25)=  0  LPL(26)=  0  LPL(27)=  0  LPL(28)=  0  LPL(29)=  0  LPL(30)=  0
 INITIAL EPOCH (COORD.TIME)      -0.50000000 IHR= 0 IMIN= 0 SEC= 0.0000 INT1=         0 INT2=  0 FRACT= 0.0000000000000  ICND=  0
         A= 2.766333786805837D+00         E= 7.862823976140470D-02       INC= 2.717884370592106D+01       ASC= 2.340674566424140D+01
       PER= 1.289601772696749D+02      ANOM= 5.845329102597486D+01    RADIUS= 0.000000000000000D+00      FLAT= 0.000000000000000D+00
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 0.000000000000000D+00
   CON( 7)= 0.000000000000000D+00   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00  CON1( 1)= 2.440000500000000D+06  CON1( 2)= 0.000000000000000D+00
  CON1( 3)= 0.000000000000000D+00  CON1( 4)= 0.000000000000000D+00  CON1( 5)= 0.000000000000000D+00  CON1( 6)= 0.000000000000000D+00
  CON1( 7)= 0.000000000000000D+00  CON1( 8)= 0.000000000000000D+00  CON1( 9)= 0.000000000000000D+00  CON1(10)= 0.000000000000000D+00
  CON1(11)= 0.000000000000000D+00  CON1(12)= 0.000000000000000D+00
   K( 1)=   0   K( 2)=   0   K( 3)=   0   K( 4)=   0   K( 5)=   0   K( 6)=   0   K( 7)=   0   K( 8)=   0   K( 9)=   0   K(10)=   0
   K(11)=   0   K(12)=   0   K(13)=   0   K(14)=   0   K(15)=   0   K(16)=   0   K(17)=   0   K(18)=   0   K(19)=   0   K(20)=   0
   K(21)=   0   K(22)=   0   K(23)=   0   K(24)=   0   K(25)=   0   K(26)=   0   K(27)=   0   K(28)=   0   K(29)=   0   K(30)=   0
   K(31)=   1   K(32)=   1   K(33)=   1   K(34)=   1   K(35)=   1   K(36)=   1   K(37)=   1   K(38)=   1   K(39)=   1   K(40)=  -1
   K(41)=  -1   K(42)=  -1   K(43)=  -1   K(44)=  -1   K(45)=  -1   K(46)=  -1   K(47)=  -1   K(48)=  -1   K(49)=  -1   K(50)=  -1
   K(51)=  -1   K(52)=  -1   K(53)=  -1   K(54)=  -1   K(55)=  -1   K(56)=  -1   K(57)=  -1   K(58)=  -1   K(59)=  -1   K(60)=  -1
   K(61)=   1   K(62)=  -1   K(63)=  -1   K(64)=  -1   K(65)=  -1   K(66)=  -1   K(67)=  -1   K(68)=  -1   K(69)=  -1   K(70)=  -1
   K(71)=  -1   K(72)=  -1   K(73)=  -1   K(74)=  -1   K(75)=  -1   K(76)=  -1   K(77)=  -1   K(78)=  -1   K(79)=  -1   K(80)=  -1
   K(81)=  -1   K(82)=  -1   K(83)=  -1   K(84)=  -1   K(85)=  -1   K(86)=  -1   K(87)=   0   K(88)=   3   K(89)=  11   K(90)=   6
   K(91)= -12   K(92)= -14   K(93)=   0   K(94)=   0   K(95)=   0   K(96)=   0   K(97)=   0   K(98)=   0   K(99)=   0   K(**)=   0
 NUMKI=  8  KI= 0 0 0 0 0 0 0   0                                                                                                    
  EPS( 1)= 0.00000E+00  EPS( 2)= 0.00000E+00  EPS( 3)= 1.00000E-16  EPS( 4)= 0.00000E+00  EPS( 5)= 0.00000E+00  EPS( 6)= 0.00000E+00
0 THERE IS NO INPUT ET-UT2, A1-UT1, OR WOBBLE TABLE
0OBSERVING SITES WITH ADJUSTABLE COORDINATES

        SITE   LSCRD   RADIUS (KM)  LONGITUDE (DEG)  LATITUDE (DEG) KSCRD
   1. HAYSTACK 0 0 0 6368.551653028   71.4886666667   42.4315183830  0 
   2. B30LNCLN 0 0 0 6368.485020842   71.2666933000   42.2678790610  0 
   3. B10LNCLN 0 0 0 6368.475524712   71.2660932800   42.2681595950  0 
   4. B04LNCLN 0 0 0 6368.478870405   71.2657733700   42.2683426048  0 
   5. B01LNCLN 0 0 0 6368.475794561   71.2657652700   42.2682041262  0 
   6. MILLSTON 0 0 0 6368.563831130   71.4913888889   42.4256609690  0 
   7. ARECIBO  0 0 0 6376.560245971   66.7530277778   18.2287613852  0 
   8. 85JPLVNS 0 0 0 6372.177000000  116.7940075000   35.0665981000  0 
   9. 11DSPION 0 0 0 5206.350322378  116.8497745987 3673.7851760000 -1*
  10. 12DSECHO 0 0 0 5212.050800000 -243.1946300000 3665.6468000000 -1*
  11. 14DSMARS 0 0 0 5203.997400000 -243.1105900000 3677.0630000000 -1*
  12. 41DSWOOM 0 0 0 5450.197800000 -136.8875900000-3302.3262000000 -1*
  13. 42DSCANB 0 0 0 5205.361028200 -148.9809579880-3674.6129890000 -1*
  14. 51DSJOHA 0 0 0 5742.938000000  -27.6854600000-2768.7193000000 -1*
  15. 61DSMADR 0 0 0 4862.604400000 -355.7510900000 4114.8518000000 -1*
  16. 62DSCEBR 0 0 0 4860.811400000 -355.6322900000 4116.9660000000 -1*
  17. AFLASER  0 0 0 5391.827000000  110.7244167000 3400.6790000000 -1*
  18. MCDONALD 0 0 0 6374.717549422  104.0222500000   30.5030810430  0 
  19. 6USNAVAL 0 0 0 6378.000000000   77.0660375000    0.0000000000  0 
  20. 8USNAVAL 0 0 0 6378.000000000   77.0655416667    0.0000000000  0 
  21. 9USNAVAL 0 0 0 6378.000000000   77.0654625000    0.0000000000  0 
  22. MUSNAVAL 0 0 0 6378.000000000   77.0655416667    0.0000000000  0 
  23. CAPETOWN 0 0 0 6378.000000000  -18.4765833333    0.0000000000  0 
1PLANETARY EPHEMERIS PROGRAM INPUT DATA   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   10
0OBSERVING SITES WITH ADJUSTABLE COORDINATES

        SITE   LSCRD   RADIUS (KM)  LONGITUDE (DEG)  LATITUDE (DEG) KSCRD
  24. GRENWICH 0 0 0 6378.000000000    0.0000000000    0.0000000000  0 
  25. CAMBRIDG 0 0 0 6378.000000000   -0.0947916652    0.0000000000  0 
  26. RADCLIFF 0 0 0 6378.000000000    1.2516666667    0.0000000000  0 
  27. OTTAWA   0 0 0 6378.000000000   75.7164583200    0.0000000000  0 
  28. PARIS    0 0 0 6378.000000000   -2.3371249980    0.0000000000  0 
  29. TOULOUSE 0 0 0 6378.000000000   -1.4624999600    0.0000000000  0 
  30. NICE     0 0 0 6378.000000000   -7.3004166650    0.0000000000  0 
  31. BESANCON 0 0 0 6378.000000000   -5.9892499600    0.0000000000  0 
  32. UCCLE    0 0 0 6378.000000000   -4.3582083333    0.0000000000  0 
  33. GTOKYO   0 0 0 6378.000000000 -139.5407500000    0.0000000000  0 
  34. STRASBRG 0 0 0 6378.000000000   -7.7683333333    0.0000000000  0 
  35. BERLIN   0 0 0 6378.000000000  -13.1066666667    0.0000000000  0 
  36. EDINBRG  0 0 0 6378.000000000    3.1833333333    0.0000000000  0 
    9 CYLINDRICAL COORDINATES(EQUAT. RADIUS(KM),LONGITUDE(DEG),Z(KM)) *
0OBSERVED SPOTS ON OTHER BODIES

      SPOT  PLANET  NO. LSPCRD  RADIUS (KM)  LONGITUDE (DEG) LATITUDE (DEG)
   1. SUR1   MOON    10  0 0 0 1735.47400000  316.6760000000  -2.5020000000     
   2. SUR3   MOON    10  0 0 0 1736.10600000  336.6825100000  -3.0550000000     
   3. SUR5   MOON    10  0 0 0 1735.11400000   23.2170000000   1.4060000000     
   4. SUR6   MOON    10  0 0 0 1736.43900000  358.6322900000   0.4590000000     
   5. SUR7   MOON    10  0 0 0 1739.05100000  348.5630000000 -40.9750000000     
   6. APL2   MOON    10  0 0 0 1735.63900000   23.4602000000   0.6707000000     
   7. AP11   MOON    10  0 0 0 1735.50100000   23.4119800000   0.6923100000     
   8. AP14   MOON    10  0 0 0 1736.35900000  -17.5387600000  -3.6244600000     
   9. AP15   MOON    10  0 0 0 1735.50870000    3.5682200000  26.1540400000     
  10. AL12   MOON    10  0 0 0 1736.01200000  -23.4850400000  -2.9905100000     
  11. AL14   MOON    10  0 0 0 1736.36300000  -17.5379800000  -3.6240900000     
  12. AL15   MOON    10  0 0 0 1735.50870000    3.5696500000  26.1547800000     
  13. AL16   MOON    10  0 0 0 1737.44700000   15.4362500000  -8.9554100000     
  14. AL17   MOON    10  0 0 0 1734.81400000   30.7084300000  20.2096600000     
  15. OLYM   MARS     4  1 1 1 3396.01537700 -133.7123400000  17.8123400000     
0THERE ARE NO INPUT PULSARS
0THERE ARE NO SKY CORRECTION STAR CATALOG ERROR MODELS
0CONSTANT BIASES FOR PLANETARY RADAR OBSERVATIONS                                                                                    
0      PLANET NO. RECV SITE NO. SEND SITE NO. SER. LRBS  BIAS1(SEC)  BIAS2(C/S)
   1. ########  1  ARECIBO    7  ARECIBO    7 6764  0 0  0.00000E+00 0.00000E+00
   2. ########  2  ARECIBO    7  ARECIBO    7 6764  0 0  0.00000E+00 0.00000E+00
      ########   2 PLANETS WHICH ARE NOT INPUT, BUT BIASES NOT ADJUSTED, SO WARNING ONLY
0EQUINOX-EQUATOR CORRECTIONS FOR OPTICAL OBSERVATIONS                                                                                

        SITE   NO. SER. LEQN   DEQUINOX    DEQUATOR    DLATITUDE          SITE   NO. SER. LEQN   DEQUINOX    DEQUATOR    DLATITUDE  
   1. 6USNAVAL  19 6956 0 0 0  0.0000E+00  0.0000E+00  0.0000E+00
0PLANETARY PHASE CORRECTIONS FOR OPTICAL OBSERVATIONS                                                                                
0      PLANET NO.   SITE   NO. SER.        LPHS              APHASE(1)   APHASE(2)   APHASE(3) APHASE(4&7) APHASE(5&8) APHASE(6&9)
   1. ########  1 6USNAVAL  19 6956 0 0 0 0 0 0 0 0 0
      ########   1 PLANETS WHICH ARE NOT INPUT, BUT PHASES NOT ADJUSTED, SO WARNING ONLY
1PLANETARY EPHEMERIS PROGRAM INPUT DATA   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   11
0PARAMETERS FOR GLOBAL INPUT DATA DELETIONS
 TDLT0=       0.000000  TDLTON=       0.000000  TDLTOF=       0.000000
-INPUT DATA FOR ALTERATIONS IN ERROR WEIGHTINGS,DELETIONS FROM OBSERVATION SERIES AND THE GENERATION OF DUMMY OBSERVATIONS        
0NTP NSEQ NCDF PLANET NPLN SITE1 NSITE SERIES SITE2 NSITE SPOT NSPT ERROR WEIGHTS  ACCTIM ITIME FDEV    FREQUENCY        CTLG

  -5 1000   1   MOON   10 HAYSTACK   1 LASR  HAYSTACK   1 AP15   9 1.0E+00 1.0E+00 1.0E-09 1   0.0000 4.317800000000D+14          
 ERRORS  1.0E-08, 0.0E+00 USED FOR DUMMY OBS FROM JD 2442326 10H  1M 30.5000S TO JD 2442327  2H  3M  0.0000S EVERY  0 DAYS  7200 SEC

  -5 1010   6   MOON   10 NICE      30 ASTG             0        0 1.0E+00 1.0E+00 1.0E-09 1   0.0000 4.317800000000D+14          
 ERRORS  1.0E-02, 1.0E-01 USED FOR DUMMY OBS FROM JD 2442326 10H  1M 30.5000S TO JD 2442327  2H  3M  0.0000S EVERY  0 DAYS  7200 SEC

  -5 1020   1   MARS    4 ARECIBO    7 RADR  ARECIBO    7 OLYM  15 1.0E+00 1.0E+00 1.0E-09 1   0.0000 2.380000000000D+09          
 ERRORS  1.0E-08, 0.0E+00 USED FOR DUMMY OBS FROM JD 2442326 20H  1M 30.5000S TO JD 2442327 12H  3M  0.0000S EVERY  0 DAYS  7200 SEC

  -5 1030   3   MARS    4 ARECIBO    7 RADR  ARECIBO    7 OLYM  15 1.0E+00 1.0E+00 1.0E-09 1   0.0000 2.380000000000D+09          
 ERRORS  1.0E-08, 0.0E+00 USED FOR DUMMY OBS FROM JD 2442326 20H  1M 30.5000S TO JD 2442327  4H  3M  0.0000S EVERY  0 DAYS  7200 SEC

  -5 1040   5   MARS    4 PARIS     28 ASTM             0        0 1.0E+00 1.0E+00 1.0E-09 1   0.0000 4.317800000000D+14          
 ERRORS  1.0E-02, 1.0E-01 USED FOR DUMMY OBS FROM JD 2442327  0H  1M 30.5000S TO JD 2442327 12H  3M  0.0000S EVERY  0 DAYS  7200 SEC
0THERE IS NO KALMAN FILTERING
0  SETTING UP, READING IN AND WRITING OUT INPUT DATA REQUIRED
           0H  0M  0.02S REAL TIME
           0H  0M  0.01S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.50000
1 PLANET  INTEGRATION (ITERAT  1)         TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   12
-
-
-
-

                       **********    **             *********    **           **   ***********   ************
                       ***********   **            ***********   ***          **   ***********   ************
                       **       **   **            **       **   ****         **   **                 **     
                       **       **   **            **       **   ** **        **   **                 **     
                       **       **   **            **       **   **  **       **   **                 **     
                       **       **   **            **       **   **   **      **   **                 **     
                       ***********   **            ***********   **    **     **   ********           **     
                       **********    **            ***********   **     **    **   ********           **     
                       **            **            **       **   **      **   **   **                 **     
                       **            **            **       **   **       **  **   **                 **     
                       **            **            **       **   **        ** **   **                 **     
                       **            **            **       **   **         ****   **                 **     
                       **            ***********   **       **   **          ***   ***********        **     
                       **            ***********   **       **   **           **   ***********        **     
-
-
                                          **           **    **           **   **       **
                                          **           **    ***          **   **      ** 
                                          **           **    ****         **   **     **  
                                          **           **    ** **        **   **    **   
                                          **           **    **  **       **   **   **    
                                          **           **    **   **      **   *** **     
                                          **           **    **    **     **   *****      
                                          **           **    **     **    **   ** **      
                                          **           **    **      **   **   **  **     
                                          **           **    **       **  **   **   **    
                                          **           **    **        ** **   **    **   
                                          **           **    **         ****   **     **  
                                          ***********  **    **          ***   **      ** 
                                          ***********  **    **           **   **       **
- PEP VERSION= 20210302 PEP.PEPLOAD.PEP790                          
1 EMBARY  INTEGRATION (ITERAT  1)         TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   13
0   INFORMATION FROM FIRST TWO RECORDS OF PERTURBING PLANET PERIPHERAL DATA SET 90
  N-BODY RUN 443-0 9-BODY INTEGRATION  442D8 INNER,  384F4 OUTER  I.C.            3/20/73 PERTURBING PLNT TAPE FROM N-BODY INTEG.
 IDIR= 1 JDT1= 2441841 JDT2= 2442400 NMOON= 0 INITIAL CONDITIONS ARE
0MERCURY RUN 443- 3/21/73 JD0=2440000.5  A = 3.870988248079130D-01  E = 2.056143416457730D-01   INC= 2.860343745748160D+01
                                        ASC= 1.085941287077680D+01 PER= 6.693261684748430D+01 ANOM0= 9.083965359314600D+01
0 VENUS  RUN 443- 3/21/73 JD0=2440000.5  A = 7.233280763362320D-01  E = 6.757538488379489D-03   INC= 2.446691434342910D+01
                                        ASC= 7.978123733702580D+00 PER= 1.237851072721730D+02 ANOM0= 2.744064790468620D+02
0 EMBARY RUN 443- 3/21/73 JD0=2440000.5  A = 9.999827653637140D-01  E = 1.675271216169400D-02   INC= 2.344335783108840D+01
                                        ASC= 6.659250140529500D-04 PER= 1.020340640186800D+02 ANOM0= 1.393658182548250D+02
0  MARS  RUN 443- 3/21/73 JD0=2440000.5  A = 1.523604534884400D+00  E = 9.345778628544589D-02   INC= 2.469301678564210D+01
                                        ASC= 3.346655760227140D+00 PER= 3.321011352392630D+02 ANOM0= 8.989975217532620D+01
0JUPITER RUN 443- 3/20/73 JD0=2440000.5  A = 5.202945960488120D+00  E = 4.819167141557811D-02   INC= 2.325294770783580D+01
                                        ASC= 3.260973418885740D+00 PER= 1.037371340811626D+01 ANOM0= 1.411862510204257D+02
0 SATURN RUN 443- 3/20/73 JD0=2440000.5  A = 9.526046869249818D+00  E = 5.461675768589707D-02   INC= 2.257317229472867D+01
                                        ASC= 5.968803885438612D+00 PER= 8.759383460863923D+01 ANOM0= 2.897400092498020D+02
0 URANUS RUN 443- 3/20/73 JD0=2440000.5  A = 1.927444351332630D+01  E = 5.124121845145173D-02   INC= 2.367141284461453D+01
                                        ASC= 1.847313091106734D+00 PER= 1.683442408256917D+02 ANOM0= 6.990628140027202D+00
0NEPTUNE RUN 443- 3/20/73 JD0=2440000.5  A = 3.011376065452675D+01  E = 6.984876998530826D-03   INC= 2.231473025943044D+01
                                        ASC= 3.519876825815031D+00 PER= 5.409517825004429D+01 ANOM0= 1.775740764613225D+02
0 PLUTO  RUN 443- 3/20/73 JD0=2440000.5  A = 3.974605251930742D+01  E = 2.522316558278015D-01   INC= 2.365507811729088D+01
                                        ASC= 4.377229536665621D+01 PER= 1.819638797966960D+02 ANOM0= 3.298440994501460D+02
0  MOON  RUN 440- 8/ 7/72 JD0=2440000.5  A = 2.571514374642509D-03  E = 5.561544735761867D-02   INC= 2.839685438634148D+01
                                        ASC= 3.312959209874043D+00 PER= 2.262712472223586D+02 ANOM0= 1.548856956849208D+02
1 EMBARY  2442304 TO  2442347 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   14
0 EMBARY     NPLNT(-3)=  3    JPLNT= 13          JDP1= 2442304        JDP0= 2442321        JDP2= 2442347     INT =  2   NCENTR= -1
 INITIAL EPOCH (COORD.TIME) JED=2442320.50000000 IHR= 0 IMIN= 0 SEC= 0.0000 INT1=         0 INT2=  0 FRACT= 0.00000000000000D+00
         A= 1.000002347957079D+00         E= 1.669041814744669D-02       INC= 2.344258335876651D+01       ASC= 8.391300825564327D-04
       PER= 1.022643964965586D+02      ANOM= 2.657488867161210D+02    RADIUS= 6.378166000000000D+03      FLAT= 3.352329869259135D-03
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 0.000000000000000D+00
   CON( 7)= 0.000000000000000D+00   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00  CON1( 1)= 2.442320500000000D+06  CON1( 2)= 0.000000000000000D+00
  CON1( 3)= 0.000000000000000D+00  CON1( 4)= 0.000000000000000D+00  CON1( 5)= 0.000000000000000D+00  CON1( 6)= 0.000000000000000D+00
  CON1( 7)= 0.000000000000000D+00  CON1( 8)= 0.000000000000000D+00  CON1( 9)= 0.000000000000000D+00  CON1(10)=-4.600000000000000D-10
  CON1(11)= 0.000000000000000D+00  CON1(12)= 0.000000000000000D+00  TCON( 1)= 0.000000000000000D+00  TCON( 2)= 0.000000000000000D+00
  TCON( 3)= 0.000000000000000D+00  TCON( 4)= 0.000000000000000D+00  TCON( 5)= 0.000000000000000D+00  TCON( 6)= 0.000000000000000D+00
  TCON( 7)= 0.000000000000000D+00  TCON( 8)= 0.000000000000000D+00  TCON( 9)= 0.000000000000000D+00  TCON(10)= 0.000000000000000D+00
  TCON(11)= 0.000000000000000D+00  TCON(12)= 0.000000000000000D+00  TCON(13)= 0.000000000000000D+00  TCON(14)= 0.000000000000000D+00
  TCON(15)= 0.000000000000000D+00  TCON(16)= 0.000000000000000D+00  TCON(17)= 0.000000000000000D+00  TCON(18)= 0.000000000000000D+00
  TCON(19)= 0.000000000000000D+00  TCON(20)= 0.000000000000000D+00  TCON(21)= 0.000000000000000D+00  TCON(22)= 0.000000000000000D+00
  TCON(23)= 0.000000000000000D+00  TCON(24)= 0.000000000000000D+00  TCON(25)= 0.000000000000000D+00  TCON(26)= 0.000000000000000D+00
  TCON(27)= 0.000000000000000D+00  TCON(28)= 0.000000000000000D+00  TCON(29)= 0.000000000000000D+00  TCON(30)= 0.000000000000000D+00
  KP( 1)=   0  KP( 2)=   0  KP( 3)=   0  KP( 4)=   0  KP( 5)=   0  KP( 6)=   0  KP( 7)=   0  KP( 8)=   0  KP( 9)=   0  KP(10)=   0
  KP(11)=   0  KP(12)=   0  KP(13)=   0  KP(14)=   0  KP(15)=   0  KP(16)=   0  KP(17)=   0  KP(18)=   0  KP(19)=   0  KP(20)=   0
  KP(21)=   0  KP(22)=   0  KP(23)=   0  KP(24)=   0  KP(25)=   0  KP(26)=   0  KP(27)=   0  KP(28)=   0  KP(29)=   0  KP(30)=   0
  KP(31)=  -1  KP(32)=  -1  KP(33)=   1  KP(34)=   1  KP(35)=   1  KP(36)=  -1  KP(37)=  -1  KP(38)=  -1  KP(39)=  -1  KP(40)=   1
  KP(41)=  -1  KP(42)=  -1  KP(43)=  -1  KP(44)=  -1  KP(45)=  -1  KP(46)=  -1  KP(47)=  -1  KP(48)=  -1  KP(49)=  -1  KP(50)=  -1
  KP(51)=  -1  KP(52)=  -1  KP(53)=  -1  KP(54)=  -1  KP(55)=  -1  KP(56)=  -1  KP(57)=  -1  KP(58)=  -1  KP(59)=  -1  KP(60)=  -1
  KP(61)=   1  KP(62)=   1  KP(63)=   1  KP(64)=  -1  KP(65)=  -1  KP(66)=  -1  KP(67)=  -1  KP(68)=  -1  KP(69)=  -1  KP(70)=   1
  KP(71)=   1  KP(72)=   1  KP(73)=  -1  KP(74)=  -1  KP(75)=  -1  KP(76)=  -1  KP(77)=  -1  KP(78)=  -1  KP(79)=   1  KP(80)=   1
  KP(81)=  -1  KP(82)=   0  KP(83)=  -1  KP(84)=  -1  KP(85)=  -1  KP(86)=  -1  KP(87)=   2  KP(88)=   2  KP(89)=   6  KP(90)=   6
  KP(91)=  -3  KP(92)=  -6  KP(93)=   0  KP(94)=   0  KP(95)=   0  KP(96)=   0  KP(97)=   0  KP(98)=   4  KP(99)=   0  KP(**)=  -1
  KKP( 1)=  0  KKP( 2)=  0  KKP( 3)=  0  KKP( 4)=  0  KKP( 5)=  0  KKP( 6)=  0  KKP( 7)=  0  KKP( 8)=  0  KKP( 9)=  0  KKP(10)=  0
  KKP(11)=  0  KKP(12)=  0  KKP(13)=  0  KKP(14)=  0  KKP(15)=  0  KKP(16)=  0  KKP(17)=  0  KKP(18)=  0  KKP(19)=  0  KKP(20)=  0
  KKP(21)=  0  KKP(22)=  0  KKP(23)=  0  KKP(24)=  0  KKP(25)=  0  KKP(26)=  0  KKP(27)=  0  KKP(28)=  0  KKP(29)=  0  KKP(30)=  0
  KKP(31)=  0  KKP(32)=  0  KKP(33)=  0  KKP(34)=  0  KKP(35)=  0  KKP(36)=  0  KKP(37)=  0  KKP(38)=  0  KKP(39)=  0  KKP(40)=  0
  KKP(41)=  0  KKP(42)=  0  KKP(43)=  0  KKP(44)=  0  KKP(45)=  0  KKP(46)=  0  KKP(47)=  0  KKP(48)=  0  KKP(49)=  0  KKP(50)=  0
  KKP(51)=  0  KKP(52)=  0  KKP(53)=  0  KKP(54)=  0  KKP(55)=  0  KKP(56)=  0  KKP(57)=  0  KKP(58)=  0  KKP(59)=  0  KKP(60)=  0
  KKP(61)=  0  KKP(62)=  0  KKP(63)=  0  KKP(64)=  0  KKP(65)=  0  KKP(66)=  0  KKP(67)=  0  KKP(68)=  0  KKP(69)=  0  KKP(70)=  0
  KKP(71)=  0  KKP(72)=  0  KKP(73)=  0  KKP(74)=  0  KKP(75)=  0  KKP(76)=  0  KKP(77)=  0  KKP(78)=  0  KKP(79)=  0  KKP(80)=  0
  KKP(81)=  0  KKP(82)=  0  KKP(83)=  0  KKP(84)=  0  KKP(85)=  0  KKP(86)=  0  KKP(87)=  0  KKP(88)=  0  KKP(89)=  0  KKP(90)=  0
  KKP(91)=  0  KKP(92)=  0  KKP(93)=  0  KKP(94)=  0  KKP(95)=  0  KKP(96)=  0  KKP(97)=  0  KKP(98)=  0  KKP(99)=  0  KKP(**)=  0
 NUMKI= 16  KI= 1 1 1 1-1-1-1   3   4  10  31  32  33  40  41  42
 EPSP( 1)= 0.00000E+00 EPSP( 2)= 0.00000E+00 EPSP( 3)= 1.00000E-09 EPSP( 4)= 0.00000E+00 EPSP( 5)= 0.00000E+00 EPSP( 6)= 0.00000E+00
    ICND =  0   IPARP = 13   NUMERICAL INTEGRATION METHOD IS  COWELL  ADAMS MOULTON       
0  SETUP FOR  EMBARY  NUMERICAL INTEGRATION  REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1 EMBARY  2442304 TO  2442347 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   15
0 STARTING PROCEDURE OF     65 STEPS REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1 EMBARY  2442304 TO  2442347 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   16
 JED-2400000         X*10**0             Y*10**0             Z*10**0           DX/DT*10**0         DY/DT*10**0         DZ/DT*10**0 
 42320.500000  0.9958267665339834  0.0978029987337707  0.0424031416132450 -0.0021118110838183  0.0156387910138846  0.0067813266678214
          2  9.9582442838097D-01 9.7802769097067D-02 4.2403042052722D-02 1.0559030626941D-03-7.8193771473804D-03-3.3906553727974D-03
          3  3.1845134765038D-01-1.8051361518714D+00-7.8274754953664D-01 1.7148272189569D-02 2.4664488443517D-04 1.0684157875727D-04
          4  6.2101872116389D-07-4.2403141608697D-02 9.7788414259925D-02 9.9316481157366D-08-6.7813266670941D-03 1.5638821940912D-02
          5  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00-1.0559023315047D-03 7.8193717326360D-03 3.3906530248441D-03
          6  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          7  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          8  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          9  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         10  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         11  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         12  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         13  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
 42312.500000  1.0033175423423898 -0.0278343963338880 -0.0120759664095351  0.0002393992541792  0.0157212704890175  0.0068170783573172
          2  1.0061885250952D+00 1.6082018723478D-01 6.9728661837647D-02-3.6424445959501D-03-7.7628481658125D-03-3.3661107278389D-03
          3  1.8298885295611D-01-1.7899389236572D+00-7.7615685740980D-01 1.6506075141750D-02-4.0149952360217D-03-1.7410898426993D-03
          4 -1.7749136991503D-07 1.2075927036823D-02-2.7849015161887D-02 9.9939589203452D-08-6.8170683299334D-03 1.5721247840694D-02
          5 -9.5776217531040D-04-6.2884814785674D-02-2.7268190327520D-02 1.2939434059307D-03 7.8280260695221D-03 3.3943908664408D-03
          6  2.1220356097675D-03 5.3031652442348D-04 1.7350292425447D-04-5.3227002456620D-04-1.2849813749785D-04-4.1454069437014D-05
          7 -4.5160511692687D-08-7.8862956243303D-08-2.2726336941020D-08-7.5770707223558D-10 1.4218437264088D-08 4.2263352779672D-09
          8  2.7941311753081D-10 9.9743780936053D-12 4.3244119654817D-12-7.0024322132247D-11-1.0335052374978D-12-4.4772202480864D-13
          9 -2.1859895530956D+01-1.2217353249823D+00-5.2964765565726D-01 5.4795208740524D+00 1.9027268019602D-01 8.2475541658808D-02
         10 -2.9198767515105D-07 1.5070430879475D-08-7.2172958307238D-08 7.3114852799241D-08-5.2914511441340D-09 1.7461548426269D-08
         11 -3.7959444958375D-07 8.4935319748213D-08 4.4970580944482D-08 9.5137214810254D-08-2.1522790458398D-08-1.1362446040554D-08
         12  1.8587867566557D-10 1.0389545240837D-11 4.5042719456242D-12-4.6602636250148D-11-1.6182404483141D-12-7.0148905439026D-13
         13  9.3390662475749D-11 2.4135252768810D-12 1.0460033121584D-12-2.3405887544037D-11-1.1648913630094D-13-5.0397907596072D-14
 42304.500000  0.9920498651655937 -0.1529541390634512 -0.0663305047033084  0.0025706026003083  0.0155102277773183  0.0067255527597789
          2  1.0537441560312D+00 2.1928970411453D-01 9.5081957038157D-02-8.2063859813337D-03-6.6879427512373D-03-2.8999770122049D-03
          3  5.7605658570667D-02-1.7418018652882D+00-7.5528283554296D-01 1.4647469266886D-02-7.9381344894151D-03-3.4422326660235D-03
          4 -9.7206030304705D-07 6.6330342130999D-02-1.5296835819852D-01 9.8326236783308D-08-6.7255318339471D-03 1.5510150292906D-02
          5 -2.0565458251395D-02-1.2408113339289D-01-5.3804095176148D-02 3.5923146163306D-03 7.3993792374548D-03 3.2085057995730D-03
          6  8.5662317834950D-03 1.9807080197522D-03 6.2868783381096D-04-1.0830518115923D-03-2.2846593612322D-04-6.9808876399283D-05
          7 -7.5580402158687D-08-1.3333829771359D-07-2.4334038744924D-08 1.8375893198744D-08 4.1701453490456D-09-2.4334031082337D-09
          8  1.1248433869829D-09-7.5690370946369D-12-3.2877923467791D-12-1.4170469746037D-10 7.0137550222029D-12 3.0418447144463D-12
          9 -8.8031028921485D+01-1.1389780451751D+00-4.9339066065678D-01 1.1092545644951D+01-3.3679530917233D-01-1.4610404916362D-01
         10 -1.1736990049438D-06 1.0934087077528D-07-2.6894805608417D-07 1.4766569882872D-07-1.9878030721855D-08 3.0991950243799D-08
         11 -1.5307870836560D-06 3.4988914477598D-07 1.8405643632583D-07 1.9349262175386D-07-4.5185352725285D-08-2.3600958973570D-08
         12  7.4884539862144D-10 9.6723778248985D-12 4.1906918378047D-12-9.4397922366832D-11 2.8698698962416D-12 1.2448750419163D-12
         13  3.7598011125614D-10-6.1723803676234D-12-2.6784145498734D-12-4.7365537983726D-11 2.7925720870920D-12 1.2111754070058D-12

 INTEGRATION COMPLETED IN ONE DIRECTION FROM EPOCH, STARTED IN OTHER

1 EMBARY  2442304 TO  2442347 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   17
 JED-2400000         X*10**0             Y*10**0             Z*10**0           DX/DT*10**0         DY/DT*10**0         DZ/DT*10**0 
 42320.500000  0.9958267665339834  0.0978029987337707  0.0424031416132450 -0.0021118110838183  0.0156387910138846  0.0067813266678214
          2  9.9582442838097D-01 9.7802769097067D-02 4.2403042052722D-02 1.0559030626941D-03-7.8193771473804D-03-3.3906553727974D-03
          3  3.1845134765038D-01-1.8051361518714D+00-7.8274754953664D-01 1.7148272189569D-02 2.4664488443517D-04 1.0684157875727D-04
          4  6.2101872116389D-07-4.2403141608697D-02 9.7788414259925D-02 9.9316481157366D-08-6.7813266670941D-03 1.5638821940912D-02
          5  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00-1.0559023315047D-03 7.8193717326360D-03 3.3906530248441D-03
          6  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          7  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          8  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          9  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         10  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         11  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         12  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         13  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
 42328.500000  0.9695920083248554  0.2215965763885594  0.0960828366725237 -0.0044392044573415  0.0152603653648169  0.0066172464046292
          2  1.0228605722435D+00 3.8471804888012D-02 1.6675688607199D-02 5.6706837288678D-03-6.8414898813745D-03-2.9666544906526D-03
          3  4.5392580939892D-01-1.7859256702639D+00-7.7441834677699D-01 1.6505401283018D-02 4.5330189746938D-03 1.9655060030569D-03
          4  1.4056865573308D-06-9.6082871860572D-02 2.2158244357386D-01 9.6493093832975D-08-6.6172549489513D-03 1.5260446742715D-02
          5 -1.7756915342071D-02 6.1041364726109D-02 2.6468950733460D-02-3.3699534735755D-03 7.3672658544376D-03 3.1946250249325D-03
          6  2.1186050169520D-03 5.9391269815005D-04 2.0399226081413D-04 5.3102286876707D-04 1.5271626117038D-04 5.3047237977448D-05
          7 -1.4635058491742D-07 3.6614810949739D-09 1.3713087203795D-08-2.4807057689699D-08 6.9344759762346D-09 5.0345998207483D-09
          8  2.7838941272082D-10 3.3550950989669D-11 1.4548190552013D-11 6.9632650481703D-11 9.9093719977274D-12 4.2965584463950D-12
          9 -2.1761795452391D+01-3.0725334210479D+00-1.3321929010689D+00-5.4424286336849D+00-8.8733615941404D-01-3.8473703561604D-01
         10 -2.9214085789654D-07-1.0033794048658D-08-8.0790627222180D-08-7.3165472239326D-08-4.1653008901315D-09-2.0711968396322D-08
         11 -3.8022793655287D-07 7.9634388113714D-08 4.2708593134169D-08-9.5344780895832D-08 1.9474211750366D-08 1.0487831901616D-08
         12  1.8490022398058D-10 2.6107472657786D-11 1.1319898165616D-11 4.6233087444321D-11 7.5380969970200D-12 3.2684643214195D-12
         13  9.3041782636146D-11 1.0260123075322D-11 4.4484213623590D-12 2.3273265755864D-11 3.0721172577104D-12 1.3320151204417D-12
 42336.500000  0.9249818945728010  0.3411874635049775  0.1479402600264554 -0.0066982012135900  0.0145891155427768  0.0063261911164779
          2  1.0857381268474D+00-8.9524824359561D-03-3.8889238557841D-03 9.9790647761062D-03-4.8491505790539D-03-2.1027636620373D-03
          3  5.7917348435647D-01-1.7334083864023D+00-7.5164656243310D-01 1.4604204500426D-02 8.5196511707826D-03 3.6942039392875D-03
          4  2.1591112996001D-06-1.4794039302857D-01 3.4117417169586D-01 9.1496055054228D-08-6.3262067319152D-03 1.4589243648365D-02
          5 -5.3586047492027D-02 1.1671292080502D-01 5.0609553107306D-02-5.5590741251247D-03 6.4793995391775D-03 2.8096414529114D-03
          6  8.5424185605799D-03 2.5245236718946D-03 8.8774345896983D-04 1.0791685382016D-03 3.3620081817320D-04 1.2087423432010D-04
          7 -3.1713816076661D-07 2.4267831547724D-08 4.4369910807247D-08-3.0591044240012D-08-6.7947361819229D-09 1.1729431204508D-09
          8  1.1160461603663D-09 1.8451664445046D-10 8.0005520119450D-11 1.3996524282824D-10 2.9597722183192D-11 1.2833706343108D-11
          9 -8.7216643986502D+01-1.6231886078110D+01-7.0379957350095D+00-1.0934862352540D+01-2.5404077287132D+00-1.1015122918819D+00
         10 -1.1742669034618D-06-9.5547893492337D-08-3.3967626402033D-07-1.4766219169227D-07-1.9207608982961D-08-4.4548412225486D-08
         11 -1.5330127059637D-06 3.0173627944955D-07 1.6346584570666D-07-1.9344235755321D-07 3.5213138591283D-08 1.9330523416957D-08
         12  7.4076301055064D-10 1.3785248228602D-10 5.9772305305145D-11 9.2839606714421D-11 2.1566002307591D-11 9.3510383997953D-12
         13  3.7300661112284D-10 5.7834383720007D-11 2.5076194196940D-11 4.6776676784496D-11 9.4105133969805D-12 4.0803331379333D-12
 42344.500000  0.8627218791672550  0.4542648103862594  0.1969733753447424 -0.0088443237237526  0.0136339580087712  0.0059120269158903
          2  1.1811189553811D+00-3.6560232530596D-02-1.5860890809174D-02 1.3760274664577D-02-1.9002138737282D-03-8.2406845680601D-04
          3  6.8456051993101D-01-1.6512423830648D+00-7.1601835763654D-01 1.1568889127992D-02 1.1895910846028D-02 5.1582386745598D-03
          4  2.8641385193381D-06-1.9697365540945D-01 4.5425271303790D-01 8.4428075017824D-08-5.9120478208588D-03 1.3634127702134D-02
          5 -1.0613286463828D-01 1.6360777580168D-01 7.0944499760217D-02-7.5348490893226D-03 5.1780349562648D-03 2.2453548197106D-03
          6  1.9469114782705D-02 6.1007413907294D-03 2.1980114918693D-03 1.6584158495024D-03 5.6784193655240D-04 2.1134234843037D-04
          7 -6.5661117906769D-07-3.5965562715235D-09 6.8002823891924D-08-4.0963621292360D-08 2.7243092672640D-09 5.5389783784191D-09
          8  2.5209121949722D-09 5.3837030714895D-10 2.3343761178865D-10 2.1141820357367D-10 6.1027698839661D-11 2.6461676283098D-11
          9 -1.9696071161227D+02-4.6157940584882D+01-2.0013937707379D+01-1.6514782773406D+01-5.1102088893968D+00-2.2157996166427D+00
         10 -2.6597528164089D-06-3.5308374659935D-07-8.0291418350057D-07-2.2399681255850D-07-4.7675951696029D-08-7.1904705199353D-08
         11 -3.4850734211647D-06 6.2542314532996D-07 3.4418283192795D-07-2.9512665501421D-07 4.4332541411489D-08 2.5237953722152D-08
         12  1.6722513343818D-09 3.9180333493530D-10 1.6988634501684D-10 1.4014026309518D-10 4.3352012665960D-11 1.8797678467737D-11
         13  8.4247456434908D-10 1.7124530485226D-10 7.4250991876006D-11 7.0644163002016D-11 1.9664208204947D-11 8.5264476064522D-12
 42352.500000  0.7838897101435031  0.5586127432389472  0.2422213278645233 -0.0108341416922834  0.0124097300538012  0.0053811873800732
          2  1.3039325203950D+00-3.7058805640120D-02-1.6077936434263D-02 1.6803337762100D-02 1.9088649345457D-03 8.2760716980159D-04
          3  7.6182046473508D-01-1.5454317312238D+00-6.7013717028748D-01 7.6154987534458D-03 1.4390044240506D-02 6.2397709777681D-03
          4  3.5044547729452D-06-2.4222179135599D-01 5.5860215227083D-01 7.5324808512352D-08-5.3812120232665D-03 1.2409935866108D-02
1 EMBARY  2442304 TO  2442347 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   18
 JED-2400000         X*10**0             Y*10**0             Z*10**0           DX/DT*10**0         DY/DT*10**0         DZ/DT*10**0 
          5 -1.7334796277439D-01 1.9855642703146D-01 8.6099414494348D-02-9.2124730725976D-03 3.5002656284956D-03 1.5178496623955D-03
          6  3.5172609964252D-02 1.1807760986215D-02 4.3604425937620D-03 2.2730132618120D-03 8.7381097521533D-04 3.3608899502850D-04
          7 -9.6838038054757D-07-6.4299963409868D-08 8.6273356949130D-08-5.0117052692542D-08-2.2627072167278D-08-2.8436257010327D-09
          8  4.4999957209448D-09 1.1994330165543D-09 5.2007817728482D-10 2.8328889978963D-10 1.0690038040876D-10 4.6353030020194D-11
          9 -3.5156470616251D+02-1.0100753211813D+02-4.3796984492797D+01-2.2132576844404D+01-8.8107183092428D+00-3.8203846358500D+00
         10 -4.7611829503167D-06-9.0311467929648D-07-1.5016673267821D-06-3.0138717381136D-07-9.2938137387719D-08-1.0358731671067D-07
         11 -6.2610833531432D-06 9.8478919388879D-07 5.5560577901189D-07-3.9899792308465D-07 4.3625791187264D-08 2.6780698276474D-08
         12  2.9838315071896D-09 8.5695106083296D-10 3.7157787993729D-10 1.8771651742297D-10 7.4697088528587D-11 3.2389317337572D-11
         13  1.5036635237689D-09 3.8529522148764D-10 1.6706392023293D-10 9.4626439590575D-11 3.4734042483737D-11 1.5060865519183D-11
-ADAMS-MOULTON NUMERICAL INTEGRATION HAS BEEN COMPLETED IN      28 STEPS FOR 5.6000000000000D+01 DAYS OR 5.00000D-01 STEPS PER DAY
                 REAL TIME ELAPSED=  0H  0M  0.01S  OR   0.0004S  PER STEP  OR   0.0002S PER DAY
                 TASK TIME ELAPSED=  0H  0M  0.00S  OR   0.0000S  PER STEP  OR   0.0000S PER DAY
           (TASK TIME)/(REAL TIME)= 0.00000
0    3 QUANTITIES TO BE   CALCULATED   FOR BODY 5:   1   7   8
                              CODES FOR THE ABOVE:   0 506   5
1  MARS   2442321 TO  2442336 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   19
0  MARS      NPLNT( 1)=  4    JPLNT= 14          JDP1= 2442321        JDP0= 2442321        JDP2= 2442336     INT =  1   NCENTR= -1
 INITIAL EPOCH (COORD.TIME) JED=2442320.50000000 IHR= 0 IMIN= 0 SEC= 0.0000 INT1=         0 INT2=  0 FRACT= 0.00000000000000D+00
         A= 1.523705402995026D+00         E= 9.330246794792663D-02       INC= 2.469315159246902D+01       ASC= 3.344609353024999D+00
       PER= 3.322243666631956D+02      ANOM= 2.255191494289617D+02    RADIUS= 3.392459000000000D+03      FLAT= 0.000000000000000D+00
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 1.514460000000000D+02
   CON( 7)= 1.025956000000000D+00   CON( 8)= 5.269510000000000D+01   CON( 9)= 3.173116000000000D+02   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00  CON1( 1)= 2.443509500000000D+06  CON1( 2)= 0.000000000000000D+00
  CON1( 3)= 0.000000000000000D+00  CON1( 4)= 0.000000000000000D+00  CON1( 5)= 0.000000000000000D+00  CON1( 6)= 0.000000000000000D+00
  CON1( 7)= 0.000000000000000D+00  CON1( 8)= 0.000000000000000D+00  CON1( 9)= 0.000000000000000D+00  CON1(10)=-8.400000000000000D-11
  CON1(11)= 0.000000000000000D+00  CON1(12)= 0.000000000000000D+00  TCON( 1)= 0.000000000000000D+00  TCON( 2)= 0.000000000000000D+00
  TCON( 3)= 0.000000000000000D+00  TCON( 4)= 0.000000000000000D+00  TCON( 5)= 0.000000000000000D+00  TCON( 6)= 0.000000000000000D+00
  TCON( 7)= 0.000000000000000D+00  TCON( 8)= 0.000000000000000D+00  TCON( 9)= 0.000000000000000D+00  TCON(10)= 0.000000000000000D+00
  TCON(11)= 0.000000000000000D+00  TCON(12)= 0.000000000000000D+00  TCON(13)= 0.000000000000000D+00  TCON(14)= 0.000000000000000D+00
  TCON(15)= 0.000000000000000D+00  TCON(16)= 0.000000000000000D+00  TCON(17)= 0.000000000000000D+00  TCON(18)= 0.000000000000000D+00
  TCON(19)= 0.000000000000000D+00  TCON(20)= 0.000000000000000D+00  TCON(21)= 0.000000000000000D+00  TCON(22)= 0.000000000000000D+00
  TCON(23)= 0.000000000000000D+00  TCON(24)= 0.000000000000000D+00  TCON(25)= 0.000000000000000D+00  TCON(26)= 0.000000000000000D+00
  TCON(27)= 0.000000000000000D+00  TCON(28)= 0.000000000000000D+00  TCON(29)= 0.000000000000000D+00  TCON(30)= 0.000000000000000D+00
  KP( 1)=   5  KP( 2)=   0  KP( 3)=   0  KP( 4)=   0  KP( 5)=   0  KP( 6)=   0  KP( 7)=   0  KP( 8)=   0  KP( 9)=   0  KP(10)=   0
  KP(11)=   0  KP(12)=   0  KP(13)=   0  KP(14)=   0  KP(15)=   0  KP(16)=   0  KP(17)=   0  KP(18)=   0  KP(19)=   0  KP(20)=   0
  KP(21)=   0  KP(22)=   0  KP(23)=   0  KP(24)=   0  KP(25)=   0  KP(26)=   0  KP(27)=   0  KP(28)=   0  KP(29)=   0  KP(30)=   0
  KP(31)=  -1  KP(32)=  -1  KP(33)=   1  KP(34)=   1  KP(35)=   1  KP(36)=  -1  KP(37)=  -1  KP(38)=  -1  KP(39)=  -1  KP(40)=  -1
  KP(41)=   1  KP(42)=  -1  KP(43)=  -1  KP(44)=  -1  KP(45)=  -1  KP(46)=  -1  KP(47)=  -1  KP(48)=  -1  KP(49)=  -1  KP(50)=  -1
  KP(51)=  -1  KP(52)=  -1  KP(53)=  -1  KP(54)=  -1  KP(55)=  -1  KP(56)=  -1  KP(57)=  -1  KP(58)=  -1  KP(59)=  -1  KP(60)=  -1
  KP(61)=   1  KP(62)=  -1  KP(63)=  -1  KP(64)=  -1  KP(65)=  -1  KP(66)=  -1  KP(67)=  -1  KP(68)=  -1  KP(69)=  -1  KP(70)=   1
  KP(71)=   1  KP(72)=  -1  KP(73)=  -1  KP(74)=  -1  KP(75)=  -1  KP(76)=  -1  KP(77)=  -1  KP(78)=  -1  KP(79)=   1  KP(80)=   1
  KP(81)=  -1  KP(82)=  -1  KP(83)=  -1  KP(84)=  -1  KP(85)=  -1  KP(86)=  -1  KP(87)=   1  KP(88)=   3  KP(89)=  11  KP(90)=   6
  KP(91)= -10  KP(92)= -14  KP(93)=   0  KP(94)=   0  KP(95)=   0  KP(96)=   0  KP(97)=   0  KP(98)=   4  KP(99)=   0  KP(**)=   0
  KKP( 1)= -1  KKP( 2)=  0  KKP( 3)=  0  KKP( 4)=  0  KKP( 5)=  0  KKP( 6)=  0  KKP( 7)=  0  KKP( 8)=  0  KKP( 9)=  0  KKP(10)=  0
  KKP(11)=  0  KKP(12)=  0  KKP(13)=  0  KKP(14)=  0  KKP(15)=  0  KKP(16)=  0  KKP(17)=  0  KKP(18)=  0  KKP(19)=  0  KKP(20)=  0
  KKP(21)=  0  KKP(22)=  0  KKP(23)=  0  KKP(24)=  0  KKP(25)=  0  KKP(26)=  0  KKP(27)=  0  KKP(28)=  0  KKP(29)=  0  KKP(30)=  0
  KKP(31)=  0  KKP(32)=  0  KKP(33)=  0  KKP(34)=  0  KKP(35)=  0  KKP(36)=  0  KKP(37)=  0  KKP(38)=  0  KKP(39)=  0  KKP(40)=  0
  KKP(41)=  0  KKP(42)=  0  KKP(43)=  0  KKP(44)=  0  KKP(45)=  0  KKP(46)=  0  KKP(47)=  0  KKP(48)=  0  KKP(49)=  0  KKP(50)=  0
  KKP(51)=  0  KKP(52)=  0  KKP(53)=  0  KKP(54)=  0  KKP(55)=  0  KKP(56)=  0  KKP(57)=  0  KKP(58)=  0  KKP(59)=  0  KKP(60)=  0
  KKP(61)=  0  KKP(62)=  0  KKP(63)=  0  KKP(64)=  0  KKP(65)=  0  KKP(66)=  0  KKP(67)=  0  KKP(68)=  0  KKP(69)=  0  KKP(70)=  0
  KKP(71)=  0  KKP(72)=  0  KKP(73)=  0  KKP(74)=  0  KKP(75)=  0  KKP(76)=  0  KKP(77)=  0  KKP(78)=  0  KKP(79)=  0  KKP(80)=  0
  KKP(81)=  0  KKP(82)=  0  KKP(83)=  0  KKP(84)=  0  KKP(85)=  0  KKP(86)=  0  KKP(87)=  0  KKP(88)=  0  KKP(89)=  0  KKP(90)=  0
  KKP(91)=  0  KKP(92)=  0  KKP(93)=  0  KKP(94)=  0  KKP(95)=  0  KKP(96)=  0  KKP(97)=  0  KKP(98)=  0  KKP(99)=  0  KKP(**)=  0
 NUMKI= 16  KI= 1-1-1-1 1 1 1   3   5  11  31  40  41  49  50 506
 EPSP( 1)= 0.00000E+00 EPSP( 2)= 0.00000E+00 EPSP( 3)= 1.00000E-16 EPSP( 4)= 0.00000E+00 EPSP( 5)= 0.00000E+00 EPSP( 6)= 0.00000E+00
    ICND =  0   IPARP = 13   NUMERICAL INTEGRATION METHOD IS  ENCKE   ROYAL ROAD (2ND SUM)
0  SETUP FOR   MARS   NUMERICAL INTEGRATION  REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
  NTRG(1)=  5   NTZONE(1)=  0   NTTESS(1)=  0  HSITB(1)=  0.000000   FOR PARTIALS ONLY:   NTZONE(1)=  0   NTTESS(1)=  0
1  MARS   2442321 TO  2442336 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   20
0 STARTING PROCEDURE OF    143 STEPS REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1  MARS   2442321 TO  2442336 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   21
 JED-2400000         X*10**0             Y*10**0             Z*10**0           DX/DT*10**0         DY/DT*10**0         DZ/DT*10**0 
 42320.500000 -1.5822790612839730 -0.3680747075174319 -0.1265083724103855  0.0038721247506113 -0.0112388165240628 -0.0052627203983654
          2  3.6807470751743D-01-1.5822790612840D+00 0.0000000000000D+00 1.1238816524063D-02 3.8721247506113D-03 0.0000000000000D+00
          3  3.8717718912453D-01-1.4345091509152D+00-6.6885551993228D-01 1.2405919262533D-02 3.6463164068696D-03 1.3409387941729D-03
          4  4.2336927849225D-01-1.2288265356502D+00-5.7541383128492D-01 1.1833061959347D-02 2.7526439085834D-03 9.4609190359474D-04
          5  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          6  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          7  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          8  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          9  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         10  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         11  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         12  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         13  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
 42324.500000 -1.5659259680807351 -0.4128199425210292 -0.1474860392320887  0.0043041138685100 -0.0111316452905353 -0.0052251152438261
          2  4.1281996741235D-01-1.5659259156735D+00-1.0601054388442D-09 1.1131657739525D-02 4.3041400504074D-03-5.2850128874780D-10
          3  4.3657953433821D-01-1.4191400882000D+00-6.6312604826968D-01 1.2292878665223D-02 4.0379107932506D-03 1.5237211249922D-03
          4  4.7060418989342D-01-1.2171094280184D+00-5.7130253969754D-01 1.1782548704480D-02 3.1061978626947D-03 1.1097355620250D-03
          5 -2.0058009348090D-03-2.1430654068959D-04-9.7166967471501D-05-1.0018575453859D-03-1.1837145510262D-04-5.3425937800554D-05
          6 -3.5763405312967D-05 1.5223263510088D-05 7.4730706718349D-06-1.7862617287381D-05 7.6143076805270D-06 3.7382645899096D-06
          7 -1.3348078230715D-04 4.4440544995716D-05 5.5582767393427D-05-6.6757535337500D-05 2.2136236680886D-05 2.7750619465632D-05
          8 -1.6489661645029D-11-2.8120851567248D-12-8.4844576956950D-13-8.2473232084968D-12-1.4439969214408D-12-4.4157795521602D-13
          9 -9.4770015620307D-08 2.0321382138315D-08 1.0857362644652D-08-4.7402194936599D-08 1.0077781625793D-08 5.3934357376586D-09
         10 -1.0490561574628D-11-2.5480598413180D-12-8.8819678102528D-13-5.2442776479034D-12-1.3008872002612D-12-4.5645336591307D-13
         11  1.6158462451757D-13 3.3454243978917D-14 2.7037055848342D-14 8.0622617230798D-14 1.7128008679896D-14 1.3634495738949D-14
         12 -1.1668953641959D-04-2.4804874676257D-05-1.8036864327964D-05-5.8235582050269D-05-1.2691779371277D-05-9.1115493102830D-06
         13 -3.3223212385595D-08-4.4645926352233D-08-1.8357949881668D-08-1.6658530847162D-08-2.2308459289955D-08-9.1719943285373D-09
 42328.500000 -1.5478488862333087 -0.4571105833789782 -0.1683012861097851  0.0047340507213492 -0.0110115094086276 -0.0051815039829011
          2  4.5711068303489D-01-1.5478486769734D+00-4.2086531178255D-09 1.1011534354811D-02 4.7341029389013D-03-1.0423745618794D-09
          3  4.8550117294931D-01-1.4022085456434D+00-6.5666647069202D-01 1.2165554160106D-02 4.4274940655067D-03 1.7059634886902D-03
          4  5.1761472218644D-01-1.2039747555195D+00-5.6653453727918D-01 1.1720825652988D-02 3.4613991913698D-03 1.2744361105147D-03
          5 -8.0043223513897D-03-1.0370527577833D-03-4.6625085263609D-04-1.9958093550258D-03-3.0429904566536D-04-1.3598262699928D-04
          6 -1.4276011722697D-04 6.0922093335853D-05 2.9914156908500D-05-3.5619437765501D-05 1.5234240961394D-05 7.4825091807944D-06
          7 -5.3425257492724D-04 1.7637266545441D-04 2.2164427995698D-04-1.3365868878007D-04 4.3734184554396D-05 5.5230532545604D-05
          8 -6.6006290099556D-11-1.1862062678492D-11-3.6742177823160D-12-1.6515398623426D-11-3.1206326940520D-12-9.8938589521616D-13
          9 -3.7939534387978D-07 7.9924728788773D-08 4.2851554529414D-08-9.4937415998422D-08 1.9632122723063D-08 1.0564777027608D-08
         10 -4.1950715870434D-11-1.0625493605328D-11-3.7519694019896D-12-1.0485914016481D-11-2.7656177276173D-12-9.8818518160751D-13
         11  6.4367800477535D-13 1.4024623928575D-13 1.0999969686396D-13 1.6026793311495D-13 3.6673231174933D-14 2.7962099007984D-14
         12 -4.6505133483311D-04-1.0386320906245D-04-7.3636506972368D-05-1.1584608862505D-04-2.7130721748027D-05-1.8781402990334D-05
         13 -1.3368632693889D-07-1.7834797455885D-07-7.3317972359250D-08-3.3630801321890D-08-4.4527414694704D-08-1.8300622166684D-08
 42332.500000 -1.5280568519317041 -0.5008946550826555 -0.1889300129912202  0.0051615200037638 -0.0108783517156513 -0.0051318492351652
          2  5.0089487968788D-01-1.5280563820160D+00-9.3626997668124D-09 1.0878389267542D-02 5.1615980832880D-03-1.5288245202608D-09
          3  5.3388483214592D-01-1.3837233210813D+00-6.4947928920471D-01 1.2023878287809D-02 4.8146879644652D-03 1.8874940660642D-03
          4  5.6435546846186D-01-1.1894162656684D+00-5.6110573593601D-01 1.1647607288905D-02 3.8180765436393D-03 1.4401225545493D-03
          5 -1.7957288610572D-02-2.7393186901603D-03-1.2240291734608D-03-2.9785189812037D-03-5.5818662447311D-04-2.4778840532074D-04
          6 -3.2059664474305D-04 1.3707463123911D-04 6.7328241489726D-05-5.3284364298337D-05 2.2837175030933D-05 1.1223034985085D-05
          7 -1.2030357053396D-03 3.9349696832254D-04 4.9699179341699D-04-2.0077324697721D-04 6.4719035468573D-05 8.2383846403273D-05
          8 -1.4865531089590D-10-2.8102194058085D-11-8.9115982602680D-12-2.4815285790279D-11-5.0411273580638D-12-1.6482774118878D-12
          9 -8.5451872959187D-07 1.7660418176875D-07 9.5047697495042D-08-1.4265932588211D-07 2.8605558219590D-08 1.5490073252497D-08
         10 -9.4383037134514D-11-2.4899773298847D-11-8.8976327724841D-12-1.5731406176964D-11-4.4004737659070D-12-1.5979048259243D-12
         11  1.4425312071970D-12 3.3010885864059D-13 2.5165823878726D-13 2.3901542517756D-13 5.8670617168084D-14 4.2983308147548D-14
         12 -1.0427026359020D-03-2.4422022179402D-04-1.6903672775249D-04-1.7288975952508D-04-4.3347278457616D-05-2.9012958804351D-05
         13 -3.0277791757308D-07-4.0074594822269D-07-1.6470411221366D-07-5.0984441889028D-08-6.6657194783314D-08-2.7385271503338D-08
 42336.500000 -1.5065605837107740 -0.5441199669367453 -0.2093479756866417  0.0055860957900605 -0.0107321216389837 -0.0050761163603122
          2  5.4412036724010D-01-1.5065597500733D+00-1.6380936575537D-08 1.0732171970975D-02 5.5861995378873D-03-1.9714072815013D-09
          3  5.8167298331530D-01-1.3636947444307D+00-6.4156770220817D-01 1.1867790725318D-02 5.1991046568763D-03 2.0681364517746D-03
          4  6.1077986030062D-01-1.1734284263923D+00-5.5501234710039D-01 1.1562598478904D-02 4.1760408442509D-03 1.6067160217970D-03
1  MARS   2442321 TO  2442336 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   22
 JED-2400000         X*10**0             Y*10**0             Z*10**0           DX/DT*10**0         DY/DT*10**0         DZ/DT*10**0 
          5 -3.1812965539905D-02-5.5935103825093D-03-2.4876450662589D-03-3.9465957360670D-03-8.8029272772988D-04-3.8890758868679D-04
          6 -5.6892345526805D-04 2.4356357854602D-04 1.1967892178520D-04-7.0865397176128D-05 3.0398114517970D-05 1.4948891936699D-05
          7 -2.1407909035214D-03 6.9319580621718D-04 8.8019832942277D-04-2.6815136493653D-04 8.5006788637239D-05 1.0914947442766D-04
          8 -2.6458445362889D-10-5.2534022771772D-11-1.7016356201416D-11-3.3157039778979D-11-7.2189009382133D-12-2.4241428796046D-12
          9 -1.5209769164724D-06 3.0790760766692D-07 1.6640695360642D-07-1.9061107209511D-07 3.6932566632486D-08 2.0141358839436D-08
         10 -1.6781512549545D-10-4.6066301610179D-11-1.6643651284163D-11-2.0986746999224D-11-6.2131402410924D-12-2.2889790900653D-12
         11  2.5547040802748D-12 6.1295010808975D-13 4.5480371627700D-13 3.1693969747888D-13 8.3172791137980D-14 5.8707997575952D-14
         12 -1.8474863656504D-03-4.5306898004397D-04-3.0650224061016D-04-2.2942097172908D-04-6.1384926805237D-05-3.9816285979111D-05
         13 -5.4216632390726D-07-7.1149849289458D-07-2.9234580408143D-07-6.8791543953908D-08-8.8707015051677D-08-3.6429310153400D-08
 42340.500000 -1.4833725256189561 -0.5867341402483472 -0.2295307974208064  0.0060073411946903 -0.0105727758686382 -0.0050142737577736
          2  5.8673476783601D-01-1.4833712259949D+00-2.5048430521648D-08 1.0572839225738D-02 6.0074704053714D-03-2.3497015011217D-09
          3  6.2880787204152D-01-1.3421347165941D+00-6.3293562316071D-01 1.1697239021098D-02 5.5803464350537D-03 2.2477094967525D-03
          4  6.5684013064245D-01-1.1560064995900D+00-5.4825091449969D-01 1.1465495159465D-02 4.5350839572909D-03 1.7741291302489D-03
          5 -4.9505969358778D-02-9.8727722816152D-03-4.3743897872529D-03-4.8966077053522D-03-1.2707276677947D-03-5.5934988107041D-04
          6 -8.8741048208600D-04 3.8016795142554D-04 1.8688394867479D-04-8.8364200819788D-05 3.7890321514850D-05 1.8648155164505D-05
          7 -3.3486377400825D-03 1.0725004450410D-03 1.3695836091970D-03-3.3582206899756D-04 1.0450624968129D-04 1.3546200509418D-04
          8 -4.1397927888839D-10-8.6217154947444D-11-2.8469897874294D-11-4.1549482060580D-11-9.6695617257125D-12-3.3239102226940D-12
          9 -2.3797563443251D-06 4.7110625929237D-07 2.5577126461513D-07-2.3882417349741D-07 4.4540532336055D-08 2.4487061698807D-08
         10 -2.6229730950363D-10-7.4854039661830D-11-2.7323283693180D-11-2.6257283656377D-11-8.2127120729812D-12-3.0654393019661D-12
         11  3.9770419958784D-12 9.9892384810207D-13 7.2228551583216D-13 3.9410858537526D-13 1.1025004591611D-13 7.5155474175801D-14
         12 -2.8774529357501D-03-7.3780250653624D-04-4.8835145678283D-04-2.8548885580015D-04-8.1300159032263D-05-5.1208168280570D-05
         13 -8.5381548983389D-07-1.1103200820603D-06-4.5609450638972D-07-8.7127638867484D-08-1.1069568700111D-07-4.5440438939631D-08
-  ROYAL ROAD  NUMERICAL INTEGRATION HAS BEEN COMPLETED IN      20 STEPS FOR 2.0000000000000D+01 DAYS OR 1.00000D+00 STEPS PER DAY
                 REAL TIME ELAPSED=  0H  0M  0.01S  OR   0.0005S  PER STEP  OR   0.0005S PER DAY
                 TASK TIME ELAPSED=  0H  0M  0.00S  OR   0.0000S  PER STEP  OR   0.0000S PER DAY
           (TASK TIME)/(REAL TIME)= 0.00000
0  ALL  2 PLANET NUMERICAL INTEGRATIONS DURING THE LAST   11 PAGES STARTING ON PAGE   12 REQUIRED
           0H  0M  0.07S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1  MOON  INTEGRATION  (ITERAT  1)         TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   23
-
-
-
-

                                   **            **    *********     *********    **           **
                                   ***          ***   ***********   ***********   ***          **
                                   ****        ****   **       **   **       **   ****         **
                                   ** **      ** **   **       **   **       **   ** **        **
                                   **  **    **  **   **       **   **       **   **  **       **
                                   **   **  **   **   **       **   **       **   **   **      **
                                   **    ****    **   **       **   **       **   **    **     **
                                   **     **     **   **       **   **       **   **     **    **
                                   **            **   **       **   **       **   **      **   **
                                   **            **   **       **   **       **   **       **  **
                                   **            **   **       **   **       **   **        ** **
                                   **            **   **       **   **       **   **         ****
                                   **            **   ***********   ***********   **          ***
                                   **            **    *********     *********    **           **
-
-
                                          **           **    **           **   **       **
                                          **           **    ***          **   **      ** 
                                          **           **    ****         **   **     **  
                                          **           **    ** **        **   **    **   
                                          **           **    **  **       **   **   **    
                                          **           **    **   **      **   *** **     
                                          **           **    **    **     **   *****      
                                          **           **    **     **    **   ** **      
                                          **           **    **      **   **   **  **     
                                          **           **    **       **  **   **   **    
                                          **           **    **        ** **   **    **   
                                          **           **    **         ****   **     **  
                                          ***********  **    **          ***   **      ** 
                                          ***********  **    **           **   **       **
- PEP VERSION= 20210302 PEP.PEPLOAD.PEP790                          
1  MOON  INTEGRATION  (ITERAT  1)         TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   24
0   INFORMATION FROM FIRST TWO RECORDS OF PERTURBING PLANET PERIPHERAL DATA SET 90
  N-BODY RUN 443-0 9-BODY INTEGRATION  442D8 INNER,  384F4 OUTER  I.C.            3/20/73 PERTURBING PLNT TAPE FROM N-BODY INTEG.
 IDIR= 1 JDT1= 2441841 JDT2= 2442400 NMOON= 0 INITIAL CONDITIONS ARE
0MERCURY RUN 443- 3/21/73 JD0=2440000.5  A = 3.870988248079130D-01  E = 2.056143416457730D-01   INC= 2.860343745748160D+01
                                        ASC= 1.085941287077680D+01 PER= 6.693261684748430D+01 ANOM0= 9.083965359314600D+01
0 VENUS  RUN 443- 3/21/73 JD0=2440000.5  A = 7.233280763362320D-01  E = 6.757538488379489D-03   INC= 2.446691434342910D+01
                                        ASC= 7.978123733702580D+00 PER= 1.237851072721730D+02 ANOM0= 2.744064790468620D+02
0 EMBARY RUN 443- 3/21/73 JD0=2440000.5  A = 9.999827653637140D-01  E = 1.675271216169400D-02   INC= 2.344335783108840D+01
                                        ASC= 6.659250140529500D-04 PER= 1.020340640186800D+02 ANOM0= 1.393658182548250D+02
0  MARS  RUN 443- 3/21/73 JD0=2440000.5  A = 1.523604534884400D+00  E = 9.345778628544589D-02   INC= 2.469301678564210D+01
                                        ASC= 3.346655760227140D+00 PER= 3.321011352392630D+02 ANOM0= 8.989975217532620D+01
0JUPITER RUN 443- 3/20/73 JD0=2440000.5  A = 5.202945960488120D+00  E = 4.819167141557811D-02   INC= 2.325294770783580D+01
                                        ASC= 3.260973418885740D+00 PER= 1.037371340811626D+01 ANOM0= 1.411862510204257D+02
0 SATURN RUN 443- 3/20/73 JD0=2440000.5  A = 9.526046869249818D+00  E = 5.461675768589707D-02   INC= 2.257317229472867D+01
                                        ASC= 5.968803885438612D+00 PER= 8.759383460863923D+01 ANOM0= 2.897400092498020D+02
0 URANUS RUN 443- 3/20/73 JD0=2440000.5  A = 1.927444351332630D+01  E = 5.124121845145173D-02   INC= 2.367141284461453D+01
                                        ASC= 1.847313091106734D+00 PER= 1.683442408256917D+02 ANOM0= 6.990628140027202D+00
0NEPTUNE RUN 443- 3/20/73 JD0=2440000.5  A = 3.011376065452675D+01  E = 6.984876998530826D-03   INC= 2.231473025943044D+01
                                        ASC= 3.519876825815031D+00 PER= 5.409517825004429D+01 ANOM0= 1.775740764613225D+02
0 PLUTO  RUN 443- 3/20/73 JD0=2440000.5  A = 3.974605251930742D+01  E = 2.522316558278015D-01   INC= 2.365507811729088D+01
                                        ASC= 4.377229536665621D+01 PER= 1.819638797966960D+02 ANOM0= 3.298440994501460D+02
0  MOON  RUN 440- 8/ 7/72 JD0=2440000.5  A = 2.571514374642509D-03  E = 5.561544735761867D-02   INC= 2.839685438634148D+01
                                        ASC= 3.312959209874043D+00 PER= 2.262712472223586D+02 ANOM0= 1.548856956849208D+02
1  MOON  INTEGRATION  (ITERAT  1)         TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   25
0DATA ON FIRST TWO RECORDS OF DATA SET 13   BODY=  EMBARY    NPLNT(-3)=  3    NCENTR= -1     IPAR= 13     IFILTR= 0
 TITLE=TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE=   14 ITERAT= 1 LEVEL= 790
     JD1=2442304          JD0=2442321          JD2=2442347          INT=  2     ICND=  0
   COND(1)= 1.000002347957079D+00   COND(2)= 1.669041814744669D-02   COND(3)= 2.344258335876651D+01   COND(4)= 8.391300825564327D-04
   COND(5)= 1.022643964965586D+02   COND(6)= 2.657488867161210D+02   CON( 1)= 6.378166000000000D+03   CON( 2)= 3.352329869259135D-03
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 0.000000000000000D+00
   CON( 7)= 0.000000000000000D+00   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00   CON(25)= 2.442320500000000D+06   CON(26)= 0.000000000000000D+00
   CON(27)= 0.000000000000000D+00   CON(28)= 0.000000000000000D+00   CON(29)= 0.000000000000000D+00   CON(30)= 0.000000000000000D+00
   CON(31)= 0.000000000000000D+00   CON(32)= 0.000000000000000D+00   CON(33)= 0.000000000000000D+00   CON(34)=-4.600000000000000D-10
   CON(35)= 0.000000000000000D+00   CON(36)= 0.000000000000000D+00   CON(
  PRM(  1)= 1.657848020429993D-07  PRM(  2)= 2.447848585877872D-06  PRM(  3)= 3.040436898620584D-06  PRM(  4)= 3.227159776680543D-07
  PRM(  5)= 9.547534692876814D-04  PRM(  6)= 2.857796067672611D-04  PRM(  7)= 4.361098996947231D-05  PRM(  8)= 5.192107995846314D-05
  PRM(  9)= 2.500000000000000D-07  PRM( 10)= 1.215052064980984D-02  PRM( 11)= 6.700000000000000D-10  PRM( 12)= 0.000000000000000D+00
  PRM( 13)= 0.000000000000000D+00  PRM( 14)= 0.000000000000000D+00  PRM( 15)= 0.000000000000000D+00  PRM( 16)= 0.000000000000000D+00
  PRM( 17)= 0.000000000000000D+00  PRM( 18)= 0.000000000000000D+00  PRM( 19)= 0.000000000000000D+00  PRM( 20)= 0.000000000000000D+00
  PRM( 21)= 0.000000000000000D+00  PRM( 22)= 0.000000000000000D+00  PRM( 23)= 0.000000000000000D+00  PRM( 24)= 0.000000000000000D+00
  PRM( 25)= 0.000000000000000D+00  PRM( 26)= 0.000000000000000D+00  PRM( 27)= 0.000000000000000D+00  PRM( 28)= 0.000000000000000D+00
  PRM( 29)= 0.000000000000000D+00  PRM( 30)= 0.000000000000000D+00  PRM( 31)= 1.000000000000000D+00  PRM( 32)= 1.000000000000000D-35
  PRM( 33)= 1.000000000000000D-08  PRM( 34)= 0.000000000000000D+00  PRM( 35)= 0.000000000000000D+00  PRM( 36)= 0.000000000000000D+00
  PRM( 37)= 0.000000000000000D+00  PRM( 38)= 0.000000000000000D+00  PRM( 39)=-7.800000000000000D-04  PRM( 40)= 0.000000000000000D+00
  PRM( 41)= 1.000000000000000D+00  PRM( 42)= 1.000000000000000D+00  PRM( 43)= 1.000000000000000D+00  PRM( 44)= 1.000000000000000D+00
  PRM( 45)= 0.000000000000000D+00  PRM( 46)= 0.000000000000000D+00  PRM( 47)= 0.000000000000000D+00  PRM( 48)= 2.344330000000000D+01
  PRM( 49)= 2.900000000000000D+00  PRM( 50)= 1.000000000000000D-09  PRM( 51)= 4.990047800000000D+02  PRM( 52)= 0.000000000000000D+00
  PRM( 53)= 1.000000000000000D+00  PRM( 54)= 0.000000000000000D+00  PRM( 55)= 0.000000000000000D+00  PRM( 56)= 0.000000000000000D+00
  PRM( 57)= 0.000000000000000D+00  PRM( 58)= 0.000000000000000D+00  PRM( 59)= 0.000000000000000D+00  PRM( 60)= 1.000000000000000D+01
  PRM( 61)= 0.000000000000000D+00  PRM( 62)= 1.000000000000000D+00  PRM( 63)= 1.000000000000000D+00  PRM( 64)= 0.000000000000000D+00
  PRM( 65)= 0.000000000000000D+00  PRM( 66)= 0.000000000000000D+00  PRM( 67)= 0.000000000000000D+00  PRM( 68)= 0.000000000000000D+00
  PRM( 69)= 0.000000000000000D+00  PRM( 70)= 0.000000000000000D+00  PRM( 71)= 0.000000000000000D+00  PRM( 72)= 0.000000000000000D+00
  PRM( 73)= 0.000000000000000D+00  PRM( 74)= 0.000000000000000D+00  PRM( 75)= 0.000000000000000D+00  PRM( 76)= 0.000000000000000D+00
  PRM( 77)= 0.000000000000000D+00  PRM( 78)= 0.000000000000000D+00  PRM( 79)= 0.000000000000000D+00  PRM( 80)= 0.000000000000000D+00
  PRM( 81)= 1.000000000000000D+00  PRM( 82)= 0.000000000000000D+00  PRM( 83)= 0.000000000000000D+00  PRM( 84)= 0.000000000000000D+00
  PRM( 85)= 0.000000000000000D+00  PRM( 86)= 0.000000000000000D+00  PRM( 87)= 0.000000000000000D+00  PRM( 88)= 0.000000000000000D+00
  PRM( 89)= 0.000000000000000D+00  PRM( 90)= 0.000000000000000D+00  PRM( 91)= 2.344578706750000D+01  PRM( 92)= 7.250000000000000D+00
  PRM( 93)= 7.506250000000000D+01  PRM( 94)= 6.960000000000000D+05  PRM( 95)= 6.960000000000000D+05  PRM( 96)= 0.000000000000000D+00
  PRM( 97)= 2.440000500000000D+06  PRM( 98)= 4.263529034000000D-05  PRM( 99)= 0.000000000000000D+00  PRM(100)= 2.997924580000000D+05
  K(  1)=   0  K(  2)=   0  K(  3)=   0  K(  4)=   0  K(  5)=   0  K(  6)=   0  K(  7)=   0  K(  8)=   0  K(  9)=   0  K( 10)=   0
  K( 11)=   0  K( 12)=   0  K( 13)=   0  K( 14)=   0  K( 15)=   0  K( 16)=   0  K( 17)=   0  K( 18)=   0  K( 19)=   0  K( 20)=   0
  K( 21)=   0  K( 22)=   0  K( 23)=   0  K( 24)=   0  K( 25)=   0  K( 26)=   0  K( 27)=   0  K( 28)=   0  K( 29)=   0  K( 30)=   0
  K( 31)=  -1  K( 32)=  -1  K( 33)=   1  K( 34)=   1  K( 35)=   1  K( 36)=  -1  K( 37)=  -1  K( 38)=  -1  K( 39)=  -1  K( 40)=   1
  K( 41)=  -1  K( 42)=  -1  K( 43)=  -1  K( 44)=  -1  K( 45)=  -1  K( 46)=  -1  K( 47)=  -1  K( 48)=  -1  K( 49)=  -1  K( 50)=  -1
  K( 51)=  -1  K( 52)=  -1  K( 53)=  -1  K( 54)=  -1  K( 55)=  -1  K( 56)=  -1  K( 57)=  -1  K( 58)=  -1  K( 59)=  -1  K( 60)=  -1
  K( 61)=   1  K( 62)=   1  K( 63)=   1  K( 64)=  -1  K( 65)=  -1  K( 66)=  -1  K( 67)=  -1  K( 68)=  -1  K( 69)=  -1  K( 70)=   1
  K( 71)=   1  K( 72)=   1  K( 73)=  -1  K( 74)=  -1  K( 75)=  -1  K( 76)=  -1  K( 77)=  -1  K( 78)=  -1  K( 79)=   1  K( 80)=   1
  K( 81)=  -1  K( 82)=   0  K( 83)=  -1  K( 84)=  -1  K( 85)=  -1  K( 86)=  -1  K( 87)=   2  K( 88)=   2  K( 89)=   6  K( 90)=   6
  K( 91)=  -3  K( 92)=  -6  K( 93)=   0  K( 94)=   0  K( 95)=   0  K( 96)=   0  K( 97)=   0  K( 98)=   4  K( 99)=   0  K(100)=  -1
 NUMKI= 16  KI= 1 1 1 1-1-1-1    3    4   10   31   32   33   40   41   42
   EPS(1)= 0.00000E+00   EPS(2)= 0.00000E+00   EPS(3)= 1.00000E-09   EPS(4)= 0.00000E+00   EPS(5)= 0.00000E+00   EPS(6)= 0.00000E+00
0    5 QUANTITIES TO BE READ FROM TAPE FOR BODY 3:   1   4   8   9  10
                              CODES FOR THE ABOVE:   0 303   3  10  32
                        RELATIVE POSITION ON TAPE:   1   4   5   7   9
1  MOON   2442321 TO  2442332 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   26
0MOON        MPLNT    = 10    IMN  = 20          JDM1= 2442321       JDMN0= 2442321        JDM2= 2442332    INTMN= -1   MCENTR=  3
       A  = 2.584844478828888D-03       E  = 4.663456366998220D-02     INC  = 2.240083560206691D+01     ASC  = 3.473239233192494D+02
     PER  = 1.425397834496920D+02    ANOM  = 2.222476213123667D+02     MRAD = 1.738090000000000D+03  MCON( 2)= 0.000000000000000D+00
  MCON( 3)= 0.000000000000000D+00  MCON( 4)= 0.000000000000000D+00  MCON( 5)= 0.000000000000000D+00  MCON( 6)= 0.000000000000000D+00
  MCON( 7)= 0.000000000000000D+00  MCON( 8)= 0.000000000000000D+00  MCON( 9)= 0.000000000000000D+00  MCON(10)= 0.000000000000000D+00
  MCON(11)= 0.000000000000000D+00  MCON(12)= 0.000000000000000D+00  MCON(13)= 0.000000000000000D+00  MCON(14)= 0.000000000000000D+00
  MCON(15)= 0.000000000000000D+00  MCON(16)= 8.936456813092100D-02  MCON(17)= 0.000000000000000D+00  MCON(18)= 0.000000000000000D+00
  MCON(19)= 0.000000000000000D+00  MCON(20)= 0.000000000000000D+00  MCON(21)= 0.000000000000000D+00  MCON(22)= 0.000000000000000D+00
  MCON(23)= 0.000000000000000D+00  MCON(24)= 0.000000000000000D+00  MCON(
 MCON1( 1)= 0.000000000000000D+00 MCON1( 2)= 0.000000000000000D+00 MCON1( 3)= 0.000000000000000D+00 MCON1( 4)= 0.000000000000000D+00
 MCON1( 5)= 0.000000000000000D+00 MCON1( 6)= 0.000000000000000D+00 MCON1( 7)= 0.000000000000000D+00 MCON1( 8)= 0.000000000000000D+00
 MCON1( 9)= 0.000000000000000D+00 MCON1(10)=-1.900000000000000D-11 MCON1(11)= 0.000000000000000D+00 MCON1(12)= 0.000000000000000D+00
  KM( 1)=   0  KM( 2)=   0  KM( 3)=   0  KM( 4)=   0  KM( 5)=   0  KM( 6)=   0  KM( 7)=   0  KM( 8)=   0  KM( 9)=   0  KM(10)=   0
  KM(11)=   0  KM(12)=   0  KM(13)=   0  KM(14)=   0  KM(15)=   0  KM(16)=   0  KM(17)=   0  KM(18)=   0  KM(19)=   0  KM(20)=   0
  KM(21)=   0  KM(22)=   0  KM(23)=   0  KM(24)=   0  KM(25)=   0  KM(26)=   0  KM(27)=   0  KM(28)=   0  KM(29)=   0  KM(30)=   0
  KM(31)=   1  KM(32)=   1  KM(33)=   1  KM(34)=   1  KM(35)=   1  KM(36)=   1  KM(37)=   1  KM(38)=   1  KM(39)=   1  KM(40)=   1
  KM(41)=  -1  KM(42)=  -1  KM(43)=  -1  KM(44)=  -1  KM(45)=  -1  KM(46)=  -1  KM(47)=  -1  KM(48)=  -1  KM(49)=  -1  KM(50)=  -1
  KM(51)=  -1  KM(52)=  -1  KM(53)=  -1  KM(54)=  -1  KM(55)=  -1  KM(56)=  -1  KM(57)=  -1  KM(58)=  -1  KM(59)=  -1  KM(60)=  -1
  KM(61)=   1  KM(62)=   1  KM(63)=  -1  KM(64)=  -1  KM(65)=  -1  KM(66)=  -1  KM(67)=  -1  KM(68)=  -1  KM(69)=  -1  KM(70)=  -1
  KM(71)=  -1  KM(72)=  -1  KM(73)=  -1  KM(74)=  -1  KM(75)=  -1  KM(76)=  -1  KM(77)=  -1  KM(78)=  -1  KM(79)=  -1  KM(80)=  -1
  KM(81)=   4  KM(82)=   3  KM(83)=   3  KM(84)=   1  KM(85)=   1  KM(86)=  -1  KM(87)=  -2  KM(88)=   2  KM(89)=   6  KM(90)=   6
  KM(91)=  -4  KM(92)=  -6  KM(93)=   0  KM(94)=   0  KM(95)=   0  KM(96)=   0  KM(97)=   0  KM(98)=  10  KM(99)=   0  KM(**)=  -1
  KKM( 1)=  3  KKM( 2)=  0  KKM( 3)=  0  KKM( 4)=  0  KKM( 5)=  0  KKM( 6)=  0  KKM( 7)=  0  KKM( 8)=  0  KKM( 9)=  0  KKM(10)=  0
  KKM(11)=  0  KKM(12)=  0  KKM(13)=  0  KKM(14)=  0  KKM(15)=  0  KKM(16)=  0  KKM(17)=  0  KKM(18)=  0  KKM(19)=  0  KKM(20)=  0
 NUMKI= 16  KI= 1 1-1-1-1 1 1  -16  -20    3   10   32  303 1031    2    3
 EPSM( 1)= 0.00000E+00 EPSM( 2)= 0.00000E+00 EPSM( 3)= 1.00000E-10 EPSM( 4)= 0.00000E+00 EPSM( 5)= 0.00000E+00 EPSM( 6)= 0.00000E+00
   IPARM = 12   NUMERICAL INTEGRATION METHOD IS  COWELL  ADAMS MOULTON       
1  MOON   2442321 TO  2442332 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   27
0MOON ROTATION        MPLNT    =-10    ILIB= 21          JDM1= 2442321       JDMR0= 2442321        JDM2= 2442332     INTMR= -1
   PSI = 6.406890562446280D-02     THETA = 4.174732965826930D-01       PHI = 4.427722151828226D+02
  DPSI =-1.116475592635000D-04    DTHETA =-3.914136528300000D-05      DPHI = 2.300859910558332D-01
     MRAD = 1738.090000000000         ALPHA= 4.039070136017979D-04     BETA = 6.317168491855110D-04     GAMMA= 2.278098937105160D-04
    MEQINC= 2.692029800000000D-02   CON( 6)= 2.415119936539230D-02   CON( 7)= 4.683189055748760D-03   CON( 8)= 0.000000000000000D+00
   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00
   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00
   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00
   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00
  CON1( 1)= 0.000000000000000D+00  CON1( 2)= 0.000000000000000D+00  CON1( 3)= 0.000000000000000D+00  CON1( 4)= 0.000000000000000D+00
  CON1( 5)= 0.000000000000000D+00  CON1( 6)= 0.000000000000000D+00  CON1( 7)= 0.000000000000000D+00  CON1( 8)= 0.000000000000000D+00
  CON1( 9)= 0.000000000000000D+00  CON1(10)= 0.000000000000000D+00  CON1(11)= 0.000000000000000D+00  CON1(12)= 0.000000000000000D+00
  KMR( 1)=  0  KMR( 2)=  0  KMR( 3)=  0  KMR( 4)=  0  KMR( 5)=  0  KMR( 6)=  0  KMR( 7)=  0  KMR( 8)=  0  KMR( 9)=  0  KMR(10)=  0
  KMR(11)=  0  KMR(12)=  0  KMR(13)=  0  KMR(14)=  0  KMR(15)=  0  KMR(16)=  0  KMR(17)=  0  KMR(18)=  0  KMR(19)=  0  KMR(20)=  0
  KMR(21)=  0  KMR(22)=  0  KMR(23)=  0  KMR(24)=  0  KMR(25)=  0  KMR(26)=  0  KMR(27)=  0  KMR(28)=  0  KMR(29)=  0  KMR(30)=  0
  KMR(31)= -1  KMR(32)= -1  KMR(33)=  1  KMR(34)= -1  KMR(35)= -1  KMR(36)= -1  KMR(37)= -1  KMR(38)= -1  KMR(39)= -1  KMR(40)=  1
  KMR(41)= -1  KMR(42)= -1  KMR(43)= -1  KMR(44)= -1  KMR(45)= -1  KMR(46)= -1  KMR(47)= -1  KMR(48)= -1  KMR(49)= -1  KMR(50)= -1
  KMR(51)= -1  KMR(52)= -1  KMR(53)= -1  KMR(54)= -1  KMR(55)= -1  KMR(56)= -1  KMR(57)= -1  KMR(58)= -1  KMR(59)= -1  KMR(60)= -1
  KMR(61)= -1  KMR(62)= -1  KMR(63)= -1  KMR(64)= -1  KMR(65)= -1  KMR(66)= -1  KMR(67)= -1  KMR(68)= -1  KMR(69)= -1  KMR(70)= -1
  KMR(71)= -1  KMR(72)= -1  KMR(73)= -1  KMR(74)= -1  KMR(75)= -1  KMR(76)= -1  KMR(77)= -1  KMR(78)= -1  KMR(79)= -1  KMR(80)= -1
  KMR(81)=  1  KMR(82)=  0  KMR(83)=  0  KMR(84)= -1  KMR(85)= -1  KMR(86)=303  KMR(87)= -2  KMR(88)=  2  KMR(89)=  6  KMR(90)=  6
  KMR(91)= -4  KMR(92)= -6  KMR(93)=  0  KMR(94)=  0  KMR(95)=  0  KMR(96)=  0  KMR(97)=  0  KMR(98)= 10  KMR(99)=  0  KMR(**)= -1
 KKMR( 1)=  0 KKMR( 2)=  0 KKMR( 3)=  0 KKMR( 4)=  0 KKMR( 5)=  0 KKMR( 6)=  0 KKMR( 7)=  0 KKMR( 8)=  0 KKMR( 9)=  0 KKMR(10)=  0
 KKMR(11)=  0 KKMR(12)=  0 KKMR(13)=  0 KKMR(14)=  0 KKMR(15)=  0 KKMR(16)=  0 KKMR(17)=  0 KKMR(18)=  0 KKMR(19)=  0 KKMR(20)=  0
 NUMKI= 18  KI= 1 1 1-1-1-1-1   -3   -6   -7 1031    2    3 1041    3    1    3    2
 EPSM( 1)= 0.00000E+00 EPSM( 2)= 0.00000E+00 EPSM( 3)= 1.00000E-10 EPSM( 4)= 0.00000E+00 EPSM( 5)= 0.00000E+00 EPSM( 6)= 0.00000E+00
   IPARMR = 10   NUMERICAL INTEGRATION METHOD IS  COWELL  ADAMS MOULTON       

-FUNDAMENTAL QUANTITIES FOR MOON INTEGRATION:
 MASS(SUN/E+M)= 328900.1000000000       MASS(E/M)= 81.30100000000000
 AU= 499.0047800000000 LT SEC       MRAD= 1738.090000000000 KM       BETA= 6.317168491855110D-04
 GAMMA= 2.278098937105160D-04      MNZ2= 2.038220000000000D-04
-DERIVED QUANTITIES:
 AU= 1.495978695499492D+08 KM   GM(EARTH)= 398601.0091671799     KM**3/SEC**2   GM(MOON)= 4902.781136359699     KM**3/SEC**2  
 MOON C/(M*R**2)=0.3938156299706623D+00 MNC22 (UNSCALED) = 2.242877420128912D-05
0  SETUP FOR MOON ORBIT AND ROTATION INTEGRATION REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
0 STARTING PROCEDURE OF     48 STEPS REQUIRED
           0H  0M  0.04S REAL TIME
           0H  0M  0.03S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.75000
1  MOON   2442321 TO  2442332 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   28
    JED              X                   Y                   Z                 DX/DT               DY/DT               DZ/DT
                    PSI                THETA                PHI               DPSI/DT            DTHETA/DT            DPHI/DT
 2442320.500  0.0026232835202694 -0.0005304792874056  0.0000239484232369  0.0000854824499752  0.0005195905088138  0.0002166804074577
              0.0640689056244628  0.4174732965826930  2.9492436802515526 -0.0001116475592635 -0.0000391413652830  0.2300859910558332
          2  1.0148709300522D+00-2.0522677157194D-01 9.2649377682290D-03-1.6535317825761D-02-1.0050711233683D-01-4.1913625603478D-02
          3  4.8154565883489D-04 2.4273344719612D-03 1.0196859167141D-03-5.6094323911689D-04 9.7151799061714D-05-1.1668569205592D-05
          4  3.7452356561803D-04 2.2764776872760D-03 9.4934011395519D-04-5.3932365166979D-04 1.0906180144394D-04-4.9235818286741D-06
          5  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          6  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          7  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 1.4057593172544D+01 8.5446685153955D+01 3.5633103840439D+01
          8  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          9  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         10  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         11  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         12  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          2  1.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          3  0.0000000000000D+00 1.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          4  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          5  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          6  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          7  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          8  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
          9  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
         10  0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00 0.0000000000000D+00
 2442325.500  0.0015443080553091  0.0018548002229826  0.0008834068655862 -0.0004831670596414  0.0003290889493570  0.0000879154586221
              0.0631480623203880  0.4174873951385355  4.0999778330108183 -0.0001928197680788  0.0000564257707753  0.2301476741857979
          2  2.0247561788871D+00-2.3202423575946D-01 8.8583559405263D-02 3.3821380783124D-01 2.2541965214609D-01 1.2061399413071D-01
          3 -2.0428915795990D-03 1.4768218699890D-03 4.0922596340194D-04-3.4024375352111D-04-4.4949965030631D-04-2.1133583254087D-04
          4 -2.1616265473102D-03 1.4267086831763D-03 3.7845370072365D-04-3.6656478158080D-04-4.3938131570894D-04-2.0959110299698D-04
          5  9.1668400753131D-16 9.0694328448134D-15 4.3614532779591D-16 1.9526041735311D-16 3.8265024010081D-15 4.7561423668535D-16
          6  5.0985761060207D-03 1.0413403821552D-03 5.9546187359314D-04 2.1851831917571D-03 8.0638953671406D-04 4.2096031599406D-04
          7 -4.0448684592601D+02 2.6909860116473D+02 7.1747370806422D+01-1.4881947452234D+02-2.7802208612400D+01-2.4541982245756D+01
          8  1.6381465155968D-07 3.1040253117822D-08 5.5015049545259D-09 6.5787494788176D-08 1.4471729253293D-08 5.0096953416771D-09
          9 -4.0546328878207D+00-9.1141366497703D-01-7.3678337347347D-01-1.6796995856273D+00-8.7097064288388D-01-5.0252309798971D-01
         10  3.8949970236256D-07-1.2680660970724D-06 3.1874987968137D-06 1.7340677549847D-07-4.4191046746610D-07 1.2247305335926D-06
         11 -8.2335262842572D-08-1.2561985620699D-08-2.2935615833372D-08-3.4743671012567D-08-1.5778651394580D-08-1.2593823506652D-08
         12  1.0020348158473D-10 9.1265903675291D-11-1.4360163628301D-10 3.7224483078429D-11 4.4954988716311D-11-4.9881923105944D-11
          2  9.9904006462375D-01 1.5620820663769D-04 5.0903921759578D-04-2.9540278082775D-04 1.1721284344968D-04 1.2032450384613D-04
          3  1.3366958776380D-03 1.0001665962688D+00-1.3312040333991D-03 5.4639042399207D-04-1.3687833125755D-05-5.1996152995843D-04
          4 -3.6849264037058D-01 8.0584206045079D-02 3.3707673271262D-01-9.6698261088775D-02 4.8303492739412D-02 8.8470223777908D-02
          5  1.6415792692045D-06-8.7624298404405D-07-1.9945076234270D-06 7.9269689368391D-08-4.7465951642930D-07-2.9823122314839D-07
          6  1.2555060643603D-06-1.0420272059921D-08-1.4142729353123D-06 6.2535918231309D-07-1.3382685098167D-07-6.6191362854902D-07
          7  2.3075134915811D-03-4.3642823258817D-04-2.7041597411540D-03 9.1973689022716D-04-3.7135191767280D-04-1.0539784260705D-03
          8 -1.9760955869745D-02 4.5582001867925D-03 1.8063674770873D-02-6.2091420226182D-03 3.0736506042384D-03 5.6754719848647D-03
          9  2.2185018471987D-02-4.9207637372362D-03-1.9639043455012D-02 5.8328398334730D-03-2.9460099975133D-03-5.0352262469376D-03
         10 -2.0018327613081D-01 4.5772523935012D-02 1.8551423235167D-01-6.2502311287447D-02 3.0814548256698D-02 5.8110710433632D-02
 2442330.500 -0.0012817240845037  0.0019935422635079  0.0006812258494335 -0.0005128861265792 -0.0002933871017466 -0.0001642966285624
              0.0627250413864410  0.4179041636144537  5.2501817691363240  0.0000304183160269  0.0000811220614314  0.2299288169210390
          2  2.5537361598637D+00 2.5198927493075D+00 1.2367065872354D+00-3.4439932475104D-01 7.6450077112157D-01 2.7442851648910D-01
          3 -2.1375608037005D-03-1.2818988530888D-03-7.0610297006229D-04 3.3146841290226D-04-5.3586312395649D-04-1.8471951521946D-04
          4 -2.3357189414515D-03-1.4087234620983D-03-7.7440557093276D-04 3.4356183401719D-04-5.6242979151349D-04-1.9419254613269D-04
          5 -6.3409059850626D-15 5.2447707931135D-14 1.4892885144896D-14-7.3557589899973D-15 1.6038628756928D-14 6.5574322892967D-15
          6  1.5858045499297D-02 1.3850100018946D-02 6.2108759783586D-03 6.1758376778156D-04 3.9794025198103D-03 1.5804286450210D-03
          7 -8.6423038163179D+02-4.9555277840834D+02-2.7578339276442D+02 4.1369149033468D+01-2.4881540484544D+02-9.5782285022332D+01
          8  3.7404991794008D-07 2.2003859104976D-07 8.6015531651271D-08-1.7663146739162D-08 4.4447690925258D-08 1.9164006744191D-08
          9 -1.0442063051460D+01-1.6046709688962D+01-7.3763561454083D+00 8.5183177295655D-01-5.4173465220052D+00-2.0885878556311D+00
         10  9.9150006602019D-07-2.5242875946475D-06 8.3055204926650D-06-4.5574170564952D-08 9.1579929941300D-08 1.1475684712862D-07
1  MOON   2442321 TO  2442332 (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   29
    JED              X                   Y                   Z                 DX/DT               DY/DT               DZ/DT
                    PSI                THETA                PHI               DPSI/DT            DTHETA/DT            DPHI/DT
         11 -2.1889809646411D-07-3.1672716317230D-07-1.6102045346680D-07 1.5672384271532D-08-1.1113888091007D-07-3.9874233366857D-08
         12  2.2176249235855D-10 4.8440125486482D-10-3.7675797587086D-10-1.4212796721325D-11 9.8660078899568D-11-3.2593230381766D-11
          2  9.9861213227681D-01 1.0719515317387D-03-2.7020875424400D-04 1.8409151545658D-04 1.7189988028124D-04-4.9185952637707D-04
          3  2.3217499340863D-03 9.9886563610648D-01-2.1118786772861D-03-6.2231758566430D-04-4.8099900447839D-04 6.4105900884605D-04
          4 -4.7016744093896D-01 3.5977469141471D-01 4.3054499784029D-01 4.2130336842617D-02 3.4191894930856D-02-3.8359023208712D-02
          5 -4.4611970914352D-06-2.6231477233594D-06 1.9062207860216D-06-2.1707639844097D-06 1.9531765186756D-07 1.5784118586031D-06
          6  2.8049996251267D-06-2.5657317470985D-06-3.1478194123068D-06-6.7851081689480D-07-8.0453849112459D-07 6.4464358547248D-07
          7  3.1993197572747D-03-6.0164456677768D-03-4.9169772941829D-03-1.9826817474238D-03-1.8781547481017D-03 1.4568414255515D-03
          8 -2.2707668043291D-02 3.4631959857890D-02 2.0756515249950D-02 9.6933427022163D-03 7.8966867151586D-03-8.8586129886865D-03
          9  2.8247877443299D-02-2.1718909922230D-02-2.2701028610657D-02-2.4186418812867D-03-1.9385376989920D-03 2.8807780738819D-03
         10 -2.3095349937046D-01 3.4586272411422D-01 2.1833184260185D-01 9.6484916869541D-02 7.8827569969586D-02-8.7560692954179D-02
 2442335.500 -0.0023852630850134 -0.0004728521123148 -0.0004059658989754  0.0001281187325473 -0.0005651276519368 -0.0002140950641390
              0.0630164348039414  0.4181572430944009  0.1165092631573015  0.0000461514728412  0.0000293992528076  0.2299170405417692
          2 -1.9662439762579D+00 4.7601781321639D+00 1.7219745747577D+00-1.2495828685420D+00-1.4854125845494D-01-1.7275142405300D-01
          3  7.4013651607672D-04-2.5931534227297D-03-9.7012330378310D-04 6.7052097068281D-04 1.0033526791212D-04 1.0107408188937D-04
          4  6.5335721670881D-04-2.7162749402255D-03-1.0267627292330D-03 6.7571162970205D-04 1.2706655712434D-04 1.1219535671243D-04
          5 -1.6510076026573D-13 1.3676328241074D-13 4.6452102077077D-14-6.2941034438348D-14 2.6752431058166D-15-1.0933096213389D-15
          6  1.4155102360226D-02 2.1922419519225D-02 9.2804012432450D-03 2.3072934999882D-03-1.5332228294306D-03-3.3619652597920D-04
          7  2.9568721048446D+02-1.4008150452862D+03-5.3249440217363D+02 3.6815766890190D+02-1.9860687965201D+01 2.5483429835976D+01
          8  2.2525374803139D-07 2.2508061821254D-07 1.0535843467449D-07 2.4541949138638D-08-3.2502467764083D-08-4.7700170031943D-09
          9  1.9316894952166D+01-3.8435256133910D+01-1.3604121153675D+01 1.0296374486243D+01-4.4745043555985D-01 7.6334671344303D-01
         10  7.9408567555813D-07 5.0765343121883D-07-1.2050047631294D-06 7.4480466014920D-08 1.2518565615330D-06-3.8914249683600D-06
         11  3.7422575816934D-07-7.8241015570245D-07-2.5583327601120D-07 2.0557236671280D-07-1.2143516176050D-08 2.3830420182480D-08
         12 -1.0261157942899D-10 7.3040975305735D-10-4.4464455128877D-10-8.3457950085695D-11-2.9381926323522D-11 4.0834274111543D-12
          2  9.9949355735563D-01 1.3831520780912D-03-3.1719682641731D-03-1.8514792622304D-05-8.7296904060621D-06-4.9687612806948D-04
          3 -5.4021431250756D-03 9.9728255722467D-01 5.4071698633503D-03-1.9587467179346D-03 3.5126912067691D-05 1.8875116983061D-03
          4 -6.2772269500026D-01 3.9217784357811D-01 5.7565781840476D-01-1.4061779445398D-01 1.0620819318760D-02 1.2883786491257D-01
          5 -8.4773139318476D-06 1.2955252662185D-06 3.5927832919662D-06 1.2412754947259D-06 8.3845785460854D-07-1.4939902903248D-06
          6 -6.3583806239957D-06-4.9340504503997D-06 6.4259610117607D-06-2.1589707221711D-06 7.0409549634177D-08 2.4338522718668D-06
          7 -3.0074854246270D-02-1.2879866004336D-02 2.2797958486995D-02-1.0903953004048D-02 3.1339842817148D-04 9.1558482681152D-03
          8  9.4723509163443D-02 5.9179179040454D-02-8.6551835847801D-02 3.4129895031892D-02-1.0067674508453D-03-3.1187042406796D-02
          9  4.2222227804697D-02-2.2970687275560D-02-3.1991086249937D-02 1.0185576458261D-02-7.0143002521617D-04-8.6913314310914D-03
         10  9.4711812940000D-01 5.9195950882438D-01-8.5560267346235D-01 3.4352461217201D-01-1.0129252263045D-02-3.1303179879693D-01
-ADAMS-MOULTON NUMERICAL INTEGRATION HAS BEEN COMPLETED IN      60 STEPS FOR 1.5000000000000D+01 DAYS OR 4.00000D+00 STEPS PER DAY
                 REAL TIME ELAPSED=  0H  0M  0.04S  OR   0.0007S  PER STEP  OR   0.0027S PER DAY
                 TASK TIME ELAPSED=  0H  0M  0.03S  OR   0.0005S  PER STEP  OR   0.0020S PER DAY
           (TASK TIME)/(REAL TIME)= 0.75000
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   30
-
-
-
-

                        **********    *********    **            **   **********     *********    ********** 
                       ***********   ***********   ***          ***   ***********   ***********   ***********
                       **            **       **   ****        ****   **       **   **       **   **       **
                       **            **       **   ** **      ** **   **       **   **       **   **       **
                       **            **       **   **  **    **  **   **       **   **       **   **       **
                       **            **       **   **   **  **   **   **       **   **       **   **       **
                       **            **       **   **    ****    **   ***********   ***********   ***********
                       **            **       **   **     **     **   **********    ***********   ********** 
                       **            **       **   **            **   **            **       **   **  **     
                       **            **       **   **            **   **            **       **   **   **    
                       **            **       **   **            **   **            **       **   **    **   
                       **            **       **   **            **   **            **       **   **     **  
                       ***********   ***********   **            **   **            **       **   **      ** 
                        **********    *********    **            **   **            **       **   **       **
-
-
                                          **           **    **           **   **       **
                                          **           **    ***          **   **      ** 
                                          **           **    ****         **   **     **  
                                          **           **    ** **        **   **    **   
                                          **           **    **  **       **   **   **    
                                          **           **    **   **      **   *** **     
                                          **           **    **    **     **   *****      
                                          **           **    **     **    **   ** **      
                                          **           **    **      **   **   **  **     
                                          **           **    **       **  **   **   **    
                                          **           **    **        ** **   **    **   
                                          **           **    **         ****   **     **  
                                          ***********  **    **          ***   **      ** 
                                          ***********  **    **           **   **       **
- PEP VERSION= 20210302 PEP.PEPLOAD.PEP790                          
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   31
0INFORMATION FROM FIRST TWO RECORDS OF 10-BODY DATA SET 90 PRODUCED ON ITERATION  1 PAGE   15 OF RUN WITH TITLE
  N-BODY RUN 443-0 9-BODY INTEGRATION  442D8 INNER,  384F4 OUTER  I.C.            3/20/73 PERTURBING PLNT TAPE FROM N-BODY INTEG.
 JDBD1=2441841  JDBD2=2442400  NMOON= 0 NBDY1= 9   INTBD= 2  JVLBD=   0  EPSBD= 1.00000E-16
  1. MERCURY RUN 443- 3/21/73   NPL= 1   NCP= 0   INT= 2   JD0=2440001      A = 3.870988248079130D-01      E = 2.056143416457730D-01
             INC= 2.860343745748160D+01     ASC= 1.085941287077680D+01     PER= 6.693261684748430D+01    ANOM= 9.083965359314600D+01
  2.  VENUS  RUN 443- 3/21/73   NPL= 2   NCP= 0   INT= 4   JD0=2440001      A = 7.233280763362320D-01      E = 6.757538488379489D-03
             INC= 2.446691434342910D+01     ASC= 7.978123733702580D+00     PER= 1.237851072721730D+02    ANOM= 2.744064790468620D+02
  3.  EMBARY RUN 443- 3/21/73   NPL= 3   NCP= 0   INT= 4   JD0=2440001      A = 9.999827653637140D-01      E = 1.675271216169400D-02
             INC= 2.344335783108840D+01     ASC= 6.659250140529500D-04     PER= 1.020340640186800D+02    ANOM= 1.393658182548250D+02
  4.   MARS  RUN 443- 3/21/73   NPL= 4   NCP= 0   INT= 4   JD0=2440001      A = 1.523604534884400D+00      E = 9.345778628544589D-02
             INC= 2.469301678564210D+01     ASC= 3.346655760227140D+00     PER= 3.321011352392630D+02    ANOM= 8.989975217532620D+01
  5. JUPITER RUN 443- 3/20/73   NPL= 5   NCP= 0   INT= 4   JD0=2440001      A = 5.202945960488120D+00      E = 4.819167141557811D-02
             INC= 2.325294770783580D+01     ASC= 3.260973418885740D+00     PER= 1.037371340811626D+01    ANOM= 1.411862510204257D+02
  6.  SATURN RUN 443- 3/20/73   NPL= 6   NCP= 0   INT= 4   JD0=2440001      A = 9.526046869249818D+00      E = 5.461675768589707D-02
             INC= 2.257317229472867D+01     ASC= 5.968803885438612D+00     PER= 8.759383460863923D+01    ANOM= 2.897400092498020D+02
  7.  URANUS RUN 443- 3/20/73   NPL= 7   NCP= 0   INT= 4   JD0=2440001      A = 1.927444351332630D+01      E = 5.124121845145173D-02
             INC= 2.367141284461453D+01     ASC= 1.847313091106734D+00     PER= 1.683442408256917D+02    ANOM= 6.990628140027202D+00
  8. NEPTUNE RUN 443- 3/20/73   NPL= 8   NCP= 0   INT= 4   JD0=2440001      A = 3.011376065452675D+01      E = 6.984876998530826D-03
             INC= 2.231473025943044D+01     ASC= 3.519876825815031D+00     PER= 5.409517825004429D+01    ANOM= 1.775740764613225D+02
  9.  PLUTO  RUN 443- 3/20/73   NPL= 9   NCP= 0   INT= 4   JD0=2440001      A = 3.974605251930742D+01      E = 2.522316558278015D-01
             INC= 2.365507811729088D+01     ASC= 4.377229536665621D+01     PER= 1.819638797966960D+02    ANOM= 3.298440994501460D+02
 10.   MOON  RUN 440- 8/ 7/72   NPL=10   NCP= 3   INT=-1   JD0=2440001      A = 2.571514374642509D-03      E = 5.561544735761867D-02
             INC= 2.839685438634148D+01     ASC= 3.312959209874043D+00     PER= 2.262712472223586D+02    ANOM= 1.548856956849208D+02
   K( 1)=  -1   K( 2)=  -1   K( 3)=   0   K( 4)=  -1   K( 5)=  -1   K( 6)=  -1   K( 7)=  -1   K( 8)=  -1   K( 9)=  -1   K(10)=  -1
   K(11)=  -1   K(12)=  -1   K(13)=  -1   K(14)=  -1   K(15)=  -1   K(16)=  -1   K(17)=  -1   K(18)=  -1   K(19)=   0   K(20)=  -1
   K(21)=   0   K(22)=  -1   K(23)=  -1   K(24)=  -1   K(25)=  -1   K(26)=  -1   K(27)=  -1   K(28)=   3   K(29)=  11   K(30)=  54
   K(31)= -12   K(32)= -14   K(33)=  -1   K(34)=  -1   K(35)=  -1   K(36)=  -1   K(37)=  -1   K(38)=  -1   K(39)=   0   K(40)=  -1
 MASS1( 1)= 1.660759315500254D-07 MASS1( 2)= 2.447848585877872D-06 MASS1( 3)= 3.040436898620584D-06 MASS1( 4)= 3.227159776680543D-07
 MASS1( 5)= 9.547041567532574D-04 MASS1( 6)= 2.857997480046462D-04 MASS1( 7)= 4.374219404144008D-05 MASS1( 8)= 5.159795806034380D-05
 MASS1( 9)= 2.500000000000000D-07 MASS1(10)= 0.000000000000000D+00
 RELFT( 1)= 1.000000000000000D+00 RELFT( 2)= 1.000000000000000D+00 RELFT( 3)= 1.000000000000000D+00 RELFT( 4)= 1.000000000000000D+00
 RELFT( 5)= 1.000000000000000D+00 RELFT( 6)= 1.000000000000000D+00 RELFT( 7)= 1.000000000000000D+00 RELFT( 8)= 1.000000000000000D+00
 RELFT( 9)= 1.000000000000000D+00 RELFT(10)= 0.000000000000000D+00
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   32
0DATA ON FIRST TWO RECORDS OF DATA SET 20   BODY=   MOON     NPLNT(-2)= 10    NCENTR=  3     IPAR= 12     IFILTR= 0
 TITLE=TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE=   26 ITERAT= 1 LEVEL= 790
     JD1=2442321          JD0=2442321          JD2=2442332          INT= -1     ICND=  0
   COND(1)= 2.584844478828888D-03   COND(2)= 4.663456366998220D-02   COND(3)= 2.240083560206691D+01   COND(4)= 3.473239233192494D+02
   COND(5)= 1.425397834496920D+02   COND(6)= 2.222476213123667D+02   CON( 1)= 1.738090000000000D+03   CON( 2)= 0.000000000000000D+00
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 0.000000000000000D+00
   CON( 7)= 0.000000000000000D+00   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 8.936456813092100D-02   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00   CON(25)= 0.000000000000000D+00   CON(26)= 0.000000000000000D+00
   CON(27)= 0.000000000000000D+00   CON(28)= 0.000000000000000D+00   CON(29)= 0.000000000000000D+00   CON(30)= 0.000000000000000D+00
   CON(31)= 0.000000000000000D+00   CON(32)= 0.000000000000000D+00   CON(33)= 0.000000000000000D+00   CON(34)=-1.900000000000000D-11
   CON(35)= 0.000000000000000D+00   CON(36)= 0.000000000000000D+00   CON(
  PRM(  1)= 1.657848020429993D-07  PRM(  2)= 2.447848585877872D-06  PRM(  3)= 3.040436898620584D-06  PRM(  4)= 3.227159776680543D-07
  PRM(  5)= 9.547534692876814D-04  PRM(  6)= 2.857796067672611D-04  PRM(  7)= 4.361098996947231D-05  PRM(  8)= 5.192107995846314D-05
  PRM(  9)= 2.500000000000000D-07  PRM( 10)= 1.215052064980984D-02  PRM( 11)= 6.700000000000000D-10  PRM( 12)= 0.000000000000000D+00
  PRM( 13)= 0.000000000000000D+00  PRM( 14)= 0.000000000000000D+00  PRM( 15)= 0.000000000000000D+00  PRM( 16)= 0.000000000000000D+00
  PRM( 17)= 0.000000000000000D+00  PRM( 18)= 0.000000000000000D+00  PRM( 19)= 0.000000000000000D+00  PRM( 20)= 0.000000000000000D+00
  PRM( 21)= 0.000000000000000D+00  PRM( 22)= 0.000000000000000D+00  PRM( 23)= 0.000000000000000D+00  PRM( 24)= 0.000000000000000D+00
  PRM( 25)= 0.000000000000000D+00  PRM( 26)= 0.000000000000000D+00  PRM( 27)= 0.000000000000000D+00  PRM( 28)= 0.000000000000000D+00
  PRM( 29)= 0.000000000000000D+00  PRM( 30)= 0.000000000000000D+00  PRM( 31)= 1.000000000000000D+00  PRM( 32)= 1.000000000000000D-35
  PRM( 33)= 1.000000000000000D-08  PRM( 34)= 0.000000000000000D+00  PRM( 35)= 0.000000000000000D+00  PRM( 36)= 0.000000000000000D+00
  PRM( 37)= 0.000000000000000D+00  PRM( 38)= 0.000000000000000D+00  PRM( 39)=-7.800000000000000D-04  PRM( 40)= 0.000000000000000D+00
  PRM( 41)= 1.000000000000000D+00  PRM( 42)= 1.000000000000000D+00  PRM( 43)= 1.000000000000000D+00  PRM( 44)= 1.000000000000000D+00
  PRM( 45)= 0.000000000000000D+00  PRM( 46)= 0.000000000000000D+00  PRM( 47)= 0.000000000000000D+00  PRM( 48)= 2.344330000000000D+01
  PRM( 49)= 2.900000000000000D+00  PRM( 50)= 1.000000000000000D-09  PRM( 51)= 4.990047800000000D+02  PRM( 52)= 0.000000000000000D+00
  PRM( 53)= 1.000000000000000D+00  PRM( 54)= 0.000000000000000D+00  PRM( 55)= 0.000000000000000D+00  PRM( 56)= 0.000000000000000D+00
  PRM( 57)= 0.000000000000000D+00  PRM( 58)= 0.000000000000000D+00  PRM( 59)= 0.000000000000000D+00  PRM( 60)= 1.000000000000000D+01
  PRM( 61)= 0.000000000000000D+00  PRM( 62)= 1.000000000000000D+00  PRM( 63)= 1.000000000000000D+00  PRM( 64)= 0.000000000000000D+00
  PRM( 65)= 0.000000000000000D+00  PRM( 66)= 0.000000000000000D+00  PRM( 67)= 0.000000000000000D+00  PRM( 68)= 0.000000000000000D+00
  PRM( 69)= 0.000000000000000D+00  PRM( 70)= 0.000000000000000D+00  PRM( 71)= 0.000000000000000D+00  PRM( 72)= 0.000000000000000D+00
  PRM( 73)= 0.000000000000000D+00  PRM( 74)= 0.000000000000000D+00  PRM( 75)= 0.000000000000000D+00  PRM( 76)= 0.000000000000000D+00
  PRM( 77)= 0.000000000000000D+00  PRM( 78)= 0.000000000000000D+00  PRM( 79)= 0.000000000000000D+00  PRM( 80)= 0.000000000000000D+00
  PRM( 81)= 1.000000000000000D+00  PRM( 82)= 0.000000000000000D+00  PRM( 83)= 0.000000000000000D+00  PRM( 84)= 0.000000000000000D+00
  PRM( 85)= 0.000000000000000D+00  PRM( 86)= 0.000000000000000D+00  PRM( 87)= 0.000000000000000D+00  PRM( 88)= 0.000000000000000D+00
  PRM( 89)= 0.000000000000000D+00  PRM( 90)= 0.000000000000000D+00  PRM( 91)= 2.344578706750000D+01  PRM( 92)= 7.250000000000000D+00
  PRM( 93)= 7.506250000000000D+01  PRM( 94)= 6.960000000000000D+05  PRM( 95)= 6.960000000000000D+05  PRM( 96)= 0.000000000000000D+00
  PRM( 97)= 2.440000500000000D+06  PRM( 98)= 4.263529034000000D-05  PRM( 99)= 0.000000000000000D+00  PRM(100)= 2.997924580000000D+05
  K(  1)=   0  K(  2)=   0  K(  3)=   0  K(  4)=   0  K(  5)=   0  K(  6)=   0  K(  7)=   0  K(  8)=   0  K(  9)=   0  K( 10)=   0
  K( 11)=   0  K( 12)=   0  K( 13)=   0  K( 14)=   0  K( 15)=   0  K( 16)=   0  K( 17)=   0  K( 18)=   0  K( 19)=   0  K( 20)=   0
  K( 21)=   0  K( 22)=   0  K( 23)=   0  K( 24)=   0  K( 25)=   0  K( 26)=   0  K( 27)=   0  K( 28)=   0  K( 29)=   0  K( 30)=   0
  K( 31)=   1  K( 32)=   1  K( 33)=   1  K( 34)=   1  K( 35)=   1  K( 36)=   1  K( 37)=   1  K( 38)=   1  K( 39)=   1  K( 40)=   1
  K( 41)=  -1  K( 42)=  -1  K( 43)=  -1  K( 44)=  -1  K( 45)=  -1  K( 46)=  -1  K( 47)=  -1  K( 48)=  -1  K( 49)=  -1  K( 50)=  -1
  K( 51)=  -1  K( 52)=  -1  K( 53)=  -1  K( 54)=  -1  K( 55)=  -1  K( 56)=  -1  K( 57)=  -1  K( 58)=  -1  K( 59)=  -1  K( 60)=  -1
  K( 61)=   1  K( 62)=   1  K( 63)=  -1  K( 64)=  -1  K( 65)=  -1  K( 66)=  -1  K( 67)=  -1  K( 68)=  -1  K( 69)=  -1  K( 70)=  -1
  K( 71)=  -1  K( 72)=  -1  K( 73)=  -1  K( 74)=  -1  K( 75)=  -1  K( 76)=  -1  K( 77)=  -1  K( 78)=  -1  K( 79)=  -1  K( 80)=  -1
  K( 81)=   4  K( 82)=   3  K( 83)=   3  K( 84)=   1  K( 85)=   1  K( 86)=  -1  K( 87)=  -2  K( 88)=   2  K( 89)=   6  K( 90)=   6
  K( 91)=  -4  K( 92)=  -6  K( 93)=   0  K( 94)=   0  K( 95)=   0  K( 96)=   0  K( 97)=   0  K( 98)=  10  K( 99)=   0  K(100)=  -1
 NUMKI= 16  KI= 1 1-1-1-1 1 1  -16  -20    3   10   32  303 1031    2    3
   EPS(1)= 0.00000E+00   EPS(2)= 0.00000E+00   EPS(3)= 1.00000E-10   EPS(4)= 0.00000E+00   EPS(5)= 0.00000E+00   EPS(6)= 0.00000E+00
 INDIRECT PARTIALS ITERATED FOR  1 BODIES:  3
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   33
0DATA ON FIRST TWO RECORDS OF DATA SET 13   BODY=  EMBARY    NPLNT(-3)=  3    NCENTR= -1     IPAR= 13     IFILTR= 0
 TITLE=TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE=   14 ITERAT= 1 LEVEL= 790
     JD1=2442304          JD0=2442321          JD2=2442347          INT=  2     ICND=  0
   COND(1)= 1.000002347957079D+00   COND(2)= 1.669041814744669D-02   COND(3)= 2.344258335876651D+01   COND(4)= 8.391300825564327D-04
   COND(5)= 1.022643964965586D+02   COND(6)= 2.657488867161210D+02   CON( 1)= 6.378166000000000D+03   CON( 2)= 3.352329869259135D-03
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 0.000000000000000D+00
   CON( 7)= 0.000000000000000D+00   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00   CON(25)= 2.442320500000000D+06   CON(26)= 0.000000000000000D+00
   CON(27)= 0.000000000000000D+00   CON(28)= 0.000000000000000D+00   CON(29)= 0.000000000000000D+00   CON(30)= 0.000000000000000D+00
   CON(31)= 0.000000000000000D+00   CON(32)= 0.000000000000000D+00   CON(33)= 0.000000000000000D+00   CON(34)=-4.600000000000000D-10
   CON(35)= 0.000000000000000D+00   CON(36)= 0.000000000000000D+00   CON(
  PRM(  1)= 1.657848020429993D-07  PRM(  2)= 2.447848585877872D-06  PRM(  3)= 3.040436898620584D-06  PRM(  4)= 3.227159776680543D-07
  PRM(  5)= 9.547534692876814D-04  PRM(  6)= 2.857796067672611D-04  PRM(  7)= 4.361098996947231D-05  PRM(  8)= 5.192107995846314D-05
  PRM(  9)= 2.500000000000000D-07  PRM( 10)= 1.215052064980984D-02  PRM( 11)= 6.700000000000000D-10  PRM( 12)= 0.000000000000000D+00
  PRM( 13)= 0.000000000000000D+00  PRM( 14)= 0.000000000000000D+00  PRM( 15)= 0.000000000000000D+00  PRM( 16)= 0.000000000000000D+00
  PRM( 17)= 0.000000000000000D+00  PRM( 18)= 0.000000000000000D+00  PRM( 19)= 0.000000000000000D+00  PRM( 20)= 0.000000000000000D+00
  PRM( 21)= 0.000000000000000D+00  PRM( 22)= 0.000000000000000D+00  PRM( 23)= 0.000000000000000D+00  PRM( 24)= 0.000000000000000D+00
  PRM( 25)= 0.000000000000000D+00  PRM( 26)= 0.000000000000000D+00  PRM( 27)= 0.000000000000000D+00  PRM( 28)= 0.000000000000000D+00
  PRM( 29)= 0.000000000000000D+00  PRM( 30)= 0.000000000000000D+00  PRM( 31)= 1.000000000000000D+00  PRM( 32)= 1.000000000000000D-35
  PRM( 33)= 1.000000000000000D-08  PRM( 34)= 0.000000000000000D+00  PRM( 35)= 0.000000000000000D+00  PRM( 36)= 0.000000000000000D+00
  PRM( 37)= 0.000000000000000D+00  PRM( 38)= 0.000000000000000D+00  PRM( 39)=-7.800000000000000D-04  PRM( 40)= 0.000000000000000D+00
  PRM( 41)= 1.000000000000000D+00  PRM( 42)= 1.000000000000000D+00  PRM( 43)= 1.000000000000000D+00  PRM( 44)= 1.000000000000000D+00
  PRM( 45)= 0.000000000000000D+00  PRM( 46)= 0.000000000000000D+00  PRM( 47)= 0.000000000000000D+00  PRM( 48)= 2.344330000000000D+01
  PRM( 49)= 2.900000000000000D+00  PRM( 50)= 1.000000000000000D-09  PRM( 51)= 4.990047800000000D+02  PRM( 52)= 0.000000000000000D+00
  PRM( 53)= 1.000000000000000D+00  PRM( 54)= 0.000000000000000D+00  PRM( 55)= 0.000000000000000D+00  PRM( 56)= 0.000000000000000D+00
  PRM( 57)= 0.000000000000000D+00  PRM( 58)= 0.000000000000000D+00  PRM( 59)= 0.000000000000000D+00  PRM( 60)= 1.000000000000000D+01
  PRM( 61)= 0.000000000000000D+00  PRM( 62)= 1.000000000000000D+00  PRM( 63)= 1.000000000000000D+00  PRM( 64)= 0.000000000000000D+00
  PRM( 65)= 0.000000000000000D+00  PRM( 66)= 0.000000000000000D+00  PRM( 67)= 0.000000000000000D+00  PRM( 68)= 0.000000000000000D+00
  PRM( 69)= 0.000000000000000D+00  PRM( 70)= 0.000000000000000D+00  PRM( 71)= 0.000000000000000D+00  PRM( 72)= 0.000000000000000D+00
  PRM( 73)= 0.000000000000000D+00  PRM( 74)= 0.000000000000000D+00  PRM( 75)= 0.000000000000000D+00  PRM( 76)= 0.000000000000000D+00
  PRM( 77)= 0.000000000000000D+00  PRM( 78)= 0.000000000000000D+00  PRM( 79)= 0.000000000000000D+00  PRM( 80)= 0.000000000000000D+00
  PRM( 81)= 1.000000000000000D+00  PRM( 82)= 0.000000000000000D+00  PRM( 83)= 0.000000000000000D+00  PRM( 84)= 0.000000000000000D+00
  PRM( 85)= 0.000000000000000D+00  PRM( 86)= 0.000000000000000D+00  PRM( 87)= 0.000000000000000D+00  PRM( 88)= 0.000000000000000D+00
  PRM( 89)= 0.000000000000000D+00  PRM( 90)= 0.000000000000000D+00  PRM( 91)= 2.344578706750000D+01  PRM( 92)= 7.250000000000000D+00
  PRM( 93)= 7.506250000000000D+01  PRM( 94)= 6.960000000000000D+05  PRM( 95)= 6.960000000000000D+05  PRM( 96)= 0.000000000000000D+00
  PRM( 97)= 2.440000500000000D+06  PRM( 98)= 4.263529034000000D-05  PRM( 99)= 0.000000000000000D+00  PRM(100)= 2.997924580000000D+05
  K(  1)=   0  K(  2)=   0  K(  3)=   0  K(  4)=   0  K(  5)=   0  K(  6)=   0  K(  7)=   0  K(  8)=   0  K(  9)=   0  K( 10)=   0
  K( 11)=   0  K( 12)=   0  K( 13)=   0  K( 14)=   0  K( 15)=   0  K( 16)=   0  K( 17)=   0  K( 18)=   0  K( 19)=   0  K( 20)=   0
  K( 21)=   0  K( 22)=   0  K( 23)=   0  K( 24)=   0  K( 25)=   0  K( 26)=   0  K( 27)=   0  K( 28)=   0  K( 29)=   0  K( 30)=   0
  K( 31)=  -1  K( 32)=  -1  K( 33)=   1  K( 34)=   1  K( 35)=   1  K( 36)=  -1  K( 37)=  -1  K( 38)=  -1  K( 39)=  -1  K( 40)=   1
  K( 41)=  -1  K( 42)=  -1  K( 43)=  -1  K( 44)=  -1  K( 45)=  -1  K( 46)=  -1  K( 47)=  -1  K( 48)=  -1  K( 49)=  -1  K( 50)=  -1
  K( 51)=  -1  K( 52)=  -1  K( 53)=  -1  K( 54)=  -1  K( 55)=  -1  K( 56)=  -1  K( 57)=  -1  K( 58)=  -1  K( 59)=  -1  K( 60)=  -1
  K( 61)=   1  K( 62)=   1  K( 63)=   1  K( 64)=  -1  K( 65)=  -1  K( 66)=  -1  K( 67)=  -1  K( 68)=  -1  K( 69)=  -1  K( 70)=   1
  K( 71)=   1  K( 72)=   1  K( 73)=  -1  K( 74)=  -1  K( 75)=  -1  K( 76)=  -1  K( 77)=  -1  K( 78)=  -1  K( 79)=   1  K( 80)=   1
  K( 81)=  -1  K( 82)=   0  K( 83)=  -1  K( 84)=  -1  K( 85)=  -1  K( 86)=  -1  K( 87)=   2  K( 88)=   2  K( 89)=   6  K( 90)=   6
  K( 91)=  -3  K( 92)=  -6  K( 93)=   0  K( 94)=   0  K( 95)=   0  K( 96)=   0  K( 97)=   0  K( 98)=   4  K( 99)=   0  K(100)=  -1
 NUMKI= 16  KI= 1 1 1 1-1-1-1    3    4   10   31   32   33   40   41   42
   EPS(1)= 0.00000E+00   EPS(2)= 0.00000E+00   EPS(3)= 1.00000E-09   EPS(4)= 0.00000E+00   EPS(5)= 0.00000E+00   EPS(6)= 0.00000E+00
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   34
0DATA ON FIRST TWO RECORDS OF DATA SET 21   BODY=  MROTAT    NPLNT( 0)=-10    NCENTR=  3     IPAR= 10     IFILTR= 0
 TITLE=TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE=   27 ITERAT= 1 LEVEL= 790
     JD1=2442321          JD0=2442321          JD2=2442332          INT= -1     ICND=  0
   COND(1)= 6.406890562446280D-02   COND(2)= 4.174732965826930D-01   COND(3)= 4.427722151828226D+02   COND(4)=-1.116475592635000D-04
   COND(5)=-3.914136528300000D-05   COND(6)= 2.300859910558332D-01   CON( 1)= 1.738090000000000D+03   CON( 2)= 4.039070136017979D-04
   CON( 3)= 6.317168491855110D-04   CON( 4)= 2.278098937105160D-04   CON( 5)= 2.692029800000000D-02   CON( 6)= 2.415119936539230D-02
   CON( 7)= 4.683189055748760D-03   CON( 8)= 0.000000000000000D+00   CON( 9)= 0.000000000000000D+00   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00   CON(25)= 0.000000000000000D+00   CON(26)= 0.000000000000000D+00
   CON(27)= 0.000000000000000D+00   CON(28)= 0.000000000000000D+00   CON(29)= 0.000000000000000D+00   CON(30)= 0.000000000000000D+00
   CON(31)= 0.000000000000000D+00   CON(32)= 0.000000000000000D+00   CON(33)= 0.000000000000000D+00   CON(34)= 0.000000000000000D+00
   CON(35)= 0.000000000000000D+00   CON(36)= 0.000000000000000D+00   CON(
  PRM(  1)= 1.657848020429993D-07  PRM(  2)= 2.447848585877872D-06  PRM(  3)= 3.040436898620584D-06  PRM(  4)= 3.227159776680543D-07
  PRM(  5)= 9.547534692876814D-04  PRM(  6)= 2.857796067672611D-04  PRM(  7)= 4.361098996947231D-05  PRM(  8)= 5.192107995846314D-05
  PRM(  9)= 2.500000000000000D-07  PRM( 10)= 1.215052064980984D-02  PRM( 11)= 6.700000000000000D-10  PRM( 12)= 0.000000000000000D+00
  PRM( 13)= 0.000000000000000D+00  PRM( 14)= 0.000000000000000D+00  PRM( 15)= 0.000000000000000D+00  PRM( 16)= 0.000000000000000D+00
  PRM( 17)= 0.000000000000000D+00  PRM( 18)= 0.000000000000000D+00  PRM( 19)= 0.000000000000000D+00  PRM( 20)= 0.000000000000000D+00
  PRM( 21)= 0.000000000000000D+00  PRM( 22)= 0.000000000000000D+00  PRM( 23)= 0.000000000000000D+00  PRM( 24)= 0.000000000000000D+00
  PRM( 25)= 0.000000000000000D+00  PRM( 26)= 0.000000000000000D+00  PRM( 27)= 0.000000000000000D+00  PRM( 28)= 0.000000000000000D+00
  PRM( 29)= 0.000000000000000D+00  PRM( 30)= 0.000000000000000D+00  PRM( 31)= 1.000000000000000D+00  PRM( 32)= 1.000000000000000D-35
  PRM( 33)= 1.000000000000000D-08  PRM( 34)= 0.000000000000000D+00  PRM( 35)= 0.000000000000000D+00  PRM( 36)= 0.000000000000000D+00
  PRM( 37)= 0.000000000000000D+00  PRM( 38)= 0.000000000000000D+00  PRM( 39)=-7.800000000000000D-04  PRM( 40)= 0.000000000000000D+00
  PRM( 41)= 1.000000000000000D+00  PRM( 42)= 1.000000000000000D+00  PRM( 43)= 1.000000000000000D+00  PRM( 44)= 1.000000000000000D+00
  PRM( 45)= 0.000000000000000D+00  PRM( 46)= 0.000000000000000D+00  PRM( 47)= 0.000000000000000D+00  PRM( 48)= 2.344330000000000D+01
  PRM( 49)= 2.900000000000000D+00  PRM( 50)= 1.000000000000000D-09  PRM( 51)= 4.990047800000000D+02  PRM( 52)= 0.000000000000000D+00
  PRM( 53)= 1.000000000000000D+00  PRM( 54)= 0.000000000000000D+00  PRM( 55)= 0.000000000000000D+00  PRM( 56)= 0.000000000000000D+00
  PRM( 57)= 0.000000000000000D+00  PRM( 58)= 0.000000000000000D+00  PRM( 59)= 0.000000000000000D+00  PRM( 60)= 1.000000000000000D+01
  PRM( 61)= 0.000000000000000D+00  PRM( 62)= 1.000000000000000D+00  PRM( 63)= 1.000000000000000D+00  PRM( 64)= 0.000000000000000D+00
  PRM( 65)= 0.000000000000000D+00  PRM( 66)= 0.000000000000000D+00  PRM( 67)= 0.000000000000000D+00  PRM( 68)= 0.000000000000000D+00
  PRM( 69)= 0.000000000000000D+00  PRM( 70)= 0.000000000000000D+00  PRM( 71)= 0.000000000000000D+00  PRM( 72)= 0.000000000000000D+00
  PRM( 73)= 0.000000000000000D+00  PRM( 74)= 0.000000000000000D+00  PRM( 75)= 0.000000000000000D+00  PRM( 76)= 0.000000000000000D+00
  PRM( 77)= 0.000000000000000D+00  PRM( 78)= 0.000000000000000D+00  PRM( 79)= 0.000000000000000D+00  PRM( 80)= 0.000000000000000D+00
  PRM( 81)= 1.000000000000000D+00  PRM( 82)= 0.000000000000000D+00  PRM( 83)= 0.000000000000000D+00  PRM( 84)= 0.000000000000000D+00
  PRM( 85)= 0.000000000000000D+00  PRM( 86)= 0.000000000000000D+00  PRM( 87)= 0.000000000000000D+00  PRM( 88)= 0.000000000000000D+00
  PRM( 89)= 0.000000000000000D+00  PRM( 90)= 0.000000000000000D+00  PRM( 91)= 2.344578706750000D+01  PRM( 92)= 7.250000000000000D+00
  PRM( 93)= 7.506250000000000D+01  PRM( 94)= 6.960000000000000D+05  PRM( 95)= 6.960000000000000D+05  PRM( 96)= 0.000000000000000D+00
  PRM( 97)= 2.440000500000000D+06  PRM( 98)= 4.263529034000000D-05  PRM( 99)= 0.000000000000000D+00  PRM(100)= 2.997924580000000D+05
  K(  1)=   0  K(  2)=   0  K(  3)=   0  K(  4)=   0  K(  5)=   0  K(  6)=   0  K(  7)=   0  K(  8)=   0  K(  9)=   0  K( 10)=   0
  K( 11)=   0  K( 12)=   0  K( 13)=   0  K( 14)=   0  K( 15)=   0  K( 16)=   0  K( 17)=   0  K( 18)=   0  K( 19)=   0  K( 20)=   0
  K( 21)=   0  K( 22)=   0  K( 23)=   0  K( 24)=   0  K( 25)=   0  K( 26)=   0  K( 27)=   0  K( 28)=   0  K( 29)=   0  K( 30)=   0
  K( 31)=  -1  K( 32)=  -1  K( 33)=   1  K( 34)=  -1  K( 35)=  -1  K( 36)=  -1  K( 37)=  -1  K( 38)=  -1  K( 39)=  -1  K( 40)=   1
  K( 41)=  -1  K( 42)=  -1  K( 43)=  -1  K( 44)=  -1  K( 45)=  -1  K( 46)=  -1  K( 47)=  -1  K( 48)=  -1  K( 49)=  -1  K( 50)=  -1
  K( 51)=  -1  K( 52)=  -1  K( 53)=  -1  K( 54)=  -1  K( 55)=  -1  K( 56)=  -1  K( 57)=  -1  K( 58)=  -1  K( 59)=  -1  K( 60)=  -1
  K( 61)=  -1  K( 62)=  -1  K( 63)=  -1  K( 64)=  -1  K( 65)=  -1  K( 66)=  -1  K( 67)=  -1  K( 68)=  -1  K( 69)=  -1  K( 70)=  -1
  K( 71)=  -1  K( 72)=  -1  K( 73)=  -1  K( 74)=  -1  K( 75)=  -1  K( 76)=  -1  K( 77)=  -1  K( 78)=  -1  K( 79)=  -1  K( 80)=  -1
  K( 81)=   1  K( 82)=   0  K( 83)=   0  K( 84)=  -1  K( 85)=  -1  K( 86)= 303  K( 87)=  -2  K( 88)=   2  K( 89)=   6  K( 90)=   6
  K( 91)=  -4  K( 92)=  -6  K( 93)=   0  K( 94)=   0  K( 95)=   0  K( 96)=   0  K( 97)=   0  K( 98)=  10  K( 99)=   0  K(100)=  -1
 NUMKI= 18  KI= 1 1 1-1-1-1-1   -3   -6   -7 1031    2    3 1041    3    1    3    2
   EPS(1)= 0.00000E+00   EPS(2)= 0.00000E+00   EPS(3)= 1.00000E-10   EPS(4)= 0.00000E+00   EPS(5)= 0.00000E+00   EPS(6)= 0.00000E+00
 INDIRECT PARTIALS ITERATED FOR  1 BODIES:  3
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   35
0DATA ON FIRST TWO RECORDS OF UT1 DATA SET 93
 HEADING=IERS-B TAI-UT1R: ABBREVIATED FOR BIGTEST                                        
 FORMAT=(5X,I5,6(I7,1X),20X,I2)           KIND= 2  JDT1=2442259  JDT2=2442524  NPR= 6  INT= 5  UNITS=   1.000000E-04 (SEC)
0DATA ON FIRST TWO RECORDS OF WOBBLE DATA SET 94
 HEADING=IERS-B POLE: ABBREVIATED FOR BIGTEST                                            
 FORMAT=(1X,I9,12I5,8X,I2)                KIND= 0  JDT1=2442259  JDT2=2442524  NPR= 6  INT= 5  UNITS=   1.000000E-03 (ARCSEC)
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   36

 NPLNT0=10  NCODF=NCODA(1)= 1  NTAPA(1)= -5  NSEQA(1)= 1000  IABS1=IOBS0( 1)= 0  IABS2=IOBS1( 1)=31  IOBS= 0  IOBCON= 2  NREWND= 0
0 DEL/DOP RADAR   OBSERVATIONS OF   MOON   MADE FROM HAYSTACK,HAYSTACK  SERIES=LASR  (SPOT=AP15)  FREQUENCY= 4.3178000000000D+14
0                 NAME  NSITE LSITE     RADIUS       LONGITUDE      LATITUDE  KSITE
   RECEIVE SITE=HAYSTACK   1  0 0 0 6368.551653028   71.48866667    42.43151838  0
      SEND SITE=HAYSTACK   1  0 0 0 6368.551653028   71.48866667    42.43151838  0
0ACCPRC= 1.000E+01  ACCDST= 1.000E-09  ACCTIM= 1.000E-09  ERWGT1= 1.000E+00  ERWGT2= 1.000E+00  FDEV=   000.0D+00
0  NSPOT =  9   NPREC =  0   NDPREC=  0   NTMDLY=  2   NDOP  = -1   NMEDIA=  0   NEATM = -1   NPATM = -1   NPSHP = -1   NPLNG = -1
   NEION = -1  KOB(11)= -1   LOPTRF=  1   NDDIFF=  0   NLIBPR=  2  KOB(15)= -1   CALCVL= -1   ITIME =  1   NTIME =  0    IWOB =  1
0  SETUP FOR  DEL/DOP RADAR   OBSERVATION SERIES REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
0 GRNWCH  JULIAN  REC UTC (WWV)  UT1-UTC AT-UTC  CT-AT             TIME DELAY (SEC)               DOPPLER SHIFT (CYCLES/SECOND) 
   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC       OBSERVED     ERROR     OBS-TH         OBSERVED     ERROR     OBS-TH 
 SUMCOR(1)=    5.5142000286E-21
 10/ 5/74 2442326 10  1 30.5000 -0.0352 13.0344 32.1483    2.5116630865 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  2.51166308654470D+00     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00 -3.10993208547776D+04 
  0.00000000000000D+00  0.00000000000000D+00  1.33384487960790D-04  0.00000000000000D+00 -1.06525735102907D+01  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  2.55546416239231D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  4.34066864539734D-04  0.00000000000000D+00  1.09835117994902D+03 
 -5.06879838631856D-02 -1.66268439461765D-01  9.05209220603499D-12  4.57415001191956D+00  2.76311054452907D-03 -4.38872635319860D-03 
 -9.51303806300444D-04  5.79091927337183D-09  1.72134643821198D-09 -7.22635693554119D-05 -5.36960730335582D-05  5.89706042301939D-05 
 -5.37262215673206D-04
 SUMCOR(1)=    5.5355479931E-21
 10/ 5/74 2442326 12  1 30.5000 -0.0354 13.0344 32.1483    2.5211421281 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  2.52114212814960D+00     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00 -3.37750868163981D+04 
  0.00000000000000D+00  0.00000000000000D+00  1.35608529578362D-04  0.00000000000000D+00 -1.09845047134931D+01  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  2.55345929228089D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  4.36242255978607D-04  0.00000000000000D+00  1.10701729641416D+03 
 -6.23483837106993D-02 -1.77581406152724D-01  9.43904226919407D-12  4.69585517089533D+00  2.76958724082310D-03 -4.40706844780768D-03 
 -9.68555772380933D-04  5.87334417630943D-09  1.78570745408883D-09 -7.45249048796160D-05 -5.50239801297144D-05  6.01084789397896D-05 
 -5.50260642369905D-04
 SUMCOR(1)=    5.5634681660E-21
 10/ 5/74 2442326 14  1 30.5000 -0.0357 13.0344 32.1483    2.5336200273 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  2.53362002726567D+00     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00 -3.49832983637453D+04 
  0.00000000000000D+00  0.00000000000000D+00  1.37348355308459D-04  0.00000000000000D+00 -1.12951190513701D+01  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  2.55150732370803D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  4.32266386016004D-04  0.00000000000000D+00  1.11052531775656D+03 
 -6.69639786467147D-02 -1.81589509828100D-01  9.85416650616273D-12  4.80539373679650D+00  2.74434531142016D-03 -4.42021891058127D-03 
 -9.82331161912244D-04  5.95468934986628D-09  1.84886780624598D-09 -7.66092613761399D-05 -5.61783565059843D-05  6.10167647439149D-05 
 -5.61576773800949D-04
 SUMCOR(1)=    5.5900070120E-21
 10/ 5/74 2442326 16  1 30.5000 -0.0359 13.0344 32.1483    2.5454829947 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  2.54548299472129D+00     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00 -3.44821327144919D+04 
  0.00000000000000D+00  0.00000000000000D+00  1.38484840374946D-04  0.00000000000000D+00 -1.15779851607773D+01  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  2.54964544693320D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  4.23833484086265D-04  0.00000000000000D+00  1.10803631877106D+03 
 -6.37808853519976D-02 -1.77475025259394D-01  1.02971961991443D-11  4.89929921965536D+00  2.68428977479117D-03 -4.43333427323831D-03 
 -9.93159294356223D-04  6.04565106751119D-09  1.91256238059180D-09 -7.84663967540942D-05 -5.71886814531354D-05  6.17221178335808D-05 
 -5.71524621196332D-04
 SUMCOR(1)=    5.6076920365E-21
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   37
0 DEL/DOP RADAR   OBSERVATIONS OF   MOON   MADE FROM HAYSTACK,HAYSTACK  SERIES=LASR  (SPOT=AP15)  FREQUENCY= 4.3178000000000D+14
0 GRNWCH  JULIAN  REC UTC (WWV)  UT1-UTC AT-UTC  CT-AT             TIME DELAY (SEC)               DOPPLER SHIFT (CYCLES/SECOND) 
   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC       OBSERVED     ERROR     OBS-TH         OBSERVED     ERROR     OBS-TH 
 10/ 5/74 2442326 18  1 30.5000 -0.0362 13.0344 32.1483    2.5533334076 1.0E-08  0.00000E+00 
 SUMCOR(1)=    5.6113174141E-21
 10/ 5/74 2442326 20  1 30.5000 -0.0364 13.0344 32.1483    2.5548004327 1.0E-08  0.00000E+00 
 SUMCOR(1)=    5.5991201351E-21
 10/ 5/74 2442326 22  1 30.5000 -0.0366 13.0344 32.1483    2.5490757660 1.0E-08  0.00000E+00 
 SUMCOR(1)=    5.5731608810E-21
 10/ 6/74 2442327  0  1 30.5000 -0.0369 13.0344 32.1483    2.5370880120 1.0E-08  0.00000E+00 
 SUMCOR(1)=    5.5388567155E-21
 10/ 6/74 2442327  2  1 30.5000 -0.0371 13.0344 32.1483    2.5212933387 1.0E-08  0.00000E+00 
0ERROR ANALYSIS FOR  DEL/DOP RADAR   OBSERVATION SERIES EXTENDING FOR   2 PAGES STARTING ON PAGE   36
 AVERAGE NUMBER OF ITERATIONS FOR     9 OBSERVATIONS WAS (  3.00000,  3.00000)
                                   TIME DELAY      DOPPLER SHIFT      
 NUMBER OF MEASUREMENTS DELETED          0             0
 NUMBER OF MEASUREMENTS INCLUDED         9             0
                AVERAGE (OBS-TH)   0.00000E+00   0.00000E+00
             AVERAGE ABS(OBS-TH)   0.00000E+00   0.00000E+00
       ROOT MEAN SQUARE (OBS-TH)   0.00000E+00   0.00000E+00
                   AVERAGE ERROR   1.00000E-08   0.00000E+00
          ROOT MEAN SQUARE ERROR   1.00000E-08   0.00000E+00
          AVERAGE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON= 2, IOBS= 0, IABS1= 0) FOR THIS OBSERVATION SERIES WERE (    0,    0,    0)
0     9 +   0 OBSERVATION CARDS PROCESSED (NUMBER OF PARTIAL DERIVATIVES= 32)
0  PROCESSING OBSERVATION SERIES REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
-ERROR ANALYSIS FOR  DEL/DOP RADAR   OBSERVATIONS OF   MOON   PROCESSED DURING THE LAST    8 PAGES STARTING ON PAGE   30
                                   TIME DELAY      DOPPLER SHIFT      
 NUMBER OF MEASUREMENTS DELETED          0             0
 NUMBER OF MEASUREMENTS INCLUDED         9             0
          AVERAGE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON,IOBS,IABS1) = (    0,    0,    0)

       9 +    0 OBSERVATION RECORDS PROCESSED (NPARAM=  44)
0   PROCESSING OBSERVATIONS  REQUIRED
           0H  0M  0.03S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   38

 NPLNT0=10  NCODF=NCODA(1)= 6  NTAPA(1)= -5  NSEQA(1)= 1010  IABS1=IOBS0( 1)= 0  IABS2=IOBS1( 1)=31  IOBS= 0  IOBCON= 2  NREWND= 0
0 ASTROGRAPHIC    OBSERVATIONS OF   MOON   MADE FROM NICE               SERIES=ASTG  (SPOT=    )  FREQUENCY= 4.3178000000000D+14
0                 NAME  NSITE LSITE     RADIUS       LONGITUDE      LATITUDE  KSITE
   RECEIVE SITE=NICE      30  0 0 0 6378.000000000   -7.30041666     0.00000000  0
0ACCPRC= 1.000E+01  ACCDST= 1.000E-09  ACCTIM= 1.000E-09  ERWGT1= 1.000E+00  ERWGT2= 1.000E+00  FDEV=   000.0D+00
0  NSPOT =  0   NPREC = -1   NDPREC= -1   NTMDLY= -1   NDOP  = -1   NMEDIA= -1   NEATM = -1   NPATM = -1   NPSHP = -1   NPLNG = -1
   NEION = -1  KOB(11)= -1   LOPTRF= -1   NDDIFF= -1   NLIBPR= -1  KOB(15)= -1   CALCVL= -1   ITIME =  1   NTIME =  0    IWOB =  0
0TOPOCENTRIC RIGHT ASCENSION-DECLINATION REFERRED TO THE MEAN EQUINOX AND EQUATOR OF 1950.0
0  SETUP FOR  ASTROGRAPHIC    OBSERVATION SERIES REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
0 GRNWCH   JULIAN   UT  TIME OF   AT-UT    CT-AT           RIGHT ASCENSION                    DECLINATION                       
   DATE     DAY     OBSERVATION                      OBSERVED     ERROR    OBS-TH      OBSERVED     ERROR   OBS-TH              
           NUMBER  HR MIN  SEC     SEC      SEC    HR MIN  SEC     SEC      SEC      DEG  '   ''     ''      ''    CLMP LIMB OBS
 10/ 5/74 2442326  10  1 30.5000 13.0344  32.1483   3 40 31.2344  0.0107   0.00000    20 48 43.662  0.100   0.0000              
 ILDT=   0   0   0   0   0   0  1.32312343834402D+04     DERIV(1)=  1.06980416752139D-02  0.00000000000000D+00  3.01327442368364D+09 
  0.00000000000000D+00  0.00000000000000D+00 -7.75157532953950D-01  0.00000000000000D+00  4.94307026298841D+04  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -6.91591782887055D+00  0.00000000000000D+00 -1.05483418095910D+07 
  1.43697443263565D+04  1.48260850549426D+04  3.04443353842894D-08 -2.36756561332345D+04  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  3.90005696662078D-01 -1.80177890819616D-04  0.00000000000000D+00 
  0.00000000000000D+00
 ILDT=   0   0   0   0   0   0  7.49236623629305D+04     DERIV(2)=  1.00000001490116D-01  0.00000000000000D+00  6.15451604240396D+09 
  0.00000000000000D+00  0.00000000000000D+00 -3.40477412347867D+00  0.00000000000000D+00  9.89299820435242D+04  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  3.03959412201458D+02  0.00000000000000D+00 -2.15092019781745D+07 
  2.93173137932349D+04  3.02463894510705D+04 -2.16905618810373D-07 -7.22784220413192D+04  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -8.84850015387290D-02 -1.68073862636349D-02  0.00000000000000D+00 
  0.00000000000000D+00
 10/ 5/74 2442326  12  1 30.5000 13.0344  32.1483   3 46 30.8449  0.0107   0.00000    20 49 48.563  0.100   0.0000              
 ILDT=   0   0   0   0   0   0  1.35908449395709D+04     DERIV(1)=  1.06993222499805D-02  0.00000000000000D+00  3.05221460780963D+09 
  0.00000000000000D+00  0.00000000000000D+00 -8.14162730549634D-01  0.00000000000000D+00  5.18457433939955D+04  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -6.92931162587807D+00  0.00000000000000D+00 -1.07160267868035D+07 
  1.43097110807360D+04  1.47927223500044D+04  2.99377298258167D-08 -2.48985637160892D+04  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  4.07804021586392D-01 -2.04421606156128D-04  0.00000000000000D+00 
  0.00000000000000D+00
 ILDT=   0   0   0   0   0   0  7.49885631513917D+04     DERIV(2)=  1.00000001490116D-01  0.00000000000000D+00  5.81743742257433D+09 
  0.00000000000000D+00  0.00000000000000D+00 -3.37005895842738D+00  0.00000000000000D+00  9.41588757763724D+04  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  3.11860583201168D+02  0.00000000000000D+00 -2.01447846567135D+07 
  2.72642863287453D+04  2.81548584625964D+04 -2.24386027201540D-07 -7.07540342670141D+04  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -1.39495341532530D-01 -1.71687491551635D-02  0.00000000000000D+00 
  0.00000000000000D+00
 10/ 5/74 2442326  14  1 30.5000 13.0344  32.1483   3 53  5.8264  0.0107   0.00000    20 54 27.052  0.100   0.0000              
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   39
0 ASTROGRAPHIC    OBSERVATIONS OF   MOON   MADE FROM NICE               SERIES=ASTG  (SPOT=    )  FREQUENCY= 4.3178000000000D+14
0 GRNWCH   JULIAN   UT  TIME OF   AT-UT    CT-AT           RIGHT ASCENSION                    DECLINATION                       
   DATE     DAY     OBSERVATION                      OBSERVED     ERROR    OBS-TH      OBSERVED     ERROR   OBS-TH              
           NUMBER  HR MIN  SEC     SEC      SEC    HR MIN  SEC     SEC      SEC      DEG  '   ''     ''      ''    CLMP LIMB OBS
 ILDT=   0   0   0   0   0   0  1.39858263862354D+04     DERIV(1)=  1.07048307823257D-02  0.00000000000000D+00  3.10216149277548D+09 
  0.00000000000000D+00  0.00000000000000D+00 -8.58475325791946D-01  0.00000000000000D+00  5.46619988742014D+04  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -6.94049207322980D+00  0.00000000000000D+00 -1.09351592807690D+07 
  1.43007749316415D+04  1.48133448131929D+04  2.93012450415056D-08 -2.63061327861532D+04  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  4.28542475315041D-01 -2.32842823692214D-04  0.00000000000000D+00 
  0.00000000000000D+00
 ILDT=   0   0   0   0   0   0  7.52670516410563D+04     DERIV(2)=  1.00000001490116D-01  0.00000000000000D+00  5.44786423664521D+09 
  0.00000000000000D+00  0.00000000000000D+00 -3.33704366012746D+00  0.00000000000000D+00  8.94067153471260D+04  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  3.20792775354667D+02  0.00000000000000D+00 -1.87576805041731D+07 
  2.51191557231151D+04  2.59712565174798D+04 -2.33718499372366D-07 -6.92534111057462D+04  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -1.93172682386117D-01 -1.75853868758516D-02  0.00000000000000D+00 
  0.00000000000000D+00
 10/ 5/74 2442326  16  1 30.5000 13.0344  32.1483   3 59 49.2336  0.0107   0.00000    21  3 31.898  0.100   0.0000              
 ILDT=   0   0   0   0   0   0  1.43892335654173D+04     DERIV(1)=  1.07156807700680D-02  0.00000000000000D+00  3.16523184664084D+09 
  0.00000000000000D+00  0.00000000000000D+00 -9.07578858006227D-01  0.00000000000000D+00  5.78358124218615D+04  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -6.96413986072210D+00  0.00000000000000D+00 -1.12033265978963D+07 
  1.43510295687996D+04  1.48954543711465D+04  2.86465167777271D-08 -2.78815430383564D+04  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  4.51970275059672D-01 -2.64212263773569D-04  0.00000000000000D+00 
  0.00000000000000D+00
 ILDT=   0   0   0   0   0   0  7.58118977898652D+04     DERIV(2)=  1.00000001490116D-01  0.00000000000000D+00  5.07411802098999D+09 
  0.00000000000000D+00  0.00000000000000D+00 -3.31716831787614D+00  0.00000000000000D+00  8.54014118487596D+04  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  3.30968715498682D+02  0.00000000000000D+00 -1.74650372267227D+07 
  2.30133600682198D+04  2.38332643333683D+04 -2.45143306771511D-07 -6.81245397808429D+04  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -2.44667005391411D-01 -1.80761944550302D-02  0.00000000000000D+00 
  0.00000000000000D+00
 10/ 5/74 2442326  18  1 30.5000 13.0344  32.1483   4  6 12.5008  0.0107   0.00000    21 16 50.892  0.100   0.0000              
 10/ 5/74 2442326  20  1 30.5000 13.0344  32.1483   4 11 51.3373  0.0108   0.00000    21 33  7.130  0.100   0.0000              
 10/ 5/74 2442326  22  1 30.5000 13.0344  32.1483   4 16 31.3327  0.0108   0.00000    21 50  9.379  0.100   0.0000              
 10/ 6/74 2442327   0  1 30.5000 13.0344  32.1483   4 20 12.3333  0.0108   0.00000    22  5 16.936  0.100   0.0000              
 10/ 6/74 2442327   2  1 30.5000 13.0344  32.1483   4 23  9.8888  0.0108   0.00000    22 16  0.162  0.100   0.0000              
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   40
0 ASTROGRAPHIC    OBSERVATIONS OF   MOON   MADE FROM NICE               SERIES=ASTG  (SPOT=    )  FREQUENCY= 4.3178000000000D+14
0ERROR ANALYSIS FOR  ASTROGRAPHIC    OBSERVATION SERIES EXTENDING FOR   3 PAGES STARTING ON PAGE   38
 AVERAGE NUMBER OF ITERATIONS FOR     9 OBSERVATIONS WAS (  2.22222)
                               RIGHT ASCENSION    DECLINATION         
 NUMBER OF MEASUREMENTS DELETED          0             0
 NUMBER OF MEASUREMENTS INCLUDED         9             9
                AVERAGE (OBS-TH)   0.00000E+00   0.00000E+00
             AVERAGE ABS(OBS-TH)   0.00000E+00   0.00000E+00
       ROOT MEAN SQUARE (OBS-TH)   0.00000E+00   0.00000E+00
                   AVERAGE ERROR   1.07413E-02   1.00000E-01
          ROOT MEAN SQUARE ERROR   1.07414E-02   1.00000E-01
          AVERAGE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON= 2, IOBS= 0, IABS1= 0) FOR THIS OBSERVATION SERIES WERE (    0,    0,    0)
0     9 +   0 OBSERVATION CARDS PROCESSED (NUMBER OF PARTIAL DERIVATIVES= 32)
0  PROCESSING OBSERVATION SERIES REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
-ERROR ANALYSIS FOR  ASTROGRAPHIC    OBSERVATIONS OF   MOON   PROCESSED DURING THE LAST    3 PAGES STARTING ON PAGE   38
                               RIGHT ASCENSION    DECLINATION         
 NUMBER OF MEASUREMENTS DELETED          0             0
 NUMBER OF MEASUREMENTS INCLUDED         9             9
          AVERAGE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON,IOBS,IABS1) = (    0,    0,    0)

       9 +    0 OBSERVATION RECORDS PROCESSED (NPARAM=  44)
0   PROCESSING OBSERVATIONS  REQUIRED
           0H  0M  0.03S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   41
0DATA ON FIRST TWO RECORDS OF DATA SET 14   BODY=   MARS     NPLNT( 1)=  4    NCENTR= -1     IPAR= 13     IFILTR= 0
 TITLE=TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE=   19 ITERAT= 1 LEVEL= 790
     JD1=2442321          JD0=2442321          JD2=2442336          INT=  1     ICND=  0
   COND(1)= 1.523705402995026D+00   COND(2)= 9.330246794792663D-02   COND(3)= 2.469315159246902D+01   COND(4)= 3.344609353024999D+00
   COND(5)= 3.322243666631956D+02   COND(6)= 2.255191494289617D+02   CON( 1)= 3.392459000000000D+03   CON( 2)= 0.000000000000000D+00
   CON( 3)= 0.000000000000000D+00   CON( 4)= 0.000000000000000D+00   CON( 5)= 0.000000000000000D+00   CON( 6)= 1.514460000000000D+02
   CON( 7)= 1.025956000000000D+00   CON( 8)= 5.269510000000000D+01   CON( 9)= 3.173116000000000D+02   CON(10)= 0.000000000000000D+00
   CON(11)= 0.000000000000000D+00   CON(12)= 0.000000000000000D+00   CON(13)= 0.000000000000000D+00   CON(14)= 0.000000000000000D+00
   CON(15)= 0.000000000000000D+00   CON(16)= 0.000000000000000D+00   CON(17)= 0.000000000000000D+00   CON(18)= 0.000000000000000D+00
   CON(19)= 0.000000000000000D+00   CON(20)= 0.000000000000000D+00   CON(21)= 0.000000000000000D+00   CON(22)= 0.000000000000000D+00
   CON(23)= 0.000000000000000D+00   CON(24)= 0.000000000000000D+00   CON(25)= 2.443509500000000D+06   CON(26)= 0.000000000000000D+00
   CON(27)= 0.000000000000000D+00   CON(28)= 0.000000000000000D+00   CON(29)= 0.000000000000000D+00   CON(30)= 0.000000000000000D+00
   CON(31)= 0.000000000000000D+00   CON(32)= 0.000000000000000D+00   CON(33)= 0.000000000000000D+00   CON(34)=-8.400000000000000D-11
   CON(35)= 0.000000000000000D+00   CON(36)= 0.000000000000000D+00   CON(
  PRM(  1)= 1.657848020429993D-07  PRM(  2)= 2.447848585877872D-06  PRM(  3)= 3.040436898620584D-06  PRM(  4)= 3.227159776680543D-07
  PRM(  5)= 9.547534692876814D-04  PRM(  6)= 2.857796067672611D-04  PRM(  7)= 4.361098996947231D-05  PRM(  8)= 5.192107995846314D-05
  PRM(  9)= 2.500000000000000D-07  PRM( 10)= 1.215052064980984D-02  PRM( 11)= 6.700000000000000D-10  PRM( 12)= 0.000000000000000D+00
  PRM( 13)= 0.000000000000000D+00  PRM( 14)= 0.000000000000000D+00  PRM( 15)= 0.000000000000000D+00  PRM( 16)= 0.000000000000000D+00
  PRM( 17)= 0.000000000000000D+00  PRM( 18)= 0.000000000000000D+00  PRM( 19)= 0.000000000000000D+00  PRM( 20)= 0.000000000000000D+00
  PRM( 21)= 0.000000000000000D+00  PRM( 22)= 0.000000000000000D+00  PRM( 23)= 0.000000000000000D+00  PRM( 24)= 0.000000000000000D+00
  PRM( 25)= 0.000000000000000D+00  PRM( 26)= 0.000000000000000D+00  PRM( 27)= 0.000000000000000D+00  PRM( 28)= 0.000000000000000D+00
  PRM( 29)= 0.000000000000000D+00  PRM( 30)= 0.000000000000000D+00  PRM( 31)= 1.000000000000000D+00  PRM( 32)= 1.000000000000000D-35
  PRM( 33)= 1.000000000000000D-08  PRM( 34)= 0.000000000000000D+00  PRM( 35)= 0.000000000000000D+00  PRM( 36)= 0.000000000000000D+00
  PRM( 37)= 0.000000000000000D+00  PRM( 38)= 0.000000000000000D+00  PRM( 39)=-7.800000000000000D-04  PRM( 40)= 0.000000000000000D+00
  PRM( 41)= 1.000000000000000D+00  PRM( 42)= 1.000000000000000D+00  PRM( 43)= 1.000000000000000D+00  PRM( 44)= 1.000000000000000D+00
  PRM( 45)= 0.000000000000000D+00  PRM( 46)= 0.000000000000000D+00  PRM( 47)= 0.000000000000000D+00  PRM( 48)= 2.344330000000000D+01
  PRM( 49)= 2.900000000000000D+00  PRM( 50)= 1.000000000000000D-09  PRM( 51)= 4.990047800000000D+02  PRM( 52)= 0.000000000000000D+00
  PRM( 53)= 1.000000000000000D+00  PRM( 54)= 0.000000000000000D+00  PRM( 55)= 0.000000000000000D+00  PRM( 56)= 0.000000000000000D+00
  PRM( 57)= 0.000000000000000D+00  PRM( 58)= 0.000000000000000D+00  PRM( 59)= 0.000000000000000D+00  PRM( 60)= 1.000000000000000D+01
  PRM( 61)= 0.000000000000000D+00  PRM( 62)= 1.000000000000000D+00  PRM( 63)= 1.000000000000000D+00  PRM( 64)= 0.000000000000000D+00
  PRM( 65)= 0.000000000000000D+00  PRM( 66)= 0.000000000000000D+00  PRM( 67)= 0.000000000000000D+00  PRM( 68)= 0.000000000000000D+00
  PRM( 69)= 0.000000000000000D+00  PRM( 70)= 0.000000000000000D+00  PRM( 71)= 0.000000000000000D+00  PRM( 72)= 0.000000000000000D+00
  PRM( 73)= 0.000000000000000D+00  PRM( 74)= 0.000000000000000D+00  PRM( 75)= 0.000000000000000D+00  PRM( 76)= 0.000000000000000D+00
  PRM( 77)= 0.000000000000000D+00  PRM( 78)= 0.000000000000000D+00  PRM( 79)= 0.000000000000000D+00  PRM( 80)= 0.000000000000000D+00
  PRM( 81)= 1.000000000000000D+00  PRM( 82)= 0.000000000000000D+00  PRM( 83)= 0.000000000000000D+00  PRM( 84)= 0.000000000000000D+00
  PRM( 85)= 0.000000000000000D+00  PRM( 86)= 0.000000000000000D+00  PRM( 87)= 0.000000000000000D+00  PRM( 88)= 0.000000000000000D+00
  PRM( 89)= 0.000000000000000D+00  PRM( 90)= 0.000000000000000D+00  PRM( 91)= 2.344578706750000D+01  PRM( 92)= 7.250000000000000D+00
  PRM( 93)= 7.506250000000000D+01  PRM( 94)= 6.960000000000000D+05  PRM( 95)= 6.960000000000000D+05  PRM( 96)= 0.000000000000000D+00
  PRM( 97)= 2.440000500000000D+06  PRM( 98)= 4.263529034000000D-05  PRM( 99)= 0.000000000000000D+00  PRM(100)= 2.997924580000000D+05
  K(  1)=   5  K(  2)=   0  K(  3)=   0  K(  4)=   0  K(  5)=   0  K(  6)=   0  K(  7)=   0  K(  8)=   0  K(  9)=   0  K( 10)=   0
  K( 11)=   0  K( 12)=   0  K( 13)=   0  K( 14)=   0  K( 15)=   0  K( 16)=   0  K( 17)=   0  K( 18)=   0  K( 19)=   0  K( 20)=   0
  K( 21)=   0  K( 22)=   0  K( 23)=   0  K( 24)=   0  K( 25)=   0  K( 26)=   0  K( 27)=   0  K( 28)=   0  K( 29)=   0  K( 30)=   0
  K( 31)=  -1  K( 32)=  -1  K( 33)=   1  K( 34)=   1  K( 35)=   1  K( 36)=  -1  K( 37)=  -1  K( 38)=  -1  K( 39)=  -1  K( 40)=  -1
  K( 41)=   1  K( 42)=  -1  K( 43)=  -1  K( 44)=  -1  K( 45)=  -1  K( 46)=  -1  K( 47)=  -1  K( 48)=  -1  K( 49)=  -1  K( 50)=  -1
  K( 51)=  -1  K( 52)=  -1  K( 53)=  -1  K( 54)=  -1  K( 55)=  -1  K( 56)=  -1  K( 57)=  -1  K( 58)=  -1  K( 59)=  -1  K( 60)=  -1
  K( 61)=   1  K( 62)=  -1  K( 63)=  -1  K( 64)=  -1  K( 65)=  -1  K( 66)=  -1  K( 67)=  -1  K( 68)=  -1  K( 69)=  -1  K( 70)=   1
  K( 71)=   1  K( 72)=  -1  K( 73)=  -1  K( 74)=  -1  K( 75)=  -1  K( 76)=  -1  K( 77)=  -1  K( 78)=  -1  K( 79)=   1  K( 80)=   1
  K( 81)=  -1  K( 82)=  -1  K( 83)=  -1  K( 84)=  -1  K( 85)=  -1  K( 86)=  -1  K( 87)=   1  K( 88)=   3  K( 89)=  11  K( 90)=   6
  K( 91)= -10  K( 92)= -14  K( 93)=   0  K( 94)=   0  K( 95)=   0  K( 96)=   0  K( 97)=   0  K( 98)=   4  K( 99)=   0  K(100)=   0
 NUMKI= 16  KI= 1-1-1-1 1 1 1    3    5   11   31   40   41   49   50  506
   EPS(1)= 0.00000E+00   EPS(2)= 0.00000E+00   EPS(3)= 1.00000E-16   EPS(4)= 0.00000E+00   EPS(5)= 0.00000E+00   EPS(6)= 0.00000E+00
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   42

 NPLNT0= 4  NCODF=NCODA(1)= 1  NTAPA(1)= -5  NSEQA(1)= 1020  IABS1=IOBS0( 1)= 0  IABS2=IOBS1( 1)=31  IOBS= 0  IOBCON= 2  NREWND= 0
0 DEL/DOP RADAR   OBSERVATIONS OF   MARS   MADE FROM ARECIBO ,ARECIBO   SERIES=RADR  (SPOT=OLYM)  FREQUENCY= 2.3800000000000D+09
0                 NAME  NSITE LSITE     RADIUS       LONGITUDE      LATITUDE  KSITE
   RECEIVE SITE=ARECIBO    7  0 0 0 6376.560245971   66.75302778    18.22876139  0
      SEND SITE=ARECIBO    7  0 0 0 6376.560245971   66.75302778    18.22876139  0
0ACCPRC= 1.000E+01  ACCDST= 1.000E-09  ACCTIM= 1.000E-09  ERWGT1= 1.000E+00  ERWGT2= 1.000E+00  FDEV=   000.0D+00
0  NSPOT = 15   NPREC =  1   NDPREC=  0   NTMDLY=  1   NDOP  = -1   NMEDIA=  0   NEATM = -1   NPATM = -1   NPSHP = -1   NPLNG = -1
   NEION = -1  KOB(11)= -1   LOPTRF=  1   NDDIFF=  0   NLIBPR=  2  KOB(15)= -1   CALCVL= -1   ITIME =  1   NTIME =  0    IWOB =  1
0  SETUP FOR  DEL/DOP RADAR   OBSERVATION SERIES REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
0 GRNWCH  JULIAN  REC UTC (WWV)  UT1-UTC AT-UTC  CT-AT             TIME DELAY (SEC)               DOPPLER SHIFT (CYCLES/SECOND) 
   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC       OBSERVED     ERROR     OBS-TH         OBSERVED     ERROR     OBS-TH 
 SUMCOR(1)=    4.3555974116E-06
 10/ 5/74 2442326 20  1 30.5000 -0.0364 13.0344 32.1483 2617.0545447399 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  2.61705454473989D+03     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00  5.48707522428892D+03 
  1.16381749555061D+00  6.39552465463359D-02 -1.66381283151416D+00  1.84315238656524D+01 -3.13645685767279D+01 -7.10335175973525D+00 
 -2.52060370055994D-07  2.45379403751976D-02  1.56624577075177D+01 -3.51915733141158D-10  2.54228623159737D-01  2.61709628210117D+03 
  4.35559721626144D-07  9.91613603406313D+02 -9.17771273537662D+01 -2.20223612383447D+00  1.48534615680467D+01 -2.71619637505828D+01 
  2.34344161616905D-02  2.52596528891836D-02 -4.98633776984807D-14 -8.89838800456253D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -4.97243401371669D+01 -4.93268385445934D+01 -1.36926819922801D+02 
 -7.10924244884674D+00  8.41151519624445D+03  2.80868704333887D+00 -6.53365498238249D+00  3.57806229760567D+00 -0.00000000000000D+00 
 -0.00000000000000D+00 -0.00000000000000D+00  9.39253667733500D-05
 SUMCOR(1)=    4.3969989747E-06
 10/ 5/74 2442326 22  1 30.5000 -0.0366 13.0344 32.1483 2617.0110268180 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  2.61701102681803D+03     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00  5.65038696671143D+03 
  1.19750920510079D+00  6.57439193678879D-02 -1.62697871748024D+00  1.89630489111033D+01 -3.22690899998327D+01 -7.30846501938760D+00 
 -2.66649825504063D-07  2.52461189739496D-02  1.56993387508726D+01 -3.62056211147657D-10  2.61556706159775D-01  2.61703945454075D+03 
  4.39699881816021D-07  9.91615913875774D+02 -9.22444821206790D+01 -2.21467286832819D+00  1.52818066037856D+01 -2.75646562636471D+01 
  2.37966538130346D-02  2.56502735487983D-02 -5.14229506050343D-14 -9.18182408129797D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -4.92948011413439D+01 -4.88645092009647D+01 -1.36598063730386D+02 
 -2.89222255413032D+00  3.42177956029176D+03  1.59019315486823D+00 -2.06742767675962D+00  3.68124230731769D+00 -0.00000000000000D+00 
 -0.00000000000000D+00 -0.00000000000000D+00  9.67213494114564D-05
 SUMCOR(1)=    4.4397820602E-06
 10/ 6/74 2442327  0  1 30.5000 -0.0369 13.0344 32.1483 2616.9740739739 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  2.61697407397387D+03     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00  5.81436354237600D+03 
  1.23168538714696D+00  6.75553936948752D-02 -1.58959968304674D+00  1.95021216437208D+01 -3.31864888961787D+01 -7.51649894898038D+00 
 -2.81863442871977D-07  2.59643808470679D-02  1.57370929275390D+01 -3.72340202532751D-10  2.68988615490542D-01  2.61698221467863D+03 
  4.43978183324389D-07  9.91619369113174D+02 -9.27182437615424D+01 -2.22759908549644D+00  1.57162358244083D+01 -2.79671290485909D+01 
  2.41511816463906D-02  2.60334465952596D-02 -5.30195358588101D-14 -9.46931812490381D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -4.88588529912904D+01 -4.83970480836334D+01 -1.36264713013701D+02 
  2.05965605269336D+00 -2.43660106850489D+03 -3.68910544434725D-01  2.47968111498534D+00  3.78588581931245D+00 -0.00000000000000D+00 
 -0.00000000000000D+00 -0.00000000000000D+00  9.95609922145573D-05
 SUMCOR(1)=    4.4838684516E-06
 10/ 6/74 2442327  2  1 30.5000 -0.0371 13.0344 32.1483 2616.9382351344 1.0E-08  0.00000E+00 
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   43
0 DEL/DOP RADAR   OBSERVATIONS OF   MARS   MADE FROM ARECIBO ,ARECIBO   SERIES=RADR  (SPOT=OLYM)  FREQUENCY= 2.3800000000000D+09
0 GRNWCH  JULIAN  REC UTC (WWV)  UT1-UTC AT-UTC  CT-AT             TIME DELAY (SEC)               DOPPLER SHIFT (CYCLES/SECOND) 
   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC       OBSERVED     ERROR     OBS-TH         OBSERVED     ERROR     OBS-TH 
 ILDT=   0   0   0   0   0   0  2.61693823513435D+03     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00  5.97889153110867D+03 
  1.26634623489243D+00  6.93895230415396D-02 -1.55168935846688D+00  2.00487391137988D+01 -3.41167645039473D+01 -7.72745321949975D+00 
 -2.97718713112586D-07  2.66927238707422D-02  1.57756405673281D+01 -3.82767688725398D-10  2.76524339994217D-01  2.61692456254038D+03 
  4.48386866222877D-07  9.91623751767654D+02 -9.32004198010346D+01 -2.24093903525829D+00  1.61567468271228D+01 -2.83690457560114D+01 
  2.44978555604845D-02  2.64089899701702D-02 -5.46534042355988D-14 -9.76074885552056D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -4.84152554851717D+01 -4.79228089135633D+01 -1.35925326520145D+02 
  6.48937280784647D+00 -7.67647533722216D+03 -2.56818492610857D+00  5.95297630898323D+00  3.89199228667604D+00 -0.00000000000000D+00 
 -0.00000000000000D+00 -0.00000000000000D+00  1.02444503799808D-04
 SUMCOR(1)=    4.5290362323E-06
 10/ 6/74 2442327  4  1 30.5000 -0.0373 13.0344 32.1483 2616.8977215209 1.0E-08  0.00000E+00 
 SUMCOR(1)=    4.5749752644E-06
 10/ 6/74 2442327  6  1 30.5000 -0.0376 13.0344 32.1483 2616.8479439992 1.0E-08  0.00000E+00 
 SUMCOR(1)=    4.6213767746E-06
 10/ 6/74 2442327  8  1 30.5000 -0.0378 13.0344 32.1483 2616.7867142540 1.0E-08  0.00000E+00 
 SUMCOR(1)=    4.6680229389E-06
 10/ 6/74 2442327 10  1 30.5000 -0.0381 13.0344 32.1483 2616.7147952432 1.0E-08  0.00000E+00 
 SUMCOR(1)=    4.7148578233E-06
 10/ 6/74 2442327 12  1 30.5000 -0.0383 13.0344 32.1483 2616.6356625950 1.0E-08  0.00000E+00 
0ERROR ANALYSIS FOR  DEL/DOP RADAR   OBSERVATION SERIES EXTENDING FOR   2 PAGES STARTING ON PAGE   42
 AVERAGE NUMBER OF ITERATIONS FOR     9 OBSERVATIONS WAS (  3.11111,  3.00000)
                                   TIME DELAY      DOPPLER SHIFT      
 NUMBER OF MEASUREMENTS DELETED          0             0
 NUMBER OF MEASUREMENTS INCLUDED         9             0
                AVERAGE (OBS-TH)   0.00000E+00   0.00000E+00
             AVERAGE ABS(OBS-TH)   0.00000E+00   0.00000E+00
       ROOT MEAN SQUARE (OBS-TH)   0.00000E+00   0.00000E+00
                   AVERAGE ERROR   1.00000E-08   0.00000E+00
          ROOT MEAN SQUARE ERROR   1.00000E-08   0.00000E+00
          AVERAGE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON= 2, IOBS= 0, IABS1= 0) FOR THIS OBSERVATION SERIES WERE (    0,    0,    0)
0     9 +   0 OBSERVATION CARDS PROCESSED (NUMBER OF PARTIAL DERIVATIVES= 40)
0  PROCESSING OBSERVATION SERIES REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   44
-ERROR ANALYSIS FOR  DEL/DOP RADAR   OBSERVATIONS OF   MARS   PROCESSED DURING THE LAST    4 PAGES STARTING ON PAGE   41
                                   TIME DELAY      DOPPLER SHIFT      
 NUMBER OF MEASUREMENTS DELETED          0             0
 NUMBER OF MEASUREMENTS INCLUDED         9             0
          AVERAGE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON,IOBS,IABS1) = (    0,    0,    0)

       9 +    0 OBSERVATION RECORDS PROCESSED (NPARAM=  44)
0   PROCESSING OBSERVATIONS  REQUIRED
           0H  0M  0.03S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   45

 NPLNT0= 4  NCODF=NCODA(1)= 3  NTAPA(1)= -5  NSEQA(1)= 1030  IABS1=IOBS0( 1)= 0  IABS2=IOBS1( 1)=31  IOBS= 0  IOBCON= 2  NREWND= 0
0 DEL/DOP RADAR   OBSERVATIONS OF   MARS   MADE FROM ARECIBO ,ARECIBO   SERIES=RADR  (SPOT=OLYM)  FREQUENCY= 2.3800000000000D+09
0                 NAME  NSITE LSITE     RADIUS       LONGITUDE      LATITUDE  KSITE
   RECEIVE SITE=ARECIBO    7  0 0 0 6376.560245971   66.75302778    18.22876139  0
      SEND SITE=ARECIBO    7  0 0 0 6376.560245971   66.75302778    18.22876139  0
0ACCPRC= 1.000E+01  ACCDST= 1.000E-09  ACCTIM= 1.000E-09  ERWGT1= 1.000E+00  ERWGT2= 1.000E+00  FDEV=   000.0D+00
0  NSPOT = 15   NPREC =  1   NDPREC=  0   NTMDLY= -1   NDOP  = -1   NMEDIA=  0   NEATM = -1   NPATM = -1   NPSHP = -1   NPLNG =  1
   NEION = -1  KOB(11)= -1   LOPTRF=  1   NDDIFF=  0   NLIBPR=  2  KOB(15)= -1   CALCVL= -1   ITIME =  1   NTIME =  0    IWOB =  1
0  SETUP FOR  DEL/DOP RADAR   OBSERVATION SERIES REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
0 GRNWCH  JULIAN  REC UTC (WWV)  UT1-UTC AT-UTC  CT-AT             TIME DELAY (SEC)               DOPPLER SHIFT (CYCLES/SECOND) 
   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC       OBSERVED     ERROR     OBS-TH         OBSERVED     ERROR     OBS-TH 
 10/ 5/74 2442326 20  1 30.5000 -0.0364 13.0344 32.1483    0.0061946482 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  6.19464817524173D-03     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00  5.48707522428891D+03 
  1.16381749555061D+00  6.39552465463359D-02 -1.66381283151416D+00  1.84315238656524D+01 -3.13645685767279D+01 -7.10335175973525D+00 
 -2.52060370055994D-07  2.45379403751976D-02  9.91993343582204D-03 -3.51915733141158D-10  2.54228623159737D-01  2.61709628210117D+03 
  0.00000000000000D+00  9.91613603406313D+02 -9.17771273537669D+01 -2.20223612383449D+00  1.48534615680467D+01 -2.71619637505828D+01 
  2.34344161616905D-02  2.52596528891836D-02 -4.98633776984808D-14 -8.89838800456253D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -4.97243401371662D+01 -4.93268385445928D+01 -1.36926819922800D+02 
 -7.10924244884674D+00  8.41151519624445D+03  2.80868704333887D+00 -6.53365498238249D+00  3.57806229760567D+00 -7.24081938180805D+02 
 -7.10924244884674D+00 -2.18488167082606D+00  9.39253667733500D-05
 10/ 5/74 2442326 22  1 30.5000 -0.0366 13.0344 32.1483    0.0009769870 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  9.76987013320052D-04     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00  5.65038696671143D+03 
  1.19750920510079D+00  6.57439193678879D-02 -1.62697871748024D+00  1.89630489111033D+01 -3.22690899998327D+01 -7.30846501938760D+00 
 -2.66649825504063D-07  2.52461189739496D-02  1.02059023092503D-02 -3.62056211147657D-10  2.61556706159775D-01  2.61703945454075D+03 
  0.00000000000000D+00  9.91615913875774D+02 -9.22444821206790D+01 -2.21467286832819D+00  1.52818066037856D+01 -2.75646562636471D+01 
  2.37966538130346D-02  2.56502735487983D-02 -5.14229506050343D-14 -9.18182408129797D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -4.92948011413439D+01 -4.88645092009647D+01 -1.36598063730386D+02 
 -2.89222255413032D+00  3.42177956029176D+03  1.59019315486823D+00 -2.06742767675962D+00  3.68124230731769D+00 -9.53926796884012D+02 
 -2.89222255413032D+00 -1.34651330752751D+00  9.67213494114564D-05
 10/ 6/74 2442327  0  1 30.5000 -0.0369 13.0344 32.1483    0.0005422088 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  5.42208762946750D-04     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00  5.81436354237600D+03 
  1.23168538714696D+00  6.75553936948752D-02 -1.58959968304674D+00  1.95021216437208D+01 -3.31864888961787D+01 -7.51649894898038D+00 
 -2.81863442871977D-07  2.59643808470679D-02  1.04959284732061D-02 -3.72340202532751D-10  2.68988615490542D-01  2.61698221467863D+03 
  0.00000000000000D+00  9.91619369113174D+02 -9.27182437615424D+01 -2.22759908549644D+00  1.57162358244083D+01 -2.79671290485909D+01 
  2.41511816463906D-02  2.60334465952596D-02 -5.30195358588101D-14 -9.46931812490381D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -4.88588529912904D+01 -4.83970480836334D+01 -1.36264713013701D+02 
  2.05965605269336D+00 -2.43660106850489D+03 -3.68910544434725D-01  2.47968111498534D+00  3.78588581931245D+00 -9.73079373170895D+02 
  2.05965605269336D+00 -1.27495865442448D+00  9.95609922145573D-05
 10/ 6/74 2442327  2  1 30.5000 -0.0371 13.0344 32.1483    0.0050015863 1.0E-08  0.00000E+00 
 ILDT=   0   0   0   0   0   0  5.00158630832037D-03     DERIV(1)=  9.99999993922529D-09  0.00000000000000D+00  5.97889153110867D+03 
  1.26634623489243D+00  6.93895230415396D-02 -1.55168935846688D+00  2.00487391137988D+01 -3.41167645039473D+01 -7.72745321949975D+00 
 -2.97718713112586D-07  2.66927238707422D-02  1.07900100492004D-02 -3.82767688725398D-10  2.76524339994217D-01  2.61692456254038D+03 
  0.00000000000000D+00  9.91623751767654D+02 -9.32004198010346D+01 -2.24093903525829D+00  1.61567468271228D+01 -2.83690457560114D+01 
  2.44978555604845D-02  2.64089899701702D-02 -5.46534042355988D-14 -9.76074885552056D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -4.84152554851717D+01 -4.79228089135633D+01 -1.35925326520145D+02 
  6.48937280784647D+00 -7.67647533722216D+03 -2.56818492610857D+00  5.95297630898323D+00  3.89199228667604D+00 -7.76637847136683D+02 
  6.48937280784647D+00 -1.98806105985636D+00  1.02444503799808D-04
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   46
0 DEL/DOP RADAR   OBSERVATIONS OF   MARS   MADE FROM ARECIBO ,ARECIBO   SERIES=RADR  (SPOT=OLYM)  FREQUENCY= 2.3800000000000D+09
0 GRNWCH  JULIAN  REC UTC (WWV)  UT1-UTC AT-UTC  CT-AT             TIME DELAY (SEC)               DOPPLER SHIFT (CYCLES/SECOND) 
   DATE   DAY NUM HR MIN  SEC     SEC     SEC     SEC       OBSERVED     ERROR     OBS-TH         OBSERVED     ERROR     OBS-TH 
 10/ 6/74 2442327  4  1 30.5000 -0.0373 13.0344 32.1483    0.0132239950 1.0E-08  0.00000E+00 
0ERROR ANALYSIS FOR  DEL/DOP RADAR   OBSERVATION SERIES EXTENDING FOR   2 PAGES STARTING ON PAGE   45
 AVERAGE NUMBER OF ITERATIONS FOR     5 OBSERVATIONS WAS (  3.20000,  3.00000)
                                   TIME DELAY      DOPPLER SHIFT      
 NUMBER OF MEASUREMENTS DELETED          0             0
 NUMBER OF MEASUREMENTS INCLUDED         5             0
                AVERAGE (OBS-TH)   0.00000E+00   0.00000E+00
             AVERAGE ABS(OBS-TH)   0.00000E+00   0.00000E+00
       ROOT MEAN SQUARE (OBS-TH)   0.00000E+00   0.00000E+00
                   AVERAGE ERROR   1.00000E-08   0.00000E+00
          ROOT MEAN SQUARE ERROR   1.00000E-08   0.00000E+00
          AVERAGE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON= 2, IOBS= 0, IABS1= 0) FOR THIS OBSERVATION SERIES WERE (    0,    0,    0)
0     5 +   0 OBSERVATION CARDS PROCESSED (NUMBER OF PARTIAL DERIVATIVES= 40)
0  PROCESSING OBSERVATION SERIES REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
-ERROR ANALYSIS FOR  DEL/DOP RADAR   OBSERVATIONS OF   MARS   PROCESSED DURING THE LAST    2 PAGES STARTING ON PAGE   45
                                   TIME DELAY      DOPPLER SHIFT      
 NUMBER OF MEASUREMENTS DELETED          0             0
 NUMBER OF MEASUREMENTS INCLUDED         5             0
          AVERAGE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON,IOBS,IABS1) = (    0,    0,    0)

       5 +    0 OBSERVATION RECORDS PROCESSED (NPARAM=  44)
0   PROCESSING OBSERVATIONS  REQUIRED
           0H  0M  0.03S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   47

 NPLNT0= 4  NCODF=NCODA(1)= 5  NTAPA(1)= -5  NSEQA(1)= 1040  IABS1=IOBS0( 1)= 0  IABS2=IOBS1( 1)=31  IOBS= 0  IOBCON= 2  NREWND= 0
0  ASTROMETRIC    OBSERVATIONS OF   MARS   MADE FROM PARIS              SERIES=ASTM  (SPOT=    )  FREQUENCY= 4.3178000000000D+14
0                 NAME  NSITE LSITE     RADIUS       LONGITUDE      LATITUDE  KSITE
   RECEIVE SITE=PARIS     28  0 0 0 6378.000000000   -2.33712500     0.00000000  0
0ACCPRC= 1.000E+01  ACCDST= 1.000E-09  ACCTIM= 1.000E-09  ERWGT1= 1.000E+00  ERWGT2= 1.000E+00  FDEV=   000.0D+00
0  NSPOT =  0   NPREC = -1   NDPREC= -1   NTMDLY= -1   NDOP  = -1   NMEDIA= -1   NEATM = -1   NPATM = -1   NPSHP = -1   NPLNG = -1
   NEION = -1  KOB(11)= -1   LOPTRF= -1   NDDIFF= -1   NLIBPR= -1  KOB(15)= -1   CALCVL= -1   ITIME =  1   NTIME =  0    IWOB =  0
0TOPOCENTRIC RIGHT ASCENSION-DECLINATION REFERRED TO THE MEAN EQUINOX AND EQUATOR OF 1950.0
0  SETUP FOR   ASTROMETRIC    OBSERVATION SERIES REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
0 GRNWCH   JULIAN   UT  TIME OF   AT-UT    CT-AT           RIGHT ASCENSION                    DECLINATION                       
   DATE     DAY     OBSERVATION                      OBSERVED     ERROR    OBS-TH      OBSERVED     ERROR   OBS-TH              
           NUMBER  HR MIN  SEC     SEC      SEC    HR MIN  SEC     SEC      SEC      DEG  '   ''     ''      ''    CLMP LIMB OBS
          2442326 IS NOT ON PLANET DATA SET (2442326-2442331)
 10/ 6/74 2442327   2  1 30.5000 13.0344  32.1483  12 55 40.4263  0.0100   0.00000   - 5 17 13.710  0.100   0.0000              
 ILDT=   0   0   0   0   0   0  4.65404262986885D+04     DERIV(1)=  1.00427273588703D-02  0.00000000000000D+00 -2.14347379221108D+04 
  1.81786488029860D-01 -2.83450744621463D-01 -9.75058666867348D+00 -1.29093910781129D+01  2.07474381009570D+01  8.94260468830511D+00 
  2.45789667732113D-06 -1.29813059926790D-02 -8.26032040148548D-03  6.67669448924047D-11 -4.06986817268139D-02  0.00000000000000D+00 
  0.00000000000000D+00 -1.01802818026053D+03 -9.70348847789654D+03 -4.26127293226741D+02 -1.14822640032189D+01  3.06893556412249D+01 
 -9.48488407307056D-02 -9.42871216687709D-02 -8.46739457668539D-13 -2.18287281375470D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  8.51014399959033D+03  7.79586177538714D+03  6.81487320540693D+03 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -1.42712707489405D+00  4.27953808075970D-04
 ILDT=   0   0   0   0   0   0 -1.90337097262179D+04     DERIV(2)=  1.00000001490116D-01  0.00000000000000D+00  7.37639792275385D+04 
  1.36242662716974D-01  1.86600969488513D+00  6.37828085584656D+01  6.89734274521359D+01 -9.95067166536462D+01  1.11853545157456D+02 
 -1.55830523308852D-05  6.44427275289637D-02  4.17002114851848D-02  2.16981111768206D-09 -1.29903080316143D+00  0.00000000000000D+00 
  0.00000000000000D+00  5.45566001708934D+03  6.07616007762132D+04 -1.52029842630261D+04  5.68868018649336D+01  3.39703544339530D+01 
  3.40888594528286D-01  3.24356553177838D-01  8.02500142821530D-13  4.41742780630996D-01  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -3.51780516458064D+02 -5.24764491893238D+04 -4.59331514088186D+04 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  1.20866255872022D+01 -2.58937784534603D-03
 10/ 6/74 2442327   4  1 30.5000 13.0344  32.1483  12 55 52.6588  0.0100   0.00000   - 5 18 32.399  0.100   0.0000              
 ILDT=   0   0   0   0   0   0  4.65526588103143D+04     DERIV(1)=  1.00430826515864D-02  0.00000000000000D+00 -2.12010853041344D+04 
  1.81673396760395D-01 -2.91646859703546D-01 -9.89135130867386D+00 -1.33210854871110D+01  2.14177772847888D+01  9.20773610223238D+00 
  2.49104847015723D-06 -1.34121973893114D-02 -8.51531671227422D-03  6.99970804898581D-11 -4.28283402394025D-02  0.00000000000000D+00 
  0.00000000000000D+00 -1.02629734735966D+03 -9.70426936890217D+03 -4.28900075454498D+02 -1.18394869633144D+01  2.92152218265270D+01 
 -9.25107076870942D-02 -9.19772726962423D-02 -8.73497001942708D-13 -2.58649498919887D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  8.51015410578128D+03  7.79635543883022D+03  6.81603454260009D+03 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -1.48159852379658D+00  4.39308499333885D-04
 ILDT=   0   0   0   0   0   0 -1.91123993054692D+04     DERIV(2)=  1.00000001490116D-01  0.00000000000000D+00  7.12967210892271D+04 
  1.72563540264664D-01  1.91943442272793D+00  6.42320492608075D+01  7.12158235446126D+01 -1.02847067694717D+02  1.14797165645976D+02 
 -1.57599208894145D-05  6.66720310989073D-02  4.30142118225648D-02  2.21999995736817D-09 -1.32797553337298D+00  0.00000000000000D+00 
  0.00000000000000D+00  5.50631608133437D+03  6.07489538008056D+04 -1.53053417219039D+04  5.87018999380525D+01  4.47143080896556D+01 
  3.25650834845933D-01  3.09082605098689D-01  8.78570974943823D-13  4.75637779226486D-01  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -3.49957309343130D+02 -5.24670334421509D+04 -4.59297473097551D+04 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  1.25139236065601D+01 -2.65716171100450D-03
 10/ 6/74 2442327   6  1 30.5000 13.0344  32.1483  12 56  4.8395  0.0100   0.00000   - 5 19 51.119  0.100   0.0000              
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   48
0  ASTROMETRIC    OBSERVATIONS OF   MARS   MADE FROM PARIS              SERIES=ASTM  (SPOT=    )  FREQUENCY= 4.3178000000000D+14
0 GRNWCH   JULIAN   UT  TIME OF   AT-UT    CT-AT           RIGHT ASCENSION                    DECLINATION                       
   DATE     DAY     OBSERVATION                      OBSERVED     ERROR    OBS-TH      OBSERVED     ERROR   OBS-TH              
           NUMBER  HR MIN  SEC     SEC      SEC    HR MIN  SEC     SEC      SEC      DEG  '   ''     ''      ''    CLMP LIMB OBS
 ILDT=   0   0   0   0   0   0  4.65648395246387D+04     DERIV(1)=  1.00434397731506D-02  0.00000000000000D+00 -2.09440787073502D+04 
  1.81375163932433D-01 -2.99968384168249D-01 -1.00287860705972D+01 -1.37402635718081D+01  2.21007627384142D+01  9.47708644151137D+00 
  2.52397780993603D-06 -1.38514030908570D-02 -8.77468304759811D-03  7.33201016904931D-11 -4.50225020309953D-02  0.00000000000000D+00 
  0.00000000000000D+00 -1.03455099359522D+03 -9.70505883170342D+03 -4.31672130973660D+02 -1.22029129010557D+01  2.76743301457919D+01 
 -9.01375119406811D-02 -8.96270513779680D-02 -9.00933724213770D-13 -3.01731846862533D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  8.51018123991933D+03  7.79686507411218D+03  6.81720877769520D+03 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -1.53735067075237D+00  4.50801086976969D-04
 ILDT=   0   0   0   0   0   0 -1.91911185359786D+04     DERIV(2)=  1.00000001490116D-01  0.00000000000000D+00  6.87046433470168D+04 
  2.10722979855720D-01  1.97366528837997D+00  6.46589648210923D+01  7.35029799751748D+01 -1.06258513707616D+02  1.17775124680299D+02 
 -1.59338831251522D-05  6.89505845711304D-02  4.43529881780288D-02  2.27047349441213D-09 -1.35703004741560D+00  0.00000000000000D+00 
  0.00000000000000D+00  5.55698313592345D+03  6.07363395237673D+04 -1.54077141316551D+04  6.05520568373359D+01  5.57938206908667D+01 
  3.10301040770205D-01  2.93667294437498D-01  9.59480306473464D-13  5.11179455662557D-01  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -3.48154179259645D+02 -5.24577247726508D+04 -4.59264414918460D+04 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  1.29509231378389D+01 -2.72571593152131D-03
 10/ 6/74 2442327   8  1 30.5000 13.0344  32.1483  12 56 16.9619  0.0100   0.00000   - 5 21  9.824  0.100   0.0000              
 ILDT=   0   0   0   0   0   0  4.65769618729472D+04     DERIV(1)=  1.00437981226714D-02  0.00000000000000D+00 -2.06634334744101D+04 
  1.80891280087422D-01 -3.08414750932371D-01 -1.01628137866062D+01 -1.41668688719821D+01  2.27963179979236D+01  9.75062215192723D+00 
  2.55666793250204D-06 -1.42988562575580D-02 -9.03838283508858D-03  7.67357077702238D-11 -4.72809716774160D-02  0.00000000000000D+00 
  0.00000000000000D+00 -1.04278434166307D+03 -9.70583537518487D+03 -4.34442565962377D+02 -1.25724890857472D+01  2.60660369467417D+01 
 -8.77295750378002D-02 -8.72367299813013D-02 -9.29071412519382D-13 -3.47619149200625D-02  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  8.51020640816779D+03  7.79737310645878D+03  6.81838030546287D+03 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -1.59437978623483D+00  4.62430495570018D-04
 ILDT=   0   0   0   0   0   0 -1.92698241193038D+04     DERIV(2)=  1.00000001490116D-01  0.00000000000000D+00  6.59874014818689D+04 
  2.50734617464111D-01  2.02869902652106D+00  6.50632690964883D+01  7.58347848902464D+01 -1.09740959352776D+02  1.20787137404905D+02 
 -1.61048095490739D-05  7.12782974240598D-02  4.57164498349474D-02  2.32122264043268D-09 -1.38618849603931D+00  0.00000000000000D+00 
  0.00000000000000D+00  5.60763350569804D+03  6.07236328546820D+04 -1.55100674252356D+04  6.24371505141732D+01  6.72087741736186D+01 
  2.94842560340198D-01  2.78114068118115D-01  1.04543368315881D-12  5.48401937821153D-01  0.00000000000000D+00  0.00000000000000D+00 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00 -3.46375782473936D+02 -5.24484125657123D+04 -4.59231355678017D+04 
  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  0.00000000000000D+00  1.33976343760733D+01 -2.79503242492770D-03
 10/ 6/74 2442327  10  1 30.5000 13.0344  32.1483  12 56 29.0355  0.0100   0.00000   - 5 22 28.473  0.100   0.0000              
 10/ 6/74 2442327  12  1 30.5000 13.0344  32.1483  12 56 41.0835  0.0100   0.00000   - 5 23 47.036  0.100   0.0000              
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   49
0  ASTROMETRIC    OBSERVATIONS OF   MARS   MADE FROM PARIS              SERIES=ASTM  (SPOT=    )  FREQUENCY= 4.3178000000000D+14
0ERROR ANALYSIS FOR   ASTROMETRIC    OBSERVATION SERIES EXTENDING FOR   3 PAGES STARTING ON PAGE   47
 AVERAGE NUMBER OF ITERATIONS FOR     6 OBSERVATIONS WAS (  3.16667)
                               RIGHT ASCENSION    DECLINATION         
 NUMBER OF MEASUREMENTS DELETED          0             0
 NUMBER OF MEASUREMENTS INCLUDED         6             6
                AVERAGE (OBS-TH)   0.00000E+00   0.00000E+00
             AVERAGE ABS(OBS-TH)   0.00000E+00   0.00000E+00
       ROOT MEAN SQUARE (OBS-TH)   0.00000E+00   0.00000E+00
                   AVERAGE ERROR   1.00436E-02   1.00000E-01
          ROOT MEAN SQUARE ERROR   1.00436E-02   1.00000E-01
          AVERAGE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON= 2, IOBS= 0, IABS1= 0) FOR THIS OBSERVATION SERIES WERE (    0,    0,    0)
0     6 +   1 OBSERVATION CARDS PROCESSED (NUMBER OF PARTIAL DERIVATIVES= 37)
0  PROCESSING OBSERVATION SERIES REQUIRED
           0H  0M  0.01S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
-ERROR ANALYSIS FOR   ASTROMETRIC    OBSERVATIONS OF   MARS   PROCESSED DURING THE LAST    3 PAGES STARTING ON PAGE   47
                               RIGHT ASCENSION    DECLINATION         
 NUMBER OF MEASUREMENTS DELETED          0             0
 NUMBER OF MEASUREMENTS INCLUDED         6             6
          AVERAGE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON,IOBS,IABS1) = (    0,    0,    0)

       6 +    1 OBSERVATION RECORDS PROCESSED (NPARAM=  44)
0   PROCESSING OBSERVATIONS  REQUIRED
           0H  0M  0.03S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1COMPARISON OF THEORY AND OBS (ITERAT  1) TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   50
-ERROR ANALYSIS FOR ALL OBSERVATIONS PROCESSED IN THE COMPAR LINK DURING THE PAST   21 PAGES STARTING ON PAGE   30

 NUMBER OF MEASUREMENTS DELETED          0
 NUMBER OF MEASUREMENTS INCLUDED        53
          AVERAGE (OBS-TH)/ERROR   0.00000E+00
       AVERAGE ABS(OBS-TH)/ERROR   0.00000E+00
 ROOT MEAN SQUARE (OBS-TH)/ERROR   0.00000E+00
     AVERAGE ((OBS-TH)/ERROR)**2   0.00000E+00
         SUM ((OBS-TH)/ERROR)**2   0.00000E+00
 NUMBER OF ERROR RECORDS SKIPPED ON (IOBCON,IOBS,IABS1) = (    0,    0,    0)

      38 +    1 OBSERVATION RECORDS PROCESSED (NPARAM=  44)
0   PROCESSING OBSERVATIONS  REQUIRED
           0H  0M  0.16S REAL TIME
           0H  0M  0.00S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.00000
1LEAST SQUARES ANALYSIS   (ITERAT  1)     TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   51













                            *********    **           **    *********    **            **   *************
                           ***********   ***          **   ***********   **            **   *************
                           **       **   ****         **   **       **   **            **             **
                           **       **   ** **        **   **       **   **            **            **
                           **       **   **  **       **   **       **   **            **           **
                           **       **   **   **      **   **       **   **            **          **
                           ***********   **    **     **   ***********   **            **         **
                           ***********   **     **    **   ***********   **            **        **
                           **       **   **      **   **   **       **   **            **       **
                           **       **   **       **  **   **       **   **            **      **
                           **       **   **        ** **   **       **   **            **     **
                           **       **   **         ****   **       **   **            **    **
                           **       **   **          ***   **       **   ***********   **   *************
                           **       **   **           **   **       **   ***********   **   *************
-
-
                                          **           **    **           **   **       **
                                          **           **    ***          **   **      ** 
                                          **           **    ****         **   **     **  
                                          **           **    ** **        **   **    **   
                                          **           **    **  **       **   **   **    
                                          **           **    **   **      **   *** **     
                                          **           **    **    **     **   *****      
                                          **           **    **     **    **   ** **      
                                          **           **    **      **   **   **  **     
                                          **           **    **       **  **   **   **    
                                          **           **    **        ** **   **    **   
                                          **           **    **         ****   **     **  
                                          ***********  **    **          ***   **      ** 
                                          ***********  **    **           **   **       **
- PEP VERSION= 20210302 PEP.PEPLOAD.PEP790                          


 ANALIZ LINK CANCELLED BY ICT(15)= -1
1COMPLETION OF LEAST SQUARES ITERATION 1  TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   52
0    LEAST SQUARES ITERATION  1 DURING THE LAST   41 PAGES STARTING ON PAGE   12 REQUIRED
           0H  0M  0.33S REAL TIME
           0H  0M  0.06S TASK TIME
  (TASK TIME)/(REAL TIME)= 0.18182
0NORMAL STOP  NO MATRIX INVERSION  ICT(15)=   -1
1PLANETARY EPHEMERIS PROGRAM TERMINATED   TPL -  TEST INTEGRATION OF EARTH, MOON, AND PLANET WITH PARTIALS          3/ 2/21 PAGE   53

  NORMAL STOP IN MAIN
0 TERMINATION DATE       =  3/ 2/21
  TERMINATION DAY OF YEAR=  61
  TERMINATION CLOCK TIME = 12H 30M 42.36S

                            REAL TIMER      TASK TIMER
         CLOCK INCREMENT =  0H  0M  0.01S   0H  0M  0.00S 
    TOTAL EXECUTION TIME =  0H  0M  0.36S   0H  0M  0.07S 
 (TASK TIME)/(REAL TIME) = 0.19444
-
-
-
-
-
                                    ***************   ***               ***    **************  
                                    ***************   ****              ***    *************** 
                                    ***************   *****             ***    ****************
                                    ***               ******            ***    ***          ***
                                    ***               *** ***           ***    ***          ***
                                    ***               ***  ***          ***    ***          ***
                                    ***               ***   ***         ***    ***          ***
                                    ***               ***    ***        ***    ***          ***
                                    *************     ***     ***       ***    ***          ***
                                    *************     ***      ***      ***    ***          ***
                                    *************     ***       ***     ***    ***          ***
                                    ***               ***        ***    ***    ***          ***
                                    ***               ***         ***   ***    ***          ***
                                    ***               ***          ***  ***    ***          ***
                                    ***               ***           *** ***    ***          ***
                                    ***               ***            ******    ***          ***
                                    ***************   ***             *****    ****************
                                    ***************   ***              ****    *************** 
                                    ***************   ***               ***    **************  
