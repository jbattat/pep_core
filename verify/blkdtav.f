C|BLKDTAV
      BLOCK DATA BLKDTAV
      IMPLICIT NONE
 
      include 'chars.inc'
      DATA BLANK/' '/,COMMA/','/,QUOTE/''''/,DASH/'-'/,DOT/'.'/,ONE/'1'/
      DATA VERSN/' 3.8'/
      include 'dflt.inc'
      DATA  JSIZD/ 80 /,   IRFMD/ 0 /,   CHAND/ '1', ' ','+' /,
     1      IPUND/ 7 /,    ISQND/ 0 /,   LCMPD/ -1 /,
     2      JCLD/ 0 /,     KFIND/ 72 /,  NOPRTD/ 1 /
      include 'ibfl.inc'
      DATA IBFL,IBFLL,ISUB1 /24000,2*0,1,24001/
      include 'pgfba.inc'
      DATA NKPG/20/
      END
