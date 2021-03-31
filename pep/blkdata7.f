      block data BLKDATA7
 
      implicit none
 
c
c|blkdata7
c m.e.ash   oct 1966    constants for comparison of theory and obs

c array dimensions
      include 'globdefs.inc'
c
      include 'comdateq.inc'
      data Secday, Secdy2/8.64E4_10, 4.32E4_10/,
     .  fact/3._10,5._10,7._10,9._10/
      include 'namobs.inc'
      data Typobs/' RADAR', ' RADIO', ' RADAR', ' MERID', ' ASTMT',
     .     ' ASTGR', ' OCCLT', ' TRNST', ' SPTRN', ' PHASE', ' INTRF',
     .     ' N-CNT', ' 2-DLY', ' 2-DLY', 3*' .. ', ' PULSR', ' DPCNT',
     .     ' AZ-EL'/
      data Typlka/' LKANG'/
      end
