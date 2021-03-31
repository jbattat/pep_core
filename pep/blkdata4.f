      block data BLKDATA4
 
      implicit none
 
c
c blkdata4  - j.f.chandler - 1984 jul
c fixed special parameter names for input link
c
 
      include 'prmnms.inc'
c
c locals
      data Qprmtr/'PRMTER'/, Prmsmx/29/
 
c prmsmx = size of prms and iprms arrays (must be in order)
      data Prms/'RELFCT  ', 'GMVARY  ', 'SUNHAR  ', 'BETA    ',
     .     'GAMMA   ', 'BETA''  ', 'GAMMA'' ', 'ASCABELT', 'INCABELT',
     .     'DSTABELT', 'MASABELT', 'AULTSC  ', 'LTVARY  ', 'RELDEL  ',
     .     'RELDOP  ', 'PLASMAC ', 'PLASMAV ', 'ATMFACT ', 'IONFACT ',
     .     'CTVARY  ', 'SFATTOCT', 'PHATTOCT', 'ECINC   ', 'SEQINC  ',
     .     'SEQASC  ', 'SUNRAD  ', 'MDSTAU  ', 'MDSTSC  ', 'LTVEL   '/
      data Iprms/31, 32, 33, 41, 42, 43, 44, 47, 48, 49, 50, 51, 52, 53,
     .     54, 60, 61, 62, 63, 72, 81, 82, 91, 92, 93, 94, 98, 99, 100/
      end
