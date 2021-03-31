      real*10 function NUTRM(a,b,c,dd,e,coff,ngo)
 
      implicit none

c
c m.e.ash   oct 1967    function nutrm
c evaluate earth nutation periodic terms
c
c arguments
      real*10 a, b, c, dd, e, coff
      integer   ngo

c common 
      include 'argfun.inc'

c local
      real*10 sum, arg8

c external functions
      real*10 ARGMNT

c evaluate argument between -twopi and twopi radians
      sum  = a*L + b*Lp + c*F + dd*D + e*Ascm
      arg8 = ARGMNT(sum)
c
      if(ngo.eq.2) then
         sum = COS(arg8)
      else
         sum = SIN(arg8)
      endif
c
c multiply by coefficient
      NUTRM = coff*sum
      return
      end
