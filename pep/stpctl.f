      subroutine STPCTL(hcx)
 
      implicit none
 
c
c        dec. 1979   rdr/kcl
c        stepsize control subroutine using formula :
c        tau= p/n0 * (r/a)**n
c        kk(51) = 0 no control
c                -1 control in transition stage
c                 1 control in locked-on stage
c argument
      real*10 hcx
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'stint.inc'

      hcx = (Dtcon(2)/Dtcon(3))*(Rsbsng/Dtcon(1))**Dtcon(4)
      return
      end
