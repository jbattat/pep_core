      block data BLKDATA2
 
      implicit none
 
c
c|blkdata2  r.reasenberg   jan 1973   block data for integration
c
c array dimensions
      include 'globdefs.inc'

c commons
      include 'incon.inc'
      data Intm/0/, Ints/0/, Iptr3/3/
      include 'stint.inc'
c        v0    = initial conditions
c        n     = number of differential equations
c        m     = number of controling equations
c        lthrst= 0 starting proceedure in subroutine nint
c                  uses zero initial guess at higher derivatives
c                  of approximating polynomial
c        lthrst= 1 starting proceedure uses current value of these
c                  derivatives (set in subroutine sbout for starting
c                  proceedure at initiation and termination of les-8/9
c                  station keeping thrusts)
      data Lthrst/0/
      end
