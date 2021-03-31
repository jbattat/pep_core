      subroutine ADJAST(old,fact)
 
      implicit none
c
c m.ash   sept 1967    subroutine adjast
c calculate adjustment to parameter and standard deviation
c increment counters
c
c arguments
      real*10 fact, old
c     old  =old value of parameter
c     fact =scale factor for relating units of solution to normal
c           equations and units of parameter
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'adjstf.inc'
      include 'estvec.inc'
      include 'fcntrl.inc'
      include 'rtsidesl.inc'
c
c     adj  =adjustment to parameter in units of parameter
c     sig  =standard deviation of parameter estimate in units of paramtr
c     nwv  =new value of parameter in units of parameter
c     fract=adj/sig (large number to printout **** if sig=0.0_10)
c     ntype=1 if abs(adj).ge.eps(9)*sig
c     ntype=2 if abs(adj).ge.eps(10)*sig
c     ntype=3 if neither of the above
c     nkind(ntype)=running total of three types of adjustments
c     n    =index for solution vector of normal equations
c     line =line counter
c     il=1 blank line not yet written as separator
c     il=2 blank line written as separator
c     fradj= fraction (between zero and one) by which adjustment
c            is to be multiplied (input as eps(14))
c     statf= statistics for fract
c     nf   = number of points included in fract statistics
c
      N     = N + 1
      Sig   = Sigma(N)*ABS(fact)
      Adj   = Solut(N)*fact
      Adj   = Adj*Fradj
      Nwv   = old + Adj
      Fract = 1.0E19_10
      if(Sig.ne.0.0_10) then
         Fract = Adj/Sig
         Nf    = Nf + 1
 
c fill xbar with new estimate (unscaled) to be used in formu
         Xbar(N)  = Nwv/fact
         Statf(1) = Statf(1) + Fract
         Statf(2) = Statf(2) + ABS(Fract)
         Statf(3) = Statf(3) + Fract**2
         call TALLY
      endif
 
      Ntype = 1
      if(ABS(Adj).lt.Eps(9)*Sig) then
         Ntype = 2
         if(ABS(Adj).lt.Eps(10)*Sig) Ntype = 3
      endif
      Nkind(Ntype) = Nkind(Ntype) + 1
      call PAGCHK(59,1,1)
 
      Il = 1
      return
      end
