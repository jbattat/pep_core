      subroutine QSET
 
      implicit none
 
c        r. goldstein   july 1977
c        sets up the q matrix, used in transforming the
c        partials wrt psi0,i0 to partials wrt alpha,
c        delta
c        reference: memo by rdr dated june 22, 1977
c
c        sept. 1977
c        memo by rdr dated june 22 1977 was found to be in error. q
c        should be a 2x2 matrix. the error caused the partials wrt
c        r.a. and dec. to be wrong. in the transform of partials wrt
c        i0,psi0 to alpha0,delta0 those terms that included partials
c        wrt phi0 should be eliminated. this routine has been
c        modified to retain q as 3x3, but set the last row and
c        column to zero's (q33=1)
c
      include 'rotcom.inc'
 
      real*10 caln, saln, cd0, sd0, si0, ci0, ti0, sps0, cps0, tps0, 
     .          sj, cj, ppsal, ppsd, pphps, ppsi, pial, pid, q(3, 3)
 
      equivalence(Trig(5), sd0), (Trig(6), cd0), (Trig(14), si0),
     .            (Trig(15), ci0), (Trig(17), sps0), (Trig(18), cps0),
     .            (Trig(11), sj), (Trig(12), cj), (Trig(23), saln),
     .            (Trig(24), caln), (q(1,1), Qmat(1,1))
c
c
      ti0   = si0/ci0
      tps0  = sps0/cps0
      ppsal = (cd0*saln)/(cps0*si0)
      ppsd  = (sd0*caln)/(cps0*si0)
      pphps = -ci0
      ppsi  = -tps0/ti0
      pial  = -sj*sps0
      pid   = -(sj*sd0*saln + cj*cd0)/si0
c
c
      q(3, 1) = 0._10
      q(3, 2) = 0._10
      q(3, 3) = 1._10
      q(1, 1) = ppsal + ppsi*pial
      q(1, 2) = pial
c
c q(1,3)=pphps*(ppsal+ppsi*pial)
      q(1, 3) = 0._10
 
      q(2, 1) = ppsd + ppsi*pid
      q(2, 2) = pid
c
c q(2,3)=pphps*(ppsd+ppsi*pid)
      q(2, 3) = 0._10
c
c
      return
      end
