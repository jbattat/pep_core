      subroutine ROCHNG(i0,psi0,delta0,alpha0,n,j,ctl)
 
      implicit none
 
      real*10 i0, psi0, delta0, alpha0, n, j
      integer*4 ctl
c
c        r. goldstein   june 1977
c        all angles are in degrees
c
c           n     euler angles that relate an intermediate inertial
c           j     frame to inertial frame of reference epoch
c
c           i0    inclination and precission phase angle at epoch
c           psi0  con1(1) referred to intermediate inertial frame
c
c           delta0 right ascention and declination of pole at
c           alpha0 epoch con1(1) referred to reference inertial frame
c
c        ctl-control integer
c           0: calculate i0,psi0 from delta0,alpha0
c           1: calculate delta0, alpha0 from i0,psi0
c
      real*10 sj, cj, si, ci, sd, cd, saln, caln, aln, spsi, cpsi
 
      include 'funcon.inc'
 
      cj = COS(j*Convd)
      sj = SIN(j*Convd)
      if(ctl.eq.0) then
         sd   = SIN(delta0*Convd)
         cd   = COS(delta0*Convd)
         aln  = alpha0 - n
         saln = SIN(aln*Convd)
         caln = COS(aln*Convd)
         i0   = ACOS(cj*sd - sj*cd*saln)/Convd
         psi0 = ATAN2(-(caln*cd),(cj*cd*saln+sj*sd))/Convd
c
c
      else if(ctl.ne.1) then
 
         write(6,50) ctl
   50    format(' ROCHNGE...CTL.NE.1 OR 2...CTL=', i20)
      else
         si     = SIN(i0*Convd)
         ci     = COS(i0*Convd)
         spsi   = SIN(psi0*Convd)
         cpsi   = COS(psi0*Convd)
         delta0 = ASIN(cj*ci + sj*si*cpsi)/Convd
         alpha0 = n + ATAN2(-(sj*ci-cj*si*cpsi),-(si*spsi))/Convd
      endif
      return
      end
