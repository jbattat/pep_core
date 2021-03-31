      subroutine ASBSET(ibelt,icm,ica,ici,ico)
 
      implicit none

c arguments
      integer*4 ibelt,icm,ica,ici,ico

c  ibelt  - selector for asteroid/comet belt
c  icm    - index into prmter of the belt mass
c  ica    - index into prmter of the belt radius
c  ici    - index into prmter of the belt inclination
c  ico    - index into prmter of the belt node
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'astroi.inc'
      include 'funcon.inc'
      include 'param.inc'
      include 'petuna.inc'

c local
      real*10 asinc,asnod,cosinc,cosnod,sininc,sinnod

c initial setup for asteroid/comet belt perturbation
c for now, assume belt is sun-centered unless lcentr=0
c
      Icasbm(ibelt) = icm
      Icasba(ibelt) = ica
      Icasbi(ibelt) = ici
      Icasbo(ibelt) = ico
      Kpasb(ibelt)  = Kp(icm+30)
      Ncnasb(ibelt) = 0
      if(lcentr.le.0 .or. Ncentr.lt.0) Ncnasb(ibelt) = Ncentr
      Asbpm(ibelt) = prmter(icm)
      Asbpa(ibelt) = prmter(ica)
      asnod  = Convd*prmter(ico)
      sinnod = SIN(asnod)
      cosnod = COS(asnod)
      asinc  = Convd*prmter(ici)
      cosinc = COS(asinc)
      sininc = SIN(asinc)
      Astnrm(1,ibelt) = sinnod*sininc
      Astnrm(2,ibelt) = -cosnod*sininc
      Astnrm(3,ibelt) = cosinc
      Astcom(1,ibelt) = sinnod*cosinc
      Astcom(2,ibelt) = -cosnod*cosinc
      Astcom(3,ibelt) = -sininc
      Astvc3(1,ibelt) = -Astnrm(2,ibelt)
      Astvc3(2,ibelt) = Astnrm(1,ibelt)
      Astvc3(3,ibelt) = 0._10
      Astapr(ibelt)   = (ici.le.50 .and. Kp(ici+30).gt.0) .or.
     .                  (ico.le.50 .and. Kp(ico+30).gt.0)
      return
      end
