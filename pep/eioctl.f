      subroutine EIOCTL(icall,active,kick,kobj,f,eiocor)
 
      implicit none
c
c     r.king      march 1978      subroutine eioctl
c     calculates the theoretical zenith delay of the earth's
c     ionosphere using a static model or input values of peak
c     electron density
c
c icall= 1  delay (group delay assumed)
c 2  delay rate
c
      logical*4 active
      integer*4 icall,kick,kobj
c          active= .true.  values of peak electron density input for
c                          each observation
c                  .false. static (rectified cosine) model assumed
c
c           kick= 1 propco called by radar link
c               = 4 propco called by fermtr link
c
c           kobj= 1 observed body is nplnt0
c               = 2 observed body is nplnt2 or nspot2
c
      real*4    eiocor(2)
c ionospheric corrections (delay or rate) for sites 1 and 2
c
      real*10 f(2)
c signal frequencies for sites 1 and 2
c

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'coord.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'number.inc'
      include 'obscrd.inc'
      include 'param.inc'
      include 'prpgat.inc'
      include 'radcrd.inc'
      include 'yvect.inc'

c external functions
      real*4    EIONDL,EIONDP
c
c quantities internal to this routine
      integer   i,nag,numsit
      real*4    em(2),emdot(2),frt4(2),ionprm(2),iondl(2),hm(2),tmdly
c
c
      eiocor(2) = 0.
      iondl(2)  = 0.
      numsit    = 2
      if(Nsite2.eq.0) numsit = 1
c
c for static ionosphere calculate peak electron densities for
c both sites using a rectified cosine with a nighttime value
c of 1.e5 elec/cm**3 and a daytime peak of 1.e6 elec/cm**3
      if(active) then
c
c for active model read in peak electron density
         call SUICID(
     .'ACTIVE IONOSPHERE MODEL NO LONGER IMPLEMENTED, STOP IN EIOCTL LAB
     .EL=100 ', 18)
      else
         frt4(1) = 0.
         if(ABS(Utrec).ge.0.29) frt4(1) = Utrec/86400.
c note:  time difference between stations ignored in static
c model, tmdly=0.
         tmdly   = 0.
         frt4(2) = tmdly/86400. + frt4(1)
         call DIURN(Imonth,frt4,em,emdot)
 
c height of peak density assumed = 300 km.
         hm(1) = 350.
         hm(2) = 350.
      endif
c        save(31-34), formerly used for peak elec. density and its
c        rate of change, have been preempted for other uses.
c
c
c           optional printout of em and emdot
c     if ict(24) <= -2,then print em and emdot for both sites
      if(Ict(24).le.-2) then
         if(Line.gt.57) call OBSPAG
         do i = 1, numsit
            Line = Line + 1
            write(Iout,50) i,em(i),emdot(i)
         end do
   50    format(1x,' SITE= ', i6, ' EM= ', e14.7, ' EMDOT= ', e14.7)
      endif
c
c calculate ionospheric group delay along observing path
c (sign changed in propco if phase delay desired)
      nag = 0
      if(icall.eq.2) then
c
c
c
c calculate rate of change of ionospheric delay
         do i = 1, numsit
            iondl(i) = EIONDP(Za(i,kobj),em(i),hm(i),f(i),
     .                 Zar(i,kobj),emdot(i))
         end do
      else
         do i = 1, numsit
            iondl(i) = EIONDL(Za(i,kobj),em(i),hm(i),f(i),nag)
         end do
c
c           store quantities for partial derivatives
c
c     this code to be changed when old routines removed from pep
c     and new ionosphere parameters added.  stored quanties should
c     represent 1-way path delays and partl1 should compute observable
c     partials.
c
         if(kobj.eq.2) then
            Radum2(9)  = iondl(1)
            Radum2(10) = iondl(2)
            if(kick.eq.4) Radum2(10) = -Radum2(10)
            if(Nsite1.eq.Nsite2) Radum2(9) = iondl(1) + iondl(2)
         else
            Raddum(9)  = iondl(1)
            Raddum(10) = iondl(2)
            if(kick.eq.4) Raddum(10) = -Raddum(10)
c for single-station radio or radar observable, adjustable parameter
c scales round-trip time delay
            if(Nsite1.eq.Nsite2) Raddum(9) = iondl(1) + iondl(2)
         endif
      endif
c
c modify delay or rate corrections by input scale factor(s)
      do i = 1, numsit
         ionprm(i) = prmter(63)
 
c default value of prmter(63) is 1
         if(Ncph.gt.2) ionprm(i) = Aphs(i + 2)*ionprm(i)
         eiocor(i) = ionprm(i)*iondl(i)
      end do
c rate partials w.r.t. ionosphere parameters not programmed
c
      return
      end
