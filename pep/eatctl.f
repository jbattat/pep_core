      subroutine EATCTL(icall,active,kick,kobj,eatcor)
 
      implicit none
c
c     r.king      march 1978   subroutine eatctl
c     calculates the effect of the earth's neutral atmosphere
c     on delay and delay rate for use with all radio and laser-
c     ranging observables
c
c
c icall= 1  delay
c = 2  delay rate
c
c arguments
      logical*4 active
      integer*4 icall,kick,kobj
c           active= .true.   meteorological data available for each
c                            observation
c                =  .false.  static atmosphere assumed
c
c           kick= 1 propco called by radar link
c               = 4 propco called by fermtr link
c
c           kobj= 1 observed body is nplnt0
c               = 2 observed body is nplnt2 or nspot2
c
      real*4    eatcor(2)
c atmospheric correction (delay or rate) for sites 1 and 2

c array dimensions
      include 'globdefs.inc'

c common
      include 'coord.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'number.inc'
      include 'param.inc'
      include 'prpgat.inc'
      include 'radcrd.inc'
c
c external functions
      real*4    EATMDL,EATMDP,EATMMP
c
c quantities internal to this routine
      integer   i,numsit,ict23
      real*4    wetz(2),dryz(2),atmprm(2),atmdl(2)
c
c
c-----------------------------------------------------------------------
c
      eatcor(2) = 0.
      atmdl(2)  = 0.
      numsit    = 2
      if(Nsite2.eq.0) numsit = 1
c
c determine wet and dry components of delay at zenith
      call EATZDL(active,wetz,dryz)
c
c
c calculate path delay using mapping function
      if(icall.eq.2) then
c
c calculate delay rate using mapping function
         do i = 1,numsit
            atmdl(i) = EATMDP(wetz(i),dryz(i),Za(i,kobj),0.,0.,
     .                 Zar(i,kobj))
         end do
 
c rates of change of wetz and dryz assumed=0.
         Raddum(4) = atmdl(1)
         Raddum(5) = atmdl(2)
         if(kick.eq.4) Raddum(5) = -Raddum(5)
         if(Nsite1.eq.Nsite2) Raddum(4) = atmdl(1) + atmdl(2)
c rate partials stored in same locations as delay partials
c this is ok because there is a suicid elsewhere is nice=0
         if(kobj.eq.2) call SUICID(
     .' PARTIALS OF DELAY RATE WRT ATM. PARM. NOT CODED FOR TWO OBJECT O
     .BSERVABLES, STOP IN EATCTL ',23)
      else
         ict23 = Ict(23)
         if (MOD(ict23,2).eq.1) then
c     use the new mapping function
            do i = 1,numsit
               atmdl(i) = EATMMP(i,wetz(i),dryz(i),Za(i,kobj))
            end do
         else
c     otherwise use the old MF
            do i = 1,numsit
               atmdl(i) = EATMDL(wetz(i),dryz(i),Za(i,kobj))
            end do
         endif
c
c     store delay quantities for partial derivatives
c
c     this code to be changed when old routines removed from pep
c     and new atmosphere parameters added.  stored quanties should
c     represent 1-way path delays and partl1 should compute observable
c     partials.
c
         if(kobj.eq.2) then
 
            Radum2(4) = atmdl(1)
            Radum2(5) = atmdl(2)
            if(kick.eq.4) Radum2(5) = -Radum2(5)
            if(Nsite1.eq.Nsite2) Radum2(4) = atmdl(1) + atmdl(2)
         else
 
            Raddum(4) = atmdl(1)
            Raddum(5) = atmdl(2)
            if(kick.eq.4) Raddum(5) = -Raddum(5)
c for single station radio or radar observable, adjustable
c parameter scales the round trip time delay
            if(Nsite1.eq.Nsite2) Raddum(4) = atmdl(1) + atmdl(2)
         endif
      endif
c returned quantity is delay rate - freq needed for doppler shift
c ----------------------------------------------------------------------
c
c modify delay or rate corrections by input scale factor(s)
      if(Ncph.le.0) then
         do i = 1,2
            atmprm(i) = 1.0
         end do
      else
         do i = 1,numsit
            atmprm(i) = Aphs(i)
         end do
         if(Nsite1.eq.Nsite2) atmprm(2) = Aphs(1)
      endif
      do i = 1,numsit
         atmprm(i) = atmprm(i)*prmter(62)
 
c default value of prmter(62) is 1
         eatcor(i) = atmprm(i)*atmdl(i)
      end do
c
c
      return
      end
