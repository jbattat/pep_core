      subroutine PHSRED(in0,nstop,init)
 
      implicit none

c
c k.m.becker      august 1967       subroutine phsred
c phase corrections for planetary optical observation series are
c initialized and read
c
c parameters
      integer*4 in0,nstop
      logical*4 init
c
c        modified for *command july 1978  r.b. goldstein
c
c        setup and read quantities into temporary storage , then
c      rearrange into labelled common, such that corrections belonging
c      to a given planet are grouped together.  the order of the
c      planets is as follows
c
c              10=moon
c               0=sun
c        nplnt(i)=planet(aplnt(1-2, i)), i=1,u_mxpl
c
c       with planet numbers not coinciding with 10,0,or
c      (nplnt(i),i=1,u_mxpl) following in the order in which they were
c      setup or read

c array dimensions
      include 'globdefs.inc'

c commons
      include 'inodta.inc'
      include 'namtim.inc'
      include 'phase.inc'
 
c shared work space in input link
      common/WRKCOM/ Phsitt,Phsert,Aphast,Ncphst,Lphst,Nplpht
      character*4 Phsitt(u_mxphs),Phsert(u_mxphs)
      real*4 Aphast(9,u_mxphs)
      integer*2 Ncphst(u_mxphs),Lphst(9,u_mxphs),Nplpht(u_mxphs)
c
c     for j=1,...,numphs we have
c           phsit(j)= first four characters of eight character optical
c                     observing site name
c           phser(j)= four characters giving observation series name
c           nplphs(j)= planet number of observed body
c                      0= sun              1= mercury          2= venus
c                      4= mars             5= jupiter          6= saturn
c                      7= uranus           8= neptune          9= pluto
c                     10= moon
c           aphase(i,j)= dimensionless phase coefficients
c           ncphs(j)= number of coefficients in phase model
c           lphs(i,j) =0 aphase(i,j) not adjusted
c           lphs(i,j) =1 aphase(i,j) adjusted in least squares analysis
c
c               viii. planetary phase corrections  (phsred)
c  for each planet and observation series
c  card 1
c  columns
c    1- 4  observing site name  first four characters        (1a4)
c    6- 9  observation series name  first four characters    (1a4)
c   10-12  planet number                                     (i3)
c   13-15  number of coefficients in phase correction mode   (i3)
c  limitation is 9 (if less than or equal 3, only one card read)
c   16-17  l1=1,0  first   coefficient is adjusted or not    (i2)
c   18-19  l2=1,0  second  coefficient is adjusted or not    (i2)
c   20-21  l3=1,0  third   coefficient is adjusted or not    (i2)
c   22-23  l4=1,0  fourth  coefficient is adjusted or not    (i2)
c   24-25  l5=1,0  fifth   coefficient is adjusted or not    (i2)
c   26-27  l6=1,0  sixth   coefficient is adjusted or not    (i2)
c   28-29  l7=1,0  seventh coefficient is adjusted or not    (i2)
c   30-31  l8=1,0  eigth   coefficient is adjusted or not    (i2)
c   32-33  l9=1,0  nineth  coefficient is adjusted or not    (i2)
c   37-48  first   phase coefficient                         (e12.5)
c   49-60  second  phase coefficient                         (e12.5)
c   61-72  third   phase coefficient                         (e12.5)
c  card 2
c    1-12  fourth  phase coefficient                         (e12.5)
c   13-24  fifth   phase coefficient                         (e12.5)
c   25-36  sixth   phase coefficient                         (e12.5)
c   37-48  seventh phase coefficient                         (e12.5)
c   49-60  eigth   phase coefficient                         (e12.5)
c   61-72  nineth  phase coefficient                         (e12.5)
c  limitation  no more than u_mxphs sites allowed
c
c
      real*4 a(9)
      integer*2 npl,nc,lp(9)
      character*4 usn6/'6USN'/
      character*4 srusn6/'6956'/
      character*4 blank/'    '/,pser,pst
      integer*4 nnphs,nplcop,i,j,l,n,ns
c
c initialize site data
      do j = 1,u_mxphs
         Phsitt(j) = blank
         Phsert(j) = blank
         Ncphst(j) = 0
         Nplpht(j) = 0
         Phsit(j)  = blank
         Phser(j)  = blank
         Ncphs(j)  = 0
         Nplphs(j) = 0
         do i = 1,9
            Aphast(i,j) = 0.0
            Lphst(i,j)  = 0.0
            Aphase(i,j) = 0.0
            Lphs(i,j)   = 0
         enddo
      enddo
c
c set up standard site names and series
      l = 0
      Numphs    = 1
      Phsitt(1) = usn6
      Phsert(1) = srusn6
      Nplpht(1) = 1
      if(.not. (init)) then
c
c
c spool phase corr.cards from in to in0 with a-format printout
         call PEPTIC(In,Iout,in0,6,'PHASE CORRECTION CARDS  ',
     .               nstop,1)
         do while( .true. )
c
c read phase data
            read(in0,20) pst,pser,npl,nc,lp,(a(i),i=1,3)
   20       format(a4,1x,a4,2I3,9I2,3x,3E12.5)
            if(pst.eq.blank) goto 100
            if(nc.gt.3) then
               read(in0,30) (a(i),i=4,nc)
   30          format(6E12.5)
            endif
c
c search to see if standard site and standard series
            if(Numphs.gt.0) then
               do i = 1,Numphs
                  if(pst.eq.Phsitt(i) .and. pser.eq.Phsert(i) .and.
     .             npl.eq.Nplpht(i)) then
                     ns = i
                     goto 40
                  endif
               enddo
            endif
            Numphs = Numphs + 1
            if(Numphs.gt.u_mxphs)
     .           call SUICID('TOO MANY INPUT SITES, STOP IN PHSRED',9)
            Phsitt(Numphs) = pst
            Phsert(Numphs) = pser
            Nplpht(Numphs) = npl
            ns   = Numphs
   40       do i = 1,nc
               Aphast(i,ns) = a(i)
               Lphst(i,ns)  = lp(i)
            enddo
            Ncphst(ns) = nc
            Nplpht(ns) = npl
         enddo
      endif
c
c rearrange data to labelled common so that corrections to a given
c planet are grouped together
  100 if(Numphs.gt.0) then
         nnphs = Numphs
 
c find moon phase corrections
         call REDCPY(Ncphst,Nplpht,Phsert,Phsitt,Aphast,Lphst,
     .              9,1,nnphs,l,10,
     .              Ncphs,Nplphs,Phser,Phsit,Aphase,Lphs)
c
c find sun phase corrections
         call REDCPY(Ncphst,Nplpht,Phsert,Phsitt,Aphast,Lphst,
     .              9,1,nnphs,l,0,
     .              Ncphs,Nplphs,Phser,Phsit,Aphase,Lphs)
c
c find planet phase corrections
         do n = 1,Numpln
            nplcop = Nplnt(n)
            if(nplcop.gt.0)
     .       call REDCPY(Ncphst,Nplpht,Phsert,Phsitt,Aphast,Lphst,
     .                  9,1,nnphs,l,nplcop,
     .                  Ncphs,Nplphs,Phser,Phsit,Aphase,Lphs)
         enddo
 
         do j = 1,30
            call REDCPY(Ncphst,Nplpht,Phsert,Phsitt,Aphast,Lphst,
     .                9,1,nnphs,l,j,
     .                Ncphs,Nplphs,Phser,Phsit,Aphase,Lphs)
         enddo
c
c next get phase corrections for stars. (spots on planet no=-4)
         call REDCPY(Ncphst,Nplpht,Phsert,Phsitt,Aphast,Lphst,
     .              9,1,nnphs,l,-4,
     .              Ncphs,Nplphs,Phser,Phsit,Aphase,Lphs)
         if(l.ne.Numphs) call SUICID(
     .       'ERROR IN NUMBER OF PHASE CORRECTIONS, STOP IN PHSRED',
     .       13)
      endif
 
      rewind in0
      return
      end
