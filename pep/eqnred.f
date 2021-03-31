      subroutine EQNRED(in0,nstop,init)
 
      implicit none
c
c K.M.Becker       August 1967       subroutine EQNRED
c Equinox-equator-latitude corrections for optical observation
c series are initialized and read.
c
c arguments
      integer*4 in0,nstop
      logical*4 init
c
c Modified for *command July 1978  R.B. Goldstein

c array dimensions
      include 'globdefs.inc'

c common
      include 'eqenox.inc'
      integer*2 l(3)
      include 'inodta.inc'
c
c     for j=1,...,numeqn we have
c           Eqnsit(j)= First four characters of eight character optical
c                      observing site name.
c                      If '@REF', then the correction is a pure rotation
c                      vector expressed in arcsec, instead of eq-eq-lat.
c           Eqnser(j)= Four characters giving observation series name.
c           Denox(j)= Correction in seconds of time to equinox.
c                     Includes correction to longitude (mu=phi+psi).
c           Dequat(j)= Correction in seconds of arc to equator (theta).
c           Dlat(j)=   Correction in seconds of arc to site latitude.
c           Leqn(i,j)=0  no adjustment to Denox,Dequat,Dlat  for i=1,2,3
c           Leqn(i,j)=1     adjustment to Denox,Dequat,Dlat  for i=1,2,3
c
c                vii. equinox-equator-latitude corrections  (EQNRED)
c  for each optical site and observation series
c  columns
c    1- 4  site name  first four characters                  (1a4)
c    6- 9  observation series name  (four characters)        (1a4)
c      17  l1=1,0  equinox  correction is adjusted or not    (i1)
c      19  l2=1,0  equator  correction is adjusted or not    (i1)
c      21  l3=1,0  latitude correction is adjusted or not    (i1)
c   37-48  equinox  correction                               (e12.5)
c   49-60  equator  correction                               (e12.5)
c   61-72  latitude correction                               (e12.5)
c  Limitation:  no more than u_mxeqx sites allowed
c
c local
      real      deq, dl, dnx
      integer   i, j, ns
      character*4 usn6/'6USN'/
      character*4 srusn6/'6956'/
      character*4 blank /'    '/,est,eser
c
c initialize site data
      do j = 1, u_mxeqx
         Eqnsit(j) = blank
         Eqnser(j) = blank
         Denox(j)  = 0.0
         Dequat(j) = 0.0
         Dlat(j)   = 0.0
         do i = 1, 3
            Leqn(i,j) = 0
         end do
      end do
c
c set up standard site names and data
c there are 15 standards but only 1 is used to provide
c sample output
      Numeqn    = 1
      Eqnsit(1) = usn6
      Eqnser(1) = srusn6
      if(.not. (init)) then
c
c
c spool equinox-equator-declination bias or interferometer
c clock bias cards from in to in0 with a-format printout
         call PEPTIC(In,Iout,in0,10,
     .               'EQ-EQ-DECL OR INTERFEROMETER BIAS CARDS ', nstop,
     .               1)
         do while( .true. )
c
c read equator - equinox data
            read(in0,20) est,eser,(l(i),i = 1,3),dnx,deq,dl
   20       format(a4,1x,a4,6x,3(1x,i1),15x,3E12.5)
            if(est.eq.blank) goto 100
c
c search to see if this is standard site and standard series
            if(Numeqn.gt.0) then
               do i = 1, Numeqn
                  if(est.eq.Eqnsit(i)) then
                     if(eser.eq.Eqnser(i)) then
                        ns = i
                        goto 40
                     endif
                  endif
               end do
            endif
            Numeqn = Numeqn + 1
            if(Numeqn.gt.u_mxeqx)
     .           call SUICID('TOO MANY INPUT SITES, STOP IN EQNRED', 9)
            Eqnsit(Numeqn) = est
            Eqnser(Numeqn) = eser
            ns   = Numeqn
   40       do i = 1, 3
               Leqn(i,ns) = l(i)
            end do
            Denox(ns)  = dnx
            Dequat(ns) = deq
            Dlat(ns)   = dl
         end do
      endif
 
  100 rewind in0
      return
      end
