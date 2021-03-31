      program PEPMAIN
 
      implicit none

 
c|main      m.e.ash   oct 1969   main program (always in core storage)
c     Main program for PEP (Planetary Ephemeris Program)
c     Designed to run on an IBM 360 model 67 computer using overlay, it
c     now runs on an IBM 370 model 168 or any plug-compatible.
c
c     MAIN first calls subroutine TIMDAT to set the interval timers and
c     determine the calendar date.  Then it calls subroutine INPUT to
c     set up and read the input control and data constants.
c
c     Then it enters the least squares analysis loop, incrementing the
c     variable ITERAT each time through the loop.  A passage through
c     the loop consists of
c           (1) numerical integrations of equations of motion and
c               equations for partials of motion,
c           (2) computation of observed minus theoretical value of
c               observations and partial derivatives of theoretical
c               value of observations with respect to parameters to
c               be adjusted, and
c           (3) formation and solution of normal equations and
c               adjustment of parameters.
c
c     The perturbing planet data sets can be updated at the completion
c     of each passage through the loop, but this feature is unused for
c     now because of the n-body integration feature.
c
c     The time for each passage through the loop is printed out
c     and tapes can be saved before starting the next loop if desired.
c
c     Upon exiting from the least squares analysis loop, one can:
c           (1) compute theoretical value of dummy observations with
c               newly adjusted parameters (skipped if no adjustments),
c           (2) predict new observed minus theory from old observed
c               minus theory, partials of observations and adjustments,
c           (3) punch adjusted parameters in format required by input
c               link so program can be restarted where this iteration
c               stopped, and
c           (4) calculate and save uncorrelated portion of observable
c               partial derivatives (mutually exclusive with predict
c               feature).
c
c     The block data subroutine for MAIN has comment cards listing
c     (1) overlay structure of PEP;
c     (2) alphabetic list of all subroutines with for each subroutine
c         its entry points, the overlay link in which it appears, the
c         subroutines which call it and the subroutines it calls;
c     (3) alphabetic list of all labelled commons with for each labelled
c         common its length (unreliable), the overlay link in which it
c         appears and the subroutines which reference it.
c
c     INPUT and its subroutines have comment cards explaining all the
c     input control and data constants.
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c array dimensions
      include 'globdefs.inc'

c        commons
      include 'bdctrl.inc'
      include 'crdbuf.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'namtim.inc'
      include 'timstf.inc'
c
c     comment cards below explain control constants used in main program
c     for performing some or all of the options noted above.
c     (duplicates of some of the comment cards in subroutine input)
c
c
c ITERAT  =iteration number (initialized and incremented in main prog.)
c
c ICT(1) =-1 no comparison of theory and observation, no least squares
c            analysis
c ICT(1) = 0 comparison of theory and observation, no parameters
c            adjusted
c ICT(1) =positive, maximum number of least squares iterations involving
c            integration,comparison of theory and observation,solution
c            of normal equations,adjustment of parameters
c
c ICT(5) =-1 no initialization or forming of normal equations from
c            subset of saved normal equations
c ICT(5) = 0 initialize normal equations from subset of saved normal
c            equations, process additional observations and solve
c            normal equations for ict(1).gt.0 and iterat=1,2,3,...
c ICT(5) = 1 for iterat=1, form normal equations from subset of saved
c            normal equations and solve them if ict(1).gt.0 (for
c            iterat.gt.1, behaves like ict(5)=-1 if ict(1).gt.1)
c ICT(5) = 2 for iterat=1, restore statistics, scale factors, solution
c            and inverse of normal equations from saved quantities on
c            data set jmat if ict(1).gt.0 (for
c            iterat.gt.1, behaves like ict(5)=-1 if ict(1).gt.1)
c ICT(5) = 3 same as ict(5)=2, iterat=1 except go to compar link
c            first to generate an observation library tape
c ICT(5) = 4 restore solution from dataset jmat plus zbar and
c            fbar-adjoint from partially prereduced normal equations
c            on imat0.  enables computation of adjusts but not sigmas
c            for "reduced away" parameters; uncertainty of predict
c            cannot be calculated.
c
c ICT(6) =   number of least squares analysis iterations during which
c            perturbing planet data set is to be updated with just
c            completed integrations before starting new integrations
c
c ICT(77)= 0 no punch output of input data at end of program run
c ICT(77)= 1 punch output of input data at end of program run so that
c            job can restarted at termination point with last adjusted
c            values of parameters if ict(1).gt.0
c ICT(77).gt.1 write punch output to data set no. ict(77)
c            after every iteration plus do same as ict(77)=1
c
c ICT(78)= 0 no tapes are saved at end of least squares iteration
c ICT(78)=positive   at end of least squares iterations 1 to ict(78),
c            inclusive, there is a pause with typewriter message
c            'remove and mount tapes as specified' for operator to
c            follow special instructions (except for last iteration)
c
c ICT(79)=-2 dummy observations not calculated each iteration nor
c            after convergence
c ICT(79)=-1 dummy observations not calculated each iteration, but are
c            after convergence
c ICT(79)= 0 dummy observations calculated each iteration, but not
c            after convergence
c ICT(79)= 1 dummy observations calculated each iteration and after
c            convergence
c no calculations after convergence if nothing adjusted (ict(1).le.0)
c
c JCT(79)= 0 no reintegrations before calculating dummy observations
c            after orbit fit convergence
c JCT(79)= 1 reintegration of motion without partials after orbit fit
c            convergence to be used in calculation of dummy observations
c            if iabs(ict(79))=1
c
c ICT(80)=-1 input observation library tapes not used
c ICT(80)= 0 input observation library tapes are used, program eats its
c            tail from iteration to iteration
c ICT(80)= 1 input observation library tapes only are used to form the
c            normal equations for first iteration with no computation
c            of observed minus theory.
c            subsequent iterations same as ict(80)=0
c
c
c NBODY  = 0 no n-body integration or n-body data set
c NBODY  = n.gt.0, n-body integration and/or n-body data set with n
c                  bodies plus moon if not one of n=9 bodies
c
c NPLBDY(j),j=1,nbody, are the planet numbers of the bodies in the
c              n-body integration in the order which they are integrated
c
c JDBDY0=initial epoch of n-body integration is midnight beginning of
c        day of day with julian day number jdbdy0 (if negative there is
c        checkpoint restart at the epoch -jdbdy0 unless
c        jdbdy0.le.-3000000) (if zero or -3000000 there is no
c        numerical integration, comparison of theory and observation
c        assumes motion is already on tape, but if initial conditions
c        are adjusted, jdbdy0 is set to positive value on second record
c        of tape so there will be reintegration on next least squares
c        iteration, whereas if no initial conditions are adjusted,
c        jdbdy0 is set to zero in the compar link so there will be no
c        reintegration on subsequent least squares iterations)
c        (jdbdy0.le.-3000000 is signal to compar link on first
c        iteration to override end time on tape by input jdbdy2)
c
c
c NUMPLN = number of input planets,asteroids,satellites and
c            artificial space probes (exclusive of earth-moon
c            barycenter,moon,earth rotation,moon rotation)
c
c NPLNT(j)= planet number, j=1,numpln
c           1 mercury
c           2 venus
c           4 mars
c           5 jupiter
c           6 saturn
c           7 uranus
c           8 neptune
c           9 pluto
c           11-30  other natural planets,asteroids,satellites
c           31,... artificial space probes
c           -1,-2,-4,...,-9,-11,... rotation for planet iabs(nplnt)
c           3,10,-3,-10 do not exist in nplnt(j),j=1,numpln, imagined
c             to be associated with earth-moon barycenter, moon,
c             earth rotation, moon rotation
c           0 does not exist in nplnt(j),j=1,numpln, imagined to be
c             associated with sun
c
c JDEM0   = initial time for earth-moon barycenter integration
c JDMN0   = initial time for moon integration
c JDER0   = initial time for earth rotation integration
c JDMR0   = initial time for moon rotation integration
c JDPL0(j)= initial time for planet nplnt(j) integration, j=1,numpln
c same logic applies for these initial times as for jdbdy0
c
c        local
      real*4 zero4/0./
      character*4 czero
      equivalence (zero4,czero)
      integer   i, ict79, j, keep, mpage, n5
      character*80 tmes/'    LEAST SQUARES ITERATION*** DURING THE LAST*
     .**** PAGES STARTING ON PAGE***** '/
      integer*4 irltsk(2)
      character*40 typmes/'COMPLETION OF LEAST SQUARES ITERATION  1'/

      integer*2 repeat
c     repeat=0  job ends at end of run with given input data (default)
c     repeat=1  job repeats at end of run with new input data
c               if the pep run has had a normal end
c     repeat=2  job repeats even on abnormal pep ending
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c*start=100
      call PLOG0
      Ieof = 0
 
c set interval timer, determine today's date and time
  100 Date(1) = czero
      call TIMDAT(Date)
c
c read input data and control constants
      Npage = 1
      call PAGSET('********************************', 8)
      call INPUT(repeat)
      Iterat = 0
c
c-----------this is start of least squares loop-------------
c
c*start=200
c           increment iteration counter, save present status of timers
  200 Iterat = Iterat + 1
 
c save times and page
      irltsk(1) = Ireal0
      irltsk(2) = Itotsk
      mpage     = Npage
c
c integration of motion of n=nbody bodies, earth-moon
c barycenter, planets and asteroids, and satellites and
c probes with results written on magnetic tape or disk
      if(Nbody.gt.0) then
         if(Jdbdy0.ne.0 .and. Jdbdy0.gt.-3000000) then
            call PLANET(0)
            goto 500
         else
            do i = 1, Nbody
               if(Nplbdy(i).eq.3) goto 300
            end do
         endif
      endif
      if(Jdem0.ne.0 .and. Jdem0.gt.-3000000) then
         call PLANET(0)
         goto 500
      else if(Numpln.le.0) then
         goto 500
      endif
  300 do j = 1, Numpln
         if(Nplnt(j).gt.0) then
            if(Nbody.gt.0) then
               do i = 1, Nbody
                  if(Nplbdy(i).eq.Nplnt(j)) goto 400
               end do
            endif
            if(Jdpl0(j).ne.0 .and. Jdpl0(j).gt.-3000000) then
               call PLANET(0)
               goto 500
            endif
         endif
  400 end do
c
c*start=400
c integration of motion of moon with results written on
c magnetic tape or disk
  500 if(Nbody.gt.0) then
         do i = 1, Nbody
            if(Nplbdy(i).eq.10) goto 600
         end do
      endif
      if((Jdmn0.ne.0 .and. Jdmn0.gt.-3000000) .or.
     .   (Jdmr0.ne.0 .and. Jdmr0.gt.-3000000)) call MOON
c
c*start=600
c           integration of motion of earth, moon and planets about their
c           centers of mass with results written on magnetic tape
c           or disk
  600 if(Jder0.ne.0 .and. Jder0.gt.-3000000) then
         call ROTAT
      else if(Numpln.gt.0) then
         do j = 1, Numpln
            if(Nplnt(j).lt.0) then
               if(Jdpl0(j).ne.0 .and. Jdpl0(j).gt.-3000000) then
                  call ROTAT
                  goto 700
               endif
            endif
         end do
      endif
c
c*start=1000
c           comparison of theory and observation for observations of
c           sun,moon,planets,asteroids,satellites,stars and artificial
c           space probes by interpolating from motions on tape or disk.
c           calculate partials of observations with respect to
c           parameters to be adjusted and write observed minus theories
c           and partials on magnetic tape or disk
  700 if(Ict(1).lt.0) goto 1100
      if(Ict(1).ne.0) then
         if(Ict(5).ne.3) then
            if(Ict(5).le.0) then
               if(Ict(80).le.0) goto 800
            endif
            if(Iterat.le.1) goto 900
         endif
      endif
  800 call COMPAR(0)
      if(Ict(1).le.0) goto 1200
c
c*start=1200
c           form normal equations from subset of saved normal equations
c           and/or from observed minus theory and partials of
c           observations on tapes or disks.
c           solve normal equations, adjust parameters
  900 call ANALIZ
c
c update perturbing planet data set
      if(Ict(6).ge.Iterat) call PRTERB
c
c decide if tapes are to be saved at end of least squares
c analysis iteration
      write(Iout, 1000) Iterat, Heding, Date, Npage
 1000 format('1COMPLETION OF LEAST SQUARES ITERATION', i2, 2x, 18A4,
     .       1x, 2A4, ' PAGE', i5)
      Npage = Npage + 1
      Line  = 1
      if(Iterat.le.Ict(78)) then
         if(Ict(1).gt.Iterat) then
            pause 'REMOVE AND MOUNT TAPES AS SPECIFIED'
            call TIMRIT('PAUSE TO REMOVE AND MOUNT TAPES ', 8)
         endif
      endif
c
c type-to-operator completion of least squares iteration
      call EBCDIX(Iterat, typmes, 38, 3)
      call OPRMSG(typmes, 10)
c
c printout time for least squares analysis iteration
      call EBCDIX(Iterat, tmes, 28, 3)
      call EBCDIX(Npage - mpage, tmes, 47, 5)
      call EBCDIX(mpage, tmes, 75, 5)
      call TIMRTC(tmes, 20, irltsk)
c
c*start=1400
c           determine if iteration is to continue (if convergence is
c           obtained before iterat=ict(1), ict(1) is set equal to iterat
c           in subroutine adjust of analiz link. if coefficient matrix
c           of normal equations could not be inverted, ict(1) was set
c           equal to iterat-1 in subroutine snorml or xnorml of analiz
c           link to prevent least squares analysis from iterating
c           further.)
      if(Iterat.lt.Ict(1)) then
c
c before reiterating, write punch data set of adjusted
c parameters onto disk for protection
c
         if(Ict(77).gt.1) then
            keep   = Ipunch
            Ipunch = Ict(77)
            call PUNCH
            if(Ipunch.eq.7) then
               endfile Ipunch
               rewind Ipunch
            endif
            Ipunch = keep
         endif
         goto 200
      else if(Iterat.ne.Ict(1)) then
         if(Ict(15).ne.0) then
            write(Iout, 1020) Ict(15)
 1020       format('0NORMAL STOP  NO MATRIX INVERSION  ICT(15)=', i5)
            goto 1300
         else
            n5 = 20
            if(repeat.gt.1) n5 = -n5
            call SUICID(
     .' NORMAL EQUATION INVERSION FAILED IN SNORML, OR SOME OTHER PROBLE
     .M. STOP IN MAIN', n5)
            goto 100
         endif
      endif
c
c-----------this is end of least squares loop---------------
c*start=2000
c
c           prediction of new observed minus theory and/or harmonic
c           analysis of residuals
 1100 if(Ict(10).ge.-1) then
         if(Ict(1).ge.0 .or. Ict(11).eq.-3) then
 
c filter prediction done in filctl
            if(Ict(42).le.0) call PRDICT(0._10, 3.E7_10, 1)
         endif
      endif
c
c computation of dummy observations at end of least squares
c analysis iteration (skipped if ict(1).le.0)
 1200 ict79 = Ict(79)
      if(Jct(79).gt.0) call PLANET(1)
      if(iabs(ict79).eq.1) call COMPAR(1)
c
c punch adjusted parameters
      if(Ict(1).gt.0) then
c
c
c*start=9000
c the order of operations here is poor
         if(Ict(77).gt.0) call PUNCH
      endif
 1300 n5 = 5
      if(repeat.gt.0) n5 = -5
      call SUICID(' NORMAL STOP IN MAIN', n5)
      goto 100
      end
