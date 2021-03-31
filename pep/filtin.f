      subroutine FILTIN(in0,nstop,init)
 
      implicit none

c arguments 
      integer*4 in0, nstop
      logical*4 init
c
c     d. white  april 1974  subroutine filtin
c    r. goldstein march 1975  modified to allow prdict to use filtered
c        solutions
c        paul macneil january, 1978 modified for iteration
c        modified for *command july 1978  r.b. goldstein
c          modified for split into filtin, filtpn, filtep
c                    z. goldberg   march 1980
c
c         filter input to pep consists of a namelist (nmlst3) followed
c         by a list of names.  quantities in the namelist are:
c
c             timez - r*8 - jd.fract at beginning of first filter epoch
c
c             delta - r*8 - length of time in one epoch, expressed as
c                 jd.fract
c        sweigh - r*8 - weight factor for individual s-matrices
c                                                (default is 1)
c
c
c            delprd - r*8 - time span expressed as jd.fract for which to
c                do individual time span predicts
c
c             nepoch - i*4 - number of filter epochs (max set in prepda)
c
c             npnp - i*4 - number of process noise parameters (dfault 6)
c
c             lfilt(nepoch) - i*2 - if lfilt(i).ne.0, do solution
c                for epoch i, if any lfilt(i).ne.0  where i<nepoch
c                 then by implication a backwards filtering operation
c                 takes place.
c                lfilt(i)>0  for which i is smallest defines the pep
c                    internal saved solutions  for iteration purposes
c
c            lprdct(nepoch) -i*2- if lprdct(i).eq.0  do not do a prdict
c                for that epoch. if lprdct(i).eq.1 do a time span prdict
c                for that epoch. the time span is centered at the epoch
c                center and of a width equal to delprd
c
c
c             fict(50) - i*2 - filter controls
c                 fict(1) = 0  nothing
c                 fict(1) = 1  print propagated info matrix each epoch
c                 fict(1) = 2  print propagated rhs  matrix each epoch
c                 fict(1) = 3  print both
c
c                 fict(2) = 0  forward only filtering
c                 fict(2) = 1  two way filtering (can be implied by
c                     values of lfilt)
c
c                 fict(3) =-1  no cleanup nor error matrix print
c        must use fict(3) = -1 because of changes for extended
c        precision     pem  june, 1979
c                 fict(3) = 0  only error matrix print
c                 fict(3) = n  n iterations of inverse cleanup in addsmr
c                     with error matrix print
c
c                 fict(4) = 0  no scaling for dpinv nor printout
c        must use fict(4) = 0 because of changes for extended
c        precision     pem  june, 1979
c                 fict(4) = 1  print unscaled h matrix and inverse
c                 fict(4) = 2  scaling by 1/sqrt(max(row))*h*
c                     1/sqrt(max(col))
c                 fict(4) = 3  print scaled and unscaled h and h inverse
c
c                fict(5)=0  no saved solution (complete) prdict
c                fict(5)=1  do a complete saved solution prdict
c
c
c                fict(6)=0  do not supress printout in analyze link when
c                           using filter
c                fict(6)=1  supress all but the last page of printout in
c                           analyze when using filter..also supresses
c                           jout data set
c
c          fict(7)=0     do not save normal equations for any epoch
c          fict(7)>0     save normal equations for this epoch
c          fict(7)<0     save normal equations counting backwards
c                        from last epoch. the usual situation would be
c                        fict(7)=-1 , filtering forward only, and
c                        saving normal equations from last epoch.
c
c          fict(8)<=0    read filter epochs from *state (if input), or
c                        calculate from timez & delta (if not)
c          fict(8)>0     read filter epochs from w data set (iconof)
c
c             feps(20) -r*4 - filter testing quantities
c                 feps(1)  convergence criterion for inverse cleanup
c                     in addsmr
c
c                 feps(2)  scale factor for smearing matrices:
c                     smat = smat*feps(2)
c
c                 feps(3) tolerance for smat versus fep epoch
c                      discrepancies in rsmat
c
c        maxe - i*4 - nepoch.le.maxe.le.maxep, used for space
c                     allocation by prepda
c
c        maxp - i*4 - nparam.le.maxp, used by prepda
c                     for space allocation
c
c             insne - i*4 - input saved normal eqn data set (created in
c                 nrmict stage)
c             filter - i*4 - filter data set (created in form stage)
c             outsne - i*4 - output saved normal eqn data set (created
c                 in use stage)
c             smat - i*4 - smearing matrix data set (input to pep)
c        iconof - i*4 - initial condition (really, process
c                     noise parameter) offset data set
c         wzero - i*4 - filter epochs & initial conditions (pnp's)
c                       before first integration
c
c         common
      include 'filnit.inc'
      include 'filtda.inc'
      include 'filtds.inc'
      include 'filtim.inc'
      include 'inodta.inc'
c
c namelist
      namelist /NMLST3/timez, delta, Nepoch, Insne, Filter, Outsne,
     .         delprd, Lprdct, Filprm, Sweigh, Wzero, Smat, Npnp, Lfilt,
     .         Fict, Feps, Iconof, Maxe, Maxp
c
c local
      integer*4 i, n
      integer*4 maxep/250/
c
c initialize
      Nml3rd = .false.
      Pnpnrd = .false.
      Feprd  = .false.
      Nepcht = 0
 
      Maxfep = 100
      do i = 1, Maxfep
         Fep(i) = 0._10
      end do
      Nepoch = 0
      Npnp   = 6
      Insne  = 0
      Filter = 0
      Outsne = 0
      Smat   = 0
      Iconof = 0
      Wzero  = 0
      Maxe   = 0
      Maxp   = 0
      do i = 1, 250
         Sweigh(i) = 1.0_10
         Lprdct(i) = 0
         Lfilt(i)  = 0
      end do
      do i = 1, 50
         Filprm(i) = 0._10
         Fict(i)   = 0
      end do
      Fict(3) = -1
      do i = 1, 20
         Filflg(i) = .false.
         Feps(i)   = 0.0
      end do
 
c mfile, lfile, kfile initialized by subroutine prepda
      Feps(2) = 1.0
 
      if(init) return
c
c spool, echo, and read &nmlst3
c spool  &nmlst3
      call PEPTIC(In,Iout,in0,6,'FILTER CONTROLS  &NMLST3', nstop,
     .            0)
      read(in0,NMLST3)
      rewind in0
      Nml3rd = .true.
c
c
c check out input parameters
      if(Nepoch.le.Maxe .and. Maxe.le.maxep) then
 
         n = Nepoch - 1
         do i = 1, n
            if(Lfilt(i).ne.0) then
               Fict(2) = 1
               goto 50
            endif
         end do
   50    return
      else
         write(Iout,100)
  100    format(' *** NOT ENOUGH EPOCH SPACE IN FILTIN AT LABEL=50 ***'
     .         )
         nstop = nstop + 1
         return
      endif
      end
