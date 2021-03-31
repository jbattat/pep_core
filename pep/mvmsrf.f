      subroutine MVMSRF(s,srfmvm)
 
      implicit none

c mckinnis  sept 1974     subroutine mvmsrf
c theory by rad m georgevic
c solar radiation force for mvm flight
c interface between plnorb and ephemeris data set
      real*10 s
      real*10 srfmvm(3)

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'petuna.inc'
 
      integer*4 jsrf
      integer*4 init/0/
      real*4    timel, forcel(3),timer,forcer(3)
      data timel/0.E0/, forcel/3*0.E0/, timer/0.E0/, forcer/3*0.E0/
c epoch for srf calculations
c 1973, november 9, 19 hrs, 27 min, 58 sec
c julian day number = 2441996.31108586520
      real*10 epoch/2441996.31108586520_10/, dt
      real*4   factor,t
      integer   j

      dt = s - epoch
      t  = dt
      if(t.lt.0._10) call SUICID('MVMSRF TIME BEFORE EPOCH', 6)
      if(init.le.0) then
         init = 1
c on initial entry to program, must fill arrays
c the data set for read will be a value of kk,   kkp(78)
         jsrf = Kkp(78)
c
c timel,timer are left and right values of time for which the
c solar radiation force on mariner10 is forcel, forcer.
         read(jsrf,50) timel,forcel
         read(jsrf,50) timer,forcer
   50    format(1x,4E16.8)
      endif
c end of special initialization section
c test whether time is within old bounds
      if(t.lt.timel) call SUICID('MVMSRF TIME REVERSED', 5)
      do while( t.ge.timer )
 
c read in a new time interval
         timel = timer
         do j = 1, 3
            forcel(j) = forcer(j)
         end do
         read(jsrf,50,end=100) timer,forcer
      end do
 
c obtain components of force by linear interpolation.
      if(timer.le.timel) call SUICID('MVMSRF DATA WRONG   ', 5)
      factor = (t - timel)/(timer - timel)
      do j = 1, 3
         srfmvm(j) = forcel(j) + (forcer(j) - forcel(j))*factor
      end do
      return
  100 call SUICID('MVMSRF END OF DATA  ', 5)
      return
      end
