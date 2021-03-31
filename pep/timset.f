      subroutine TIMSET(ihr, imin, sec, itot, date)
 
      implicit none
 
 
c
c m.e.ash   nov 1966    subroutine timset to set and read timers
c
      character*4 date(2), czero
      real*4    sec(2), zero4/0./
      equivalence (czero,zero4)
      integer*4 ihr(2), imin(2), itot(2)
 
c if the value of date(1) sent to this subroutine is 0, then
c the calendar date in a format is returned in date(1-2).
c if the value of date(1) sent to this subroutine is not 0,
c then the date array is not altered.
c in the following, i=1 denotes real time and i=2 denotes
c task time elapsed since timers were last read and set.
c   itot(i)= binary integer number of hundreths of seconds elapsed.
c   ihr(i) = hours elapsed.
c   imin(i)= minutes elapsed.
c   sec(i) = seconds elapsed.
 
      include 'timstf.inc'
c dat0  = packed date for previous reading and setting.
c dat1  = packed date for this reading and setting.
c ireal0= real time (binary integer hundreths of seconds) for
c         previous reading and setting.
c ireal1= real time (binary integer hundreths of seconds) for
c         this reading and setting.
c itask = increment of task interval timer counted from when TIMDAT was
c         previously called (units of centiseconds).
c itotsk= total task time elapsed from start of program run in
c         binary integer number of hundreths of seconds.
c idayr = day of year in a format (only calculated if date(1)
c         has value 0 when sent to subroutine timdat).
c
      call TIMDAT(date)
      itot(2) = Itask
      Itotsk  = Itotsk + itot(2)
      itot(1) = Ireal1 - Ireal0
      if(Dat1.ne.Dat0) itot(1) = 8640000 + itot(1)
      call CONVRT(ihr(1), imin(1), sec(1), itot(1))
      if(itot(1).le.0) itot(1) = 1
      call CONVRT(ihr(2), imin(2), sec(2), itot(2))
      Dat0   = Dat1
      Ireal0 = Ireal1
      return
      end
