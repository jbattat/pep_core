      subroutine CONVRT(ihr, imin, sec, itime)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   ihr, imin, increm, itime
      real      sec
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   nov 1966    subroutine convrt to obtain hours,minutes,
c seconds from binary integer number of hundreths of seconds
c
      ihr    = itime/360000
      increm = itime - ihr*360000
      imin   = increm/6000
      sec    = increm - imin*6000
      sec    = sec/1.0E2
      return
      end
