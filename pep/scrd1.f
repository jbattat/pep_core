      subroutine SCRD1(lice)
 
      implicit none

c
c m.e.ash    may 1970        subroutine scrd1
c reference to s-body tape added 1977 jul - j.f.chandler
c first five records of observing body tape are read
c
c parameters 
      integer*2 lice
c     lice =0 printout of data on first two records of satellite-probe
c     tape
c     lice =1 no such printout

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'namtim.inc'
      include 'number.inc'
      include 'scdta.inc'
      include 'trpcom.inc'
      include 'yvectrd1.inc'

c local variables
      integer   intsc

c
c read body constants from disk
c test if coordinates available from n-body tape and,
c if so, set up controls as if from individual tape
      call XXRDBD(Nplnt(Klans1),Jsc,Klans1,Jdsc,intsc,Idirsc,Ksc,
     . Scintx,Nkisc,Kisc,2,Sccom)
      if(Jtest.eq.0) then
c
c read first two records of satellite-probe peripheral data set
         call XXRD1(lice,Nplnt(Klans1),Jsc,Klans1,Jdsc1,Jdsc2,
     .    Iparsc,i_mxeqn,intsc,Idirsc,Ksc,Scintx,Frsc1,Nkisc,Kisc)
 
c check if input central body matches that on tape
         if(Npcent(Klans1).ne.L2) call SUICID(
     .'WARNING IN SCRD1: INPUT NCENTR DOES NOT MATCH NCENTR ON PROBE TAP
     .E  ', -17)
      endif
 
      Scint = Scintx
      if(Idirsc.lt.0) Scint = -Scintx
      T0svsc = Jdxx0
c
c count partials and reset jd0 for next iteration
      call XXRDCK(Lparsc,Kisc,Klans1,1)
c
c read first three data records of satellite-probe data set
      if(Jtest.eq.0) then
         if(Ksc(88).gt.0) then
            call SCRED1(0)
         else
            call SCRAD1(0)
         endif
      endif
      return
      end
