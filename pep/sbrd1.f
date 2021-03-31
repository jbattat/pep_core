      subroutine SBRD1(lice)
 
      implicit none

c
c m.e.ash   november 1967    subroutine sbrd1
c reference to s-body tape added 1977 jul - j.f.chandler
c first five records of satellite-probe tape are read
c modified for kalman filter by paul macneil dec., 1977
c
c parameters
      integer*2 lice
c lice =0 printout of data on first two records of satellite-probe tape
c lice =1 no such printout
c
c array dimensions
      include 'globdefs.inc'

c        common
      include 'comdateq.inc'
      include 'emcke.inc'
      include 'empcnd.inc'
      include 'flcomp.inc'
      include 'funcon.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'param.inc'
      include 'sbdta.inc'
      include 'trpcom.inc'
      include 'yvectrd1.inc'

c local variables
      integer   imass, intsb, kind
 
      Ifiltr = -1
c read body constants from disk
c test if coordinates available from n-body tape and,
c if so, set up controls as if from individual tape
      call XXRDBD(Nplnt(Klanb),Jsb,Klanb,Jdsb,intsb,Idirsb,Ksb,
     .            Sbintx,Nkisb,Kisb,1,Sbcom)
      if(Jtest.eq.0 .and. Ksb(88).ge.0) then
c
c read first two records of satellite-probe peripheral data set
         Iepoch = -2
         rewind Jsb
         call XXRD1(lice,Nplnt(Klanb),Jsb,Klanb,Jdsb1,Jdsb2,
     .    Iparsb,i_mxeqn,intsb,Idirsb,Ksb,Sbintx,Frsb1,Nkisb,Kisb)
         Ifiltr = Ifltrx
         if(Ifiltr.gt.0) then
 
c integration tape is in filter span format; initialize
            Tprev  = 0._10
            Tthis  = 0._10
            Tnext  = 0._10
            Iepoch = -1
            Nlflag = -1
         endif
 
c check if input central body matches that on tape
         if(Npcent(Klanb).ne.L2) call SUICID(
     .'WARNING IN SBRD1: INPUT NCENTR DOES NOT MATCH NCENTR ON PROBE TAP
     .E  ', -17)
      endif
 
      Sbint = Sbintx
      if(Idirsb.lt.0) Sbint = -Sbintx
      T0svsb = Jdxx0
 
c count partials and reset jd0 for next iteration
      call XXRDCK(Lparsb, Kisb, Klanb, 1)
c
      if(Jtest.eq.0) then
         if(Ksb(88).eq.-9) then
c
c set up for quasi-elliptic orbit instead of interpolation
            call PNITL(Icnd(Klanb), Pcond(1,Klanb), Tconx, Elptg(1,4))
            return
 
         else if(Ksb(88).eq.-8) then
c set up for elliptic orbit with partials
            Kisb(1)=-2
c case 1: satellite of major planet
            imass=Npcent(Klanb)
            kind=0
            if(imass.le.0) then
c case 2: asteroid, etc.
               imass=Nplnt(Klanb)
               kind=1
            end if
            call IMITL(Gauss,Mass(imass),Pcond(1,Klanb),kind,4)
 
c read first three data records of satellite-probe data set
         else if(Ksb(88).gt.0) then
            call SBRED1(0)
         else
            call SBRAD1
         endif
      endif
      return
      end
