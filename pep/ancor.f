      real*10 function ANCOR(ind)
 
      implicit none
 
c           a.t.y. ng 1976
c     revised and commented by j.f.chandler 1979
c     compute the correction to the delay referred to station location
c     to make it referred to the tracking point on the 2ndary axis.
c       tmdly(s.l.) + ancor = tmdly(t.p.)
c
c     see moyer's report, jpl t.r. 32-1527, p.81
c
c     revised feb. 1991  marc murison
c
c arguments
      integer ind, n
c     ind - number (1 or 2) of site
c     n   - number of sites to set up in ANCORS

c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'coord.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'nutprc.inc'
      include 'param.inc'
      include 'sitcrd.inc'
      character*4 sitf4(2,2)
      equivalence (sitf4(1,1),sitf(1))

c local variables
      real*10 unit(3),stht,ctht
      integer i, index, j
      character*4 names(10)/'11DS', '12DS', '13DS', '41DS', '42DS',
     .   '51DS', '61DS', '62DS', 2*'    '/
      integer*4 kind(10)/    1,      1,      2,      1,      1,
     .    1,      1,      1,     2*0/
c
c antenna axis offsets in km.
      real*10 minusb(3)/ -6.706E-3_10, -.9144E-3_10, -1.2192E-3_10/

c external functions
      real*10 DOT

c alternate function entry point
      real*10 ANCORS
c
c
c entry ancor--no set up
c
c
c compute antenna correction (if any)
      index = Jsite(ind)
      if(index.eq.1) then
c
c index=1: ha-dec type mounting
c unit is true pole direction
         do j = 1, 3
            unit(j) = Nutpr(3, j)
         end do
      else if(index.eq.2) then
c
c index=2: az-elev type mounting
c unit is vector normal to geoid
         do j = 1, 3
            unit(j) = Sitnrm(j, ind)
         end do
      else
c
c index=0: no correction
         ANCOR = 0._10
         return
      endif
c
c calculate the correction
      stht  = DOT(unit, Xsitp0(1,ind))
      ctht  = SQRT(1._10 - stht**2)
      ANCOR = (minusb(index)/Ltvel)*ctht
      if(mod(Jct(6)/4096,2).eq.1) then
         if(Line.gt.56) call OBSPAG
         write(Iout,50) ind,ANCOR
   50    format(1x,'ANCOR(',i1,'): ',1pd19.11)
         Line = Line + 1
      endif
      return
 
c
c
c entry ancors
c
c set up station types in advance (called from radar)
c
      entry ANCORS(n)
      ANCORS = 0._10
      do i = 1, n
c
c search for match in station name
         do j = 1, 10
            if(sitf4(1,i).eq.names(j)) then
               Jsite(i) = kind(j)
               goto 100
            endif
         end do
c
c station not in table, tmdly=0.
         Jsite(i) = 0
  100 end do
c
c
      return
      end
