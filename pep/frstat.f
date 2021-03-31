      subroutine FRSTAT(plevel,ulevel,title)
 
      implicit none
c
c r.reasenberg,r.white   oct 1972   subroutine frstat
c
c         with entry tally.
c         collects and writes statistics on parameter adjustments
c         three levels of statistics can be compiled.  level 3 stats
c         will include several groups of level 2 stats which, in turn,
c         will include several groups of level 1 stats.  an example of
c         level 1 stats would be the s series of harmonic coefficients
c         which could be subsumed at level 2 with other harmonic coeff.
c         and at level 3 with initial conditions and other parameters
c         for a planet.
c
c         at entry tally the level 1 statistics are compiled on call
c         from adjast.
c
c         at the main entry to the routine there are three steps in
c         processing.  first, a summary is written for the level indi-
c         cated by plevel.  then the level indicated by ulevel is
c         updated by adding the totals of level plevel to those of
c         ulevel, where ulevel.gt.plevel.  finally, level plevel is
c         reinitialized (zeroed).  if ulevel.le.plevel then plevel
c         is summarized and reinitialized but no higher level is
c         updated.  the special case where plevel = 0 results in
c         (re)initialization of levels 1 through ulevel.  title is an
c         8 character string which can serve as a label on the written
c         summary.
c         to avoid redundant summarys, nu(ulevel) is incremented upon
c         updating and, therefore, nu(plevel) must be.gt.1 for
c         summarizing.
c
c         parameters
      integer*4 plevel,ulevel
      character*8 title
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'adjstf.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
c
c         local variables
c             nu(level) = number of updates at level
c             m(level) = number of adjustments in group for which stats
c                        included (sig.ne.0)
c
c             s((level-1)*3 + 1) = sum of fracts
c             s((level-1)*3 + 2) = sum of abs of fracts
c             s((level-1)*3 + 3) = sum of squares of fracts
c
c             t(1) = average of fracts
c             t(2) = average of abs of fracts
c             t(3) = average of squares of fracts
c             t(4) = root mean square of fracts
      real*10 s(9),t(4),x
      integer*4 nu(3),m(3),i,j,k
      integer   linexp
 
c summary section
c
      if(plevel.ge.1) then
         if(plevel.gt.3) return
         if(m(plevel).le.0) return
         if(nu(plevel).gt.1) then
            x = m(plevel)
            j = (plevel - 1)*3
            do i = 1,3
               t(i) = s(j + i)/x
            end do
            t(4) = SQRT(t(3))
            j    = j + 1
            k    = j + 2
            if(Ict(45).gt.0) then
               linexp = 3
               if(plevel.lt.2) linexp = 4
               call PAGCHK(60,linexp,0)
               write(Iout,100) plevel,m(plevel),title
               if(Jout.gt.0) write(Jout,100) plevel,m(plevel),
     .            title
  100          format(' LEVEL',i2,' STATISTICS FOR',i3,
     .          ' ADJUSTMENTS OF ',a8,' PARAMETERS')
               if(plevel.lt.2) then
                  write(Iout,200)
                  if(Jout.gt.0) write(Jout,200)
  200             format('   SUM FRACT',6x,'SUM ABS FRACT',2x,
     .             'SUM FRACT**2',3x,'AVG FRACT',6x,'AVG ABS FRACT',2x,
     .             'AVG FRACT**2',3x,'ROOT MEAN SQUARE FRACT')
               endif
               write(Iout,210) (s(i),i = j,k),t
               if(Jout.gt.0) write(Jout,210) (s(i),i = j,k),t
  210          format(1p,7E15.5)
               if(Line.lt.56) call PAGHED(0)
            endif
         endif
c
c update section
c
         if(ulevel.gt.plevel) then
            if(ulevel.ge.2) then
               if(ulevel.le.3) then
                  nu(ulevel) = nu(ulevel) + 1
                  m(ulevel)  = m(ulevel) + m(plevel)
                  j = (ulevel - 1)*3
                  k = (plevel - 1)*3
                  do i = 1,3
                     s(j + i) = s(j + i) + s(k + i)
                  end do
               endif
            endif
         endif
      endif
c
c reinitialization section
c
      j = max0(plevel,1)
      k = 0
      if(plevel.le.0) k = ulevel
      k = max0(k,plevel)
      do i = j,k
         nu(i) = 0
         m(i)  = 0
      end do
      j = (j - 1)*3 + 1
      k = k*3
      do i = j,k
         s(i) = 0.0_10
      end do
      return
c
c tally section
c
      entry TALLY
      nu(1) = nu(1) + 1
      m(1)  = m(1) + 1
      s(1)  = s(1) + Fract
      s(2)  = s(2) + ABS(Fract)
      s(3)  = s(3) + Fract**2
      return
      end
