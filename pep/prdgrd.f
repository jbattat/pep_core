      subroutine PRDGRD(l,ldim,lind,m,length)
 
      implicit none
 
c
c        r.b. goldstein  sept. 1978
c        called by prdeqs for setting up iptr and ipts vectors
c        for grid parameters for prdict or nrmict
c
c arguments
      integer*4 ldim
      integer*2 l(ldim,1),lind,length,m(1)

c array dimensions
      include 'globdefs.inc'

c        common
      include 'bernum.inc'
      include 'ktrap.inc'
      include 'wrkcompm.inc'

c local variables
      integer*4 i
c
c        continue to set up ipts and iptr vectors. ipts contains
c        complete topog. parameter set (all adjustable grid points),
c        but iptr is hardwired to contain only 4. therefore, jfk
c        is bumped for all l that are set, but lbj is bumped by
c        4 if any grid partials at all are anywhere in the complete
c        series.
c
      if(length.gt.0) then
         Lbj = Lbj + 4
 
c this line is only for partials on tape, but no l's input
         if(lind.gt.0) then
            do i = 1, length
               Lgrx(i) = 0
               if(l(lind,i).gt.0) then
                  Jfk = Jfk + 1
                  Ipts(Jfk) = 0
                  Lgrx(i)   = Jfk
               endif
            end do
         endif
      endif
      return
      end
c
c
c GRDSID was originally an entry in PRDGRD, sharing the local array of
c grid pointers. However, a bug in gfortran caused an addressing problem
c until the two were separated, 2011 Feb 8.

      subroutine GRDSID
      implicit none

c array dimensions
      include 'globdefs.inc'

c        common
      include 'bernum.inc'
      include 'ktrap.inc'
      include 'wrkcompm.inc'

      integer i,ipr
c
c change iptr vector for every observation
c
      do i = 1, 4
         ipr = Mixshp + i - 1
         Iptr(ipr) = 0
         if(i.le.Mnshob) Iptr(ipr) = Lgrx(Mshobs(i))
      end do
c
c
      return
      end
