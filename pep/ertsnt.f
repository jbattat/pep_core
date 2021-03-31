      subroutine ERTSNT(jd,fract,x,nvel,ntype)
 
      implicit none
c
c m.e.ash   dec 1971    subroutine ertsnt
c
c arguments
      real*10 fract, x(6)
      integer jd,nvel,ntype
c     x      = output coordinates of earth or earth-moon barycenter
c              relative to sun determined from n-body tape
c     output units are au and au/day
c     input quantities are
c        jd     = julian day number
c        fract  = coordinate time fraction of day
c        nvel   = 0 x(1-3) position output
c        nvel   = 1 x(4-6) velocity output
c        ntype  =-1 moon relative to earth
c        ntype  = 0 earth-moon barycenter relative to sun
c        ntype  = 1 earth relative to sun

c array dimensions
      include 'globdefs.inc'

c commons
      include 'bdctrl.inc'
      include 'bddta.inc'
      include 'bdydta.inc'
      include 'comdateq.inc'
      include 'tabval.inc'

c external functions
      real*10 TERPF,DTRPF
c
c local variables
      real*10 yent(7,3,6),ymnt(7,3,3)
      real*10 pm(4),xm(6)
      real*10 pe(4)
      integer i,j,jo,ngo,nm,nr
c
c read n-body tape
      if(Nbody.le.0) call SUICID(' NO N-BODY TAPE, STOP IN ERTCNT ', 8)
      call BDREED(jd,0,0)
      if(jd.le.0) return
c
c perform embary interpolation
      if(ntype.ge.0) then
         call EVTRP(jd,fract,nvel,0,9,yent,x,Body(1,1,2),1,
     .              Jdbd, Frct, pe)
         if(ntype.le.0) return
      endif
c
c calculate moon interpolation coefficients
      ngo = 2
      if(Ibdsgn.le.0) ngo = 3
      nm    = 2*(jd - Jdbd(ngo)) + 35
      pm(1) = (fract - 0.5_10)*2.0_10
      if(pm(1).ge.0.0_10) nm    = nm + 1
      if(pm(1).lt.0.0_10) pm(1) = fract*2.0_10
      pm(2) = pm(1)**2
      pm(3) = 1.0_10 - pm(1)
      pm(4) = pm(3)**2
      if(Ibdsgn.le.0) nm = 122 - nm
c
c calculate moon y-vector if necessary
c use 10-point interpolator for the moon in this one exception
      if(nm.ne.Nmmm) then
         Nmmm  = nm
         Ntab1 = 2
         Ntab2 = 3
         do j = 1, 3
            do i = Ntab1, 11
               nr = nm + ISIGN(i,Ibdsgn)
               Tabvl(i) = Mon(j,nr)
            end do
            call YCOFF(ymnt(1,1,j))
         end do
      endif
      if(nvel.gt.0) then
c
c get velocity
         do j = 4, 6
            xm(j) = DTRPF(pm,ymnt(1,2,j-3),0.5_10*Ibdsgn)
         end do
      else
c
c determine moon coordinates
         do j = 1, 3
            xm(j) = TERPF(pm,ymnt(1,2,j))
         end do
      endif
 
      jo = 3*nvel
c
c is moon output desired
      if(ntype.ge.0) then
c
c determine earth relative to sun
         do i = 1, 3
            x(i + jo) = x(i + jo) - Mnfct*xm(i + jo)
         end do
      else
         do i = 1, 3
            x(i + jo) = xm(i + jo)*Mnau
         end do
      endif
      call TRPLST(jd,fract,nvel,nvel,'ERTSNT',x)
 
      return
      end
