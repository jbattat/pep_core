      subroutine SOTRP(jd,fract,x,nvel1)
 
      implicit none

c subroutine SOTRP - J.F.Chandler - 1983 Aug
c based on 'SOLCNT' - M.E.Ash - 1971 Dec

c arguments
      integer*4 jd,nvel1
      real*10 fract,x(6)
c     x     = output coordinates of center of mass of solar system
c              relative to sun determined from n-body tape
c              plus asteroid center-of-mass tape (if any)
c     output units are au and au/day
c     input quantities are
c        jd    = julian day number
c        fract = coordinate time fraction of day
c        nvel1 = 0 x(1-3) position output
c        nvel1 = 1 x(4-6) velocity output

c local
      real*10 xcnt(6,15),ycnt(7,3,6),ymrc(7,3,6),pqbd(4),pqmrc(4),
     . sum,xoff(6)
      integer*4 i,in1,index,j,k,nvel,nvelsv

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'bdctrl.inc'
      include 'bddta.inc'
      include 'bdydta.inc'
      include 'cmcke.inc'
      include 'number.inc'
      include 'param.inc'
      include 'trpcom.inc'
 
      nvel = nvel1
c
c read n-body tape
  100 if(Nbody.le.0) call SUICID(' NO N-BODY TAPE, STOP IN SOTRP  ',8)
      call BDREED(jd,0,0)
      if(jd.le.0) return
c
  200 in1 = 1
      index = 3
      if(nvel.gt.0) index = 6
      if(Jdcnts.ne.Jdbd(1)) then
         Ntb1s(8) = 9999
         Ntb1s(10)= 9999
      endif
      if(Nplnt0.eq.-4) then

c do individual interpolation for each body, so that partial derivatives
c with respect to mass can be calculated later
         if(nvel.gt.0) in1 = 4
         do i=in1,index
            xoff(i)=0._10
         end do
         call PLCMC(jd,fract,xoff,nvel)
         if(Sumcom.ne.0._10) then
            call EVTRP(jd,fract,-nvel,-1,8,ycnt,x,Comcrd,1,Jdbd,Frct,
     .       pqbd)
            do i=in1,index
               xoff(i)=xoff(i)-x(i)*Sumcom
            end do
         endif
c PLCMC returns position of Sun relative to SSBC, but we want SSBC
c relative to Sun
         do i=in1,index
            x(i) = -xoff(i)/Mascnt
         end do

      else

c interpolate only Mercury and, separately, the composite of the other
c planets and asteroids
c setup center of mass of solar system coordinate tabulation
         if(Jdcnts.eq.Jdbd(1)) then
            if(nvel.le.nvelsv) goto 300
            in1 = 4
         endif
         nvelsv = nvel
         do i = in1, index
            do j = 1, 15
               sum = Comcrd(i,j)*Sumcom
               do k = 1, 8
                  sum = sum + Body(i,j,k)*Mass(k + 1)
               end do
               xcnt(i,j) = sum
            end do
         end do
c
c perform interpolation
  300    call EVTRP(jd,fract,nvel,-1,8,ycnt,x,xcnt,1,Jdbd,Frct,
     .    pqbd)
         call EVTRP(jd,fract,nvel,-1,10,ymrc,xoff,Merc,1,Jdbd,Frct,
     .    pqmrc)
         do i=in1,index
            x(i)=(x(i)+xoff(i)*Mass(1))/Mascnt
         end do
      endif
      Jdcnts = Jdbd(1)
      call TRPLST(jd,fract,nvel,0,'SOTRP(CM)',x)

      if(nvel.ge.nvel1) return
      nvel = nvel1
      goto 200
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry SOLCNT(jd,fract,x,nvel1)
c
c same as main entry point except that nvel1=1 means do both
c position and velocity
c
      nvel = 0
      goto 100
      end
