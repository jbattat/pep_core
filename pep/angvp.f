      subroutine ANGVP
 
      implicit none

c j.f.chandler - 1979 february - subroutine angvp
c calculate relative velocity of observed body in au/sec

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'coord.inc'
      include 'coordxeq.inc'
      include 'number.inc'
      include 'obscrd.inc'
      include 'jdfrequv.inc'
      include 'param.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'

c local
      integer   i
c
c decide if central body is moon or planet
      if(Klan.gt.0 .and. Klan.le.u_mxpl) then
 
c get planet velocity
         call PLTRP(1,Jdy,fr(3),-1,0)
      else if(Klan.ne.u_mxpl+1) then
 
c must be earth satellite or asteroid
         do i = 4,6
            Xp(i) = 0._10
         end do
      else
 
c get moon velocity in units of au/day
         call MNTRP(1,Jdx,fr(2),-1,0,1)
         do i = 4,6
            Xp(i) = Xm(i,1)*Mnau
         end do
      endif
      if(Klanb.gt.0) then
 
c get satellite velocity
         call SBTRP(1,Jdy,fr(3),-1,0)
         do i = 4,6
            Xp(i) = Xp(i) + Xsb(i)*Cmfct
         end do
      endif

c convert to au/sec
c check if cis-lunar object -- if not, then get rel. to earth
      do i=4,6
         Xsitep(i,1)=Xsite(i,1)/Aultsc
         if(Kst1.ne.2) Xsitep(i,1)=Xsitep(i,1)+x(i,2)/Secday
         Xsitep(i,1)=Xsitep(i,1)-Xp(i)/Secday
      end do
      return
      end
