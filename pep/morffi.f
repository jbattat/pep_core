      subroutine MORFFI(n0)
 
      implicit none
 
c
c rj cappallo  january 1978  sr morffi
c
c morffi calculates the torque due to the earth-figure moon-figure
c interaction and adds it to the rigid body torque, n0. formulae
c can be found in rjc's thesis.

c arguments
      real*10 n0(3)
      integer*2 iparm

c array dimensions
      include 'globdefs.inc'
 
c common
      include 'empcnd.inc'
      include 'harmor.inc'
      real*10 alpha,beta,gamma
      equivalence (Mirat(1),alpha),(Mirat(2),beta),(Mirat(3),gamma)
      include 'intstf.inc'
      include 'mnrtlb.inc'
      include 'morcrd.inc'
      include 'morstf.inc'
      real*10 am,bm,cm
      equivalence (am,I0(1,1)),(bm,I0(2,2)),(cm,I0(3,3))
      include 'nutces.inc'
      include 'param.inc'

c local variables
      real*10 g(3,3),unp(3),ume(3),unps(3),umes(3),udec(3),udecs(3),
     . ncr(3),temp(3),uphis(3),dumesdy(3),dudecsdy(3),duphisdy(3),
     . dgdy(3,3),dumesdx(3),dudecsdx(3),duphisdx(3),dgdx(3,3)
      real*10 sind,cosd,cos2d,sin2d,sinsqd,udecmg,
     . cofrr,cofdr,cofdd,cofpp,dudotd,dsindx,dsin2dx,dsinsqdx,
     . dcofrrdx,dcofdrdx,dcofdddx,dcofppdx,mfact
      integer i,j,j1,l


c external function
      real*10 DOT

c
      do i = 1,3
 
c get north pole of date components in reference and selenodetic systems
         unp(i) = Nutprc(3,i)
         unps(i) = Mrotlb(i,1)*Nutprc(3,1)+Mrotlb(i,2)*Nutprc(3,2)
     .    +Mrotlb(i,3)*Nutprc(3,3)
 
c also unit vector from earth's com to moon's
         ume(i) = Mecor(i)/Rem
         umes(i) = -Smecor(i)/Rem
      end do
 
c udec=unit vector in direction of increasing declination
      call CROSS(unps,umes,ncr)
      call CROSS(umes,ncr,udecs)
      udecmg = SQRT(DOT(udecs,udecs))
      do i = 1,3
         udecs(i) = udecs(i)/udecmg
      end do
 
c uphi transverse to np-lunar com plane
      call CROSS(udecs,umes,uphis)
 
c calculate declination of lunar com
      sind = DOT(unps,umes)
c      dec = ASIN(sind)
 
c calculate and store the three intermediate vectors in g
      sinsqd = sind*sind
      cosd   = SQRT(1._10-sinsqd)
      sin2d  = 2._10*sind*cosd
      cos2d  = 1._10-2._10*sinsqd
      cofrr  = 6._10-18._10*sinsqd
      cofdr  = 6._10*sin2d
      cofdd  = -4.5_10+10.5_10*sinsqd
      cofpp  = 7.5_10*sinsqd-1.5_10

c j is index to vector components; i tells which g we have
      do j = 2,3
         j1 = j - 1
         do i = 1,j1
            g(i,j) = Gamat3*Erzhar(1)*Erad**2/Rem5*(
     .       (cofrr*umes(i)+cofdr*udecs(i))*umes(j)
     .       +(cofdr*umes(i)+cofdd*udecs(i))*udecs(j)
     .       +cofpp*uphis(i)*uphis(j))
         end do
      end do
 
c finally: calculate the torque components and add to n0
      n0(1) = n0(1) + (cm - bm)*g(2,3)
      n0(2) = n0(2) + (am - cm)*g(1,3)
      n0(3) = n0(3) + (bm - am)*g(1,2)
 
      return
c
c now calculate correction to partial derivatives

c Note: although entry MORFFI calculates a correction to the torque,
c the partial derivatives are taken of the angular accelerations.

c Entry MORFFI2 calculates only corrections to the indirect terms in
c the partials due to the dependence of the acceleration on the state
c vector.  The direct contributions due to the dependence of the
c figure-figure component of acceleration on parameter values is
c done with one call to MORFFIP for each parameter.

      entry MORFFI2(iparm)

      call PRODCT(Mrotlb,udecs,udec,-3,3,1)

c partials w.r.t. rotation state
      do l=1,3
         call PRODCT(Drtdy(1,1,l),ume,dumesdy,-3,3,1)
         call PRODCT(Drtdy(1,1,l),udec,dudecsdy,-3,3,1)
         call CROSS(dudecsdy,umes,duphisdy)
         call CROSS(udecs,dumesdy,temp)
         do i=1,3
            duphisdy(i)=duphisdy(i)+temp(i)
         end do

         do j = 2,3
            j1 = j - 1
            do i = 1,j1
               dgdy(i,j) = Gamat3*Erzhar(1)*Erad**2/Rem5*(
     .          (cofrr*umes(i)+cofdr*udecs(i))*dumesdy(j)
     .          +(cofrr*dumesdy(i)+cofdr*dudecsdy(i))*umes(j)
     .          +(cofdr*umes(i)+cofdd*udecs(i))*dudecsdy(j)
     .          +(cofdr*dumesdy(i)+cofdd*dudecsdy(i))*udecs(j)
     .          +cofpp*duphisdy(i)*uphis(j)+cofpp*uphis(i)*duphisdy(j))
            end do
         end do
         ddwdy(1,l) = ddwdy(1,l) + alpha*dgdy(2,3)
         ddwdy(2,l) = ddwdy(2,l) - beta*dgdy(1,3)
         ddwdy(3,l) = ddwdy(3,l) + gamma*dgdy(1,2)
      end do

      if(iparm.le.1) return
c compute partials w.r.t. orbit state
      do l=1,3
         do i=1,3
            dumesdx(i)=(Mrotlb(i,l)-umes(i)*ume(l))/Rem
         end do
         call CROSS(unps,dumesdx,temp)
         call CROSS(umes,temp,dudecsdx)
         call CROSS(dumesdx,ncr,temp)
         do i=1,3
            dudecsdx(i)=(dudecsdx(i)+temp(i))/udecmg
         end do
         dudotd=DOT(dudecsdx,udecs)
         do i=1,3
            dudecsdx(i)=dudecsdx(i)-udecs(i)*dudotd
         end do
         call CROSS(dudecsdx,umes,duphisdx)
         call CROSS(udecs,dumesdx,temp)
         do i=1,3
            duphisdx(i)=duphisdx(i)+temp(i)
         end do
         dsindx=(unp(l)-ume(l)*sind)/Rem
         dsin2dx=2._10*cos2d*dsindx/cosd
         dsinsqdx=2._10*sind*dsindx
         dcofrrdx=-18._10*dsinsqdx
         dcofdrdx=6._10*dsin2dx
         dcofdddx=10.5_10*dsinsqdx
         dcofppdx=7.5_10*dsinsqdx

         do j = 2,3
            j1 = j - 1
            do i = 1,j1
               dgdx(i,j) = Gamat3*Erzhar(1)*Erad**2/Rem5*(
     .          (cofrr*umes(i)+cofdr*udecs(i))*dumesdx(j)
     .          +(cofrr*dumesdx(i)+cofdr*dudecsdx(i))*umes(j)
     .          +(cofdr*umes(i)+cofdd*udecs(i))*dudecsdx(j)
     .          +(cofdr*dumesdx(i)+cofdd*dudecsdx(i))*udecs(j)
     .          +cofpp*duphisdx(i)*uphis(j)+cofpp*uphis(i)*duphisdx(j)
     .          +(dcofrrdx*umes(i)+dcofdrdx*udecs(i))*umes(j)
     .          +(dcofdrdx*umes(i)+dcofdddx*udecs(i))*udecs(j)
     .          +dcofppdx*uphis(i)*uphis(j))
     .          -5._10*g(i,j)*Mecor(l)/Rem2
            end do
         end do
         Ddwdx(1,l) = Ddwdx(1,l) + alpha*dgdx(2,3)
         Ddwdx(2,l) = Ddwdx(2,l) - beta*dgdx(1,3)
         Ddwdx(3,l) = Ddwdx(3,l) + gamma*dgdx(1,2)
      end do

      return

c corrections to sensitivity to parameters
      entry MORFFIP

      if(Icrtrl(kkk).eq.3 .or. Icrtrl(kkk).eq.10) then
         if(Icrtrl(kkk).eq.3) then
            mfact=1._10/Mass(3)
         else
            mfact=-1._10/Masse
         endif
         Ddwdp(1) = Ddwdp(1) + mfact*alpha*g(2,3)
         Ddwdp(2) = Ddwdp(2) - mfact*beta*g(1,3)
         Ddwdp(3) = Ddwdp(3) + mfact*gamma*g(1,2)
      else
         Ddwdp(1) = Ddwdp(1) + Dalpha(Kkk)*g(2,3)
         Ddwdp(2) = Ddwdp(2) - Dbeta(Kkk)*g(1,3)
         Ddwdp(3) = Ddwdp(3) + Dgamma(Kkk)*g(1,2)
      endif
      return

      end
