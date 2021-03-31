      subroutine MOREDT(dw,didk,didt,diddk,diddt)
 
      implicit none
 
c
c rj cappallo  january 1980   sr moredt
c
c moredt calculates useful inertia tensor perturbation quantities
c for an elastic moon with constant-time-lag response.
c didk,didt   partials of delta i wrt k2 or k2t
c diddk,diddt similar partials for delta i dot

c array dimensions
      include 'globdefs.inc'

c commons
      include 'empcnd.inc'
      include 'harmor.inc'
      real*10 alpha,beta,gamma
      equivalence (Mirat(1),alpha),(Mirat(2),beta),(Mirat(3),gamma)
      include 'intstf.inc'
      include 'mnrtlb.inc'
      real*10 w(3)
      equivalence (w,W1)
      include 'morcrd.inc'
      include 'morstf.inc'
      include 'param.inc'
      include 'prtcod.inc'

      real*10 del(3,3)/1._10,3*0._10,1._10,3*0._10,1._10/
c plastic lunar deformation diagonal elements
c (+1/3,+1/3,-2/3) * (mean motion)**2
      real*10 plas(3)/2*1.76296E-2_10,-3.52591E-2_10/
 
c parameters
      real*10 dw(3),didk(3,3),didt(3,3),diddk(3,3),diddt(3,3)
      real*10 k,t,didy(3,3,6),diddy(3,3,6),didx(3,3,6),diddx(3,3,6)
      integer*2 iparm
      real*10 ddidp(3,3),ddiddp(3,3)
c local variables
      real*10 temp(3),temq(3),temr(3),tems(3),velspa(3),wcsx(3),
     . accspa(3),acc195(3),ndot(3),wdd(3),emass,mrad5,didddk(3,3),
     . wddw,wddwd,term,termr,term0,term1,term2,g3,rijd,wdotw,
     . wdotwd,wdotwj,dsvdy(3),dsady(3),
     . rdotv,rdotvd,rfct1,rfct2,rfct3,wddwdd,iw(3),idw(3),wi(3,3),
     . dwi(3,3),ir(3),iv(3),dnddy(3),dddwdy(3),wcr(3,3),dwcr(3,3),
     . wdr(3,3),wr(3,3),dddrdx(3,3),drdv,drdvd,rijx,rijxd,x3s(3),rij,
     . rijy,rijyd,dnddx(3),dddwdx(3),rir(3),x10s(3),mfact,
     . didkm(3,3),diddkm(3,3),didddkm(3,3),daccdm(3),daccdp(3),
     . didpk(3,3),diddpk(3,3),didddpk(3,3),dwdddi(3),dwdddp(3),
     . g3r5,rdaccdm,rdaccdp,wdotwjm,wddn(3),dmdp
      integer i,j,l
c external functions
      real*10 DOT

c Note: local arrays for symmetric tensors are only half filled.  Arrays
c returned as arguments are fully filled.
 
c precalculate useful variables and expressions
      wdotw  = DOT(w,w)/3._10
      wdotwd = DOT(w,dw)*2._10/3._10
      emass  = Mass(3)*Masse
      mrad5  = Mrad**5
      g3     = 3._10*Gamat
      g3r5   = 3._10*Gamat3/Rem5
c
c calculate earth velocity and accel. in spa frame
      call PRODCT(Mrotlb,Mecor(4),temp,3,3,1)
      call CROSS(w,Smecor,wcsx)
      do i = 1, 3
         velspa(i) = -temp(i) - wcsx(i)
      end do
      rdotv = DOT(Smecor,velspa)
 
c need accel. of earth in equatorial frame first
      do i = 1, 3
         acc195(i) = Gamtem*Mecor(i)/Rem3
     .                + Gamat*(Mcor(i)/Rm3 - Ecor(i)/Re3)
      end do
 
c resolve components of acc195 in spa system
      call PRODCT(Mrotlb,acc195,tems,3,3,1)
 
c then combine with kinematics to get spa acceleration
c also get partial of accel. wrt Mass(3), at constant position
c note that partial of accel. wrt Mass(10) is 2nd order in Rem/Re
      call CROSS(w,wcsx,temr)
      call CROSS(dw,Smecor,temq)
      call CROSS(w,velspa,temp)
      do i = 1, 3
         accspa(i) = tems(i) - 2._10*temp(i) - temq(i) - temr(i)
         daccdm(i) = -Gamat*Smecor(i)/Rem3
      end do
      rdotvd = DOT(Smecor,accspa) + DOT(velspa,velspa)
      rdaccdm= -Gamat/Rem
      rfct1=rdotv/Rem2
      rfct2=rfct1*rfct1
      rfct3=rdotvd/Rem2
c
c must also calculate angular jerk
c -first calculate time-derivative of torque (2nd order rigid body
      call PRODCT(I0,Smecor,ir,3,3,1)
      call CROSS(Smecor,ir,rir)
      call CROSS(velspa,ir,tems)
      call PRODCT(I0,velspa,iv,3,3,1)
      call CROSS(Smecor,iv,temq)
      do i = 1, 3
         ndot(i)=3._10*Gamat3*(temq(i)+tems(i)-5._10*rfct1*rir(i))/Rem5
      end do
 
c -then calculate the ang. jerk
      call PRODCT(I0,w,iw,3,3,1)
      call CROSS(dw,iw,temp)
      call PRODCT(I0,dw,idw,3,3,1)
      call CROSS(w,idw,temq)
      do i = 1, 3
         temr(i) = ndot(i) - temp(i) - temq(i)
      end do
      call PRODCT(I0i,temr,wdd,3,3,1)
      wdotwj = (DOT(w,wdd) + DOT(dw,dw))*2._10/3._10
      do i=1,3
         wddn(i)=I0i(i,i)*ndot(i)/emass
      end do
      wdotwjm=DOT(w,wddn)*2._10/3._10
 
c calculate perturbations to inertia tensor & time derivs.
c also save partials of these wrt earth mass for later use in calculating
c partials of the angular acceleration
c note that the angular jerk has a component proportional to the earth
c mass and so gives a mass dependence of the rotational part of the
c second time derivative
      do i = 1, 3
         do j = 1, i
            rij=Smecor(i)*Smecor(j)
            mfact=(rij-Rem2*del(i,j)/3._10)/Rem5
            didk(i,j)=mrad5*((w(i)*w(j)+(plas(i)-wdotw)*del(i,j))/g3
     .                - emass*mfact)
            didkm(i,j)=-mrad5*mfact

            rijd = velspa(i)*Smecor(j) + Smecor(i)*velspa(j)
            mfact=(rijd+del(i,j)*rdotv-5._10*rij*rfct1)/Rem5
            diddk(i,j)=mrad5*((w(i)*dw(j)+dw(i)*w(j)-del(i,j)*wdotwd)/g3
     .                 -emass*mfact)
            diddkm(i,j)=-mrad5*mfact
            mfact=(accspa(i)*Smecor(j)+2._10*velspa(i)*velspa(j)
     .       +Smecor(i)*accspa(j)+del(i,j)*rdotvd
     .       -5._10*((2._10*rijd+del(i,j)*rdotv)*rfct1+rij*rfct3)
     .       + 35._10*rij*rfct2)/Rem5
            didddk(i,j)=mrad5*((wdd(i)*w(j)+2._10*dw(i)*dw(j)
     .       +w(i)*wdd(j)-del(i,j)*wdotwj)/g3 - emass*mfact)
            didddkm(i,j)=mrad5*((wddn(i)*w(j)+w(i)*wddn(j)
     .       - del(i,j)*wdotwjm)/g3 - mfact)
            if(i.ne.j) then
               didk(j,i)   = didk(i,j)
               diddk(j,i)  = diddk(i,j)
               didddk(j,i) = didddk(i,j)
            endif
         end do
      end do

      do i = 1, 3
         do j = 1, 3
            didt(i,j)  = -diddk(i,j)
            diddt(i,j) = -didddk(i,j)
         end do
      end do
      return

      entry MOREDT2(dw,k,t,didy,diddy,iparm,didx,diddx)
c
c calculate sensitivity of delta I and delta I dot to change in the
c rotation state, to be used in calculating the partial derivatives
c of the angular acceleration, and set up for partials wrt parameters
c if iparm.gt.1 then also calculate sensitivity to orbit state


c dw = angular acc, as supplied in the main entry
c k  = love number
c t  = scaled time lag
c didy, diddy = 3x3x6 arrays, partials of delta I and delta I dot
c     wrt rotation state vector
c additional information passed through common:
c dwdy = partials of w wrt rotation state
c ddwdy = partials (to leading order) of dw wrt rotation state
c dsxdy = partials of selenodetic earth coords wrt euler angles
c drtdy = partials of selenodetic transformation matrix wrt euler angles
c dwddx = partials (to leading order) of dw wrt orbit state

c set up cross products as matrix multiplications, making use of the
c fact that I0 is diagonal
      dwi(1,1)=0._10
      dwi(2,1)=dw(3)*I0(1,1)
      dwi(3,1)=-dw(2)*I0(1,1)
      dwi(1,2)=-dw(3)*I0(2,2)
      dwi(2,2)=0._10
      dwi(3,2)=dw(1)*I0(2,2)
      dwi(1,3)=dw(2)*I0(3,3)
      dwi(2,3)=-dw(1)*I0(3,3)
      dwi(3,3)=0._10
      wi(1,1)=0._10
      wi(2,1)=w(3)*I0(1,1)
      wi(3,1)=-w(2)*I0(1,1)
      wi(1,2)=-w(3)*I0(2,2)
      wi(2,2)=0._10
      wi(3,2)=w(1)*I0(2,2)
      wi(1,3)=w(2)*I0(3,3)
      wi(2,3)=-w(1)*I0(3,3)
      wi(3,3)=0._10

      do l=1,6
         wddw=DOT(w,Dwdy(1,l))*(2._10/3._10)
         wddwd=(DOT(dw,Dwdy(1,l))+DOT(w,ddwdy(1,l)))*(2._10/3._10)
         call CROSS(Dwdy(1,l),Smecor,temr)
         do i=1,3
            dsvdy(i)=-temr(i)
            tems(i)=-2._10*velspa(i)-wcsx(i)
         end do
         if(l.le.3) then
            call PRODCT(Drtdy(1,1,l),Mecor(4),temp,-3,3,1)
            call CROSS(w,Dsxdy(1,l),temq)
            do i=1,3
               dsvdy(i)=dsvdy(i)-temp(i)-temq(i)
               temr(i)=temr(i)+temq(i)+2._10*temp(i)
            end do
         endif

         call CROSS(Dwdy(1,l),tems,temp)
         call CROSS(w,temr,tems)
         call CROSS(ddwdy(1,l),Smecor,temq)
         do i=1,3
            dsady(i)=temp(i)+tems(i)-temq(i)
         end do
         if(l.le.3) then
            call PRODCT(Drtdy(1,1,l),acc195,temr,-3,3,1)
            call CROSS(dw,Dsxdy(1,l),temq)
            do i=1,3
               dsady(i)=dsady(i)+temr(i)-temq(i)
            end do
         endif
c get partial of ndot
         call CROSS(dsvdy,ir,dnddy)
         if(l.le.3) then
            call PRODCT(I0,Dsxdy(1,l),temp,3,3,1)
            call CROSS(velspa,temp,temq)
            call CROSS(Smecor,temp,temr)
            call CROSS(Dsxdy(1,l),ir,tems)
            call CROSS(Dsxdy(1,l),iv,temp)
            do i=1,3
               dnddy(i)=dnddy(i)+temq(i)+temp(i)
     .          -5._10*(temr(i)+tems(i))*rfct1
            end do
         endif
         call PRODCT(I0,dsvdy,temp,3,3,1)
         call CROSS(Smecor,temp,tems)
         do i=1,3
            dnddy(i)=(dnddy(i)+tems(i))*3._10*Gamat3/Rem5
         end do
c get partial of angular jerk
         call PRODCT(dwi,Dwdy(1,l),temp,3,3,1)
         call CROSS(ddwdy(1,l),iw,temq)
         call PRODCT(wi,ddwdy(1,l),temr,3,3,1)
         call CROSS(Dwdy(1,l),idw,tems)
         do i=1,3
            temp(i)=dnddy(i)-temp(i)-temq(i)-temr(i)-tems(i)
         end do
         call PRODCT(I0i,temp,dddwdy,3,3,1)
         wddwdd=(2._10*DOT(dw,ddwdy(1,l))+DOT(w,dddwdy)
     .    +DOT(Dwdy(1,l),wdd))*(2._10/3._10)
         do i=1,3
            do j=1,i
               term0=Dwdy(i,l)*w(j)+w(i)*Dwdy(j,l)
               term1=ddwdy(i,l)*w(j)+dw(i)*Dwdy(j,l)+Dwdy(i,l)*dw(j)
     .          +w(i)*ddwdy(j,l)
               term2=dddwdy(i)*w(j)+wdd(i)*Dwdy(j,l)
     .          +2._10*(ddwdy(i,l)*dw(j)+dw(i)*ddwdy(j,l))
     .          +Dwdy(i,l)*wdd(j)+w(i)*dddwdy(j)
               if(i.eq.j) then
                  term0=term0-wddw
                  term1=term1-wddwd
                  term2=term2-wddwdd
               endif
               term0=term0/g3
               term1=term1/g3
               term2=term2/g3
               rijyd=dsvdy(i)*Smecor(j)+Smecor(i)*dsvdy(j)
               if(l.le.3) then
                  rijyd=rijyd+Dsxdy(i,l)*velspa(j)+velspa(i)*Dsxdy(j,l)
                  rijy=Dsxdy(i,l)*Smecor(j)+Smecor(i)*Dsxdy(j,l)
               endif
               term1=term1-emass*rijyd/Rem5
               term2=term2-emass*(dsady(i)*Smecor(j)+2._10*(dsvdy(i)
     .          *velspa(j)+velspa(i)*dsvdy(j))+Smecor(i)*dsady(j)
     .          -10._10*rijyd*rfct1)/Rem5
               if(l.le.3) then
                  termr=rijy/Rem5
                  term0=term0-emass*termr
                  term1=term1+emass*5._10*termr*rfct1
                  term2=term2-emass*((Dsxdy(i,l)*accspa(j)
     .             +accspa(i)*Dsxdy(j,l))/Rem5
     .             +termr*(35._10*rfct2-5._10*rfct3))
               endif
               term=(k*term0-t*term1)*mrad5
               didy(i,j,l)=term
               didy(j,i,l)=term
               term=(k*term1-t*term2)*mrad5
               diddy(i,j,l)=term
               diddy(j,i,l)=term
            end do
         end do
      end do

c set up for partials with respect to parameters

c dependence of angular jerk on moment-of-inertia ratios
c (it would be more efficient to calculate the jerk this way, by
c  multiplying these components respectively by alpha,beta,gamma --
c  but the intermediate quantities in the more laborious matrix
c  and cross product method are saved and reused in other
c  calculations)
      dwdddi(1)=g3r5*(Smecor(2)*velspa(3)+Smecor(3)*velspa(2)
     . - 5._10*rfct1*Smecor(2)*Smecor(3)) - (W2*dw(3)+W3*dw(2))
      dwdddi(2)=-g3r5*(Smecor(3)*velspa(1)+Smecor(1)*velspa(3)
     . - 5._10*rfct1*Smecor(3)*Smecor(1)) + (W3*dw(1)+W1*dw(3))
      dwdddi(3)=g3r5*(Smecor(1)*velspa(2)+Smecor(2)*velspa(1)
     . - 5._10*rfct1*Smecor(1)*Smecor(2)) - (W1*dw(2)+W2*dw(1))

      if(iparm.le.1) return
c compute sensitivity to orbit state
      dwcr(1,1)=0._10
      dwcr(2,1)=dw(3)
      dwcr(3,1)=-dw(2)
      dwcr(1,2)=-dw(3)
      dwcr(2,2)=0._10
      dwcr(3,2)=dw(1)
      dwcr(1,3)=dw(2)
      dwcr(2,3)=-dw(1)
      dwcr(3,3)=0._10
      wcr(1,1)=0._10
      wcr(2,1)=w(3)
      wcr(3,1)=-w(2)
      wcr(1,2)=-w(3)
      wcr(2,2)=0._10
      wcr(3,2)=w(1)
      wcr(1,3)=w(2)
      wcr(2,3)=-w(1)
      wcr(3,3)=0._10
      call PRODCT(dwcr,Mrotlb,wdr,3,3,3)
      call PRODCT(wcr,Mrotlb,wr,3,3,3)
      call PRODCT(wcr,wr,dddrdx,3,3,3)
      call PRODCT(Mrotlb,Ecor,x3s,3,3,1)
      call PRODCT(Mrotlb,Mcor,x10s,3,3,1)
      do l=1,3
         call CROSS(Ddwdx(1,l),Smecor,temp)
         do i=1,3
            dddrdx(i,l)=-dddrdx(i,l)+wdr(i,l)-temp(i)+
     .       Gamtem*(Mrotlb(i,l)+3._10*Smecor(i)*Mecor(l)/Rem2)/Rem3
     .       -Gamat*(Mass(10)*(3._10*x3s(i)*Ecor(l)/Re2-Mrotlb(i,l))/Re3
     .       +Masse*(3._10*x10s(i)*Mcor(l)/Rm2-Mrotlb(i,l))/Rm3)
c     .       -Gamat*(3._10*x3s(i)*Xpert(l,3)/Rpert2(3)-Mrotlb(i,l))
c     .       /Rpert3(3)
         end do
      end do
      do l=1,6
c get partial of rdotv and rdotv dot
         if(l.le.3) then
            drdv=Mecor(l+3)
            drdvd= Gamtem*Mecor(l)/Rem3
     .                - 2._10*Gamat*(Mcor(l)/Rm3 - Ecor(l)/Re3)
         else
            drdv=Mecor(l-3)
            drdvd=2._10*Mecor(l)
         endif
c get partial of ndot
         if(l.le.3) then
            call CROSS(wr(1,l),ir,dnddx)
            call PRODCT(I0,Mrotlb(1,l),temp,3,3,1)
            call CROSS(temp,velspa,temq)
            call CROSS(temp,Smecor,temr)
            call CROSS(Mrotlb(1,l),ir,tems)
            call CROSS(Mrotlb(1,l),iv,temp)
            do i=1,3
               dnddx(i)=dnddx(i)+temq(i)-temp(i)
     .          -5._10*(temr(i)-tems(i)-2._10*rir(i)*Mecor(l)/Rem2)
     .          *rfct1
            end do
            call PRODCT(I0,wr(1,l),temp,3,3,1)
         else
            call CROSS(Mrotlb(1,l-3),ir,dnddx)
            call PRODCT(I0,Mrotlb(1,l-3),temp,3,3,1)
            do i=1,3
               dnddx(i)=-dnddx(i)
               temp(i)=-temp(i)
            end do
         endif
         call CROSS(Smecor,temp,tems)
         do i=1,3
            dnddx(i)=(dnddx(i)+tems(i)-5._10*rir(i)*drdv/Rem2)*3._10
     .       *Gamat3/Rem5
            if(l.le.3) dnddx(i)=dnddx(i)-5._10*ndot(i)*Mecor(l)/Rem2
         end do
c get partial of angular jerk
         call CROSS(Ddwdx(1,l),iw,temq)
         call PRODCT(wi,Ddwdx(1,l),temr,3,3,1)
         do i=1,3
            temp(i)=dnddx(i)-temq(i)-temr(i)
         end do
         call PRODCT(I0i,temp,dddwdx,3,3,1)
         wddwdd=(2._10*DOT(dw,Ddwdx(1,l))+DOT(w,dddwdx))*(2._10/3._10)
         wddwd=DOT(w,Ddwdx(1,l))*(2._10/3._10)         
         do i=1,3
            do j=1,i
               term1=(Ddwdx(i,l)*w(j)+w(i)*Ddwdx(j,l)-del(i,j)*wddwd)/g3
               term2=(dddwdx(i)*w(j)+2._10*(Ddwdx(i,l)*dw(j)
     .          +dw(i)*Ddwdx(j,l))+w(i)*dddwdx(j)-del(i,j)*wddwdd)/g3
               rij=Smecor(i)*Smecor(j)
               rijd=velspa(i)*Smecor(j)+Smecor(i)*velspa(j)
               if(l.le.3) then
                  rijx=-Mrotlb(i,l)*Smecor(j)-Smecor(i)*Mrotlb(j,l)+
     .             del(i,j)*Mecor(l)
                  rijxd=wr(i,l)*Smecor(j)-velspa(i)*Mrotlb(j,l)
     .             -Mrotlb(i,l)*velspa(j)+Smecor(i)*wr(j,l)+
     .             del(i,j)*drdv
                  term0=-(rijx-5._10*Mecor(l)*rij/Rem2)*emass/Rem5
                  term1=term1-(rijxd-5._10*(rijx*rfct1+rij*drdv/Rem2)
     .             - 5._10*Mecor(l)*(rijd-7._10*rij*rfct1)/Rem2)
     .             *emass/Rem5
                  term2=term2-(dddrdx(i,l)*Smecor(j)
     .             -accspa(i)*Mrotlb(j,l)+2._10*(wr(i,l)*velspa(j)
     .             +velspa(i)*wr(j,l)) - Mrotlb(i,l)*accspa(j)
     .             +Smecor(i)*dddrdx(j,l) + del(i,j)*drdvd
     .             - 10._10*rijxd*rfct1
     .             - 5._10*rijx*rfct3 + 35._10*rijx*rfct2
     .             + (-10._10*rijd*drdv
     .             + Smecor(i)*Smecor(j)*(-5._10*drdvd+70._10*rfct1*drdv
     .             + (35._10*rfct3-315._10*rfct2)*Mecor(l))
     .             + (-5._10*(accspa(i)*Smecor(j)
     .             + 2._10*velspa(i)*velspa(j)
     .             + Smecor(i)*accspa(j))
     .             + 70._10*rijd*rfct1)*Mecor(l))/Rem2)*emass/Rem5

               else
                  rijxd=-Mrotlb(i,l-3)*Smecor(j)-Smecor(i)*Mrotlb(j,l-3)
     .             +del(i,j)*drdv
                  term0=0._10
                  term1=term1-(rijxd-5._10*rij*drdv/Rem2)*emass/Rem5
                  term2=term2-(2._10*(wr(i,l-3)*Smecor(j)
     .             -Mrotlb(i,l-3)*velspa(j)-velspa(i)*Mrotlb(j,l-3)
     .             +Smecor(i)*wr(j,l-3)) + del(i,j)*drdvd
     .             - 10._10*rijxd*rfct1 + (-10._10*rijd*drdv
     .             +rij*(-5._10*drdvd+70._10*rfct1*drdv))/Rem2)
     .             *emass/Rem5
               endif
               term=(k*term0-t*term1)*mrad5
               didx(i,j,l)=term
               didx(j,i,l)=term
               term=(k*term1-t*term2)*mrad5
               diddx(i,j,l)=term
               diddx(j,i,l)=term
            end do
         end do
      end do
      return

      entry MOREDTP(dw,k,t,ddidp,ddiddp)
c Calculate partials of elasticity/dissipation corrections wrt parameters
c  Input:
c dw = angular acc, as supplied in the main entry, i.e., before correction
c k  = love number
c t  = scaled time lag
c  Other information in common:
c Ddwdp = partial of rigid-body angular acceleration, not yet corrected
c W1,W2,W3= angular velocity
c I0,I0i = mean inertia tensor and its inverse
c alpha,beta,gamma = moment-of-inertia ratios
c Dalpha,Dbeta,Dgamma = partials of moment-of-inertia ratios wrt all parms
c  Output:
c ddidp  = partial of inertia tensor correction
c ddiddp = partial of inertia tensor rate correction

c partial of angular jerk, excluding direct dependence of Ndot on mass
      dwdddp(1)=Dalpha(Kkk)*dwdddi(1) - alpha*(W2*Ddwdp(3)+W3*Ddwdp(2))
      dwdddp(2)=Dbeta(Kkk)*dwdddi(2) + beta*(W3*Ddwdp(1)+W1*Ddwdp(3))
      dwdddp(3)=Dgamma(Kkk)*dwdddi(3) - gamma*(W1*Ddwdp(2)+W2*Ddwdp(1))
c The only parameter dependence of the gravitational acceleration is upon
c Mass(3), aside from a second-order dependence upon Mass(10), which we
c ignore here.  This dependence is purely radial, while the Coriolis
c acceleration dependence is purely tangential.
      call CROSS(Smecor,Ddwdp,daccdp)
      if(Icrtrl(Kkk).eq.3) then
         do i=1,3
            daccdp(i)=daccdp(i)+daccdm(i)
         end do
         rdaccdp=rdaccdm
      else
         rdaccdp=0._10
      endif

c calculate perturbation dependence, excluding direct mass effect
      wdotwd = DOT(w,Ddwdp)*2._10/3._10
      wdotwj = (DOT(w,dwdddp) + 2._10*DOT(Ddwdp,dw))*2._10/3._10
      do i = 1, 3
         do j = 1, i
            rij=Smecor(i)*Smecor(j)
            diddpk(i,j)=mrad5*((w(i)*Ddwdp(j)+Ddwdp(i)*w(j)
     .       -del(i,j)*wdotwd)/g3)
            mfact=(daccdp(i)*Smecor(j)+Smecor(i)*daccdp(j)
     .            +del(i,j)*rdaccdp-5._10*rij*rdaccdp/Rem2)/Rem5
            didddpk(i,j)=mrad5*((dwdddp(i)*w(j)
     .       +2._10*(Ddwdp(i)*dw(j)+dw(i)*Ddwdp(j))+w(i)*dwdddp(j)
     .       -del(i,j)*wdotwj)/g3 - emass*mfact)
         end do
      end do

c add direct mass dependence if any
      if(Icrtrl(Kkk).eq.3 .or. Icrtrl(Kkk).eq.10) then
         dmdp=Masse
         if(Icrtrl(Kkk).eq.10) dmdp=-Mass(3)
         do i=1,3
            do j=1,i
               didpk(i,j)=didkm(i,j)*dmdp
               diddpk(i,j)=diddpk(i,j)+diddkm(i,j)*dmdp
               didddpk(i,j)=didddpk(i,j)+didddkm(i,j)*dmdp
            end do
         end do
      else
         do i=1,3
            do j=1,i
               didpk(i,j)=0._10
            end do
         end do
      endif
c combine time series into perturbation and rate
      do i=1,3
         do j=1,i
            ddidp(i,j) = k*didpk(i,j) - t*diddpk(i,j)
            ddiddp(i,j)= k*diddpk(i,j) - t*didddpk(i,j)
            if(i.ne.j) then
               ddidp(j,i)  = ddidp(i,j)
               ddiddp(j,i) = ddiddp(i,j)
            endif
         end do
      end do

      return
      end
