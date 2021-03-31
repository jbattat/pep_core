      subroutine MORED(s,dw,kmr83)
 
      implicit none
c
c rj cappallo   january 1979   sr mored, entry moredp
c
c     mored calculates and applys corrections to the lunar rotational
c     equations of motion to model elasticity and dissipation. entry
c     moredp calculates the rhs of the variational eqs. wrt:
c             k2(moon) = mrcond(12), and
c             k2 * lunar elastic response lag = mrcond(13).
c     see rjc's thesis for relevant formulae.

c        calling parameters
      real*10 s, dw(3)
      integer*2 kmr83,iparm

c array dimensions
      include 'globdefs.inc'

c commons
      include 'empcnd.inc'
      real*10 k2, reslag
      equivalence (Mrcond(12), k2), (Mrcond(13), reslag)
      include 'harmor.inc'
      real*10 mmabc(3),dmmabc(i_mxplprt,3)
      equivalence (mmabc, Mma),(dmmabc,Dmma)
      include 'intstf.inc'
      include 'mnrtlb.inc'
      real*10 w(3)
      equivalence (w, W1)
      include 'morcrd.inc'
      include 'morstf.inc'
      include 'param.inc'

c external functions
      real*10 DOT
c
c                  l o c a l   v a r i a b l e s
      real*10 deli(3,3),delid(3,3),deln(3),temp(3),temq(3),temr(3),
     . ddidk2(3,3), ddiidk(3,3), ddiddk(3,3), ddidp(3,3),
     . ddiidp(3,3), ddiidl(3,3), ddiddp(3,3), ddidl(3,3),
     . ddiddl(3,3),deli0ii(3,3),didy(3,3,6),diddy(3,3,6),
     . didx(3,3,6),diddx(3,3,6),ddndx(3),mfact
      real*10 fact,factj,diw(3),delir(3),ddndy(3),dwin(3),delnk(3)
      integer i,j
      real*10 deli0i(3,3)

c save unperturbed angular acceleration for use in calculations
      do i=1,3
         dwin(i)=dw(i)
      end do
c
c      call subroutines to calculate perturbations to inertia tensor
c     and partials wrt k2 and a dissipation parameter.
c     MOREDT does calculation for kt (constant-time-lag, kmr83=0 or 1)
c     MOREDQ does calculation for k/q (constant-q model, kmr83=2)
c     if partial derivatives are needed, these subroutines are called
c     again to return the partials wrt the rotational state

      if(kmr83.le.1) call MOREDT(dwin,ddidk2,ddidl,ddiddk,ddiddl)
      if(kmr83.eq.2) call MOREDQ(s,I0(3,3),ddidk2,ddidl,ddiddk,ddiddl)

      do i = 1,3
         do j = 1,3
 
c form deli and delid, corrected for dissipation
            deli(i,j)  = k2*ddidk2(i,j) + reslag*ddidl(i,j)
            delid(i,j) = k2*ddiddk(i,j) + reslag*ddiddl(i,j)
 
c then calculate perturbation to i0 inverse (and also multiply this
c perturbation by i0)
            factj = I0i(i,i)
            fact = factj*I0i(j,j)
            deli0ii(i,j) = -deli(i,j)*factj
            deli0i(i,j) = -deli(i,j)*fact
 
c also compute partial of i0 inverse perturbation wrt k2
            ddiidk(i,j) = -ddidk2(i,j)*fact
 
c and wrt tidal response lag
            ddiidl(i,j) = -ddidl(i,j)*fact
         end do
      end do
c
c find perturbation to torque from new inertia tensor
      call PRODCT(deli, Smecor, delir, 3, 3, 1)
      call CROSS(Smecor, delir, deln)
      do i = 1, 3
         deln(i) = deln(i)*3._10*Gamat3/Rem5
      end do

c now apply effect of perturbations to angular accelerations
      call PRODCT(deli0ii, dw, temp, 3, 3, 1)
      call PRODCT(deli, w, diw, 3, 3, 1)
      call CROSS(w, diw, temr)
      call PRODCT(delid, w, temq, 3, 3, 1)
      do i = 1, 3
         delnk(i) = deln(i) - temr(i) - temq(i)
      end do
      call PRODCT(I0i, delnk, temr, 3, 3, 1)
      do i = 1, 3
         dw(i) = dw(i) + temp(i) + temr(i)
      end do

      return
c
c              v a r i a t i o n a l   e q u a t i o n s
c calculate corrections to dependence of angular acceleration on state
      entry MORED2(s,kmr83,iparm)
      if(kmr83.le.1) call MOREDT2(dwin,k2,reslag,didy,diddy,
     . iparm,didx,diddx)
      if(kmr83.eq.2) call MOREDQ2(s,I0(3,3),k2,reslag,didy,diddy,
     . iparm,didx,diddx)

      do i=1,6
         call PRODCT(deli0ii,ddwdy(1,i),temp,3,3,1)
         do j=1,3
            ddwdy(j,i)=ddwdy(j,i)+temp(j)
         end do
         call PRODCT(didy(1,1,i),dwin,temq,3,3,1)
         do j=1,3
            ddwdy(j,i)=ddwdy(j,i)-temq(j)/I0(j,j)
         end do
         call PRODCT(didy(1,1,i),Smecor,temp, 3,3,1)
         if(i.le.3) then
            call PRODCT(deli,Dsxdy(1,i),temq,3,3,1)
            do j=1,3
               temp(j)=temp(j)+temq(j)
            end do
         endif
         call CROSS(Smecor, temp, ddndy)
         if(i.le.3) then
            call CROSS(Dsxdy(1,i),delir,temq)
            do j=1,3
               ddndy(j)=ddndy(j)+temq(j)
            end do
         endif
         do j = 1, 3
            ddndy(j) = ddndy(j)*3._10*Gamat3/Rem5
         end do
         call CROSS(Dwdy(1,i),diw,temp)
         call PRODCT(didy(1,1,i),w,temq,3,3,1)
         call PRODCT(deli,Dwdy(1,i),temr,3,3,1)
         do j=1,3
            temq(j)=temq(j)+temr(j)
         end do
         call CROSS(w,temq,temr)
         do j=1,3
            temp(j)=temp(j)+temr(j)
         end do
         call PRODCT(delid,Dwdy(1,i),temq,3,3,1)
         call PRODCT(diddy(1,1,i),w,temr,3,3,1)
         do j=1,3
            temp(j)=ddndy(j)-(temp(j)+temq(j)+temr(j))
         end do
         call PRODCT(I0i,temp,temq,3,3,1)
         do j=1,3
            ddwdy(j,i)=ddwdy(j,i)+temq(j)
         end do
      end do

      if(iparm.le.1) return
c calculate sensitivity of angular acceleration to orbit state
      do i=1,6
         if(i.le.3) then
            call PRODCT(deli0ii,Ddwdx(1,i),temp,3,3,1)
            do j=1,3
               Ddwdx(j,i)=Ddwdx(j,i)+temp(j)
            end do
         endif
         call PRODCT(didx(1,1,i),dwin,temq,3,3,1)
         do j=1,3
            Ddwdx(j,i)=Ddwdx(j,i)-temq(j)/I0(j,j)
         end do
         call PRODCT(didx(1,1,i),Smecor,temp, 3,3,1)
         if(i.le.3) then
            call PRODCT(deli,Mrotlb(1,i),temq,3,3,1)
            do j=1,3
               temp(j)=temp(j)-temq(j)
            end do
         endif
         call CROSS(Smecor, temp, ddndx)
         if(i.le.3) then
            call CROSS(Mrotlb(1,i),delir,temq)
            do j=1,3
               ddndx(j)=ddndx(j)-temq(j)
            end do
         endif
         do j = 1, 3
            ddndx(j) = ddndx(j)*3._10*Gamat3/Rem5
         end do
         if(i.le.3) then
            do j=1,3
               ddndx(j)=ddndx(j)-5._10*deln(j)*Mecor(i)/Rem2
            end do
         endif
         call PRODCT(didx(1,1,i),w,temq,3,3,1)
         call CROSS(w,temq,temp)
         call PRODCT(diddx(1,1,i),w,temr,3,3,1)
         do j=1,3
            temp(j)=ddndx(j)-(temp(j)+temr(j))
         end do
         call PRODCT(I0i,temp,temq,3,3,1)
         do j=1,3
            Ddwdx(j,i)=Ddwdx(j,i)+temq(j)
         end do
      end do

      return
      
c
c calculate direct partial derivative wrt Love number or lag
      entry MOREDP
      if(Icrtrl(Kkk).eq.-7) then
c
c time lag of lunar solid body tidal and rot'l. deformation
c or partial wrt k2/q
         do i = 1, 3
            do j = 1, 3
               ddidp(i, j)  = ddidl(i, j)
               ddiidp(i, j) = ddiidl(i, j)
               ddiddp(i, j) = ddiddl(i, j)
            end do
         end do
      else
 
c lunar k2 (by process of elimination)
         do i = 1, 3
            do j = 1, 3
               ddidp(i, j)  = ddidk2(i, j)
               ddiidp(i, j) = ddiidk(i, j)
               ddiddp(i, j) = ddiddk(i, j)
            end do
         end do
      endif
c
c form rhs of var'l eqns for parameter p(=k2 or reslag if entered at MOREDP,
c or any other parameter if entered at MOREDPC)
  200 call PRODCT(ddiidp, N0mkin, temp, 3, 3, 1)
      do i=1,3
         Ddwdp(i)=Ddwdp(i)+temp(i)
      end do

      call PRODCT(ddidp, Smecor, temq, 3, 3, 1)
      call CROSS(Smecor, temq, temr)
      call PRODCT(ddidp, w, temp, 3, 3, 1)
      call CROSS(w, temp, temq)
      call PRODCT(ddiddp, w, temp, 3, 3, 1)
      do i = 1, 3
         temr(i) = temr(i)*3._10*Gamat3/Rem5 - temq(i) - temp(i)
      end do
      if(Icrtrl(Kkk).eq.3 .or. Icrtrl(Kkk).eq.10) then
         if(Icrtrl(Kkk).eq.3) then
            mfact=1._10/Mass(3)
         else
            mfact=-1._10/Masse
         endif
         do i=1,3
            temr(i)=temr(i)+deln(i)*mfact
         end do
      endif
      call PRODCT(I0i, temr, temp, 3, 3, 1)
      do i = 1, 3
         Ddwdp(i) = Ddwdp(i) + temp(i)
      end do
      return

c calculate effect of elasticity on partial derivatives

      entry MOREDPC(kmr83)

c get dependence of perturbations on parameters
      if(kmr83.le.1) call MOREDTP(dwin,k2,reslag,ddidp,ddiddp)
      if(kmr83.eq.2) call MOREDQP(ddidp,ddiddp)

c first calculate ddiidp treating I0 as a constant -
      do i=1,3
         do j=1,3
            ddiidp(i,j)= -ddidp(i,j)*I0i(i,i)*I0i(j,j)
         end do
      end do
c now add the dependence of I0
c also calculate I0^-1 dI0/dp term for simple correction
      if(Icrtrl(Kkk).eq.3) then
         do i=1,3
            do j=1,3
               ddiidp(i,j)=ddiidp(i,j)-deli0i(i,j)*2._10/Mass(3)
            end do
            temq(i)=Ddwdp(i)+dwin(i)/Mass(3)
         end do
      else if(Icrtrl(Kkk).eq.10) then
         do i=1,3
            do j=1,3
               ddiidp(i,j)=ddiidp(i,j)-deli0i(i,j)*2._10/Mass(10)
            end do
            temq(i)=Ddwdp(i)+dwin(i)/Mass(10)
         end do
      else
         do i=1,3
            mfact=dmmabc(Kkk,i)/mmabc(i)
            do j=1,3
               ddiidp(i,j)=ddiidp(i,j)+deli0i(i,j)*(mfact
     .          +dmmabc(Kkk,j)/mmabc(j))
            end do
            temq(i)=Ddwdp(i)-dwin(i)*mfact
         end do
      endif

      call PRODCT(deli0ii,temq,temp,3,3,1)
      do i=1,3
         Ddwdp(i)=Ddwdp(i) + temp(i) + delnk(i)*dmmabc(Kkk,i)/Mmoon
      end do

c add contribution of dependence of elastic perturbation on parms

      if(Icrtrl(Kkk).eq.3) then
c
c partial wrt embary mass
         do i = 1,3
            Ddwdp(i)=Ddwdp(i)-delnk(i)*I0i(i,i)/Mass(3)
         end do
      else if(Icrtrl(Kkk).eq.10) then
c
c partial wrt Moon mass fraction
         do i = 1,3
            Ddwdp(i)=Ddwdp(i)-delnk(i)*I0i(i,i)/Mass(10)
         end do
      endif
      goto 200

      end
