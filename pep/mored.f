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
      real*10 s,dw(3),fno
      integer*2 kmr83,iparmc

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
      include 'metuna.inc'
      include 'mnrtlb.inc'
      real*10 w(3)
      equivalence (w, W1)
      include 'morcrd.inc'
      include 'morstf.inc'
      include 'param.inc'

c external functions
      real*10 DOT,DOTN
c
c                  l o c a l   v a r i a b l e s
      real*10 deli(3,3),delid(3,3),deln(3),temp(3),temq(3),temr(3),
     . ddidk2(3,3), ddiidk(3,3), ddiddk(3,3), ddidp(3,3),
     . ddiidp(3,3), ddiidl(3,3), ddiddp(3,3), ddidl(3,3),
     . ddiddl(3,3),deli0ii(3,3),didy(3,3,6),diddy(3,3,6),
     . didx(3,3,6),diddx(3,3,6),ddndx(3),mfact,delirr(6),dadyed(6,3),
     . dadxed(6,3),didxrr,dadxed1(6,3),tempfn(3,6),tempda(6,3),
     . ddidkm(3,3),ddidlm(3,3),ddidm(3,3),didmr(3),dfndm(3),
     . dtridx,masfct,dtridm,deltri,dadxedt(6,3),delnp(3,6),gamorb,
     . dadxedp(6,3,6),ddid2l(3,3),ddidd2l(3,3),ddid2lm(3,3),t2k

      real*10 fact,factj,diw(3),delir(3,6),ddndy(3),dwin(3),delnk(3)
      integer i,is,j,k,kkkrot
      real*10 deli0i(3,3),dit1(3,3),dit2(3,3),didyt(3,3,6),dd,tems(3),
     . ddndyp(3),ddndxp(3),delns(3)

      t2k=0._10
ct      if(k2.ne.0._10 .and. kmr83.le.1) t2k=0.5_10*reslag*reslag/k2
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

      if(kmr83.le.1) call MOREDT(dwin,ddidk2,ddidl,ddid2l,ddiddk,ddiddl,
     . ddidd2l,ddidkm,ddidlm,ddid2lm)
      if(kmr83.eq.2) call MOREDQ(s,I0(3,3),ddidk2,ddidl,ddid2l,ddiddk,
     . ddiddl,ddidd2l,ddidkm,ddidlm,ddid2lm)

      do i = 1,3
         do j = 1,3
 
c form deli and delid, corrected for dissipation
            deli(i,j)  = k2*ddidk2(i,j) + reslag*ddidl(i,j)
ct     .       + t2k*ddid2l(i,j)
            delid(i,j) = k2*ddiddk(i,j) + reslag*ddiddl(i,j)
ct     .       + t2k*ddidd2l(i,j)

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
      call PRODCT(deli,Smecor,delir,3,3,1)
      call CROSS(Smecor,delir,deln)
      do i = 1,3
         delnp(i,1)= deln(i)*3._10*Gamat3/Rem5
         deln(i)   = delnp(i,1)
      end do

c include the sun and planet effects if so indicated
      if(Kmr(81).ge.0) then
         call PRODCT(deli,Smcor,delir(1,2),3,3,1)
         call CROSS(Smcor,delir(1,2),delns)
         do i = 1,3
            delnp(i,2)= delns(i)*3._10*Gamat/Rm5
            deln(i)   = deln(i) + delnp(i,2)
         end do
      endif
      if(Kmr(84).ge.0) then
         do is=3,6
            if(Npmhar(is).gt.0) then
               call PRODCT(deli,Spcor(1,is),delir(1,is),3,3,1)
               call CROSS(Spcor(1,is),delir(1,is),delns)
               do i = 1,3
                  delnp(i,is)=delns(i)*3._10*Gamtrq(is)/Spc5(is)
                  deln(i)    = deln(i) + delnp(i,is)
               end do
            endif
         end do
      endif

c calculate and save orbit perturbation from deformed inertia tensor
c and also set up for partials w.r.t. masses
      if(Orbint) then
c first set up for mass partials
         do i=1,3
            do j=1,3
               ddidm(i,j) = k2*ddidkm(i,j) + reslag*ddidlm(i,j)
ct     .          + t2k*ddid2lm(i,j)
            end do
         end do
         dtridm = 1.5_10*(ddidm(1,1)+ddidm(2,2)+ddidm(3,3))

         deltri = 1.5_10*(deli(1,1)+deli(2,2)+deli(3,3))
         do i=1,3
            Delfn(i)=0._10
            dfndm(i)=0._10
         end do
c the relevant mass for torque due to earth is just the earth mass,
c but for the orbital effect it's earth+moon,
c since orbit is moon relative to earth
         gamorb=Gamtem
         do is=1,6
            if(is.eq.1 .or. (is.eq.2.and.Kmr(81).ge.0) .or.
     .       (is.gt.2.and.Kmr(84).ge.0.and.Npmhar(is).gt.0)) then
               if(is.gt.1) gamorb=Gamtrq(is)
               delirr(is) = 7.5_10*DOT(delir(1,is),Spcor(1,is))/Spc2(is)
               do i=1,3
                  tempfn(i,is)=3._10*delir(i,is)-(delirr(is)-deltri)*
     .             Spcor(i,is)
               end do
               call PRODCT(Mrotlb,tempfn(1,is),temp,-3,3,1)
               do i=1,3
                  Delfn(i)=Delfn(i)+temp(i)*gamorb/Mmoon/Spc5(is)
               end do
               call PRODCT(ddidm,Spcor(1,is),didmr,3,3,1)
               didxrr = 7.5_10*DOT(didmr,Spcor(1,is))/Spc2(is)
               do i=1,3
                  temp(i)=3._10*didmr(i)-(didxrr-dtridm)*Spcor(i,is)
               end do
               call PRODCT(Mrotlb,temp,temq,-3,3,1)
               do i=1,3
                  dfndm(i)=dfndm(i) + temq(i)*gamorb/Mmoon/Spc5(is)
               end do
            endif
         end do
      endif

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
      entry MORED2(s,kmr83,iparmc)
      if(iparmc.ne.Iparm) call SUICID('IPARM MISMATCH IN MORED2   ',7)
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
            ddndy(j)=0._10
         end do
         masfct=1._10
         do is=1,6
            if(is.eq.2) masfct=Masse
            if(is.eq.1 .or. (is.eq.2.and.Kmr(81).gt.0) .or.
     .       (is.gt.2.and.Kmr(84).gt.0.and.Npmhar(is).gt.0)) then
               call PRODCT(didy(1,1,i),Spcor(1,is),temp, 3,3,1)
               if(i.le.3) then
                  call PRODCT(deli,Dsxdy(1,i),temq,3,3,1)
                  do j=1,3
                     temp(j)=temp(j)+temq(j)*masfct
                  end do
               endif
               call CROSS(Spcor(1,is),temp,ddndyp)
               if(i.le.3) then
                  call CROSS(Dsxdy(1,i),delir(1,is),temq)
                  do j=1,3
                     ddndyp(j)=ddndyp(j)+temq(j)*masfct
                  end do
               endif
               do j = 1, 3
                  ddndy(j) = ddndy(j) +
     .             ddndyp(j)*3._10*Gamtrq(is)/Spc5(is)
               end do
            endif
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
c set up for dependence of orbit perturbation on rotation state
      do i=1,Index2
         call PRODCT(Mrotlb,didy(1,1,i),dit1,-3,3,3)
         call PRODCT(dit1,Mrotlb,didyt(1,1,i),3,3,3)
         if(i.le.3) then
            call PRODCT(Mrotlb,deli,dit1,-3,3,3)
            call PRODCT(dit1,Drtdy(1,1,i),dit2,3,-3,3)
            do j=1,3
               do k=1,j
                  dd=dit2(j,k)+dit2(k,j)
                  didyt(j,k,i)=didyt(j,k,i)+dd
                  if(k.lt.j) didyt(k,j,i)=didyt(k,j,i)+dd
               end do
            end do
         endif
      end do
c set up for partials of orbit wrt orbit state
      do i=1,3
         do j=1,3
            dadxed(i,j)=0._10
            dadxed(i+3,j)=0._10
            dadxedt(i+3,j)=0._10
            dadyed(i,j)=0._10
            dadyed(i+3,j)=0._10
         end do
      end do
      gamorb=Gamtem
      do is=1,6
         if(is.eq.1 .or. (is.eq.2.and.Kmr(81).gt.0) .or.
     .    (is.gt.2.and.Kmr(84).gt.0.and.Npmhar(is).gt.0)) then
            if(is.gt.1) gamorb=Gamtrq(is)
            do i=1,3
               do j=1,3
                  dadxedt(i,j)=3._10*deli(i,j) + ((2._10*delirr(is)*
     .             Spcor(i,is) - 15._10*delir(i,is))*Spcor(j,is)
     .             - 5._10*tempfn(j,is)*Spcor(i,is))/Spc2(is)
                  if(j.eq.i) dadxedt(i,j)=dadxedt(i,j)-delirr(is)+deltri
                  if(is.gt.1) dadxedt(i,j)=dadxedt(i,j)*Masse
               end do
            end do
            do i=1,Index2
               call PRODCT(didx(1,1,i),Spcor(1,is),temq,3,3,1)
               didxrr=7.5_10*DOT(temq,Spcor(1,is))/Spc2(is)
               dtridx=1.5_10*(didx(1,1,i)+didx(2,2,i)+didx(3,3,i))
               do j=1,3
                  dadxed1(i,j)=3._10*temq(j)-(didxrr-dtridx)*Spcor(j,is)
               end do
            end do
c since positions and velocities transform the same way, dadxed can be
c pictured as (3,2,3) and thus as (3,6) for premultiplication
c (the velocity dependence of the geometric term is zero anyway)
            call PRODCT(Mrotlb,dadxedt,tempda,-3,3,6)
            do i=1,Index2
               do j=1,3
                  tempda(i,j)=tempda(i,j)-dadxed1(i,j)
                  tempda(i,j)=-tempda(i,j)*gamorb/Mmoon/Spc5(is)
               end do
            end do
            call PRODCT(tempda,Mrotlb,dadxedp(1,1,is),6,3,3)
            do i=1,Index2
               do j=1,3
                  dadxed(i,j)=dadxed(i,j)+dadxedp(i,j,is)
               end do
            end do
c
c get dependence of orbit perturbation on rotation state
            do i=1,Index2
               call PRODCT(didyt(1,1,i),Mpcor(1,is),temq,3,3,1)
               didxrr=7.5_10*DOT(temq,Mpcor(1,is))/Spc2(is)
               dtridx=1.5_10*(didyt(1,1,i)+didyt(2,2,i)+didyt(3,3,i))
               do j=1,3
                  dadyed(i,j)=dadyed(i,j)+(3._10*temq(j)-
     .            (didxrr-dtridx)*Mpcor(j,is))*gamorb/Mmoon/Spc5(is)
               end do
            end do
         endif
      end do

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
            ddndx(j)=0._10
         end do
         do is=1,6
            if(is.eq.1 .or. (is.eq.2.and.Kmr(81).gt.0) .or.
     .       (is.gt.2.and.Kmr(84).gt.0.and.Npmhar(is).gt.0)) then
               call PRODCT(didx(1,1,i),Spcor(1,is),temp, 3,3,1)
               if(i.le.3) then
                  call PRODCT(deli,Mrotlb(1,i),temq,3,3,1)
                  do j=1,3
                     temp(j)=temp(j)-temq(j)
                  end do
               endif
               call CROSS(Spcor(1,is),temp,ddndxp)
               if(i.le.3) then
                  call CROSS(Mrotlb(1,i),delir(1,is),temq)
                  do j=1,3
                     ddndxp(j)=ddndxp(j)-temq(j)
                  end do
               endif
               do j = 1, 3
                  ddndx(j) = ddndx(j) +
     .             ddndxp(j)*3._10*Gamtrq(is)/Spc5(is)
               end do
               if(i.le.3) then
                  do j=1,3
                     ddndx(j)=ddndx(j)+5._10*delnp(j,is)*Mpcor(i,is)/
     .                Spc2(is)
                  end do
               endif
            endif
         end do
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
c calculate direct libration partial derivative wrt Love number or lag
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
      if(Kmr(81).gt.0) then
         call PRODCT(ddidp,Smcor,temq,3,3,1)
         call CROSS(Smcor,temq,temp)
         do i = 1,3
            temr(i) = temr(i) + temp(i)*3._10*Gamat/Rm5
         end do
      endif
      if(Kmr(84).gt.0) then
         do is=3,6
            if(Npmhar(is).gt.0) then
               call PRODCT(ddidp,Spcor(1,is),temq,3,3,1)
               call CROSS(Spcor(1,is),temq,temp)
               do i = 1,3
                  temr(i) = temr(i)+temp(i)*3._10*Gamtrq(is)/Spc5(is)
               end do
            endif
         end do
      endif
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
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

c calculate effect of elasticity on orbit partial derivatives

      entry MOREDPO(fno,kmr83)

c dependence on love number
      if(Icmtrl(Kkk).eq.-1012) then
         if(Kk.eq.1) then
            call PRODCT(ddidk2,Smecor,temq,3,3,1)
            if(Kmr(81).gt.0) call PRODCT(ddidk2,Smcor,temr,3,3,1)
            dtridx=1.5_10*(ddidk2(1,1)+ddidk2(2,2)+ddidk2(3,3))
            goto 390
         endif
         goto 400
c dependence on lag
      else if(Icmtrl(Kkk).eq.-1013) then
         if(Kk.eq.1) then
            call PRODCT(ddidl,Smecor,temq,3,3,1)
            if(Kmr(81).gt.0) call PRODCT(ddidl,Smcor,temr,3,3,1)
            dtridx=1.5_10*(ddidl(1,1)+ddidl(2,2)+ddidl(3,3))
            goto 390
         endif
         goto 400
      else if(Icmtrl(Kkk).eq.3) then
         if(kmr83.le.1) then
            fno=fno+dfndm(Kk)*Masse
         else
            fno=fno+dfndm(Kk)/Mass(3)
         endif
         goto 410
      else if(Icmtrl(Kkk).eq.10) then
         if(kmr83.le.1) then
            fno=fno-Delfn(Kk)/Mass(10)-dfndm(Kk)*Mass(3)
         else
            fno=fno-Delfn(Kk)/Mass(10)+dfndm(Kk)/Mass(10)
         endif
         goto 410
      endif
      return
c contributions directly to deformation
  390 didxrr=7.5_10*DOT(temq,Smecor)/Rem2-dtridx
      do i=1,3
         temp(i)=(3._10*temq(i)-didxrr*Smecor(i))*Gamtem/Mmoon/Rem5
      end do
      if(Kmr(81).gt.0) then
         didxrr=7.5_10*DOT(temr,Smcor)/Rm2-dtridx
         do i=1,3
            temp(i)=temp(i)+(3._10*temr(i)-didxrr*Smcor(i))*Gamat/
     .       Mmoon/Rm5
         end do
      endif
      call PRODCT(Mrotlb,temp,temr,-3,3,1)
  400 fno=fno+temr(Kk)
c indirect contributions
  410 fno=fno+DOTN(Dmcor(1,Kkk),dadxed(1,Kk),Index2)
      kkkrot=Icmref(Kkk)
      if(kkkrot.gt.0) fno=fno+DOTN(Dmlib(1,kkkrot),dadyed(1,Kk),Index2)
      return
      end
