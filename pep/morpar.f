      subroutine MORPAR
 
      implicit none
 
c     king/cappallo may 1976   subroutine morpar
c     calculate partials of various constants in lunar librraion
c     equations with respect to parameters for which partials are
c     integrated
c
c     morpar is libpar modified for mor routines

c array dimensions
      include 'globdefs.inc'

c commons
      include 'empcnd.inc'
      real*10 alpha,beta,gamma
      equivalence (Mrcond(8),alpha),(Mrcond(9),beta),(Mrcond(10),gamma)
      include 'harmor.inc'
      real*10 mz2, mc22
      equivalence (Zhar,mz2),(Char(2),mc22)
      include 'intstf.inc'
      include 'morstf.inc'
      real*10 acore,bcore,ccore
      equivalence (I0c0(1,1),acore),(I0c0(2,2),bcore),(I0c0(3,3),ccore)
      include 'param.inc'
c
c quantities local to this routine
      real*10 a2,b1,bg1,bg2,fcc,g1
      integer i,msub
c
c derivatives of contants w.r.t. parameters integrated

c In principle, these derivatives should be set up also in parallel for
c the orbit integration, but there would be no point in integrating the
c orbit partial if the rotation partial were not done simultaneously.
c Therefore, a single setup suffices, and the orbit uses pointers to
c the derivatives organized by rotation order.

c The only exception is the occasionally desired artificial case with
c integration of orbit without rotation.  In this case, ROTINT is
c false, and the rotation pointers would all be null, but it is
c sufficient to fake the setup for the lunar J2 partial if the orbit
c partial w.r.t. J2 is desired.

      if(Imzone(1).gt.0 .and. .not.Rotint) then
         Izone(1) = 1
         Icrtrl(1)=Jzone
         Icrref(1)=Imzone(1)
         Icmref(Imzone(1))=1
      endif

      a2  = Mrad**2
      b1  = beta + 1.0_10
      g1  = gamma + 1.0_10
      bg1 = 1.0_10 - beta*gamma
      bg2 = (beta*(2.0_10+gamma) - gamma)**2
      fcc = Mrcond(25)*Mrcond(24)
      do i = 1, i_mxplprt
         Dalpha(i) = 0._10
         Dbeta(i)  = 0._10
         Dgamma(i) = 0._10
         Dmma(i)   = 0._10
         Dmmb(i)   = 0._10
         Dmmc(i)   = 0._10
         Dmmam(i)  = 0._10
         Dmmbm(i)  = 0._10
         Dmmcm(i)  = 0._10
         Dc22(i)   = 0._10
         if(Icrtrl(i).eq.0) goto 900
         if(Icrtrl(i).eq.-3) then
 
c partials w.r.t. beta
            Dmmc(i)   = g1/a2/mz2/b1**2
            Dmma(i)   = (Dmmc(i)*b1 + Mmcw*g1/bg1)/bg1
            Dmmb(i)   = (Dmmc(i)*b1 + Mmcw)/g1
            Dalpha(i) = (1.0_10 - gamma**2)/bg1**2
            Dbeta(i)  = 1.0_10
            Dc22(i)   = -mz2*gamma*g1/bg2
            if(Corint) then
               Dmma(i)   = Dmma(i)*(Awhole/Amantl)**2
               Dmmb(i)   = Dmmb(i)*(Bwhole/Bmantl)**2
               Dmmc(i)   = Dmmc(i)*(Cwhole/Cmantl)**2
               Dmmam(i)  = Dmma(i)/Mmoon
               Dmmbm(i)  = Dmmb(i)/Mmoon
               Dmmcm(i)  = Dmmc(i)/Mmoon
               Dalpha(i) = Dalpha(i)*Awhole/Amantl +
     .          Dmmam(i)*(alpha*acore-fcc)
               Dbeta(i)  = Bwhole/Bmantl +
     .          Dmmbm(i)*(beta*bcore-fcc)
               Dgamma(i) = Dmmcm(i)*gamma*ccore
            endif
         else if(Icrtrl(i).eq.-4) then
 
c partials wrt gamma
            Dmmc(i)   = (beta - 1.0_10)/a2/mz2/b1/2.0_10
            Dmma(i)   = (Dmmc(i)*b1 + Mmcw*beta*b1/bg1)/bg1
            Dmmb(i)   = (Dmmc(i) - Mmcw/g1)*b1/g1
            Dalpha(i) = (beta**2 - 1.0_10)/bg1**2
            Dgamma(i) = 1.0_10
            Dc22(i)   = mz2*beta*b1/bg2
            if(Corint) then
               Dmma(i)   = Dmma(i)*(Awhole/Amantl)**2
               Dmmb(i)   = Dmmb(i)*(Bwhole/Bmantl)**2
               Dmmc(i)   = Dmmc(i)*(Cwhole/Cmantl)**2
               Dmmam(i)  = Dmma(i)/Mmoon
               Dmmbm(i)  = Dmmb(i)/Mmoon
               Dmmcm(i)  = Dmmc(i)/Mmoon
               Dalpha(i) = Dalpha(i)*Awhole/Amantl +
     .          Dmmam(i)*(alpha*acore-fcc)
               Dbeta(i)  = Dmmbm(i)*(beta*bcore-fcc)
               Dgamma(i) = Cwhole/Cmantl +
     .          Dmmcm(i)*gamma*ccore
            endif

c partials w.r.t. j2
         else if(Icrtrl(i).eq.Jzone) then
            if(i.eq.Izone(1)) then
               Dmmc(i) = -Mmcw/mz2
               Dmma(i) = Dmmc(i)*b1/bg1
               Dmmb(i) = Dmmc(i)*b1/g1
               Dc22(i) = mc22/mz2
               if(Corint) then
                  Dmma(i)   = Dmma(i)*(Awhole/Amantl)**2
                  Dmmb(i)   = Dmmb(i)*(Bwhole/Bmantl)**2
                  Dmmc(i)   = Dmmc(i)*(Cwhole/Cmantl)**2
                  Dmmam(i)  = Dmma(i)/Mmoon
                  Dmmbm(i)  = Dmmb(i)/Mmoon
                  Dmmcm(i)  = Dmmc(i)/Mmoon
                  Dalpha(i) = Dmmam(i)*(alpha*acore-fcc)
                  Dbeta(i)  = Dmmbm(i)*(beta*bcore-fcc)
                  Dgamma(i) = Dmmcm(i)*gamma*ccore
               endif
            endif

c partials w.r.t. core moment
         else if(Corint .and. Icrtrl(i).eq.-18) then
            Dmmam(i)  = (1._10-Mrcond(25))/Amantl**2
            Dmmbm(i)  = (1._10-Mrcond(25))/Bmantl**2
            Dmmcm(i)  = 1._10/Cmantl**2
            Dmma(i)   = Mmoon*Dmmam(i)
            Dmmb(i)   = Mmoon*Dmmbm(i)
            Dmmc(i)   = Mmoon*Dmmcm(i)
            Dalpha(i) = (Cmantl-Bmantl)*Dmmam(i) - Mrcond(25)/Amantl
            Dbeta(i)  = (Cmantl-Amantl)*Dmmbm(i) - Mrcond(25)/Bmantl
            Dgamma(i) = gamma*ccore*Dmmcm(i) + gamma/Cmantl

c partials w.r.t. core flattening
         else if(Corint .and. Icrtrl(i).eq.-19) then
            Dmmam(i)  = -Mrcond(24)/Amantl**2
            Dmmbm(i)  = -Mrcond(24)/Bmantl**2
            Dmma(i)   = Mmoon*Dmmam(i)
            Dmmb(i)   = Mmoon*Dmmbm(i)
            Dalpha(i) = (Cmantl-Bmantl)*Dmmam(i) - ccore/Amantl
            Dbeta(i)  = (Cmantl-Amantl)*Dmmbm(i) - ccore/Bmantl

c partials w.r.t. earth or moon mass
         else if(Corint .and. (Icrtrl(i).eq.3.or.Icrtrl(i).eq.10)) then
            msub=13-Icrtrl(i)
            Dmma(i) = -Mass(msub)*acore/Amantl**2
            Dmmb(i) = -Mass(msub)*bcore/Bmantl**2
            Dmmc(i) = -Mass(msub)*ccore/Cmantl**2
            msub=Icrtrl(i)
            Dmmam(i)  = -Awhole/Mass(msub)/Amantl**2
            Dmmbm(i)  = -Bwhole/Mass(msub)/Bmantl**2
            Dmmcm(i)  = -Cwhole/Mass(msub)/Cmantl**2
            Dalpha(i) = (alpha*acore-fcc)*Dmmam(i)
            Dbeta(i)  = (beta*bcore-fcc)*Dmmbm(i)
            Dgamma(i) = Mrcond(24)*gamma*Dmmcm(i)
         endif
      end do
  900 return
      end
