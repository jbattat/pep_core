      subroutine PLNSAT(ncall)
 
      implicit none

c        computation of forces acting on   natural  planetary orbiter
c arguments
      integer*4 ncall
c        ncall=-2 setup once per iteration of a given step for partials
c        ncall=-1 setup once per iteration of a given step for motion
c        ncall= 0 setup once per step
c        ncall= k evaluate acceleration for equation k.gt.0
c relativity added 1978 sep - j.f.chandler
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'intstf.inc'
      include 'param.inc'
      include 'petuna.inc'
      include 'sbstuf.inc'
      include 'sbthng.inc'
      real*10 vccor(3),vsbcor(3)
      equivalence (vccor,Ccor(4)),(vsbcor,Sbcor(4))

c external functions
      real*10 DOT
c
c local variables
      real*10 c2r3,rcnt,rf0,rf1,rf2,rsun,rv2r,tf2,tf3,
     .          tf4,tg,tgm,tvr,vc2,vcvsb,vcxsb,vsb2,xc3xsb
      integer   i,j
      real*10 u4v(3)
      integer*4 iout/6/
 
      call HRTGCN(ncall)
      if(ncall.lt.0) then
c
c set-up once per iteration of a given step for motion
         if(ncall.lt.-1) then
c
c set-up once per iteration for partials
c set up for relativity
            if(Kp(61).ge.1) then
 
c compute dadx, dadv matrices
               tg = ((-2._10*rv2r-rcnt)*rf0 - 3._10*rf1)/Rsb2
               do i = 1,3
                  do j = 1,3
                     Dadxr(i,j) = tg*Sbcor(i)*Sbcor(j)
                     if(j.eq.i) goto 10
                     Dadxr(j,i) = Dadxr(i,j)
                  end do
   10             Dadxr(i,i) = Dadxr(i,i) + rf1
               end do
 
c finish up matrix with non-symmetric part, do dadvr
               tvr = 3._10*vcxsb/Rsb2*rf0
               tgm = 0.5_10*Gamat/Rc3*rf0
               tf3 = 3._10*rf2/Rsb2
               tf2 = rf0 + rf0
               tf4 = tf2 + tf2
               do i = 1,3
                  do j = 1,3
                     Dadxr(i,j) = Dadxr(i,j)
     .                             + (tvr*vccor(i) - tgm*Ccor(i))
     .                             *Sbcor(j)
     .                             + (rf0*u4v(i) - tf3*Sbcor(i))
     .                             *vsbcor(j)
                     Dadvr(i,j) = tf2*(vccor(i) - vsbcor(i))
     .                             *Sbcor(j) + tf4*Sbcor(i)
     .                             *vsbcor(j)
                  end do
                  Dadvr(i,i) = Dadvr(i,i) + rf2
               end do
               if(Pntflg) then
                  if(Line.ge.53) then
                     call NEWPGT(iout,Npage,24000)
                     Line = 2
                  endif
                  Line = Line + 7
                  write(iout,20) Dadxr,Dadvr
   20             format('0  IN PLNSAT  DADXR:',3(t26,1p,3D25.15/),
     .                   14x,'DADVR:',(t26,1p,3D25.15))
               endif
            endif
 
c set up for relativity
         else if(Kp(61).ge.0) then
            vsb2   = DOT(vsbcor,vsbcor)
            vcvsb  = DOT(vccor,vsbcor)
            vcxsb  = DOT(vccor,Sbcor)
            rv2r   = 1.5_10*vcxsb**2/Rsb2
            xc3xsb = DOT(Ccor3,Sbcor)*Gamat*0.5_10
            rcnt   = -4._10*Gamat3/Rsb
            c2r3   = Cvel2*Rsb3
            do i = 1,3
               u4v(i) = 4._10*vsbcor(i) + vccor(i)
            end do
            rf0 = -Gamat3/c2r3
            rf1 = (rcnt+rsun-vsb2+2._10*vcvsb+vc2+rv2r-xc3xsb)*rf0
            rf2 = DOT(Sbcor,u4v)*rf0
         endif
      else if(ncall.eq.0) then
c
c set-up once per step
c set up for relativity
         if(Kp(61).ge.0) then
            vc2  = DOT(vccor,vccor)
            rsun = 5._10*Gamat/Rc
         endif
c
c
c
c compute acceleration for motion
      else if(Kkk.gt.0) then
c
c
c compute accelerations for partials
c effect of relativity on partials
         if(Kp(61).ge.1) then
 
c check if alpha is relativity motion factor
            if(Icntrl(Kkk).eq.31) then
               Fn(1) = Fn(1) + Sumr(Kk)
 
c check if alpha is central body mass
            else if(Icntrl(Kkk).eq.Ncentr) then
               Fn(1) = Fn(1) + Relfct*(Sumr(Kk)/rf0 + rcnt*Sbcor(Kk))
     .                 *Gamat/c2r3
 
c check if alpha is time variable gravitation constant
            else if(Icntrl(Kkk).eq.32) then
               Fn(1) = Fn(1) - Relfct*(Sumr(Kk)/rf0 + (rcnt+rsun-xc3xsb)
     .                 *Sbcor(Kk))*Gama3*Tvary
            endif
            Fn(1) = Fn(1)
     .              + (DOT(Dsbcor(1,Kkk),Dadxr(1,Kk)) + DOT(Dsbcor(4,
     .              Kkk),Dadvr(1,Kk)))*Relfct
         endif
 
c effect of relativity on motion
      else if(Kp(61).ge.0) then
         Sumr(Kk) = rf1*Sbcor(Kk) + rf2*vsbcor(Kk)
         Fn(1)    = Fn(1) + Relfct*Sumr(Kk)
      endif
      return
      end
