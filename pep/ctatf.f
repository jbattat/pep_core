      real*10 function CTATF(jda,fra,nterm,n)
 
      implicit none
c
c rj cappallo   november 1973   sf ctatf, entry dctatf
c
c
c     this function computes ct - at for a given epoch using moyer's
c     formulation for annual and monthly terms (see tr 32-1527, p.14
c     from jpl), and an accurate( approx. 1 picosecond), but simple
c     diurnal term from d. robertson at mit 24-404.
c     the diurnal term suffers from a misconception.  asg 1977
c     the nominal variable part (rate) is saved in (d)ctatv
c     nterm=1  annual term only, site info not required
c     nterm=2  monthly and diurnal terms included
c     nterm=3  annual and monthly terms
c     nterm=4  integrated annual term, analytical monthly and diurnal
c     nterm=5  integrated annual term, analytical monthly
c     n=1      receive site
c     n=2      send site
c
c arguments
      integer jda,nterm,n
      real*10 fra

c array dimensions
      include 'globdefs.inc'

c commons
      include 'comdateq.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'param.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'
 
c local variables
      real*10 annl,DCTATF,frc
      integer*4 j,jdc,jdi
      real*10 t,jd1950/2433283._10/,manom,eanom,eandot,e/1.67527E-2_10/,
     . dluna,x(9),m1/6.248291_10/,m2/1.99096871E-7_10/,piece(4),
     . r3,xm(3),rm3,djd,d1/2.518410_10/,d2/2.462600818E-6_10/
      logical*4 ctmon

c external functions
      real*10 DOT

      CTATF    = 32.15_10
      Ctatv(n) = 0._10
      do j=1,4
         piece(j)=0._10
      end do
 
c input time is a.1--we want ct for orbital quantities
      call TIMINC(jda,fra,jdc,frc,32.15_10/864E2_10)
      djd = jdc
c
c scale variation term from dirac gdot theory
      t= Ctvary*(djd - Prm97 - 0.5_10 + frc)**2*Secday
      piece(1)=t
      CTATF= CTATF + t
      if(Atctsc.eq.0._10) goto 400
      t = (djd - jd1950 + frc)*Secday
      if(nterm.gt.3) then
 
c interpolate integrated annual term + monthly term
         jdi = jdc
         call CTREED(jdi,frc,annl,ctmon)
         if(jdi.gt.0) then
            Ctatv(n) = Ctatv(n) + annl
            piece(2) = annl
            if(ctmon) goto 200
            goto 100
         endif
      endif
 
c analytical annual term
      manom    = m1 + m2*t
      eanom    = manom + e*SIN(manom)
      piece(2) = 1.658E-3_10*SIN(eanom)
      Ctatv(n) = Ctatv(n) + piece(2)
      if(nterm.le.1) goto 300
 
c diurnal and monthly terms have been requested
  100 dluna    = d1 + d2*t
      piece(3) = 1.55E-6_10*SIN(dluna)
      Ctatv(n) = Ctatv(n) + piece(3)
  200 if(mod(nterm,2).eq.0) then
         call ERTCNT(jdc,frc,x,1,1)
         piece(4) = Aultvl*DOT(Xsite(1,n),x(4))
         Ctatv(n) = Ctatv(n) + piece(4)
      endif
  300 CTATF = CTATF + Atctsc*Ctatv(n)
  400 if(mod(Jct(6)/8192,2).eq.1) then
         if(Line.gt.56) call OBSPAG
         write(Iout,450) jda,fra,nterm,n,piece,CTATF
  450    format(' CTATF: JD.F=',i7,f13.12,' NT,N=',2i2,1p4d15.7,
     .    0pf18.14)
         Line = Line + 1
      endif
      return
 
 
 
c
c entry dctatf
c
      entry DCTATF(jda,fra,nterm,n)
 
      Dctatv(n) = 0._10
      if(Atctsc.ne.0._10) then
         call TIMINC(jda,fra,jdc,frc,32.15_10/864E2_10)
         t     = (jdc - jd1950 + frc)*Secday
         manom = m1 + m2*t
         eanom = manom + e*SIN(manom)
         eandot    = m2 + e*COS(manom)*m2
         Dctatv(n) = Dctatv(n) + 1.658E-3_10*COS(eanom)*eandot
         if(nterm.gt.1) then
            dluna     = d1 + d2*t
            Dctatv(n) = Dctatv(n) + 1.55E-6_10*COS(dluna)*d2
            if(mod(nterm,2).eq.0) then
               call ERTCNT(jdc,frc,x,0,1)
               call ERTCNT(jdc,frc,x,1,1)
 
c get moon position vector for accel. computation
               call ERTSNT(jdc,frc,xm,0,-1)
 
c need acceleration of earth for rate of ct - at
               r3  = DOT(x,x)**1.5_10
               rm3 = DOT(xm,xm)**1.5_10
               do j = 1,3
                  x(j + 6) = -(x(j)/r3 - Mass(10)*Mass(3)*xm(j)/rm3)
               end do
               Dctatv(n) = Dctatv(n) + Aultvl*DOT(Xsite(4,n),x(4))
     .                     + Auacc*Gauss**2*DOT(Xsite(1,n),x(7))
            endif
         endif
      endif
 
      DCTATF = Atctsc*Dctatv(n)
      if(mod(Jct(6)/8192,2).eq.1) then
         if(Line.gt.56) call OBSPAG
         write(Iout,550) jda,fra,nterm,n,DCTATF
  550    format(' DCTATF: JD.F=',i7,f13.12,' NT,N=',2i2,1pd23.15)
         Line = Line + 1
      endif
 
      return
      end
