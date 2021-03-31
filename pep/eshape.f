      subroutine ESHAPE(n)
 
      implicit none
 
c     m.e.ash   jan 1967    subroutine eshape
c     calculate quantities relating to sites on a non-spherical earth
c     or to observing sites and observed spots on other planets

c parameters
      integer*4 n
c     n =1 subroutine optic is calling program (receiving site only)
c     n =2 subroutine radar is calling program (receiving and sending
c          sites)
c parameter to entry point PLATMO
      integer*4 jd

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'empcnd.inc'
      real*10 eflat,erad,mrad
      equivalence (Econd(7),erad), (Econd(8),eflat), (Mcond(7),mrad)
      include 'eqnphs.inc'
      include 'funcon.inc'
      include 'mnsprt.inc'
      include 'number.inc'
      include 'param.inc'
      include 'sitcrd.inc'

c local
      real*10 crds1,h1,latr,latr0,q2,q3,q4,qq,rc1,rs1,fltc,flts,velfct
      integer*4 i,jdlast(2),k,ngo,niter,npspot,nlast
c
c           ksite(n)   =-2  for geodetic coordinates on earth
c           ksite(n)   =-1  for cylindrical polar coords.on earth
c           ksite(n)   = 0  for spherical polar coords.on earth
c           ksite(n)=positive   spherical polar coords.on body ksite(n)
c
c           setup for observing site coordinates
      velfct=1E-6_10/Ltvel/365.25_10
      nlast=n
      do i = 1,nlast
         jdlast(i)=-999
         ngo = 1
         Longr(i) = Coords(2,i)*Convd
         if(Ksite(i).lt.-1) then
c
c geodetic site coordinates (explanatory supplement, pp 57-58)
            latr = Convd*Coords(3,i)
         else if(Ksite(i).eq.-1) then
c
c cylindrical site coordinates
            Rc(i) = Coords(1,i)/Ltvel
            Rs(i) = Coords(3,i)/Ltvel
            latr  = ATAN(Rs(i)/Rc(i))
            goto 150
         else
            goto 100
         endif
   50    Cnrm(i) = COS(latr)
         Snrm(i) = SIN(latr)
         q2 = COS(2.0_10*latr)
         q3 = COS(4.0_10*latr)
         q4 = COS(6.0_10*latr)
         fltc=Cc(1)+Cc(2)*q2+Cc(3)*q3+Cc(4)*q4
         flts=Ss(1)+Ss(2)*q2+Ss(3)*q3+Ss(4)*q4
         Shgt(i) = Coords(1,i)/1000._10
         Rc(i)   = (Shgt(i) + erad*fltc)*Cnrm(i)/Ltvel
         Rs(i)   = (Shgt(i) + erad*flts)*Snrm(i)/Ltvel
         if(ngo.eq.1) goto 250
         latr = ATAN((rs1*Rc(i)*Snrm(i))/(rc1*Rs(i)*Cnrm(i)))
         if(niter.le.0) goto 200
         h1 = ((rc1-Rc(i))*Cnrm(i)+(rs1-Rs(i))*Snrm(i))*Ltvel*1000._10
         Coords(1,i) = Coords(1,i) + h1
         if(h1.gt.1.E-4_10) goto 200
         if(ABS(latr-latr0).ge.1.E-14_10) goto 200
         Rc(i) = rc1
         Rs(i) = rs1
         Coords(1,i) = crds1
         Cnrm(i) = COS(latr)
         Snrm(i) = SIN(latr)
         goto 250
c
c spherical site coordinates
  100    Rsite(i) = Coords(1,i)/Ltvel
         latr     = Coords(3,i)*Convd
         Rc(i)    = Rsite(i)*COS(latr)
         Rs(i)    = Rsite(i)*SIN(latr)
         if(Ksite(i).gt.0 .and. Ksite(i).ne.3) goto 250
c
c iteration to determine site normal quantities
  150    Cnrm(i) = 1.0_10
         Snrm(i) = 0.0_10
         Shgt(i) = Coords(1,i) - erad
         if(ABS(Rs(i)).le.1.E-12_10) goto 250
         rc1   = Rc(i)
         rs1   = Rs(i)
         crds1 = Coords(1,i)
         Coords(1,i) = 0.0_10
         latr0 = latr
         niter = 0
         ngo   = 2
         goto 50
  200    latr0 = latr
         niter = niter + 1
         if(niter.gt.100) call SUICID(
     .       ' MORE THAN 100 SITE NORMAL ITERATIONS, STOP IN ESHAPE   '
     .       , 14)
         goto 50
c
c site coordinates fixed in earth
  250    Clong(i)  = COS(Longr(i))
         Slong(i)  = SIN(Longr(i))
         Xb0(1,i) = Rc(i)*Clong(i)
         Xb0(2,i) = Rc(i)*Slong(i)
         if(Ksite(i).lt.0) Rsite(i) = SQRT(Rc(i)**2 + Rs(i)**2)

c some useful quantities
         Cphi(i) = Rc(i)/Rsite(i)

         if(Ksite(i).eq.-1) then

c cylindrical coordinates, must form spherical partials
            Cyl1(i) = Cphi(i)
            Cyl2(i) = Rs(i)/Rsite(i)
            Cyl3(i) = -Rs(i)
            Cyl4(i) = Rc(i)
         else
 
c spherical coordinates for this site
            Cyl1(i) = 1._10
            Cyl2(i) = 0._10
            Cyl3(i) = 0._10
            Cyl4(i) = 1._10
         endif

         Rc0(i)=Rc(i)
         Rs0(i)=Rs(i)
         Cnrm0(i)=Cnrm(i)
         Snrm0(i)=Snrm(i)
         Rsite0(i)=Rsite(i)
         Longr0(i)=Longr(i)
         Xb00(1,i)=Xb0(1,i)
         Xb00(2,i)=Xb0(2,i)
         Clong0(i)=Clong(i)
         Slong0(i)=Slong(i)
         Shgt0(i)=Shgt(i)

         Tecton(i)= Coords(4,i).ne.0._10 .or. Coords(5,i).ne.0._10 .or.
     .    Coords(6,i).ne.0._10
c prepare velocities of site quantities
         if(Tecton(i)) then
            Drc(i)=(Coords(4,i)*Rc(i)-Coords(6,i)*Rs(i))/Rsite(i)
            Drs(i)=(Coords(4,i)*Rs(i)+Coords(6,i)*Rc(i))/Rsite(i)
            Drsite(i)=Coords(4,i)
            Dlongr(i)=Coords(5,i)/Rc(i)
            Dxb0(1,i)=((Coords(4,i)*Rc(i)-Coords(6,i)*Rs(i))*
     .       Clong(i)/Rsite(i) - Coords(5,i)*Slong(i))
            Dxb0(2,i)=((Coords(4,i)*Rc(i)-Coords(6,i)*Rs(i))*
     .       Slong(i)/Rsite(i) + Coords(5,i)*Clong(i))
            Dclong(i)=-Coords(5,i)*Slong(i)/Rc(i)
            Dslong(i)=Coords(5,i)*Clong(i)/Rc(i)
            Dshgt(i)=(Coords(4,i)*(Cnrm(i)*Rc(i)+Snrm(i)*Rs(i))+
     .       Coords(6,i)*(Cnrm(i)*Rs(i)-Snrm(i)*Rc(i)))*Ltvel/Rsite(i)
            Dcnrm(i)=(Drc(i)*Ltvel-Cnrm(i)*Dshgt(i))/(Shgt(i)+erad*fltc)
            Dsnrm(i)=(Drs(i)*Ltvel-Snrm(i)*Dshgt(i))/(Shgt(i)+erad*flts)
         endif

      end do
c
c setup for spot coordinates on observed body
      if(Nspot.le.0) goto 400
      k = 1
      npspot = Nplnt0
      if(Nplsr.gt.0) then
         Deltc(1)=Spcdx(3,k)*Convd
         Cdeltc=COS(Deltc(1))
         Sdeltc=SIN(Deltc(1))
         Deltc(2)=Psrprm(3)*Convds/365.25_10
         Alphc(1)=Spcdx(2,k)*Convd
         Calphc=COS(Alphc(1))
         Salphc=SIN(Alphc(1))
         Alphc(2)=Psrprm(2)*Convds/365.25_10/Cdeltc
      endif
         
  300 qq = Spcdx(1,k)/Ltvel
      if(npspot.lt.0) qq = 1.0_10
      latr = Spcdx(3,k)*Convd
      Yspcd(3,k) = qq*SIN(latr)
      qq   = qq*COS(latr)
      latr = Spcdx(2,k)*Convd
      Yspcd(2,k)  = qq*SIN(latr)
      Yspcd(1,k)  = qq*COS(latr)
      Rspot(k)    = SQRT(Yspcd(1,k)**2 + Yspcd(2,k)**2 + Yspcd(3,k)**2)
      Dydphi(1,k) = -Yspcd(3,k)*COS(latr)
      Dydphi(2,k) = -Yspcd(3,k)*SIN(latr)
      Dydphi(3,k) = qq
      if(k.eq.2) return
  400 if(Nspot2.gt.0) then
         k = 2
         npspot = Nplnt2
         goto 300
      endif
 
      return

      entry PLATMO(jd)
c update plate motion to a specified epoch

c arguments
c jd - Julian day number

      do i=1,nlast
         if(jd.ne.jdlast(i)) then
            jdlast(i)=jd
            Dtsit(i)=(jd-T0sit(i))*velfct
            if(Tecton(i)) then
               Rc(i)   =Rc0(i)+Dtsit(i)*Drc(i)
               Rs(i)   =Rs0(i)+Dtsit(i)*Drs(i)
               Cnrm(i) =Cnrm0(i)+Dtsit(i)*Dcnrm(i)
               Snrm(i) =Snrm0(i)+Dtsit(i)*Dsnrm(i)
               Rsite(i)=Rsite0(i)+Dtsit(i)*Drsite(i)
               Longr(i)=Longr0(i)+Dtsit(i)*Dlongr(i)
               Xb0(1,i)=Xb00(1,i)+Dtsit(i)*Dxb0(1,i)
               Xb0(2,i)=Xb00(2,i)+Dtsit(i)*Dxb0(2,i)
               Clong(i)=Clong0(i)+Dtsit(i)*Dclong(i)
               Slong(i)=Slong0(i)+Dtsit(i)*Dslong(i)
               Shgt(i)=Shgt0(i)+Dtsit(i)*Dshgt(i)
            endif
         endif
      end do
      return
      end
