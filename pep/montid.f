      subroutine MONTID(ncall, s)
      implicit none

c
c     r.king and r.cappallo   july 1977   subroutine montid
c     evaluation of tidal friction terms in right side of differential
c     equations for the moon's orbit.  called by morfn, this routine
c     incorporates m.slade's coding from monfn.
c Added indirect terms in partials due to dependence of tidal force on
c position and velocity - J.F.Chandler 2007 Nov

c arguments
      integer ncall
c ncall=  0 setup once per step
c      =  k evaluate acceleration for equation k.gt.0
c
      real*10 s
c s    = time argument

c array dimensions
      include 'globdefs.inc'

c commons
      include 'empcnd.inc'
      real*10 mcon(24)
      equivalence(mcon, Mcond(7))
      include 'intstf.inc'
      include 'morcrd.inc'
c
c     ecor  = earth relative to sun (re is distance)
c     mcor  = moon relative to sun (rm is distance)
c     mecor = moon relative to earth (rem is distance)
c     mcor etc. are equatorial coordinates
c     smcor etc. are selenodetic coordinates
      include 'morstf.inc'
      include 'nutces.inc'

c external functions
      real*10 DOT

c quantities internal to this routine
      real*10 aprm/9.37E-16_10/, eomg/6.300388099_10/
      real*10 eomga(3),uhat(3),nhat(3),that(3),ahat(3),tidac(3)
      real*10 amag,bhat12,bhat3,comn1,comn2,comn3,cos2d,cosd2,
     . coseps,cr,drdc1,drdc2,drdc3,dwdc1,dwdc2,dwdc3,
     . dzdc1,dzdc2,dzdc3,eomt,eomr,omgzro,omgzro2,r,ratlg,sin2d,
     . tcen,vrad,vtrn,vtrnv(3),w,xmagn,z
      integer i,j

      real*10 db12dx(3),db3dx(3),odoidx(3),
     . db12dv(3),db3dv(3),odoidv(3)
c
c setup once per step for partials
      if(ncall-Knbd.eq.10 .and. ncall.lt.Korb) then
         call CROSS(Mecor(4),eomga,odoidx)
         do i=1,3
            db12dx(i)=eomg*nhat(i)*(vrad*eomt-vtrn*eomr)/vtrn/Rem +
     .       (that(i)*vrad + uhat(i)*vtrn)/Rem2
            db3dx(i)=eomg*(nhat(i)*coseps*vrad/vtrn+that(i)*eomr)/Rem
            odoidx(i)=(eomg*(odoidx(i)-2._10*coseps*vtrn*uhat(i)+
     .       eomr*eomg*Rem*(eomga(i)-eomr*uhat(i)))
     .       + vtrn/Rem*(that(i)*vrad+uhat(i)*vtrn))/Rem2/omgzro2
            db12dv(i)=-eomg*eomt*nhat(i)/vtrn - that(i)/Rem
            db3dv(i)=-eomg*coseps*nhat(i)/vtrn
            odoidv(i)=(eomg*(coseps*that(i)-eomt*nhat(i))/Rem - 
     .       vtrnv(i)/Rem2)/omgzro2
         end do
         do i=1,3
            do j=1,3
               Dadxt(j,i)=-(7._10*uhat(j)*tidac(i)+tidac(j)*uhat(i))/Rem
     .          + aprm*cr*sin2d/omgzro* (that(i)*
     .          (db12dx(j)+bhat12*odoidx(j)+bhat3*nhat(j)*vrad/vtrn/Rem)
     .          +nhat(i)*(db3dx(j)+bhat3*odoidx(j)
     .          -bhat12*nhat(j)*vrad/vtrn/Rem))
               Dadvt(j,i)=aprm*cr*sin2d/omgzro *
     .          (that(i)*(db12dv(j) + odoidv(j)*bhat12
     .          - bhat3*nhat(j)/vtrn)
     .          + nhat(i)*(-nhat(j)/Rem + odoidv(j)*bhat3))
            end do
            Dadxt(i,i)=Dadxt(i,i)+r/Rem
         end do
      endif
c determine tidal friction quantities for motion
      if(ncall.le.0) then
         cr   = (2.5668E-3/Rem)**7
         tcen = (s - 2415020.813_10)/36524.21988_10
 
c 1900.0 = jed 2415020.313
         sin2d = mcon(16) - (mcon(17) + mcon(18)*tcen)*tcen
         cos2d = sqrt(1. - sin2d**2)
         ratlg = sin2d/cos2d
         cosd2 = 0.5*(cos2d+1.)
 
c eomg is Earth spin in radians/day, eomga is unit vector
         eomga(1)= Nutprc(1,3)
         eomga(2)= Nutprc(2,3)
         eomga(3)= Nutprc(3,3)
         uhat(1)= Mecor(1)/Rem
         uhat(2)= Mecor(2)/Rem
         uhat(3)= Mecor(3)/Rem
         comn1= uhat(2)*Mecor(6) - uhat(3)*Mecor(5)
         comn2= uhat(3)*Mecor(4) - uhat(1)*Mecor(6)
         comn3= uhat(1)*Mecor(5) - uhat(2)*Mecor(4)
         vtrn= sqrt(comn1**2 + comn2**2 + comn3**2)
         nhat(1)= comn1/vtrn
         nhat(2)= comn2/vtrn
         nhat(3)= comn3/vtrn
         coseps=DOT(eomga,nhat)
         vrad=DOT(Mecor(4),uhat)
         do i=1,3
            vtrnv(i)= Mecor(i+3) - vrad*uhat(i)
            that(i)= vtrnv(i)/vtrn
         end do
         xmagn = vtrn/Rem
         eomr=DOT(eomga,uhat)
         eomt=DOT(eomga,that)
         ahat(1)= eomga(2)*nhat(3) - eomga(3)*nhat(2)
         ahat(2)= eomga(3)*nhat(1) - eomga(1)*nhat(3)
         ahat(3)= eomga(1)*nhat(2) - eomga(2)*nhat(1)
         amag   = sqrt(ahat(1)**2 + ahat(2)**2 + ahat(3)**2)
         ahat(1)= ahat(1)/amag
         ahat(2)= ahat(2)/amag
         ahat(3)= ahat(3)/amag
c         cosphi = ahat(1)*uhat(1) + ahat(2)*uhat(2) + ahat(3)*uhat(3)
c         if(abs(cosphi).ge.1.) then
c            sinphi=0.
c         else
c            sinphi=sign(sqrt(1. - cosphi**2),
c     .       uhat(1)*eomga(1) + uhat(2)*eomga(2) + uhat(3)*eomga(3))
c         endif
c         bhat1  = -(eomg*coseps - xmagn)*sinphi
c         bhat2  = (eomg*coseps - xmagn)*cosphi
         bhat12  = eomg*coseps - xmagn
         bhat3=-eomg*eomt
c         omgzro = sqrt(bhat1**2 + bhat2**2 + bhat3**2)
         omgzro2 = bhat12**2 + bhat3**2
         omgzro = SQRT(omgzro2)
         r      = -aprm*cr*(3.*cosd2 + 0.33)
c         z      = aprm*cr*sin2d*bhat12/omgzro
         dzdc1 = aprm*cr*bhat12/omgzro
         z = dzdc1*sin2d
c         w      = aprm*cr*sin2d*bhat3/omgzro
         dwdc1 = aprm*cr*bhat3/omgzro
         w = dwdc1*sin2d

         drdc1  = aprm*cr*1.5*ratlg
         drdc2  = -drdc1*tcen
         drdc3  = drdc2*tcen
         dzdc2  = -dzdc1*tcen
         dzdc3  = dzdc2*tcen
         dwdc2  = -dwdc1*tcen
         dwdc3  = dwdc2*tcen
c
c tidal effects on the motion of the moon
      else if(Kkk.le.0) then
         tidac(Kk) = r*uhat(Kk) + w*nhat(Kk) + z*that(Kk)
         Fn(1)     = Fn(1) + tidac(Kk)
c
c tidal effects on partial derivatives
      else if(Icmtrl(Kkk).eq.-16) then
         Fn(1) = Fn(1)+drdc1*uhat(Kk)+dwdc1*nhat(Kk)+dzdc1*that(Kk)
      else if(Icmtrl(Kkk).eq.-17) then
         Fn(1) = Fn(1)+drdc2*uhat(Kk)+dwdc2*nhat(Kk)+dzdc2*that(Kk)
      else if(Icmtrl(Kkk).eq.-18) then
         Fn(1) = Fn(1)+drdc3*uhat(Kk)+dwdc3*nhat(Kk)+dzdc3*that(Kk)
      endif
 
      return
      end
