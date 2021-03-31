      subroutine SITCOR(sidtm,n,nvel,norm)
 
      implicit none
 
c  m.e.ash   nov 1966    subroutine sitcor
c  coordinates referred to mean equinox and equator of reference epoch
c  determined for site n.

c arguments
      real*10 sidtm
      integer*4 n,nvel,norm

c         sidtm= true sidereal time
c         n    = 1 receiving site coordinates determined
c         n    = 2 sending site coordinates determined
c         nvel =-1 position and nothing else determined
c         nvel = 0 position determined
c         nvel = 1 position and velocity determined
c         nvel = 2 position, velocity and acceleration determined
c         nvel = 3 position,velocity,acceleration and jerk determined
c         norm = 0 normal to site not determined
c         norm = 1 normal to site determined (referred to mean
c                 equinox and equator of reference epoch)
c         norm = 2 normal to site determined and partial
c                  derivatives w.r.t. site coordinates determined
c                  even if lsite(1-3)=0
c                  (for solid body tides or fluid displacements)

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdateq.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'nutprc.inc'
      include 'param.inc'
      include 'sitcrd.inc'
c
c          kindnp = 1 derivative of precession-nutation used in
c                     computing velocity
c          kindnp = 0 derivative of precession-nutation not used
c
c           ksite(n)   =0   for spherical polar coords.
c           ksite(n)   =-1  for cylindrical polar coords.
c           xsite(.,n) =position,velocity,acceleration of site n
c                       (distance and time units are light-seconds and
c                       seconds)
c           sitnrm(.,n)=unit normal to spheroid at site n
c           pstlat(.,n)=partial derivatives of site position and
c                       velocity with respect to center of earth
c                       latitude in radians (if ksite(n)=0)
c                       or else with respect to rs=z (if ksite(n)=-1)
c           pstlon(.,n)=partial derivatives of site position and
c                       velocity with respect to longitude in radians
      real*10 rsv,sidtp,sidvl2,sidvlt,tfct,z(3),zm(3)
      integer*4 i,j
      character*4 type(4)/'COR', 'VEL', 'ACC', 'JRK'/
c
c*  start=100
c main effect of rotation
c         Cw&Sw  = cos&sin (sidtm - Longr(n))
      Sv(n) = SIN(sidtm)
      Cv(n) = COS(sidtm)
      Sw(n) = Sv(n)*Clong(n) - Cv(n)*Slong(n)
      Cw(n) = Cv(n)*Clong(n) + Sv(n)*Slong(n)
c include effect of wobble on site position
      if(Iwob.le.0) then
c
c determine site position assuming wobble is zero
         Y(1,n) = Rc(n)*Cw(n)
         Y(2,n) = Rc(n)*Sw(n)
         Y(3,n) = Rs(n)
      else
c xb0(1,n),-xb0(2,n),rs(n)    right hand coordinates
c xb1(1,n),-xb1(2,n),rs(n)    right hand coordinates
c ysite(1,n),ysite(2,n),ysite(3,n)   right hand coordinates
         Xb1(1,n) = Xb0(1,n) - Xwob*Rs(n)
         Xb1(2,n) = Xb0(2,n) - Ywob*Rs(n)
         Y(1,n) = Cv(n)*Xb1(1,n) + Sv(n)*Xb1(2,n)
         Y(2,n) = Sv(n)*Xb1(1,n) - Cv(n)*Xb1(2,n)
         Y(3,n) = Rs(n) + (Xwob*Xb0(1,n) + Ywob*Xb0(2,n))
      endif
c
c apply precession-nutation matrix to position vector
      call CHNCOR(Xsite(1,n),Y(1,n))
 
      if(mod(Jct(6)/4,2).ne.0) then
         if(Line.gt.56) call OBSPAG
         do i = 1, 3
            z(i) = Xsite(i,n)/Aultsc
            zm(i)= Xsite(i,n)*Ltvel
         end do
         sidtp = MOD(sidtm,Twopi)
         if(sidtp.lt.0._10) sidtp=sidtp+Twopi
         write(Iout,50) type(1),n,sidtp,(Xsite(i,n),i = 1,3),z,zm
   50    format(' SIT',a3,': N=',i1,'  SIDTM=',f13.11,'  X=', 1p,
     .    3D22.14,' (S)'/ 37x,3D22.14,' (AU)'/ 37x,3D22.14,' (KM)')
         Line = Line + 3
      endif
c
c*  start=300
c
c determine if nothing except position is to be calculated
      if(nvel.lt.0) return
      entry SITVEL(n,nvel,norm)

      tfct=Dtsit(n)/Aultsc
c
c determine site normal
      if(norm.gt.0) then
         z(1) = Cw(n)*Cnrm(n)
         z(2) = Sw(n)*Cnrm(n)
         z(3) = Snrm(n)
         call CHNCOR(Sitnrm(1,n),z)
      endif
c
c determine partial derivatives of position
c note: wobble terms not included in partials
      if(Lsite(1,n).gt.0 .or. norm.eq.2) then
         if(Ksite(n).lt.0) then
            z(1) = Cw(n)
            z(2) = Sw(n)
            z(3) = 0._10
            call CHNCOR(Pstrad(1,n),z)
         else
            do i = 1, 3
               Pstrad(i,n) = Xsite(i,n)/Rsite(n)
            end do
         endif
      endif
      if(Lsite(2,n).gt.0 .or. norm.eq.2) then
         do i = 1, 3
            Pstlon(i,n) = Nutpr(1,i)*Y(2,n) - Nutpr(2,i)*Y(1,n)
         end do
      endif
      if(Lsite(3,n).gt.0 .or. norm.eq.2) then
         if(Ksite(n).ge.0) then
            z(1) = -Cw(n)*Rs(n)
            z(2) = -Sw(n)*Rs(n)
            z(3) = Rc(n)
         else
            z(1) = 0.0_10
            z(2) = 0.0_10
            z(3) = 1.0_10
         endif
         call CHNCOR(Pstlat(1,n),z)
      endif
      if(Lsite(4,n).gt.0 .or. norm.eq.2) then
         z(1) = (Cv(n)*Xb00(1,n) + Sv(n)*Xb00(2,n))*tfct/Rsite0(n)
         z(2) = (Sv(n)*Xb00(1,n) - Cv(n)*Xb00(2,n))*tfct/Rsite0(n)
         z(3) = Rs(n)*tfct/Rsite0(n)
         call CHNCOR(Pstup(1,n),z)
      endif
      if(Lsite(5,n).gt.0 .or. norm.eq.2) then
         z(1) = (-Cv(n)*Slong0(n) + Sv(n)*Clong0(n))*tfct
         z(2) = -(Sv(n)*Slong0(n) + Cv(n)*Clong0(n))*tfct
         z(3) = 0._10
         call CHNCOR(Pstwes(1,n),z)
      endif
      if(Lsite(6,n).gt.0 .or. norm.eq.2) then
         z(1) = -(Cv(n)*Clong0(n) + Sv(n)*Slong0(n))*Rs0(n)*tfct
     .    /Rsite0(n)
         z(2) = -(Sv(n)*Clong0(n) - Cv(n)*Slong0(n))*Rs0(n)*tfct
     .    /Rsite0(n)
         z(3) = Rc0(n)*tfct/Rsite0(n)
         call CHNCOR(Pstnor(1,n),z)
      endif

c always get spherical partials for coordinate transformations
      if(norm.eq.2) then
         do i=1,3
            Dxdrad(i,n) = Pstrad(i,n)*Cyl1(n) + Pstlat(i,n)*Cyl2(n)
            Dxdlat(i,n) = Pstrad(i,n)*Cyl3(n) + Pstlat(i,n)*Cyl4(n)
            Dxdlon(i,n) = Pstlon(i,n)
         end do
      endif

      if(nvel.le.0) return
c
c*  start=500
c determine site velocity
      sidvlt = Sidvel
      if(Kindnp.gt.0) sidvlt = Sidvel + Pc(2)
      do i = 1, 3
         j = i + 3
         Xsite(j,n) = (-Nutpr(1,i)*Y(2,n) + Nutpr(2,i)*Y(1,n))*sidvlt
      end do
      if(Kindnp.gt.0) then
         do i = 1, 3
            j = i + 3
            Xsite(j,n) = Xsite(j,n)
     .       + (Dnutpr(1,i)*Y(1,n) + Dnutpr(2,i)*Y(2,n)
     .       + Dnutpr(3,i)*Y(3,n))
         end do
      endif
      if(mod(Jct(6)/8,2).ne.0) then
         if(Line.gt.56) call OBSPAG
         do i = 1, 3
            z(i) = Xsite(i + 3,n)/Aultvl
         end do
         write(Iout,50) type(2),n,sidtp,(Xsite(i,n),i = 4,6),z
         Line = Line + 2
      endif
c
c*  start=700
c determine partial derivatives of velocity
      if(Lsite(1,n).gt.0 .or. norm.eq.2) then
         if(Ksite(n).lt.0) then
            z(1) = -Sw(n)*sidvlt
            z(2) = Cw(n)*sidvlt
            z(3) = 0._10
            call CHNCOR(Pstrad(4,n),z)
         else
            do i = 4, 6
               Pstrad(i,n) = Xsite(i,n)/Rsite(n)
            end do
         endif
      endif
      if(Lsite(2,n).gt.0 .or. norm.eq.2) then
         do i = 1, 3
            j = i + 3
            Pstlon(j,n) = (Nutpr(1,i)*Y(1,n) + Nutpr(2,i)*Y(2,n))
     .       *sidvlt
         end do
      endif
      if(Lsite(3,n).gt.0 .or. norm.eq.2) then
         if(Ksite(n).lt.0) then
            do j = 4, 6
               Pstlat(j,n) = 0.0_10
            end do
         else
            rsv = Rs(n)*sidvlt
            do i = 1, 3
               j = i + 3
               Pstlat(j,n) = (Nutpr(1,i)*Sw(n)-Nutpr(2,i)*Cw(n))*rsv
            end do
         endif
      endif
      if(Lsite(4,n).gt.0 .or. norm.eq.2) then
         z(1) = (-Sv(n)*Xb00(1,n) + Cv(n)*Xb00(2,n))*sidvlt*tfct
     .    /Rsite0(n)
         z(2) = (Cv(n)*Xb00(1,n) + Sv(n)*Xb00(2,n))*sidvlt*tfct
     .    /Rsite0(n)
         z(3) = 0._10
         call CHNCOR(Pstup(4,n),z)
      endif
      if(Lsite(5,n).gt.0 .or. norm.eq.2) then
         z(1) = (Sv(n)*Slong0(n) + Cv(n)*Clong0(n))*sidvlt*tfct
         z(2) = (-Cv(n)*Slong0(n) + Sv(n)*Clong0(n))*sidvlt*tfct
         z(3) = 0._10
         call CHNCOR(Pstwes(4,n),z)
      endif
      if(Lsite(6,n).gt.0 .or. norm.eq.2) then
         z(1) = (Sv(n)*Clong0(n) - Cv(n)*Slong0(n))*Rs0(n)*sidvlt
     .    *tfct/Rsite0(n)
         z(2) = -(Cv(n)*Clong0(n) + Sv(n)*Slong0(n))*Rs0(n)*sidvlt
     .    *tfct/Rsite0(n)
         z(3) = 0._10
         call CHNCOR(Pstnor(4,n),z)
      endif

c always get spherical partials for coordinate transformations
      if(norm.eq.2) then
         do i=4,6
            Dxdrad(i,n) = Pstrad(i,n)*Cyl1(n) + Pstlat(i,n)*Cyl2(n)
            Dxdlat(i,n) = Pstrad(i,n)*Cyl3(n) + Pstlat(i,n)*Cyl4(n)
            Dxdlon(i,n) = Pstlon(i,n)
         end do
      endif

      if(nvel.lt.2) return
c
c*  start=1000
c determine site acceleration
      sidvl2 = sidvlt**2
      do i = 1, 3
         j = i + 6
         Xsite(j,n) = (-Nutpr(1,i)*Y(1,n) - Nutpr(2,i)*Y(2,n))*sidvl2
      end do

      if(mod(Jct(6)/16,2).ne.0) then
         if(Line.gt.56) call OBSPAG
         do i = 1, 3
            z(i) = Xsite(i + 6,n)/Auacc
         end do
         write(Iout,50) type(3),n,sidtp,(Xsite(i,n),i = 7,9),z
         Line = Line + 2
      endif

      if(nvel.lt.3) return
c
c*  start=1200
c determine site jerk
      sidvl2 = sidvl2*sidvlt
      do i = 1, 3
         Sitjrk(i,n) = (Nutpr(1,i)*Y(2,n) - Nutpr(2,i)*Y(1,n))*sidvl2
      end do
      return
      end
