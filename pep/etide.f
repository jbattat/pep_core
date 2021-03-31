      subroutine ETIDE(etide1,nvel,ns,kobj)
 
      implicit none

c
c          rw king and rj cappallo    august 1978   subroutine etide
c           interface between vetide and pep.  etide sets up
c           various vectors, calls vetide to calculate the solid body
c           tide raised by the sun and moon on the and alters the
c           appropriate site coordinates .
c
c
c parameters
      integer*4 etide1,kobj,ns,nvel
c
c           etide1 = 0  first call of etide for a given site:
c                       calculate tides and correct site coordinates
c                  = 1  calls of etide on successive iterations for a
c                       given site:  correct site coordinates but do
c                       not recalculate tides
c
c           nvel = 0  do not calculate velocities
c                  1  calculate velocities
c
c           ns = 1  first site
c                2  second site
c
c           kobj = 1  observed object is nplnt0
c                = 2  observed object is nplnt2
c
c
c         etide requires that the relative coordinates of the
c         earth, sun, and moon be available in /coord/ arrays
c         xem and xm

c array dimensions
      include 'globdefs.inc'
c
c        commons
      include 'comdateq.inc'
      include 'coord.inc'
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'obscrd.inc'
      include 'param.inc'
      include 'sitcrd.inc'
      character*4 sitf4(2,2)
      equivalence (Sitf(1),sitf4(1,1))
      include 'tidal.inc'
c
c
c variables internal to this routine
      real*10 dellat,dellon,delrad,dlat2,dlatm,dlats,dlon2,dlonm,dlons,
     . drad2,dradm,drads,gmsun,r
      integer   i,icem,ices,iers,index,ipt,iptln,j,jct10
      logical*1 init/.false./
      real*10 g(2),hlov(2),llov(2),tlag(2),xemtid(6),
     . xestid(6),dxdhem(6),dxdlem(6),dxdhes(6),dxdles(6),
     . xetid(6),xml(6),xsl(6),xiers(6),xsbod(6,2),menrad,fract,tt
c
c table of site-dependent love numbers and lag angles
c default love numbers are theoretical values from alsop
c and kuo, who used bullen's earth model.
      character*4 sittab(5) /'HAYS', 'OVRO', 'NRAO', 'TEXL', 'HAWL'/
      real*4 hlttab(3,5)/
     1    0.618,    0.088,        0.,
     2    0.618,    0.088,        0.,
     3    0.618,    0.088,        0.,
     4    0.618,    0.088,        0.,
     5    0.618,    0.088,        0./
      integer*4 maxshl/5/

c nominal second degree and third degree love numbers and shida numbers
c for IERS tide model
      real*10 h20/0.6078_10/,l20/0.0847_10/,h3/0.292_10/,l3/0.015_10/
c
c----------------------------------------------------------------------
c
c
c           check etide1
      if(etide1.le.0) then
c
c           calculations performed only once
c
c     these calculations might be moved to eshape in the future
c     with values stored in /tidal/
c
         if(.not. (init)) then
 
            gmsun  = Gauss**2*Aultsc**3/Secday**2
            Gmerth = gmsun*Mass(3)*(1._10 - Mass(10))
            Gmmoon = gmsun*Mass(3)*Mass(10)
            menrad = 6378.1366_10/Ltvel
            jct10  = Jct(10)
            icem   = mod(jct10,2)
            ices   = mod(jct10/2,2)
            iptln  = mod(jct10/4,2)
            ipt    = mod(jct10/8,2)
            iers   = mod(jct10/64,2)
            init   = .true.
         endif
c
c
c determine earth site quantities once per series
         xsbod(1,ns)=Xb0(1,ns)
         xsbod(2,ns)=-Xb0(2,ns)
         xsbod(3,ns)=Rs(ns)
         g(ns) = Gmerth/Rsite(ns)**2
 
c determine love numbers and time lag
         hlov(ns) = 0.618
         llov(ns) = 0.088
         tlag(ns) = 0.
         do j = 1, maxshl
            if(sitf4(1,ns).eq.sittab(j)) then
               hlov(ns)= hlttab(1,j)
               llov(ns)= hlttab(2,j)
               tlag(ns)= hlttab(3,j)
               goto 50
            endif
         end do
c
c set love numbers equal to input parameters if non-zero
   50    if(iers.ne.0) then
            hlov(ns)=h20
            llov(ns)=l20
         endif
         if(Econd(9).ne.0._10) hlov(ns)  = Econd(9)
         if(Econd(10).ne.0._10) llov(ns) = Econd(10)
         if(Econd(11).ne.0._10) tlag(ns) = Econd(11)
         if(iers.ne.0) then
            hlov(ns)=hlov(ns)-0.0006_10*(1._10-1.5_10*Cphi(ns)**2)
            llov(ns)=llov(ns)+0.0002_10*(1._10-1.5_10*Cphi(ns)**2)  
         endif
c
c initialization once per observation
         index = 3
         if(nvel.gt.0) index = 6
         do i = 1, index
            xemtid(i) = 0._10
            xestid(i) = 0._10
            dxdhem(i) = 0._10
            dxdlem(i) = 0._10
            dxdhes(i) = 0._10
            dxdles(i) = 0._10
         end do
         if(ipt.ne.0) then
            dradm = 0._10
            dlatm = 0._10
            dlonm = 0._10
            drads = 0._10
            dlats = 0._10
            dlons = 0._10
         endif
c
c calculate vector tide raised by moon
         if(icem.ne.0) then
            do i = 1, 3
               xml(i) = Xm(i,ns)*Mnltsc
               if(nvel.gt.0) xml(i + 3) = Xm(i + 3,ns)
     .             *Mnltsc/Secday
            end do
            if(iers.eq.0) then
               call VETIDE(nvel,Gmmoon,g(ns),Cphi(ns),hlov(ns),
     .          llov(ns),xml,Xsite(1,ns),Dxdrad(1,ns),Dxdlat(1,ns),
     .          Dxdlon(1,ns),xemtid,dxdhem,dxdlem,dradm,dlatm,dlonm)
            else
               call IERSTIDE(nvel,Gmmoon,Gmerth,Cphi(ns),hlov(ns),h3,
     .          llov(ns),l3,xml,xsbod(1,ns),menrad,
     .          Dxdrad(1,ns),Dxdlat(1,ns),Dxdlon(1,ns),
     .          xemtid,dxdhem,dxdlem,dradm,dlatm,dlonm)
            endif
         endif
c
c calculate vector tide raised  by sun
         if(ices.ne.0) then
            do i = 1,3
               xsl(i) = -(Xem(i,ns) - Xm(i,ns)*Mnfct)*Aultsc
               if(nvel.gt.0) xsl(i + 3)
     .             = -(Xem(i+3,ns) - Xm(i+3,ns)*Mnfct)*Aultvl
            end do
            if(iers.eq.0) then
               call VETIDE(nvel,gmsun,g(ns),Cphi(ns),hlov(ns),llov(ns)
     .          ,xsl,Xsite(1,ns),Dxdrad(1,ns),Dxdlat(1,ns),Dxdlon(1,ns),
     .          xestid,dxdhes,dxdles,drads,dlats,dlons)
            else
               call IERSTIDE(nvel,gmsun,Gmerth,Cphi(ns),hlov(ns),h3,
     .          llov(ns),l3,xsl,xsbod(1,ns),menrad,
     .          Dxdrad(1,ns),Dxdlat(1,ns),Dxdlon(1,ns),
     .          xestid,dxdhes,dxdles,drads,dlats,dlons)
            endif
         endif
c
c add together lunar and solar components of tide
         do i = 1,index
            xetid(i) = xemtid(i) + xestid(i)
         end do
         if(iers.ne.0) then
            tt=((Jdsvem(ns)-2451545)+Frsvem(ns)-0.5_10)/36525._10
            fract=Frsvem(ns)-(Ctat+Atuts)/Secday
            call IERSTID2(xsbod(1,ns),Cphi(ns),fract,tt,
     .       Dxdrad(1,ns),Dxdlat(1,ns),Dxdlon(1,ns),
     .       xiers,drad2,dlat2,dlon2)
            do i=1,3
               xetid(i)=xetid(i)+xiers(i)
            end do
         endif
      endif
c
c change cartesian coordinates of site vectors
      do i = 1,index
         Xsite(i,ns) = Xsite(i,ns) + xetid(i)
 
c lag angle coded only for delays
         if(i.le.3) Xsite(i,ns) = Xsite(i,ns) + tlag(ns)
     .                                 *xetid(i + 3)
      end do
      if(etide1.le.0) then
c
c
c partials w.r.t. love-number scale factors and time lag
         do i = 1,index
            Dxdhe(i,ns) = (dxdhem(i) + dxdhes(i))
            Dxdle(i,ns) = (dxdlem(i) + dxdles(i))
            Dxdte(i,ns) = 0._10
            if(i.le.3) then
               Dxdhe(i,ns) = Dxdhe(i,ns) + tlag(ns)
     .                        *(dxdhem(i+3) + dxdhes(i+3))
               Dxdle(i,ns) = Dxdle(i,ns) + tlag(ns)
     .                        *(dxdlem(i+3) + dxdles(i+3))
               Dxdte(i,ns) = xetid(i + 3)
            endif
         end do
c
c
c optional printout of love no.s and tidal components
         if(Nk1.le.0) then
            if(iptln.ne.0) then
               if(Line.gt.57) call OBSPAG
               write(Iout,60) ns,hlov(ns),llov(ns),tlag(ns)
               Line = Line + 1
   60          format(' LOVE NUMBERS FOR SITE', i2, ':  H=', f8.4,
     .                ' L=', f8.4, ' TLAG=', e10.4, ' SEC')
            endif
         endif
         if(ipt.ne.0) then
            delrad = (dradm + drads)*hlov(ns)
            dellat = (dlatm + dlats)*llov(ns)
            dellon = (dlonm + dlons)*llov(ns)
            if(iers.ne.0) then
               delrad= delrad+drad2
               dellat= dellat+dlat2
               dellon= dellon+dlon2
            endif
            if(Line.gt.58) call OBSPAG
            write(Iout,80) ns,delrad,dellat,dellon
            Line = Line + 1
   80       format(' EARTH TIDE COMPONENTS FOR SITE', i2, ':  RAD=',
     .             1pd14.6, '  LAT=', 1pd14.6, '  LONG=', 1pd14.6,
     .             '   (SEC)')
         endif
      endif
c
c
      return
      end
