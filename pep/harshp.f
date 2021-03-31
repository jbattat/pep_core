      subroutine HARSHP(lat,long,tmdly,topsec,lslope,slangl,
     .                  lscale)
 
      implicit none

c
c r.j. cappallo   august 1970   subroutine harshp
c
c     harshp accepts as arguments latitude and longitude in degrees,
c     time delay or a dummy variable, depending on whether lslope=1, or
c     lslope=2,respectively.if it is a slope measurement then slangl is
c     the angle form a line of longitued to the doppler axis, in radians
c     with the direction of increasing longitude being positive.
c     lscale=1 designates the legendre functions to be normalized.
c     tmdly is either incremented and returned, or contains the slope at
c     that point, i.e.- the tangent to the spherical harmonic expansion.

c parameters
      real*10 lat, long, tmdly, topsec, slangl
      integer*4 lscale, lslope

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'empcnd.inc'
      include 'funcon.inc'
      include 'number.inc'
      include 'shpcom.inc'
 
c fourier trigonometric quantities
      real*10 fa(20),fb(20),fc(20),fd(20),fap(20),fcp(20),fapp,
     .          fbpp
      equivalence (fa,Gleg),(fb,Gleg(21)),(fc,Gleg(41))
     .            , (fd,Gleg(61)),(fap,Gleg(81)),
     .            (fcp,Gleg(101)),(fapp,Gleg1),
     .            (fbpp,Gleg1(2))
      include 'shphar.inc'
      real*4    lolat, hilat, tlat(2)
      real*10 ovely(122)
      equivalence (tlat,lolat),(tlat(2),hilat)
      equivalence (Zone,ovely),(tlat,Scntrl)
 
c fourier expansion coefficients
      real*10 a(20),b(20),c(20),d(20),ap(20),cp(20),app,bpp
      equivalence (ovely,a),(ovely(21),b),(ovely(41),c),
     .            (ovely(61),d),(ovely(81),ap),
     .            (ovely(101),cp),(ovely(121),app),
     .            (ovely(122),bpp)

c local variables
      real*10 flat,height,pc11,ps11,pz1,z,zz
      integer   m, mn, n, nt4, nz4
c
c 2-d fourier work areas
      real*10 smlong(20),cmlong(20),snlat,cnlat
c branch for 2-dimensional fourier series shape (nshp=1)
c as opposed to spherical harmonic expansion (nshp=0)
c
      if(Nshp.eq.1) goto 200
 
      Sangle = slangl
      Lslop  = lslope
      z   = COS((9.E1_10-lat)*Convd)
      zz  = SQRT(1.0_10 - z*z)
      nz4 = Nz
      nt4 = Nt
      call LEGNDR(z,zz,nz4,nt4,Leg(2),Leg1(2),Gleg(2),Gleg1(2))
 
c low order polynomials not evaluated in legndr, must do this here
      Leg(1)   = z
      Leg1(1)  = 1.0_10
      Gleg(1)  = zz
      Gleg1(1) = -z/zz
      if(lscale.eq.1) then
 
c normalize legendre polynomials and functions
         call SCALE(nz4,nt4,Leg,Gleg)
         call SCALE(nz4,nt4,Leg1,Gleg1)
      endif
 
c compute longitudinal trig. quantities and multiples
      Sinmln(1) = SIN(long*Convd)
      Cosmln(1) = COS(long*Convd)
      do m = 2, Nt
         Cosmln(m) = Cosmln(m - 1)*Cosmln(1) - Sinmln(m - 1)*Sinmln(1)
         Sinmln(m) = Sinmln(m - 1)*Cosmln(1) + Cosmln(m - 1)*Sinmln(1)
      end do
      pz1  = Pcond(8,Klan)
      pc11 = Pcond(9,Klan)
      ps11 = Pcond(10,Klan)
 
c see if time delay or slope
      if(lslope.eq.2) then
 
c slope coding to go here
         call SUICID(
     .' ABILITY TO USE SLOPE DATA NOT YET PROGRAMMED,  STOP IN HARSHP 10
     .0  ', 17)
         goto 200
      else
 
c delay data
         height = 0.0_10
 
c zonal terms
         if(Nz.ge.1) then
            height = height + pz1*Leg(1)
            if(Nz.ge.2) then
               do n = 2, Nz
                  height = height + Zone(n - 1)*Leg(n)
               end do
            endif
         endif
 
c tesseral and sectoral terms
         if(Nt.ge.1) then
            height = height + Gleg(1)*(pc11*Cosmln(1) + ps11*Sinmln(1))
            if(Nt.ge.2) then
               mn = 0
               do n = 2, Nt
                  do m = 1, n
                     mn     = mn + 1
                     height = height + Gleg(mn + 1)
     .                        *(Ctess(mn)*Cosmln(m) + Stess(mn)
     .                        *Sinmln(m))
                  end do
               end do
            endif
         endif
      endif
  100 topsec = height + height
      tmdly  = tmdly - topsec
      return
c ad hoc coding to enable pep to use a 2-dimensional fourier series
c to represent topography near a planetary equator
c tlat(1) and tlat(2)=lowlat,hilat are limits of latitudinal expa
  200 z = (lat - lolat)/(hilat - lolat)*Twopi
 
c compute fourier trig. quantities
      snlat     = SIN(z)
      cnlat     = COS(z)
      smlong(1) = SIN(long*Convd)
      cmlong(1) = COS(long*Convd)
      do m = 2, 20
         cmlong(m) = cmlong(m - 1)*cmlong(1) - smlong(m - 1)*smlong(1)
         smlong(m) = smlong(m - 1)*cmlong(1) + cmlong(m - 1)*smlong(1)
      end do
 
c calculate value of theoretical and store fourier quantities
      fapp   = cnlat
      fbpp   = snlat
      height = app*fapp + bpp*fbpp
 
      call FLATC(flat,lat)
      height = height + flat
 
      do m = 1, 20
         fa(m)  = cmlong(m)*cnlat
         fb(m)  = cmlong(m)*snlat
         fc(m)  = smlong(m)*cnlat
         fd(m)  = smlong(m)*snlat
         fap(m) = cmlong(m)
         fcp(m) = smlong(m)
         height = height + a(m)*fa(m) + b(m)*fb(m) + c(m)*fc(m) + d(m)
     .            *fd(m) + ap(m)*fap(m) + cp(m)*fcp(m)
      end do
      goto 100
      end
