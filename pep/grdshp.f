      subroutine GRDSHP(lat,long,tmdly,topsec)
 
      implicit none
c
c        subroutine grdshp  r.b. goldstein and s. brody  feb 78
c        planet altitude grid for topography - local shape model
c        grdshp called from radpl if nshp=2
c        altitude at subradar point is determined by a 4 point int-
c        erpolation scheme using the grid of altitudes stored in
c        /shphar/. planetary flattening at the observation latitude
c        is also taken into account. the radar time delay (tmdly)
c        is then increased or decreased accordingly and returned. (see
c        bodred for shape model explanation).
c
c        also compute partials at entry point grdpar which is
c        called from shppar

c arguments
      real*10 lat,long,tmdly,topsec

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdat.inc'
      include 'inodta.inc'
      include 'param.inc'
      include 'shpcom.inc'
      include 'ltrapx.inc'
      integer*2 kind
      equivalence (kind,Numpar)
      integer*2 lshape(1000)
      equivalence (lshape,Lszhar)
      include 'mtrapx.inc'
      include 'shphar.inc'
      real*4    grid(1000)
      equivalence (grid,Zone)
      real*4    tlat(2),tlon(2),tlatin,tlonin
      equivalence (Scntrl,tlat),(Scntrl(3),tlon),
     .            (Scntrl(5),tlatin),(Scntrl(6),tlonin),
     .            (Scntrl(7),latdim),(Scntrl(8),londim)
c        tlat,tlon - latitude and longitude ranges of grid
c        tlatin,tlonin - latitude and long. grid spacing
c        latdim,londim - grid dimensions
c
c        locals
      real*10 alt,flat,topog
      integer   i,iflag,ilat,ilon,ipa,ipb,ita,itb,latdim,
     .          lim,londim,lst,n
      real*10 lat4,long4,h(4),p(4),tb0,t0a,pb0,p0a,ta,pb
      integer*2 lobstr(4),lps(4)
c
c arithmetic statement function igrid computes 1-dim pointer
      integer*4 IGRID
      IGRID(ilat,ilon) = (ilat-1)*londim + ilon
c
c
c
c        the flag "zero" is set to true only if both grid and lgrid
c        are all zero. even if grid is all zero, the partials are
c        non zero and values computed in grdshp are needed when
c        grdpar is called
c
      Shpnit = .true.
      call FLATC(flat,lat)
c
c determine topography contribution - skip if all coef. are zero
c
      topog = 0._10
      if(.not.Shpzer) then
c
c first demod. lat. and long. and check if they are in
c range of topography grid.
c
         lat4  = lat
         long4 = long
         do while (long4.lt.0._10)
            long4 = long4+360._10
         end do
         do while (long4.gt.360._10)
            long4 = long4-360._10
         end do
         if(lat4.lt.-90._10 .or. lat4.gt.90._10) lat4 = MOD(lat4,90._10)
         if(long4.lt.tlon(1) .or. long4.gt.tlon(2)) goto 200
         if(lat4.lt.tlat(1) .or. lat4.gt.tlat(2)) goto 200
c
c
c obtain long. indicies for nearest grid points and
c quantities for interpolation
         ta  = (long4 - tlon(1))/tlonin + 1._10
         ita = ta
         itb = ita + 1
 
c check for rollover in longitude
         if(itb.gt.londim) itb = 1
         t0a = (ta - ita)*tlonin
         tb0 = tlonin - t0a
c
c        obtain lat. indices for nearest grid points and quantities for
c        interpolation. this differs from the long. section because
c        the maximum latitude strip corresponds to ipb=1 (i.e. the
c        top row is maximum latitude as on a map).
c        also there is no rollover permitted.
         pb  = (tlat(2) - lat4)/tlatin + 1._10
         ipb = pb
         ipa = ipb + 1
         pb0 = (pb - ipb)*tlatin
         p0a = tlatin - pb0
c
c determine pointers into grid array
c
         lobstr(1) = IGRID(ipb,ita)
         lobstr(2) = IGRID(ipb,itb)
         lobstr(3) = IGRID(ipa,ita)
         lobstr(4) = IGRID(ipa,itb)
         do i = 1,4
            h(i) = grid(lobstr(i))
         end do
 
         p(1)  = tb0*pb0
         p(2)  = tb0*p0a
         p(3)  = t0a*p0a
         p(4)  = t0a*pb0
         topog = Grdf1*(h(1)*p(1) + h(2)*p(2) + h(3)*p(3) + h(4)*p(4))
      endif
c
c convert altitude in kilometers to light seconds
      alt = topog/Ltvel + flat
 
c subtract 2-way additional time delay from tmdly
      topsec = alt + alt
      tmdly  = tmdly - topsec
      return
c
c subradar point outside grid - stop
c
  200 write(Iout,300) lat,long,tlat,tlon
  300 format(
     .'  GRDSHP..ERROR..PROBABLY DUE TO SUBRADAR PT. OUTSIDE OF TOPOGRAP
     .HY GRID RANGE'/' LAT,LON,TLAT,TLON:'/6E20.5)
      call SUICID(' GRDSHP ERROR   ',4)
c
c entry point for partials calculation
c
      entry GRDPAR
c
c set lshobs from lobstr and lshape. (see if those
c points used for this obs. are being adjusted)
c
      Lnshob = 0
      Lixshp = kind + 1
      if(Ncode.ne.3 .and. Shpnit) then
         do i = 1,4
            Lshobs(i) = 0
            if(lshape(lobstr(i)).ne.0) then
               Lnshob = Lnshob + 1
               Lshobs(Lnshob) = lobstr(i)
               lps(Lnshob)    = i
            endif
         end do
         if(Lnshob.ne.0) then
c
c now compare lshobs to mshobs and compute or copy
c partials.
            call PCOPS(1,'SHOB',Iabs1)
            lim = Lnshob
            n   = 0
            do while( .true. )
               iflag = 1
               call PCOPY(n,lim,iflag,1,Lshobs,Mshobs)
               if(iflag.gt.0) go to 100
               Deriv(kind,1) = Grdf2*p(lps(n))
               Deriv(kind,2) = 0._10
               if(n.ge.lim) go to 100
            end do
         endif
      endif
c
c        make sure there are 4 partials on tape so that same no. of
c        partials for every point in series as required by analyze
c        link.
c        do not insert if there are no grid partials at all
c
  100 if(Lngd.le.0) return
      if(Lnshob.ne.4) then
         lst = Lnshob+1
         do i = lst,4
            kind = kind+1
            Deriv(kind,1) = 0._10
            Deriv(kind,2) = 0._10
         end do
      endif
 
      return
      end
