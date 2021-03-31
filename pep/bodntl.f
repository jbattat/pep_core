      subroutine BODNTL
 
      implicit none
 
c
c m.e.ash    aug 1971    subroutine bodntl
c updated feb. 1980   kcl: tcon added
c initialize quantities for specific body

c array dimensions
      include 'globdefs.inc'
c
      include 'fcntrlx.inc'
      common/WRKCOM/ Cond(6),Con(u_nmbod-6),Con1(12),Tcon(30),Sec,
     . Int1,Int2,Eps(6),Name,Jd1,Jd0,Jd2,Itape,L(u_nmbod),I4fill,
     . Nshp,Scntrl(9),Zhar(u_mxzon-1),Char((u_mxtes*(u_mxtes+1))/2-1),
     . Shar((u_mxtes*(u_mxtes+1))/2-1),Lzh(u_mxzon-1),
     . Lch((u_mxtes*(u_mxtes+1))/2-1),Lsh((u_mxtes*(u_mxtes+1))/2-1),
     . Denptr,K(u_nmprm),Int,Nplnt,Ncentr,Nzone,Ntess,Ihr,Imin,
     . Kk(100),Icnd
      real*10 Cond,Con,Con1,Tcon,Sec,Zhar,Char,Shar
      real*4    Eps,Scntrl
      integer   Int1,Int2,Jd1,Jd0,Jd2,Itape,I4fill,Nshp
      character*8 Name
      integer*2 L,Lzh,Lch,Lsh,I2fill,Denptr,K,Int,Nplnt,Ncentr,
     .          Nzone,Ntess,Ihr,Imin,Kk,Icnd
      real*4    shape(1000),tlat(2),tlon(2),tlatin,tlonin
      integer*2 lshape(1000)
      equivalence (Zhar,shape),(Lzh,lshape)
      equivalence (tlat,Scntrl),(tlon,Scntrl(3)),
     .            (tlatin,Scntrl(5)),(tlonin,Scntrl(6))
 
      real*10 prad(10), ratpoe(10)
      character*8 pname(10)/'MERCURY ',' VENUS  ',' EMBARY ','  MARS  ',
     .     'JUPITER ', ' SATURN ', ' URANUS ', 'NEPTUNE ', ' PLUTO  ',
     .     '  MOON  '/
c RATPOE values (ratio of gravitational to total energy) are
c uniform density values 3GM/5Rc**2 for Mercury, Venus, Mars,
c Pluto, and Moon.  Earth is from PRL 36, March 15, 1976, p555.
c The gas giants are based on exponential density profiles with
c scale height set to give the right moment of inertia (from Allen).
      data prad(1)/ 2438.0_10/,  ratpoe(1)/-0.60E-10_10/
      data prad(2)/ 6050.0_10/,  ratpoe(2)/-3.6E-10_10/
      data prad(3)/ 6378.166_10/,ratpoe(3)/-4.6E-10_10/
      data prad(4)/ 3394.0_10/,  ratpoe(4)/-0.84E-10_10/
      data prad(5)/71350.0_10/,  ratpoe(5)/-5.6E-7_10/
      data prad(6)/60400.0_10/,  ratpoe(6)/-2.6E-7_10/
      data prad(7)/23800.0_10/,  ratpoe(7)/-9.1E-8_10/
      data prad(8)/22200.0_10/,  ratpoe(8)/-7.1E-8_10/
      data prad(9)/ 3000.0_10/,  ratpoe(9)/-0.022E-10_10/
      data prad(10)/1738.0_10/,  ratpoe(10)/-0.19E-10_10/
      integer*2 intpl(10)/2, 4, 4, 4, 20, 20, 40, 40, 40, -1/
      integer*2 k87(10)/-1, 1, -1, 1, 4, 4, 4, 4, 4, -3/
      integer*2 lowlat(10)/-15, -15, 0, -30, 6*0/
      integer*2 hilat(10)/15, 15, 0, 30, 6*0/
c
c initialize specific body harmonics
      call HARNTL(Nplnt, Zhar, Char, Shar, Nzone, Ntess)
c
c standard setup for every body
      if(Nplnt.gt.0) then
         Itape = Nplnt + 10
         if(Nplnt.gt.10) goto 100
         Name     = pname(Nplnt)
         Int      = intpl(Nplnt)
         K(87)    = k87(Nplnt)
         Con(1)   = prad(Nplnt)
         Con1(10) = ratpoe(Nplnt)
c
c initialize quantities for embary and earth
         if(Nplnt.eq.3) then
            Con(2) = 1.0_10/298.3_10
            K(40)  = 1
            goto 200
c
c initialize quantities for moon
         else if(Nplnt.eq.10) then
            Ncentr = 3
            Kk(84) = 2
            goto 200
         endif
      endif
c
c initialize quantities for earth rotation
      if(Nplnt.eq.-3) then
         if(Jct(29).gt.0) then
            Con(19)=0.94_10
            Con(20)=0.94_10
         endif
         goto 200
c
c initialize quantities for moon rotation
      else if(Nplnt.eq.-10) then
         Con(3) = .630E-3_10
         Con(4) = .228E-3_10
         goto 200
c
c initialize quantities for planet rotation or shape
      else if(Nplnt.lt.0) then
         if(Nshp.ge.1) then
 
c additional planet shape models
            if(Nshp.le.1) then
c          for nshp=1,shape model is two dimensional
c          fourier series.(see comments in bodred)
c          internal controls necessary are:
c               tlat(1)=lower latitude on planet
c               tlat(2)=upper latitude on planet
c               nzone=20
c               ntess=10
c        nzone and ntess are needed for shppar to calc. four. partials
c
               tlat(1) = lowlat(-Nplnt)
               tlat(2) = hilat(-Nplnt)
               Nzone   = 20
               Ntess   = 10
            else if(Nshp.le.2) then
c        nshp=2, local shape model - altitude grid
c        grid points for appropriate range on planet represent
c        height in kilometers of the surface at the
c        lat and longitude corresponding to the grid point
c        (see bodred)
c        set defaults for lat & long range, grid spacing,
c        grid dimensions and total # of grid points
               tlat(1) = lowlat(-Nplnt)
               tlat(2) = hilat(-Nplnt)
               tlon(1) = 0.0
               tlon(2) = 360.0
               tlatin  = 10.0
 
c initialization for additional shape models goes here
               tlonin = 10.0
            endif
c planet rotation or spherical harmonic shape model
c no a priori initialization as yet
         endif
         goto 200
      endif
c
c initialize quantities for planets
  100 if(Nplnt.gt.30) then
c
c initialize quantities for artificial space probes
         Con1(12) = .03_10
         if(Ncentr.eq.0) Icnd = -1
      endif
c
c final initialization
  200 Kk(1) = Ncentr
      return
      end
