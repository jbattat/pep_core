      subroutine ANGOUT(iangbf, angbuf, bessel, sidtm1)
 
      implicit none

c
c m.ash   may 1975   subroutine angout
c calculate angles and range and their rates of change
c
c arguments
      integer*2 angbuf(257), iangbf
      real*10 bessel(3, 3)
      real*10 sidtm1
c     sidtm1 = greenwich true sideral time in radians at receive time
c              epoch. satellite at receive time if delit=1.e6, otherwise
c              at retarded time.
c
c jct(30) = 0 no e,f,g (+edot,fdot,gdot) tape output
c             earth centered rotating coordinates. e in true equator of
c             date towards greenwich meridian, g to north, f completes
c             right hand system
c jct(30).ne.0 such tape output on data set iabs(jct(30))
c              unless jct(40) says there is angle output instead
c        .lt.0 no frintout of e,f,g
c        .gt.0 printout of e,f,g at end of dummy observation series
c jct(31) = space defense center satellite number for e,f,g output
c jct(32) = 0 no binary e,f,g output in subroutine efgout for fitting
c             polynomials in subroutine efgfit
c jct(32).ne.0 such binary output on data set iabs(jct(32)) if
c              jct(30).ne.0
c        .gt.0 jct(30) output with jct(32) output
c        .lt.0 jct(30) output is supressed
c
c jct(40) = 0 no angle output for optical observatories
c jct(40).ne.0 such tape output on data set iabs(jct(30))
c              in format compatable with data general nova computer
c jct(40).gt.0 printout of angles, data set jct(40) used as intermediate
c              ebcdic buffer
c jct(40).le.0 no printout of angles
c jct(41) = 0 angle output referred to true equinox and equator of date
c             for jct(40) output
c jct(41).gt.0 angle output referred to mean equinox and equator of
c              jct(41). for eample, jct(41)=76 denotes output referred
c              to mean equinox and equator of 1976.0
c jct(42) = 0 no refraction correction made in computing jct(40) angles
c jct(42) = 1 refraction correction made in computing jct(42) angles
c jct(43) =   file number for jct(40) angle output (0 is file 1, 1 is
c             file 2, etc, for data general nova control word at end of
c             514 character record)
c jct(47) = 0 usual jct(40) print and tape output
c jct(47) > 0 jct(40) tape output has satellite radius instead
c             of doppler rate, print output has special les-8/9
c             format with 6 frequency doppler
c jct(47) = 8 les-8 output
c jct(47) = 9 les-9 output
c
c array dimensions
      include 'globdefs.inc'

c commons
      include 'coord.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'namtim.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 ctrecf, fdev, reflct
      equivalence (ctrecf, Dstf), (fdev, Dstf(10)), (reflct, Dstf(4))
      real*10 tmdly0, dop0
      equivalence (Result(1), tmdly0), (Result(2), dop0)
      include 'param.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'

c external functions 
      real*10 DOT

c local variables
      real*10 qxo1, qxo2, az, el, pxn, cot, coe, soe, q1, q2, ysat(6),
     .          xsat(6),drange,ddrang,d, qxn,xmerid(3),xop(3),xmerop(3)
      real*10 rzsat, zsat(3)
      real*10 sidtml, rsat2, rsat, decl, ra, ha, ddecl, dra, dha
      character*1 minus(2)/' ', '-'/
      real*10 xsun(3), psun, rsun2, de(3), dde, angsun
      real*4    eloff, elref, eltrue
      integer*4 i, i1, i2, iaz, iddecl, iddrng, idecld, ideclm, idecls,
     .          idha, idra, idrang, iel, ihah, iham, ihas, irah, iram,
     .          irang, iras
      integer*4 isec,isgndc
c
c determine range (m), range-rate (cm/min), and
c range-rate-rate (cm/min/min)
      irang  = Rsitp(1)*Ltvel*1E3_10 + 0.5_10
      drange = DOT(Xsitep(1,1),Xsitep(4,1))/Rsitp(1)
      idrang = drange*Ltvel*6E6_10 + SIGN(0.5_10, drange)
      if(Jct(47).gt.0) then
c
c geocentric radius for special les-8/9 output
         iddrng = SQRT(Xm(1,1)**2 + Xm(2,1)**2 + Xm(3,1)**2)
     .            *Ltvel*1E3_10 + 0.5_10
      else
         ddrang = ( -drange**2 + DOT(Xsitep(4,1),Xsitep(4,1))
     .            - (Xsitep(1,1)*(Accp(1)-Xsite(7,1))+Xsitep(2,1)
     .            *(Accp(2)-Xsite(8,1))+Xsitep(3,1)*(Accp(3)-Xsite(9,1))
     .            ))/Rsitp(1)
         iddrng = ddrang*Ltvel*3.6E8_10 + SIGN(0.5_10, ddrang)
      endif
c
c make range negative if satellite is in earth's shadow
      call MENSUN(Jd, ctrecf, xsun, 0)
      do i = 1, 3
         xsun(i) = xsun(i)*Aultsc - Xm(i, 1)
      end do
      rsun2 = DOT(xsun,xsun)
      psun  = DOT(xsun,Xm(1,1))
      if(psun.le.0._10) then
         do i = 1, 3
            de(i)   = xsun(i)*psun/rsun2
            xsun(i) = Xm(i, 1) - de(i)
         end do
         psun = SQRT(DOT(xsun,xsun))*Ltvel
         if(psun.le.6378.16_10) irang = -irang
c
c distinction between umbra and penumbra for les-8/9
         if(Jct(47).gt.0) then
            angsun = 6.95E5_10/SQRT(rsun2)
            dde    = SQRT(DOT(de,de))*angsun
            if(psun.le.6378.16_10 + dde) then
               if(psun .gt. 6378.16_10 - dde) iddrng = -iddrng
            endif
         endif
      endif
c
c aberration correction due to rotation of earth
      do i = 1, 3
         zsat(i) = -Xsitep(i, 1) + Rsitp(1)*Xsite(i + 3, 1)
      end do
      rzsat = SQRT(zsat(1)**2 + zsat(2)**2 + zsat(3)**2)
c
c determine azimuth,elevation (deg)
      d   = Sitnrm(1, 1)*Nutpr(3, 1) + Sitnrm(2, 1)*Nutpr(3, 2)
     .      + Sitnrm(3, 1)*Nutpr(3, 3)
      qxn = DOT(zsat,Sitnrm(1,1))
      do i = 1, 3
         xmerid(i) = Nutpr(3, i) - d*Sitnrm(i, 1)
         xop(i)    = zsat(i) - qxn*Sitnrm(i, 1)
      end do
      qxn = qxn/rzsat
      el  = ASIN(qxn)/Convd
      call CROSS(xop, xmerid, xmerop)
      qxo2 = DOT(xmerop,Sitnrm(1,1))
      qxo1 = DOT(xmerid,xop)
      az   = ATAN2(qxo2, qxo1)/Convd
      if(az.lt.0._10) az = az + 360._10
c
c           calculate azimuth,elevation rate of change
c           could be put here
c
c           refraction correction
      if(Jct(42).gt.0) then
         pxn = SQRT(1._10 - qxn**2)
c
c optical frequency refraction (smart, 'spherical and
c practical astronomy')
         if(Jct(42).gt.1) then
c
c radio frequency refraction (millstone subroutine dell)
            eltrue = el
            if(eltrue.lt.0.) eltrue = 0.
            call DELL(eltrue, eloff)
            elref = eltrue + eloff
            call DELL(elref, eloff)
            elref = eltrue + eloff
c (only one call to dell would be needed to remove refraction
c from observed value of elevation)
            call DELL(elref, eloff)
            d = eloff
         else
            cot = pxn/qxn
            d   = cot*(1.61928E-2_10 - 1.856E-5_10*cot**2)
         endif
         el = el + d
         if(Jct(47).le.0) then
            coe = d*Convd
            soe = SIN(coe)
            coe = COS(coe)
            q1  = coe - soe*qxn
            q2  = soe*qxn*(qxn + pxn)
            do i = 1, 3
               ysat(i)     = zsat(i)*q1 + Sitnrm(i, 1)*q2*rzsat
               ysat(i + 3) = -Xsitep(i + 3, 1)
            end do
            goto 100
         endif
c
c no refraction correction
      else if(Jct(47).le.0) then
         do i = 1, 3
            ysat(i)     = zsat(i)
            ysat(i + 3) = -Xsitep(i + 3, 1)
         end do
         goto 100
      endif
c
c geocentric ra,decl,ha for special les-8/9
c output to calculate subsatellite longitude,latitude
      do i = 1, 3
         ysat(i)     = Xm(i, 1)
         ysat(i + 3) = 0._10
      end do
c
c transform to true equinox and equator of date
  100 sidtml = sidtm1 - Longr(1)
      if(Jct(41).gt.0) then
c
c transform to mean equinox and equator of besselian
c beginning of year
         do i = 1, 3
            xsat(i)     = bessel(i, 1)*ysat(1) + bessel(i, 2)*ysat(2)
     .                    + bessel(i, 3)*ysat(3)
            xsat(i + 3) = bessel(i, 1)*ysat(4) + bessel(i, 2)*ysat(5)
     .                    + bessel(i, 3)*ysat(6)
         end do
         sidtml = sidtml - Pc(1)
      else
         call CORCHN(xsat(1), ysat(1))
         call CORCHN(xsat(4), ysat(4))
      endif
c
c compute right ascension,declination,hour angle
      rsat2 = DOT(xsat,xsat)
      rsat  = SQRT(rsat2)
      decl  = ASIN(xsat(3)/rsat)/Convds
      ra    = ATAN2(xsat(2), xsat(1))
      ha    = sidtml - ra
      i     = ha/Twopi
      if(ha.lt.0._10) i = i - 1
      ha = ha - i*Twopi
      if(ra.lt.0._10) ra = ra + Twopi
      ra = ra/Convhs
      ha = ha/Convhs
c
c compute rates of change
      q2    = xsat(1)**2 + xsat(2)**2
      q1    = SQRT(q2)
      ddecl = (xsat(6) - xsat(3)/rsat2*DOT(xsat(1),xsat(4)))/q1/Convds
      dra   = (xsat(1)*xsat(5) - xsat(2)*xsat(4))/q2
      dha   = (Sidvel - dra)/Convhs
      dra   = dra/Convhs
c
c convert to integer parts
      iaz    = az*1E5_10 + 0.5_10
      iel    = el*1E5_10 + SIGN(0.5_10, el)
      irah   = ra/3600._10
      ra     = ra - irah*3600
      iram   = ra/60._10
      ra     = ra - iram*60
      iras   = ra*1E3_10 + 0.5_10
      isgndc = 1
      if(decl.lt.0._10) then
         isgndc = 2
         decl   = -decl
      endif
      idecld = decl/3600._10
      decl   = decl - idecld*3600
      ideclm = decl/60._10
      decl   = decl - ideclm*60
      idecls = decl*1E2_10 + 0.5_10
      ihah   = ha/3600._10
      ha     = ha - ihah*3600
      iham   = ha/60._10
      ha     = ha - iham*60
      ihas   = ha*1E3_10 + 0.5_10
      idra   = dra*1E4_10 + SIGN(0.5_10, dra)
      iddecl = ddecl*1E3_10 + SIGN(0.5_10, ddecl)
      idha   = dha*1E4_10 + SIGN(0.5_10, dha)
      isec   = Sec + 0.5
c
c make sure rates of change fit their fields
      if(idra.ge.1000000 .or. idra.le.-100000) idra = 0
      if(iddecl.ge.1000000 .or. iddecl.le.-100000) iddecl = 0
      if(idha.ge.1000000 .or. idha.le.-100000) idha = 0
      if(iddrng.ge.1000000000 .or. iddrng.le.-100000000) iddrng = 0
      if(idrang.le.-1000000000) idrang = 0
c
c convert to characters and store in buffer
      write(Intern, 200) Imonth, Iday, Iyear, Ihr, Imin, isec, irah,
     .                   iram, iras, idra, minus(isgndc), idecld,
     .                   ideclm, idecls, iddecl, ihah, iham, ihas,
     .                   idha, iaz, iel, irang, idrang, iddrng
  200 format(6I2, i2, i2, i5, i6, a1, i2, i2, i4, i6, i2, i2, i5, i6,
     .       2I8, i10, i10, i9)
      rewind Intern
      i1     = iangbf*51 + 1
      iangbf = iangbf + 1
      i2     = iangbf*51
      read(Intern, 300) (angbuf(i), i = i1, i2)
  300 format(51A2)
      rewind Intern
 
      return
      end
