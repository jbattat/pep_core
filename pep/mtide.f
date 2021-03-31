      subroutine MTIDE(mtide1, nvel, kobj)
 
      implicit none

c
c          rw king and rj cappallo    august 1978   subroutine mtide
c           interface between vetide and pep.  mtide sets up
c           various vectors, calls vetide to calculate the solid body
c           tide raised by the earth and sun on the moon and alters the
c           appropriate spot coordinates .
c
c parameters
      integer*4 mtide1, nvel, kobj
c
c           mtide1 = 0  first call of mtide for a given spot:
c                       calculate tides and correct spot coordinates
c                  = 1  calls of mtide on successive iterations for a
c                       given spot:  correct spot coordinates but do
c                       not recalculate tides
c
c           nvel = 0  do not calculate velocities
c                  1  calculate velocities
c
c           kobj = 1  observed object is nplnt0
c                = 2  observed object is nplnt2
c
c
c         etide requires that the relative coordinates of the
c         earth, sun, and moon be available in /coord/ arrays
c         xem and xm, and that the unit vectors from site to
c         source be available in xsitp0 or ysitp0 of /coord/
c         or /radcrd/, respectively, for kobj= 1,2.

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
      include 'mnsprt.inc'
      include 'obscrd.inc'
      include 'param.inc'
      include 'radcrd.inc'
      include 'tidal.inc'
c
c variables internal to this routine
      logical*1 init/.false./
      real*10 llov, xmetid(6), xmstid(6), dxdhme(6), dxdlme(6),
     .          dxdhms(6), dxdlms(6), xel(6), xsl(6), xmtid(6)
      real*10 cphi, dellat, dellon, delrad, dlate, dlats, dlone,
     .          dlons, drade, drads, g, gmsun, hlov, tlag
      integer   i, icme, icms, index, ipt, iptln, j, jct11
c
c
c----------------------------------------------------------------------
c
      index = 3
      if(nvel.gt.0) index = 6
c
c check mtide1
      if(mtide1.le.0) then
c
c           calculations performed only once
c
c     these calculations might be moved to eshape in the future
c     with values stored in /tidal/
c
         if(.not.(init)) then
 
            gmsun  = Gauss**2*Aultsc**3/Secday**2
            Gmerth = gmsun*Mass(3)*(1._10 - Mass(10))
            Gmmoon = gmsun*Mass(3)*Mass(10)
 
c effects of rotation and inhomogeneities ignored for lunar tide
            g    = Gmmoon/(Mcond(7)/Ltvel)**2
            cphi = 1._10
            hlov = .018_10
            llov = 0.0026_10
            tlag = 0._10
c
c set love numbers equal to input parameters if non-zero
            if(Mcond(9).ne.0._10) hlov  = Mcond(9)
            if(Mcond(10).ne.0._10) llov = Mcond(10)
            if(Mcond(11).ne.0._10) tlag = Mcond(11)
 
            jct11 = Jct(11)
            icme  = mod(jct11, 2)
            icms  = mod(jct11/2, 2)
            iptln = mod(jct11/4, 2)
            ipt   = mod(jct11/8, 2)
            init  = .true.
         endif
c
c
c initialization once per observation
         do i = 1, index
            xmetid(i) = 0._10
            xmstid(i) = 0._10
            dxdhme(i) = 0._10
            dxdlme(i) = 0._10
            dxdhms(i) = 0._10
            dxdlms(i) = 0._10
            end do
         if(ipt.ne.0) then
            drade = 0._10
            dlate = 0._10
            dlone = 0._10
            drads = 0._10
            dlats = 0._10
            dlons = 0._10
         endif
c
c
c calculate vector tide raised on moon by the earth
         if(icme.ne.0) then
            do i = 1, 3
               xel(i) = -Xm(i, 1)*Mnltsc
               if(nvel.gt.0) xel(i + 3) = -Xm(i + 3, 1)
     .             *Mnltsc/Secday
               end do
            call VETIDE(nvel, Gmerth, g, cphi, hlov, llov, xel,
     .                  Xspcd(1,kobj), Dxspcd(1,1,kobj),
     .                  Dxspcd(1,3,kobj), Dxspcd(1,2,kobj), xmetid,
     .                  dxdhme, dxdlme, drade, dlate, dlone)
         endif
c
c calculate vector tide raised on moon by sun
         if(icms.ne.0) then
            do i = 1, 3
               xsl(i) = -(Xem(i,1) + (1._10-Mnfct*Xm(i,1)))*Aultsc
               if(nvel.gt.0) xsl(i + 3)
     .             = -(Xem(i+3,1) + (1._10-Mnfct*Xm(i+3,1)))*Aultvl
               end do
            call VETIDE(nvel, gmsun, g, cphi, hlov, llov, xsl,
     .                  Xspcd(1,kobj), Dxspcd(1,1,kobj),
     .                  Dxspcd(1,3,kobj), Dxspcd(1,2,kobj), xmstid,
     .                  dxdhms, dxdlms, drads, dlats, dlons)
         endif
c
c
c add together earth and solar components of tide
         do i = 1, index
            xmtid(i) = xmetid(i) + xmstid(i)
            end do
      endif
c
c change cartesian coordinates of spot vector
      do i = 1, index
         Xspcd(i, kobj) = Xspcd(i, kobj) + xmtid(i)
 
c lag angle coded only for delays
         if(i.le.3) Xspcd(i, kobj) = Xspcd(i, kobj)
     .       + tlag*xmtid(i + 3)
         end do
      if(mtide1.le.0) then
c
c partials w.r.t. moon love-number scale factors and time lag
         do j = 1, 2
            do i = 1, index
               Dxdhm(i, j) = (dxdhme(i) + dxdhms(i))
               Dxdlm(i, j) = (dxdlme(i) + dxdlms(i))
               Dxdtm(i, j) = 0._10
               if(i.le.3) then
                  Dxdhm(i, j) = Dxdhm(i, j)
     .                          + tlag*(dxdhme(i+3) + dxdhms(i+3))
                  Dxdlm(i, j) = Dxdlm(i, j)
     .                          + tlag*(dxdhme(i+3) + dxdhms(i+3))
                  Dxdtm(i, j) = xmtid(i + 3)
               endif
               end do
            end do
c
c
c optional printout of love no.s and tidal components
         if(Nk1.le.0) then
            if(iptln.ne.0) then
               if(Line.gt.57) call OBSPAG
               write(Iout, 10) hlov, llov, tlag
               Line = Line + 1
   10          format(' MOON LOVE NUMBERS:  H=', f8.4, ' L=', f8.4,
     .                ' TLAG=', e10.4, ' SEC')
            endif
         endif
         if(ipt.ne.0) then
            delrad = (drade + drads)*hlov
            dellat = (dlate + dlats)*llov
            dellon = (dlone + dlons)*llov
            if(Line.gt.58) call OBSPAG
            write(Iout, 20) kobj, delrad, dellat, dellon
   20       format(' MOON TIDE COMPONENTS FOR OBJECT', i2, ':  RAD=',
     .             1pd14.6, '  LAT=', 1pd14.6, '  LONG=', 1pd14.6,
     .             '  (SEC)')
            Line = Line + 1
         endif
      endif
c
c
c
      return
      end
