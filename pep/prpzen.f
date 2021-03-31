      subroutine PRPZEN(kobj)
 
      implicit none
c
c
c     r.w. king and r.b. goldstein   june 1978
c     routine to determine and print out zenith angles used for earth
c     atmosphere and ionosphere corrections

c arguments
      integer*4 kobj

c common
      include 'coord.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'number.inc'
      include 'obscrd.inc'
      include 'prpgat.inc'
      include 'radcrd.inc'
      include 'sitcrd.inc'
c
c quantities internal to this routine
      real*4    deg1, deg2, cutoff, deg(2)
      equivalence (deg, deg1), (deg(2), deg2)
      integer   i,izc,n,ns1
c
c
c           calculation and printout of zenith angles
c
c           izctl set in prplog
c           = 1 calculate zenith angles
c           = 3 calculate zenith angles and rates
c

c portable single-precision arithmetic statement function
      real*4 SNG10
      real*10 xx
      SNG10(xx)=xx
c
      if(Nsite2.eq.0) then
         Za(2, kobj)  = 0.
         Zar(2, kobj) = 0.
      endif
 
      if(Numsav.lt.51) then
         ns1 = Numsav + 1
         do i=ns1,51
            Save(i)=0._10
            if(i.ge.41 .and. i.le.46) Save(i)=-1._10
         end do
      endif
      if(kobj.eq.2) then
c
c           ict(24)= -1 print out all elevation angles
c                  =  0 do not print out all elevation angles
c                  = n  print out only elevation angles < n degrees
c                  = negative, ne.-1, test uses abs(n)
c
c           second object
c
         call ZENANG(Izctl, Ysitep(1,1), Ysitp0(1,1), Rsitp2(1),
     .        Xsite(1,1), Sitnrm(1,1), Rsite(1), Za(1,2), Zar(1,2))
         Save(54) = Za(1, 2)
         if(Numsav.lt.54) Numsav = 54
         if(Nsite2.ne.0) then
            call ZENANG(Izctl, Ysitep(1,2), Ysitp0(1,2), Rsitp2(2),
     .           Xsite(1,2), Sitnrm(1,2), Rsite(2), Za(2,2), Zar(2,2))
            Save(55) = Za(2, 2)
            if(Numsav.lt.55) Numsav = 55
         endif
      else
c
c first object
         call ZENANG(Izctl, Xsitep(1,1), Xsitp0(1,1), Rsitp(1),
     .        Xsite(1,1), Sitnrm(1,1), Rsite(1), Za(1,1), Zar(1,1))
         Save(52) = Za(1, 1)
         if(Numsav.lt.52) Numsav = 52
         if(Nsite2.ne.0) then
            call ZENANG(Izctl, Xsitep(1,2), Xsitp0(1,2), Rsitp(2),
     .           Xsite(1,2), Sitnrm(1,2), Rsite(2), Za(2,1), Zar(2,1))
            Save(53) = Za(2, 1)
            if(Numsav.lt.53) Numsav = 53
         endif
      endif
c
c printout
      if(Ict(24).ne.0) then
         do n = 1, 2
            deg(n) = 90. - Za(n, kobj)/SNG10(Convd)
         end do
         if(Ict(24).ne.-1) then
            cutoff = Ict(24)
            cutoff = abs(cutoff)
            if(deg1.gt.cutoff .and. deg2.gt.cutoff) return
         endif
         if(Line.gt.55) call OBSPAG
         write(Iout, 50) kobj, deg
         Line = Line + 1
   50    format(' ELEV. ANG. FOR OBJ.', i2, ', SITES 1 & 2:  ', 2F12.6,
     .          '  (DEGREES)')
         izc = Izctl
         if(mod(izc/2,2).eq.1) then
            do n = 1, 2
               deg(n) = -Zar(n, kobj)/SNG10(Convd)
            end do
            write(Iout, 60) kobj, deg
            Line = Line + 1
   60       format(' ELEV. RATE FOR OBJ.', i2, ', SITES 1 & 2:  ',
     .             2E15.6, '  (DEG/SEC)')
         endif
      endif
 
      return
      end
