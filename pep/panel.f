      subroutine PANEL(s, ap, icode)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      convd, cos2a, cos3a, cosa, dt, phi, sin2a, sin3a, sina
      integer   i, ipa, ipb, is, isx, isy, nerr, ntphi
 
c*** end of declarations inserted by spag
 
 
c
c r.reasenberg  subroutine panel 2/76
c
c cleaned up: m. murison 5/91
c
      real*10 s, ap(3)
      integer*4 icode
 
c
c        this subroutine uses equation 3 of jpl
c        tm 391-240 of c christensen to find the
c        effective area of the panels of the mvm
c        s/c.  the panel angle is stored in
c        phit(i) for the time span tphi(i)
c        to tphi(i+1), where time is j.date
c        past t0.  for times when phi is
c        changing quickly, a linear interpolation
c        is provided.  its use is indicated by ipa==
c        iphi(bit2)=1.  for times when the
c        s/c was in free roll, ap(1)=ap(2)=0,
c        and ipb==iphi(bit1)=1.
c
c        icode
c        0         error
c        1         ok
c
c        assume time is ut at s/c
c
      real*10 time, t0/2441990./, az/27640./
      data ntphi/22/, nerr/0/, is/-1/
      integer   iphi(25)/25*1/
      data convd/0.01745329/
 
c epochs
c  1 -  2 start to tcm1
c  3 -  5 tcm1 to tcm2
c  6 -  7 tcm2 to tcm3 venus encounter
c  8 - 11 tcm3 to tcm5 superior conjunction
c 12 - 18 tcm5 to tcm8
c 19 - 22 tcm8 to eom mercury (3) encount
      real      tphi(22)
      data tphi/0.0, 0.25, 41.8, 42.74, 79.88, 104.07, 117.92, 133.8,
     .     140.12, 172.30, 200., 290., 297.5, 300.72, 304.8, 307.89,
     .     311.51, 340., 480., 489.38, 497.37, 510./
 
      real      phit(22)
      data phit/90., 0., 25., 26., 45., 58., 66., 68., 71., 66., 200.,
     .     56., 66., 68., 68.5, 66., 70., 200., 66., 70., 71., 200./
 
      real      pnu/.035/, pmu/.075/
 
 
      time = s - t0
      if(is.le.0) then
         write(6, 50) t0, ntphi, (tphi(i), phit(i), i = 1, ntphi)
   50    format(' FIRST CALL TO PANEL  T0,NTPHI=', f20.10, i5,
     .          '   TPHI,PHIT VECTORS', 20(/,5x,5('(',f10.4,f8.2,') '))
     .          )
      else if(time.ge.tphi(is)) then
         if(time.le.tphi(is+1)) then
 
            if(ipb.eq.0) go to 300
            go to 100
         endif
      endif
 
c find the panel angle phi
      if(time.ge.tphi(1)) then
         if(time.le.tphi(ntphi)) then
 
            do i = 2, ntphi
               if(time.lt.tphi(i)) then
 
                  is  = i - 1
                  phi = phit(is)*convd
                  if(phi.gt.3.2) go to 400
                  ipb = iphi(is)/2
                  ipa = iphi(is) - 2*ipb
                  write(6, 60) ipb, ipa, is, phi, s, time
   60             format(' PANEL', 3I5, 3D20.10)
                  if(ipb.ne.0) go to 100
                  go to 200
               endif
 
c can not find time in table
               end do
         endif
      endif
      go to 400
  100 dt = tphi(is + 1) - tphi(is)
      if(dt.lt..001) go to 400
      if(phit(is+1).gt.190) go to 400
      phi = (phit(is)*(tphi(is+1)-time) + phit(is+1)*(time-tphi(is)))/dt
      phi = phi*convd
 
c find the trig functions
  200 cosa  = cos(phi)
      sina  = sin(phi)
      cos2a = cosa*cosa - sina*sina
      sin2a = 2*sina*cosa
      cos3a = cosa*cos2a - sina*sin2a
      sin3a = sina*cos2a + cosa*sin2a
 
c find effective area
  300 ap(3) = az*(pnu + (1.+pmu)*cos(phi) + pnu*cos(2.*phi)
     .        + pmu*cos(3.*phi))*2.
      ap(1) = 0.0
      if(ipa.eq.0) then
         ap(2) = -az*(pmu*sin(phi) + pnu*sin(2.*phi) + pmu*sin(3*phi))
     .           *2.
      else
         ap(2) = 0.0
      endif
      icode = 1
      return
  400 nerr  = nerr + 1
      ap(1) = 0.0
      ap(2) = 0.0
      ap(3) = 0.0
      icode = 0
      isx   = is
      if(is.le.0) isx = 1
      isy = isx + 1
      if(isy.gt.ntphi) isy = ntphi
      write(6, 500) time, is, (tphi(i), phit(i), i = isx, isy)
  500 format(' ERROR IN PANEL', f20.10, i5, 4F20.10)
      return
      end
