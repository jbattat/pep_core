      subroutine PATCTL(icall,kobj,patcor)
 
      implicit none

c
c
c     r.w.king   sept 1978
c     routine to determine the effect of planetary atmosphere
c     on delay and delay rate for a probe within the
c     atmosphere of a planet
c
c parameters
      integer*4 icall, kobj
      real*4 patcor(4)
c      icall= 1  calculate correction to delay
c           = 2  calculate correction to delay rate (coded = 0.)
c
c       kobj= 1  observed object is nplnt0 or nspot
c           = 2  observed object is nplnt2 or nspot2
c
c     patcor(1-2)=  delay (or delay rate) corrections for sites 1 and 2
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'atmos.inc'
      include 'comdat.inc'
      real*10 sbcom(12)
      equivalence (Comcon(61),sbcom)
      include 'coord.inc'
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'mnsprt.inc'
      include 'number.inc'
      include 'obscrd.inc'
      include 'param.inc'
      include 'radcrd.inc'
      include 'scdta.inc'

c external functions
      real*10 DOT

c variables internal to this routine
      real*10 cutoff, deg1, deg2, h, h0
      integer   i, jdpvm, numsit
      real*4    theta(2)
      real*10 costhe(2)
 
      numsit = 2
      if(Nsite2.eq.0) numsit = 1
c
c zero correction and partials
      do i = 1, 2
         patcor(i) = 0.
         Pdlpat(i,kobj,2) = 0._10
         Pdlpat(i,kobj,1) = 0._10
      end do
      if(kobj.eq.2) then
         if(Nplnt2.le.0) return
         if(Nspot2.le.0 .and. Klans1.le.0) return
         jdpvm = Sccom(1)
         if(jdpvm.ne.2443852) return
      else
         if(Nplnt0.le.0) return
         if(Nspot.le.0 .and. Klanb.le.0) return
 
c if probe, must be pioneer-venus entry vehicle
         jdpvm = sbcom(1)
 
c not coded for planetary radar
         if(jdpvm.ne.2443852) return
      endif
 
      if(icall.eq.1) then
c
c determine zenith angle of receiving site as viewed from probe
         do i = 1, numsit
            if(kobj.eq.1) costhe(i) = DOT(Xspcd(1,1),Xsitep(1,i))
     .          /Rspot(1)/Rsitp(i)
            if(kobj.eq.2) costhe(i) = DOT(Xspcd(1,2),Ysitep(1,i))
     .          /Rspot(2)/Rsitp2(i)
            theta(i) = ACOS(costhe(i))
         end do
c
c optional printout of elevation angle
         if(Jct(24).ne.0) then
            deg1 = 90. - theta(1)*57.295779
            if(Nsite2.eq.0) theta(2) = 0.
            deg2 = 90. - theta(2)*57.295779
            if(Jct(24).ne.-1) then
               cutoff = Jct(24)
               if(deg1.gt.cutoff .and. deg2.gt.cutoff) goto 50
            endif
            if(Line.gt.57) call OBSPAG
            write(Iout,20) kobj,deg1,deg2
   20       format(' PLANTOCENTRIC ELEV. ANG. FOR OBJ.', i2,
     .             '  SITES 1 AND 2:', 2F12.6, '   (DEGREES)')
            Line = Line + 1
         endif
c
c calculate total delay and partial derivatives
   50    do i = 1, numsit
 
c height above surface (radar radius) in km
            h = Rspot(kobj)*Ltvel - Pcond(7,Klan)
 
c atm. scale height
            h0 = Pcond(24,Klan)
            Pdlpat(i,kobj,1) = EXP(-h/h0)/costhe(i)
 
c total delay
            patcor(i) = Pdlpat(i,kobj,1)*Pcond(25,Klan)
 
c partial w.r.t. scale height  (pcond(24,klan), units= km)
            Pdlpat(i,kobj,2) = patcor(i)*h/h0**2
         end do
      else if(Nk1.le.0) then
         write(Iout,100)
  100    format(i17,
     .' WARNING:  DELAY RATE CORRECTION FOR PLANETARY ATMOSPHERE SET = 0
     .. ***')
         Line = Line + 1
      endif
c
c
c
      return
      end
