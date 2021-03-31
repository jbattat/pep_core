      subroutine SBRSET(w, npnpx2, natpnp)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   natpnp, npnp, npnpx2
 
c*** end of declarations inserted by spag
 
 
c paul e. macneil april, 1980
c reset (offset) process noise paramters other than s/c ics
c
c common
      include 'lothrf.inc'
c
c process noise parameter offsets for kalman filter facilty
      real*10 w(npnpx2)
c
c*       start=1000
c
c reset process noise parameter(s) on call from sbfout
      if( natpnp .gt. 0 ) then
 
c offset atmospheric process noise parameter(s)
         npnp = npnpx2/2
         Rhoz = Rhoz + w(npnp + 1)
      endif
c
c*  start=9900
      return
      end
