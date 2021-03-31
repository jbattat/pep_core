      subroutine PVCRD(jd,fract,nsite,k,nvel)
 
      implicit none

c
c rw king    december 1978  subroutine pvcrd
c calculate venus-centered coordinates of pioneer-venus entry
c probes during the atmospheric phase of their trajectories
c
c parameters
      real*10 fract
      integer*4 jd,nsite,k,nvel
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 tpvmsb(2),tpvmsc(2)
      equivalence (Angdum,tpvmsb),(Angdum(3),tpvmsc)
      include 'empcnd.inc'
      include 'funcon.inc'
      include 'mnsprt.inc'
      include 'number.inc'
      include 'param.inc'
      include 'sbcom.inc'
      include 'scdta.inc'
c
c variables internal to this routine
      real*10 dlat,dlong,drad,dum1,dum2,qq,rad,t,tsb0,tsc0
      integer   i,klb,mnspt1
      integer*2 np2/2/
      real*10 lat,long,dy(3),cl(2),sl(2)
 
      mnspt1 = 0
 
      if(jd.ne.2443852) call SUICID(
     .    ' PVCRD WILL WORK ONLY FOR JD 244385, STOP IN PVCRD  ', 13)
 
      if(k.eq.2) then
 
         klb  = Klans1
         tsc0 = (Sccom(2) - 2443852._10)*Secday
         tpvmsc(nsite) = fract*Secday - tsc0
         t = tpvmsc(nsite)
      else
         klb  = Klanb
         tsb0 = (Sbcom(2) - 2443852._10)*Secday
         tpvmsb(nsite) = fract*Secday - tsb0
         t = tpvmsb(nsite)
      endif
 
      rad  = Pcond(7,klb) + t*(Pcond(8,klb)
     . + t*(Pcond(9,klb)/2._10
     . + t*(Pcond(10,klb)/6._10
     . + t*Pcond(11,klb)/24._10)))
      long = Pcond(12,klb) + t*(Pcond(13,klb)
     . + t*(Pcond(14,klb)/2._10
     . + t*(Pcond(15,klb)/6._10
     . + t*Pcond(16,klb)/24._10)))
      lat  = Pcond(17,klb) + t*(Pcond(18,klb)
     . + t*(Pcond(19,klb)/2._10
     . + t*(Pcond(20,klb)/6._10
     . + t*Pcond(21,klb)/24._10)))
c
c
      qq  = rad/Ltvel
      lat = lat*Convd
      Yspcd(3,k) = qq*SIN(lat)
      qq   = qq*COS(lat)
      long = long*Convd
      cl(nsite)   = COS(long)
      sl(nsite)   = SIN(long)
      Yspcd(2,k)  = qq*sl(nsite)
      Yspcd(1,k)  = qq*cl(nsite)
      Rspot(k)    = SQRT(Yspcd(1,k)**2 + Yspcd(2,k)**2 + Yspcd(3,k)**2)
      Dydphi(1,k) = -Yspcd(3,k)*cl(nsite)
      Dydphi(2,k) = -Yspcd(3,k)*sl(nsite)
      Dydphi(3,k) = qq
 
      call SPOTCD(jd,fract,mnspt1,nvel,np2,dum1,dum2,k)
 
      if(k.eq.2) then
         do i = 1, 3
            Xsc(i,nsite) = Xspcd(i,2)/Aultsc
         end do
      else
         do i = 1, 3
            Xsb(i) = Xspcd(i,1)/Aultsc
         end do
      endif
      return
c
c-----------------------------------------------------------------------
c
c           calculate velocities and partial derivatives
c
      entry PVCRDV(nsite,nvel,k)
c
c           spot partials calculated for use by partl1
c     note that these partials assume a common transmit time for signals
c     received at each site.  this is an approximation for differenced
c     n-count observables.
      do i = 1, 3
         Dxspcd(i,1,k) = Xspcd(i,k)/Rspot(k)
         Dxspcd(i,2,k) = -Rot(1,i)*Yspcd(2,k) + Rot(2,i)*Yspcd(1,k)
         Dxspcd(i,3,k) = Rot(1,i)*Dydphi(1,k) + Rot(2,i)
     .                     *Dydphi(2,k) + Rot(3,i)*Dydphi(3,k)
      end do
 
      if(nvel.gt.0) then
         t = tpvmsb(nsite)
         if(k.eq.2) t = tpvmsc(nsite)
         drad  = (Pcond(8,klb) + t*(Pcond(9,klb)
     .    + t*(Pcond(10,klb)/2._10
     .    + t*Pcond(11,klb)/6._10)))/Ltvel
         dlong = (Pcond(13,klb) + t*(Pcond(14,klb)
     .    + t*(Pcond(15,klb)/2._10
     .    + t*Pcond(16,klb)/6._10)))*Convd
         dlat  = (Pcond(18,klb) + t*(Pcond(19,klb)
     .    + t*(Pcond(20,klb)/2._10
     .    + t*Pcond(21,klb)/6._10)))*Convd
         dy(1) = drad*Yspcd(1,k)/Rspot(k) - dlong*Yspcd(2,k)
     .           - dlat*cl(nsite)*Yspcd(3,k)
         dy(2) = drad*Yspcd(1,k)/Rspot(k) + dlong*Yspcd(1,k)
     .           - dlat*sl(nsite)*Yspcd(3,k)
         dy(3) = drad*Yspcd(3,k)/Rspot(k) + dlat*Yspcd(1,k)/cl(nsite)
 
         call SPOTCV(nvel,np2,k)
 
         if(k.eq.2) then
 
            do i = 1, 3
               Xsc(i + 3,nsite) = Xspcd(i + 3,2) + Rot(1,i)*dy(1)
     .                             + Rot(2,i)*dy(2) + Rot(3,i)*dy(3)
               Xsc(i + 3,nsite) = Xsc(i + 3,nsite)/Aultvl
            end do
         else
 
            do i = 1, 3
               Xsb(i + 3) = Xspcd(i + 3,1) + Rot(1,i)*dy(1)
     .                      + Rot(2,i)*dy(2) + Rot(3,i)*dy(3)
               Xsb(i + 3) = Xsb(i + 3)/Aultvl
            end do
         endif
      endif
      return
      end
