      subroutine RADNRM
 
      implicit none
c
c        r.b. goldstein  may 1978
c        routine to compute the jpl normal point pseudo observable.
c        this observable is the round trip time delay from the
c        center of mass of the earth to the center of mass
c        of a planet, where all object positions are to be
c        evaluated at time=ctrec
c        the time tag is given as ctrec
c
c        this is called from radctl when itime=5

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'coord.inc'
      include 'eqnphs.inc'
      include 'ltrapobs.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'number.inc'
      include 'obscrd.inc'
      include 'kobequiv.inc'
      real*10 ctrecf,fdev,reflct
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      include 'param.inc'
      include 'pqind.inc'
      include 'prpgat.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'
c
c local variables
      real*10 ctv30,dum,tau3t0,tmp
      integer   i,j,lemctl,lplctl
c
c setups
      Atuts  = 0._10
      Ututs  = 0.
      Ctat   = 0._10
      Utrec  = 0._10
      Ctrec  = Sec
      Ctrec  = Ihr*3600 + Imin*60 + Ctrec
      Fract  = Ctrec/Secday
      ctv30  = Ctvary*(Jds - Prm97 - 0.5_10 + Fract)**2
      Ctrec  = Ctrec + ctv30*Secday
      Fract  = Ctrec/Secday
      ctrecf = Fract
      lemctl = 1
      lplctl = 1
      Jd     = Jds
c
c
c c.m. of earth relative to sun
c
      call ETRP(1,Jd,Fract,0,lemctl,1,2)
      if(Jd.gt.0) then
c
c c.m. of planet relative to sun
c
         call PLTRP(1,Jd,Fract,0,lplctl)
         if(Jd.gt.0) then
            do j = 1,3
               Xsbsun(j) = Xp(j)*Aultsc
            end do
            if(Klanb.gt.0) then
               call SBTRP(1,Jd,Fract,0,lplctl)
               if(Jd.le.0) return
               do j = 1,3
                  Xsbsun(j) = Xsbsun(j) + Cmfct*Aultsc*Xsb(j)
               end do
            endif
 
c xplsc needed for media and maybe other propco routines
            do j = 1,3
               Xplsc(j) = Xsbsun(j)
            end do
c
c calculate earth-planet distance  (xsitep)
c
            do j = 1,3
               Xsitep(j,1) = Xemlsc(j,1) - Xsbsun(j)
            end do
 
            call UVECTR(3,Xsitep(1,1),Rsitp(1),Xsitp0(1,1),dum)
            reflct = Rsitp(1)
c
c observable is just twice the earth-planet distance. no
c send time calculations are needed. fill send arrays
c from the receive arrays.
            Rsitp(2) = Rsitp(1)
            tmdly    = Rsitp(1) + Rsitp(2)
c
c correction needed if ctvary is non-zero
c
            if(Ctvary.ne.0.0_10) then
               tmp    = tmdly/Secday
               tau3t0 = Jds - Prm97 + Fract - 0.5_10
               tmdly  = tmdly -
     .                  Ctvary*((tmp-ctv30)*(2._10*tau3t0-tmp-ctv30)
     .                  + 2._10*Ctvary*(tau3t0-tmp)**3)*Secday
            endif
            do j = 1,3
               Xem(j,2)    = Xem(j,1)
               Xm(j,2)     = Xm(j,1)
               Xemlsc(j,2) = Xemlsc(j,1)
               Xsitep(j,2) = Xsitep(j,1)
               Xsitp0(j,2) = Xsitp0(j,1)
            end do
            do i = 1,4
               Pem(i,2) = Pem(i,1)
               Pm(i,2)  = Pm(i,1)
            end do
            if(calcvl.gt.0) then
c
c---------------- need velocities for ctvary partial -------------------
c
c zero out site coordinates (used in ctvary partial)
               call ZFILL(Xsite,16*18)
               call ETRP(1,Jd,Fract,1,lemctl,1,2)
               if(Jd.le.0) return
               call PLTRP(1,Jd,Fract,1,lplctl)
               if(Jd.le.0) return
               do j = 4,6
                  Xsbsun(j) = Xp(j)*Aultvl
               end do
               if(Klanb.gt.0) then
                  call SBTRP(1,Jd,Fract,1,lplctl)
                  if(Jd.le.0) return
                  do j = 4,6
                     Xsbsun(j) = Xsbsun(j) + Cmfct*Aultvl*Xsb(j)
                  end do
               endif
               do j = 4,6
                  Xemlsc(j,2) = Xemlsc(j,1)
                  Xsitep(j,1) = Xemlsc(j,1) - Xsbsun(j)
                  Xsitep(j,2) = Xsitep(j,1)
               end do
c
c relative velocity
c
               call VLRTRD(Xemlsc(4,1),Xsbsun(4),Xemlsc(4,2),0.0_10,
     .                     7,0)
            endif
c
c corrections
c
            call RADREL(1)
            call PROPCO(1,1)
            tmdly = tmdly + Sumcor(1)
c
c bias
c constant bias in time delay
            if(Nrbias.gt.0) tmdly = tmdly + Rbsx(1)
         endif
      endif
c
c*  start=9900
      return
      end
