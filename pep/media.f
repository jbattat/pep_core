      subroutine MEDIA(ngo,rtrn)
 
      implicit none

c     m.e.ash   sept 1967     subroutine media
c     effect of interplanetary media (plasma) on time delay and doppler
c
c arguments
      integer*4 ngo
      real*4    rtrn
c     ngo =1 effect on time delay
c     ngo =2 effect on doppler shift
c     ngo =3 effect on doppler shift phase delay
c     changes to media made by t.m. eubanks ---nov 1976
c     rtrn = value of correction calculated by media
c     rtrn (ngo = 3) = - rtrn (ngo = 1)
c
c     prmter(60) and prmter(61) /param/ (on &nmlst1)
c     are defaulted to zero in prmred
c     prmter(60) = static plasma density in atoms/cc at 1 a.u.
c                  good estimate = 10  atom/cc at 1 a.u.
c     prmter(61) = time varying part of plasma -same units
c  ----------------------------------------------------------------
c     note : prmter(60) must be set in the &nmlst1
c            the default of 0 makes the media correction = 0
c  ----------------------------------------------------------------

c array dimensions
      include 'globdefs.inc'

c commons 
      include 'comdateq.inc'
      include 'coord.inc'
      include 'funcon.inc'
      include 'ltrapx.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'obscrd.inc'
      include 'kobequiv.inc'
      real*10 ctrecf
      equivalence (ctrecf,Dstf)
      include 'param.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'

c external functions
      real*10 DOT
 
c local variables
      real*10 aarct,close,close2,fofres,he,hp,
     .          rerp,resq,rexrp,rpsq,tfcons,tplas,trm3,trm4,xex
      real*10 te(3),tp(3),dxsit(3)
      integer   i
c
c the formulation here is just twice the reflect-receive
c rather than the full send-receive
c
c effect of interplanetary plasma on time delay
      Raddum(2) = 0.0_10
      Raddum(3) = 0.0_10
      rtrn = 0.
      if(Freq.gt.0._10) then
 
c do delay setup even if just doppler
         if(ngo.ne.2 .or. Nice.gt.0) then
            xex    = DOT(Xsitp0,Xemlsc)
            resq   = DOT(Xemlsc,Xemlsc)
            close2 = resq - xex**2
            close  = SQRT(close2)
            rpsq   = DOT(Xsbsun,Xsbsun)
            rerp   = DOT(Xemlsc,Xsbsun)
            rexrp  = SQRT(resq*rpsq - rerp**2)
            aarct  = ATAN2(rexrp,rerp)
c solar minimum c. 1964.0
            tplas  = Jd - 2438396
            tplas  = tplas + ctrecf - tmdly/172800.0_10
            tfcons = Twopi/(11.0_10*365.25_10)
            trm3   = tfcons*tplas
            trm4   = -COS(trm3)
 
c new fofres so that prmter60 &61 are elec/cc at 1 a.u.
            fofres    = 8.2E7_10*(Aultsc/Freq)**2
            Raddum(2) = fofres/close*aarct
            if(nddiff.lt.0) Raddum(2) = Raddum(2)*0.5_10
            if(ngo.eq.1) then
            else if(ngo.eq.2) then
               goto 50
            else
 
c phase delay
               Raddum(2) = -Raddum(2)
            endif
            Raddum(3) = Raddum(2)*trm4
            rtrn = prmter(60)*Raddum(2) + prmter(61)*Raddum(3)
            return
         endif
c
c effect of interplanetary plasma on doppler shift
   50    call UVECTR(4,Xsitep,Rsitp,Xsitp0,dxsit)
         he = rerp/resq
         hp = rerp/rpsq
         do i = 1,3
            te(i) = he*Xemlsc(i,1) - Xsbsun(i)
            tp(i) = hp*Xsbsun(i) - Xemlsc(i,1)
         end do
         Raddum(7) = fofres*((DOT(Xemlsc(4,1),te)+DOT(Xsbsun(4),tp))
     .               /rexrp -
     .               aarct*(DOT(Xemlsc,Xemlsc(4,1))-xex*(DOT(Xsitp0,
     .               Xemlsc(4,1))+DOT(dxsit,Xemlsc)))/close2)/close
         Raddum(8) = Raddum(7)*trm4 + Raddum(2)*COS(trm3)*tfcons/Secday
         rtrn = Raddum(7)*prmter(60) + Raddum(8)*prmter(61)
      endif
 
      return
      end
