      subroutine RADNS(nvel,norm)
 
      implicit none

c          r. goldstein     april, 1976
c          interpolations are now done by calls to mntrp,emtrp,
c          pltrp, and sbtrp.
c          polynomial extrapolations are now done by dguess
c
c       sept. 1976...derivative of sbdldp (sbbg)
c         this routine handles radar to a natural satellite of a
c         planet other than the earth or to a spacecraft treated as
c         a natural satellite
c         required flags:
c            ncodf=1, nspot=0, klanb.gt.0
c       oct. 1976...large sections of code replaced by calls to
c         subroutines propco and radrel
c       may 1978... converted to implicit real*8  and unit vector code
c           made compatible with rest of 'new code'
c
c     this subroutine is called by radctl
c
c arguments
      integer*4 norm,nvel

c array dimensions
      include 'globdefs.inc'

c commons
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 dutrec
      equivalence (dutrec,Angdum(10))
      include 'eqnphs.inc'
      include 'ltrapobs.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      include 'kobequiv.inc'
      real*10 ctrecf,fdev,reflct,tmdly0,dop0
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (Result(1),tmdly0),(Result(2),dop0)
      real*4    acctim
      equivalence (acctim,Estf)
      include 'param.inc'
      include 'prpgat.inc'
      include 'sbcom.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'
c
c external functions
      real*10 ANCOR,CTATF

c local
      real*10 dum,frcns,sendtm,zsb(3)
      integer   itype,j,jdcns,lemctl,loop1,loop2,lplctl,lsbctl,nn

c     xslcns(1-6,j) = position and velocity (not used) of the center
c                     of mass of the solar system relative to the sun
c                     at receive time (j=1), reflection time (j=2)
c                     and send time (j=3)
c
c          receive time site and earth coordinates determined in radctl
c
c
      if(Nspot.ne.0) call SUICID(' NSPOT NO GOOD IN RADNS ', 6)
c
c------------------setup for interpolation (reflection)-----------------
      lsbctl = -1
      lplctl = -1
      lemctl = 0
      loop1  = 0
      loop2  = 0
      do while(.true.)
c
c          reflection time iteration loop
c       tmdly1, the first guess, was generated in radctl
c*  start=1100
c
         call TIMINC(Jd, ctrecf, Jdx, Fract, -Tmdly1/Secday)
         call SBTRP(itype, Jdx, Fract, 0, lsbctl)
         if(Jdx.le.0) then
            Jd = 0
            return
         endif
         do j = 1, 3
            Xsbsun(j) = Xsb(j)*Aultsc
            zsb(j)    = Xsbsun(j)
         end do
c
c------------------obtain refl. pt. wrt receiv.site at refl. time-------
         if(Klan.ne.0) then
c
c a planet is involved so obtain refl. pt. wrt sun
c
            call PLTRP(1, Jdx, Fract, 0, lplctl)
            if(Jdx.le.0) then
               Jd = 0
               return
            endif
            do j = 1, 3
               Xplsc(j)  = Xp(j)*Aultsc
               Xsbsun(j) = Xsbsun(j)*Cmfct + Xplsc(j)
            end do
         endif
c
c------------------calculate receiving site wrt refl. pt. at refl.time--
c*  start=1400
         do j = 1, 3
            Xsitep(j, 1) = Xemlsc(j, 1) - Xsbsun(j)
         end do
         if(Nswcns.gt.0) then
            call SOTRP(Jdx, Fract, Xslcns(1,2), 0)
            do j = 1, 3
               Xsitep(j, 1) = Xsitep(j, 1)
     .                        + (Xslcns(j,2) - Xslcns(j,1))*Aultsc
            end do
         endif
c
c------------------calculate delay time for reflection------------------
         call UVECTR(3, Xsitep(1,1), Rsitp(1), Xsitp0(1,1), dum)
         Tmdly2 = Rsitp(1) - Xpdly
         loop1  = loop1 + 1
         if(loop1.gt.20)
     .        call SUICID('MORE THAN 20 REFLECT ITERATIONS ', 8)
         if(ABS(Tmdly2-Tmdly1).le.acctim) goto 2000
         Tmdly1 = Tmdly2
      end do
 
c
c*  start=2000
c iteration to determine send site position
 2000 Nit(20) = Nit(20) + loop1
      Tguess  = Tmdly1
      Dstf(4) = Tmdly1
      call TSAV(Tmdly2, 1, Jd, Utrec)
      call DGUESS(tmdly, 3, Jd, Utrec)
      tmdly = tmdly + Tmdly2
 
      do while(.true.)
         Tmdly1 = tmdly
         sendtm = Utrec - Tmdly1*fdev
         Sidtm  = Sidtm0 + Sidvel*sendtm
c in computing sidereal time, tmdly1 converted to wwv seconds as
c approximation to universal time (if tmdly1 left in coordinate
c time,sending site position would be off by 3 centimeters)
         call SITCOR(Sidtm, 2, -1, 0)
         call TIMINC(Jd, ctrecf, jdcns, frcns, -Tmdly1/Secday)
         call ETRP(1, jdcns, frcns, 0, lemctl, 2, 3)
         if(jdcns.le.0) then
            Jd = 0
            return
         endif
c
c------------------determine  sending  site wrt refl. pt. at  send  time
         do j = 1, 3
            Xsitep(j, 2) = Xemlsc(j, 2) - Xsbsun(j)
         end do
         if(Nswcns.gt.0) then
            call SOTRP(jdcns, frcns, Xslcns(1,3), 0)
            do j = 1, 3
               Xsitep(j, 2) = Xsitep(j, 2)
     .          + (Xslcns(j,2) - Xslcns(j,3))*Aultsc
            end do
         endif
c
c calculate delay time for send
         call UVECTR(3, Xsitep(1,2), Rsitp(2), Xsitp0(1,2), dum)
         tmdly = Tmdly2 + Rsitp(2) - Xpdly
         loop2 = loop2 + 1
         if(loop2.gt.20)
     .        call SUICID('MORE THAN 20 SEND ITERATIONS', 7)
         if(ABS(tmdly-Tmdly1).le.acctim) goto 2100
      end do
 
 2100 Nit(19) = Nit(19) + loop2
      call TSAV(tmdly - Tmdly2, 3, Jd, Utrec)
c
c - - - - - - - - - - end of sending iteration  - - - - - - - - - - - -
c
c test for dummy observation below horizon
      if(Idumob.eq.1) then
         nn = -2
         call HORIZN(Xsite, Xsitep, nn, Jd)
c
c test for planetary moon occulted by central body
         call CNTOCC(zsb, Xsitep, -2, Jd)
         if(Jd.le.0) return
      endif
c
c*  start=2500
c
c velocity calculations
      call SITVEL(2, nvel, norm)
 
c velocity of site
      if(nvel.gt.0) call ETRP(1,jdcns,frcns,1,0,2,3)
c
c------------------correct time delay for various effects---------------
c*  start=3000
c
c     general relativistic correction to time delay
      call RADREL(1)
c
c we can skip the corrections if we are only doing doppler
      if(Nice.le.0) then
c
c
         call PROPCO(1, 1)
         tmdly = tmdly + Sumcor(1)
c
c correct for antenna offsets
         tmdly = tmdly + ANCOR(1) + ANCOR(2)
c
c coordinate time delay changed to atomic time delay
         tmdly = tmdly +
     .           (CTATF(Jd,(Ctrec-tmdly-Ctat)/Secday,
     .           Ntrmct,2) - Ctat)
 
c time delay in atomic time changed to time delay in ut time
         if(ntime.gt.0) tmdly = tmdly*fdev
 
c constant bias in time delay
         if(Nrbias.gt.0) tmdly = tmdly + Rbsx(1)
c
c power series for clock errors
         if(Neqnox.gt.0) tmdly = tmdly + Pnox +
     .       (Pequat + Plat*dutrec/2._10)*dutrec
c
c------------------doppler---------------------------------------------
c*  start=4000
c
c       nice.lt.0  no doppler;            nice.gt.0  no time delay
         if(Nice.lt.0 .and. calcvl.le.0) return
      endif
c
c
c       no phase delay doppler for radns
c
c
c      velocity of refl. pt.
      call SBTRP(itype, Jdx, Fract, 1, lsbctl)
      do j = 4, 6
         Xsbsun(j) = Xsb(j)*Aultvl
      end do
      if(Klan.gt.0) then
         call PLTRP(1, Jdx, Fract, 1, lplctl)
         do j = 4, 6
            Xplsc(j)  = Xp(j)*Aultvl
            Xsbsun(j) = Xsbsun(j)*Cmfct + Xplsc(j)
         end do
      endif
c
c
c relative velocity
      call VLRTRD(Xemlsc(4,1), Xsbsun(4), Xemlsc(4,2), 1._10, 7, 0)
      if(Nice.ge.0) then
 
c compute doppler shift
         call RADREL(2)
c
c
c------------------correct doppler shift for various effects------------
c
c
         call PROPCO(2, 1)
         dop = dop + Sumcor(2)
 
         dop = Freq*dop
c
c constant bias in doppler shift
         if(Nrbias.gt.0) dop = dop + Rbsx(2)
c
c power series for clock errors
         if(Neqnox.gt.0) dop = dop - (Pequat + Plat*dutrec)*Freq
      endif
c
c*  start=9990
      return
      end
