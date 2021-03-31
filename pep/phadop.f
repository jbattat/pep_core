      subroutine PHADOP(npath,idopob)
 
      implicit none
c
c       sept. 1976..a derivative of sbdldp (sbbg)..r. goldstein
c        calculate doppler observable by phase delay difference method.
c
c arguments
      integer*4 idopob,npath

c        note....it is not possible at present to calculate phase delay
c          doppler and time delay observables in the same series. if
c          this feature is to be added in the future, the npath logic,
c          idopob logic, utrsav logic, and logic in calling propco all
c          must be checked and changed. it will be necessary to call
c          radctl three times: in the middle (for tmdly), and 2 ends
c          (for phase delay doppler) of the count interval.

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'comdat.inc'
      real*10 secday,gmc2
      equivalence (Comcon(73),secday)
      equivalence (Comcon(90),gmc2)
      include 'coord.inc'
      include 'difnct.inc'
      real*10 delay(2)
      equivalence (delay(1),Difdly(1,1))
      include 'ltrapobs.inc'
      real*10 tmdly,dop
      equivalence (Deriv(2,1),tmdly),(Deriv(2,2),dop)
      include 'number.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct,tmdly0,dop0,tc
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (Result,tmdly0),(Result(2),dop0)
      equivalence (Save(28),tc)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'sitcrd.inc'
      include 'spqind.inc'
      include 'yvect.inc'

c external functions
      real*10 DOT

c local 
      real*10 radmsv(2),pos(2,2),vel2(2,2)
      real*10 dum,favg,favg1,favg2,t2s 
      integer   i,idum,numsit
c
c setup
      if(Nk1.lt.0) then
         if(Atctsc.gt.0._10) then
            favg  = gmc2*1.5_10/Aultsc
            favg1 = favg + Vswrte(1)**2/2._10
            favg2 = favg + Vswrte(2)**2/2._10
         endif
         Fava   = 0._10
         Favb   = 0._10
         numsit = 2
         if(Nsite2.eq.0) numsit = 1
      endif
c
c
c
c       the first time through, npath is 1. pos, delay, and vel2 are
c       saved, utrec is updated, idopob is updated, npath is incremented
c       and then it exits. the next time it enters phadop, npath is 2.
c       the pos, vel2, and delay are saved, and then phase delay doppler
c       is calculated. npath is set to 1 in anticipation of the next
c       observation.
c
      delay(npath) = tmdly
      do i = 1, numsit
         vel2(i,npath) = DOT(Xemlsc(4,i),Xemlsc(4,i))
      end do
      pos(1,npath) = Rs1m
      pos(2,npath) = Rs2m
      radmsv(npath) = Raddum(1)
      Jd1s(npath)   = Jd
      Ctrcs(npath)  = ctrecf
      Rsave(npath,1,1) = Dstf(4)
      if(npath.eq.2) then
c
c
         npath = 1
c
c     compute average derivative of atomic time w/r/t coordinate time
c     at receive (a) and send (b) sites
c           maybe check atctsc here for consistency
c
         if(Atctsc.gt.0._10) then
            Fava = favg1 - (gmc2*(1._10/pos(1,1)+1._10/pos(1,2))/2._10
     .              + (vel2(1,1)+vel2(1,2))/4._10)
            if(Nsite2.ne.0) Favb = favg2 -
     .       (gmc2*(1._10/pos(2,1)+1._10/pos(2,2))/2._10
     .       + (vel2(2,1)+vel2(2,2))/4._10)
         endif
c
c compute doppler by phase delay difference method
c
         Tcs = tc
         if(ntime.gt.0) Tcs = Tcs/fdev
         Tcs = Tcs/(1._10 + Favb)
         dop = (delay(1) - delay(2))/Tcs + (Favb - Fava)
         if(ntmdly.gt.0) Raddum(6) = (radmsv(1) - radmsv(2))/tc
         Cnttm2 = ((Jd1s(2)-Jd1s(1)) + (Ctrcs(2)-Ctrcs(1)))/2._10
         t2s    = Jdx + Fract
      else
c
c update utrec to end of count interval
c
         npath = npath + 1
         Utrec = Utrsav + tc
         Jd    = Jdsav
         if(Utrec.ge.secday) then
            Utrec = Utrec - secday
            Jd    = Jd + 1
         endif
         idopob = 1
         Jdsav  = Jd
         Utrsav = Utrec
      endif
      return
 
      entry PHDMVE(npath,idopob)
c
c     beginning time of this observation is the same as final time of
c     previous observation. simply update quantities and move to final
c     time.
c
c           check 1st if filter epoch crossed (sat. obs.)
      if(Klanb.gt.0) call SBRED2(t2s + tc/1728E2_10)
      if(.not. Recalc) then
c
c transfer quantities
         Contig     = .true.
         pos(1,1)  = pos(1,2)
         pos(2,1)  = pos(2,2)
         vel2(1,1) = vel2(1,2)
         vel2(2,1) = vel2(2,2)
         delay(1)   = delay(2)
         radmsv(1)  = radmsv(2)
         Jd1s(1)    = Jd1s(2)
         Ctrcs(1)   = Ctrcs(2)
         Rsave(1,1,1) = Rsave(2,1,1)
         npath  = 2
         idopob = -1
         Utrec  = Utrec + tc
         call PRPSTR(3,idum,1,1,dum)
         if(Utrec.ge.secday) then
            Utrec = Utrec - secday
            Jd    = Jd + 1
         endif
      else
 
c crossed epoch, must calculate path 1 afresh
         Recalc = .false.
      endif
c should also put deriv1 into common
c and shift deriv -> deriv1 here
c
      return
      end
