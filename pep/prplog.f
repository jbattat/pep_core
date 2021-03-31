      subroutine PRPLOG(ipct,kick)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, i1, ip, ipart, ipc, ipct, ipu, itp1, itp2, jct3, 
     .          kick, lxp
 
c*** end of declarations inserted by spag
 
 
c
c        r.b. goldstein  march 1978
c        routine to set up uscor and calcor for propco
c        must only be called for the following values of ipct:
c             -2, 1, 2
c
c           kick= 1  prplog called by radar
c               = 4  prplog called by fermtr
c
c
c        this works via a series of pointer arrays, that point
c        to the following:
c             pxnudy:   external, neutral, time delay
c             pxchdy:   external, charged particle, time delay
c             pxnudo:   external, neutral, doppler
c             pxchdo:   external, charged particle, doppler
c             pinudy:   internal, neutral particle, time delay
c             pichdy:   internal, charged part., time delay
c             pinudo:   internal, neutral part., doppler
c             pichdo:   internal, charged part., doppler
c
c        the main loop for the "standard" logic is used for both
c        time delay and doppler and neutral and charged particle
c        corrections. it is made specific by changing the indices on
c        the arrays pint and pext which are equivalenced to the above
c        arrays.
c
c        each of the above pointer arrays is ordered so that the
c        "most desirable" correction type comes first. this order is
c        alterable by prplog in the case of active v.s. static
c        corrections.
c
c        common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obscrd.inc'
      include 'prpgat.inc'
c
c all doppler and time delay corections
      integer*2 paldly(25)/1,3,5,7,9,11,13,15,17,19,21,23,25,12*0/,
     .          paldop(25)/2,4,6,8,10,12,14,22,17*0/
      integer*2 pall(25,2)
      equivalence (pall(1,1),paldly),(pall(1,2),paldop)
c
c
c
c internal doppler and time delay corrections
      integer*2 pindly(10)/1,3,5,7,9,19,4*0/,
     .          pindop(10)/2,4,6,8,10,5*0/
      integer*2 pin(10,2)
      equivalence (pin(1,1),pindly),(pin(1,2),pindop)
c
c external corrections for doppler, delay, neutral and ch. part.
c order determines desirability
c Note: no planetary atmosphere or interstellar corrections here.
c
      integer*2 pxnudy(10)/10*0/,
     .          pxchdy(10)/21,11,15,13,6*0/,
     .          pxnudo(10)/10*0/,
     .          pxchdo(10)/22,12,14,7*0/
      integer*2 pext(10,2,2)
      equivalence (pext(1,1,1),pxnudy),(pext(1,2,1),pxchdy),
     .            (pext(1,1,2),pxnudo),(pext(1,2,2),pxchdo)
c
c
c        internal corrections for doppler, tmdly, neutral and charged
c        particles. order determines desirability. the order is
c        alterable by prplog in the case of the active vs. the static
c        corrections.
      integer*2 pinudy(10)/1,3,8*0/,
     .          pichdy(10)/9,5,7,7*0/,
     .          pinudo(10)/2,4,8*0/,
     .          pichdo(10)/10,6,8,7*0/
      integer*2 pint(10,3,2)
      equivalence (pint(1,1,1),pinudy),(pint(1,2,1),pichdy),
     .            (pint(1,1,2),pinudo),(pint(1,2,2),pichdo)
c
c internal planetary atmosphere corrections
      integer*2 pipady(10)/19,9*0/,
     .          pipado(10)/10*0/
      equivalence (pint(1,3,1),pipady),(pint(1,3,2),pipado)
c
c
c special corrections
c
      integer*2 pspecl(10,2),pspdy(10)/17,23,8*0/,pspdo(10)/10*0/
      equivalence (pspecl(1,1),pspdy),(pspecl(1,2),pspdo)
c
c
c corrections that require zenith angle and/or rate
      integer*2 pzen(15)/1,2,3,4,5,6,7,8,7*0/
c
c ranks
      integer*2 rnk(100)/4*2,4*3,2*1,2*20,2*17,18,0,5,0,2,0,2*25,5,0,
     . 1,75*0/
c
c
c
c
      logical*1 ovride, init/.false./
      integer*2 typ
c
c
c*  start=100
c
      if(.not. (init)) then
         jct3 = Jct(3)
         ipu  = mod(jct3/2,2)
         ipc  = mod(jct3/4,2)
         init = .true.
      endif
c
c
c        initial logic
c
c        set all uscor and calcor to .false.
      do i = 1, 100
         Uscor(i)  = .false.
         Calcor(i) = .false.
      end do
      Izctl = 0
c
c get typ: 1=tmdly    2=dop
      itp1 = 1
      itp2 = 2
      if(Nice.gt.0) itp1 = 2
      if(Nice.lt.0) itp2 = 1
c
c loop to set uscor from jclten
      do typ = itp1, itp2
         ovride = .false.
         do i = 1, 100
            ip = pall(i,typ)
            if(ip.eq.0) goto 50
            if(Jclten(ip).eq.3) then
               Uscor(ip) = .true.
               ovride    = .true.
            endif
         end do
c
c loop to set calcor from jclone
c loop only over internally calculated corrections
   50    do i  = 1, 100
            ip = pin(i,typ)
            if(ip.eq.0) goto 100
            if(Jclone(ip).ne.1) then
               if(Jclone(ip).eq.3) Calcor(ip) = .true.
               if(Uscor(ip) .and. (Ical(ip).eq.0)) Calcor(ip) = .true.
            endif
         end do
c
c test to see if normal logic was overridden
  100    if(.not. (ovride)) then
c
c
c*  start=500
c-----------------------------------------------------------------------
c        "normal" logic
c        this code is only entered if there were no 3's set in the
c        10's digit of any jcal
c        it applies to both tmdly and dop according to "typ"
c        it applies to  charged particles, neutral atmosphere, and
c        planetary atmosphere, depending on ipart:
c                  ipart:   1=neutral   2=charged   3=planet atm.
c
c
c         test to see if active neutral atm. is more
c        desirable than passive.
c
            pinudy(1) = 1
            pinudy(2) = 3
            pinudo(1) = 2
            pinudo(2) = 4
            if(Numsav.ge.46 .and. Save(45).gt.1._10) then
               pinudy(1) = 3
               pinudy(2) = 1
               pinudo(1) = 4
               pinudo(2) = 2
            endif
c
c
            do ipart = 1, 3
c No external corrections for planetary atmospheres.
               if(ipart.eq.3) goto 110
c
c external calcs. if an external cal. exists, and its not
c unconditionally shut, then use it and dont bother with internal
c calibrations
               do i = 1, 100
                  ip = pext(i,ipart,typ)
                  if(ip.eq.0) goto 110
                  if((Ical(ip).ne.0) .and. (Jclten(ip).ne.1))
     .               then
                     Uscor(ip) = .true.
                     goto 120
                  endif
               end do
c
c        internal cals. once in this loop, there are no external
c        corrections. choose the most desirable internal correction
c        that is not shut off with a 1 in the 10's place
c        if it exists on tape, do not re-calculate it.
  110          do i  = 1, 100
                  ip = pint(i,ipart,typ)
                  if(ip.eq.0) goto 120
                  if(Jclten(ip).ne.1) then
                     Uscor(ip) = .true.
                     if(Ical(ip).eq.0) Calcor(ip) = .true.
c
c this line is to make it ok for plasma to exist with ionosphere
                     if((ip.ne.9) .and. (ip.ne.10)) goto 120
                  endif
               end do
c
c
  120       end do
c
c
c*  start=600
c        special corrections - as of now, only calibrations
c        external only - can use all in this class
            do i = 1, 100
               ip = pspecl(i,typ)
               if(ip.eq.0) goto 150
               if((Jclten(ip).ne.1) .and. (Ical(ip).ne.0))
     .            Uscor(ip) = .true.
c
c
c        here put a patch for internal special corrections
c        for example antenna corrections
c
c-----------------------------------------------------------------------
c*  start=1000
c
            end do
         endif
c
c
c        out of "normal" logic section
c        test if calcor must be shifted for phase delay doppler
c        radar link only, internal only, calcor only, ipct=-2 only
  150    if(kick.eq.1) then
            if(ipct.eq.-2) then
               do i = 1, 100
                  ip = pindop(i)
                  if(ip.eq.0) goto 200
                  if(Calcor(ip)) then
                     i1 = ip - 1
                     Calcor(i1) = .true.
                     Calcor(ip) = .false.
                  endif
               end do
            endif
         endif
c
c see if zenith angle and/or rate needs to be computed
  200    do i  = 1, 100
            ip = pzen(i)
            if(ip.eq.0) goto 250
            if(Calcor(ip)) then
               Izctl = 1
               if(Nice.ge.0 .and. ipct.ge.0) Izctl = 3
               goto 250
            endif
         end do
c
c
c get ranks of all internal cals
  250    do i  = 1, 100
            ip = pin(i,typ)
            if(ip.eq.0) goto 300
            if(Uscor(ip)) Ical(ip) = rnk(ip)
         end do
 
c*  start=1500
  300 end do
c
c*  start=2000
c        set ncal by searching backwards for first non zero
c        value of ical
c
      do i = 1, 100
         Ncal = 101 - i
         if(Ical(Ncal).ne.0) goto 400
      end do
      Ncal = 0
c
c
  400 lxp = Ncal/50 + 1
      if(ipu.ne.0) then
         if(Line + lxp.gt.56) call OBSPAG
         write(Iout,450) (Uscor(i),i = 1,Ncal)
  450    format(' USCOR : ', 50L2/ 8x, 50L2)
         Line = Line + lxp
      endif
      if(ipu.ne.0) then
         if(Line + lxp.gt.58) call OBSPAG
         write(Iout,500) (Calcor(i),i = 1,Ncal)
         Line = Line + lxp
  500    format(' CALCOR: ', 50L2/ 8x, 50L2)
      endif
 
      return
      end
