      subroutine FERCTL(idopob)
 
      implicit none
 
c*** start of declarations inserted by spag
 
c*** end of declarations inserted by spag
c
c           r. king   june 1977    subroutine ferctl
c           program control for computation of interferometry
c           observables, including both conventional vlbi delay and
c           delay rate (incoherent noise source) and differential
c           cycle count (coherent source).  ferctl calls ferstr, fermn,
c           or fersb depending on whether observed object is a star,
c           s/c in cislunar space, or s/c elsewhere in the solar system.
c           this is a revised version of old pep subroutine interf
c           written by m.ash/r.preston/r.king  feb 70 - jun 75
c
c
c              idopob= 0 when ferctl first called for observation
c              idopob= 1 in ferctl for differential n-count
c                        observable on first call of ferctl
c              idopob= -1 set in ferctl on second call (ferctl called
c                        second time only if idopob was set = 1 on
c                        first call
c

c array dimensions
      include 'globdefs.inc'
c common 
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 ctat2,dutrec
      equivalence (Angdum(9),ctat2),(Angdum(10),dutrec)
      include 'difnct.inc'
      include 'eqnphs.inc'
      real*4    eqnx(3)
      equivalence (eqnx,Pnox)
      include 'fcntrl.inc'
      include 'ltrapobs.inc'
      real*10 dfdly,difnct,ddr
      equivalence (Deriv(2,1),dfdly,difnct),(Deriv(2,2),ddr)
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 tc,freql(2),xfactr,dbias
      equivalence (Save(28),tc),(Save(29),freql),
     .            (Save(33),xfactr),(Save(34),dbias)
      real*10 ctrecf,freq2,freqtr(2)
      equivalence (Dstf,ctrecf),(Dstf(5),freq2),
     .            (Dstf(7),freqtr)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'prpgat.inc'
      include 'scdta.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'
c
c local
      real*10 epsilo/1.E-6_10/,clock,dsec,fruct,fxsb,fxsc,tcsave,xffact
      integer*4 i,idopob,ipct,j,jda,jdxsb,jdxsc,jdxx,k,klb,
     . kspt,norm,npath,npc,nrvfrq,nterm,nvel
      integer*2 npl

c external functions
      real*10 A1UT1,A1WWV,CTATF,UT2UT1

      if(idopob.le.0) then
c
c----------------------setup nrvfrq, nvel, norm-------------------------
         norm   = 0
         nvel   = 0
         nrvfrq = 1
         if(nintrf.ge.0) then
            if(nddiff.gt.0 .and.
     .          (Numsav.lt.39 .or. Save(38).eq.1._10)) nrvfrq = 0
c only for uncondensed alsep ddr observations is beta coded
c using differenced delays
            if(nrvfrq.gt.0) nvel = 1
         end if
         if(Jct(10).gt.0 .or. Jct(11).gt.0) nvel = 1
         if(Nice.ge.0) nvel    = 1
         if(Jct(59).gt.0) nvel = 1
         if(Izctl.ne.0) norm   = 1
         if(Jct(10).gt.0) norm = 2
c           nvel =0 positions determined
c           nvel =1 positions and velocities determined
c           norm =0 normals to sites not determined
c           norm =1 normals to sites determined
c           norm =2 normals to site determined, and partials
c                   determined for solid body tides even if
c                   lsite(1-3)=0
c
c
c----------setup counting interval and reference frequency--------------
c----------for n-count observable---------------------------------------
c
c           idumob=-2  observations on tape
c
         if(Idumob.ne.-2) then
c set save(28-34)
c count interval set in fermtr for dummy observations
            dbias    = 0._10
            freql(1) = Freq
            freql(2) = Freq
 
c save(31-32) are received frequencies set in fercnt
            xfactr = 360._10
            if(nintrf.eq.1) xfactr = 1._10
            if(nrvfrq.eq.0) nrvfrq = 1
 
c nrvfrq = 1 but xmtr freqs are from freq & freq2, not con(16-19)
            if(Numsav.lt.34) Numsav = 34
            freqtr(1) = Freq
            freqtr(2) = 0._10
            if(nddiff.eq.1) then
 
c set save(35-39)
               do i = 35,37
                  Save(i) = 1._10
               end do
               if(Numsav.le.37) Numsav = 37
 
c save(38-39) are received frequencies set in fercnt
               freqtr(2) = freq2
            end if
         end if
 
         Cnttim = Save(28)
c
c if last observation deleted, increase count interval
         if(nast1(1).eq.2) Cnttim = tcsave + Cnttim
         tcsave = Cnttim
c
c save count interval and doppler multiplication factor for
c calculation of partials
         Dstf(9) = Cnttim
         Dstf(6) = Save(33)
c
c------------------setup time quantities--------------------------------
c
         Utrec = Ihr*3600 + Imin*60
         dsec  = Sec
         if(Numsav.ge.40 .and. ABS(Save(40)-dsec).lt.1.E-3_10)
     .       dsec = Save(40)
         Utrec = Utrec + dsec
         npath = 1
         if(nintrf.ge.0) then
            if(Nk1.eq.-1 .or. nrvfrq.eq.0) then
c           npath=1  observation at beginning of count interval
c           npath=2  observation at end of count interval
c
c
c           calculate observable at beginning of interval if not
c           saved from last observation
               Utrec = Utrec - Cnttim
               if(Utrec.lt.0.0_10) then
                  Utrec = Utrec + Secday
                  Jds   = Jds - 1
               end if
c
c beginning time of this observation is the same as final
c time of previous observation - simply update quantities
c and move to final time.
               if(nast1(1).ne.2 .and. Nk1.ne.-1 .and.
     .          Jds.eq.Jdsav .and. ABS(Utrec-Utrsav).le.epsilo) then
                  do k = 1,2
                     do j = 1,2
                        Rsave(1,j,k) = Rsave(2,j,k)
                     end do
                     Difdly(1,k) = Difdly(2,k)
                  end do
                  npath  = 2
                  idopob = -1
                  Nk2    = 2
                  Utrec  = Utrec + Cnttim
                  if(Utrec.ge.Secday) then
                     Utrec = Utrec - Secday
                     Jds   = Jds + 1
                  end if
               end if
               Utrsav = Utrec
               Jdsav  = Jds
            else
               npath = 2
            end if
         end if
      else
         idopob = -1
      end if
c
c calculate time since first observation for clock and
c transmitter frequency terms
      if(Nk1.lt.0) then
         Utrec0 = Utrec
         if(Klanb.gt.0) then
            jdxsb = Sbcom(3)
            fxsb  = jdxsb
            fxsb  = Sbcom(3) - fxsb
         end if
         if(Klans1.gt.0) then
            jdxsc = Sccom(3)
            fxsc  = jdxsc
            fxsc  = Sccom(3) - fxsc
         end if
      end if
      dutrec = Utrec - Utrec0
      if(dutrec.lt.0._10) dutrec = dutrec + Secday
      clock = eqnx(1) + eqnx(2)*dutrec + eqnx(3)*dutrec**2/2._10
      if(Klanb.gt.0) then
         jdxx   = Jds - jdxsb
         Tfrqsb = jdxx*Secday + Utrec - fxsb*Secday
      end if
 
c ct-utc ignored for now in freq. calculations
      if(Klans1.gt.0) then
         jdxx   = Jds - jdxsc
         Tfrqsc = jdxx*Secday + Utrec - fxsc*Secday
      end if
c
c determine time quantities
      if(Jct(59).gt.0) Utrec = Utrec + clock
      fruct = Utrec/Secday
      if(itime.eq.2) then
c  atuts, ututs read in

      else if(itime.eq.0) then
c  observation time is UT2
         Ututs = -UT2UT1(Jds,fruct)
         Atuts = A1UT1(Jds,fruct) + Ututs

      else
c  observation time is UTC
         Atuts = A1WWV(Jds,fruct)
         Ututs = Atuts - A1UT1(Jds,fruct)
      end if
c ututs = ut1 - given observation time
c atuts = a.1 - given observation time
c save a.1 jd and fraction of day for ctatf calls
      call TIMINC(Jds,fruct,jda,Frect,Atuts/Secday)
 
c ctat=32.15 for determining precession nutation and site coords
      Ctat  = 32.15_10
      ctat2 = Ctat
      call TIMINC(jda,Frect,Jd,Fract,Ctat/Secday)
 
c jd, fract are now ct, not utc
      Utrec = Utrec + Ututs
      if(nddiff.ge.0) then
c
c           if save(48).ne.1  (default value), offset sample time
c           at second site by save(48) seconds.  for n-count obser-
c           vables, this compensates for the difference in propagation
c           time to the two sites. this difference is at most 20 ms, so
c           a common utrec is used to calculate at-utc, ct-at, and
c           precession-nutation at both sites.
         Utrec2 = Utrec
         if(Jct(59).gt.0) Utrec2 = Utrec2 - clock
 
c clock terms apply only to first site
         if(Numsav.ge.48 .and. Save(48).ne.1._10)
     .       Utrec2 = Utrec2 + Save(48)
      end if
c utrec stored in /yvect/, utrec2 in /difnct/
c
c obtain mean sideral time at midnight and the
c (sideral time)/ut ratio for day of observation
      call SIDTIM(Jds,Ctat+Atuts-Ututs,Sidtm0,Sidvel,Dera)
      call PLATMO(Jds)
c
c read moon tape or nbody tape to get nutation angles
      call MNREED(Jd)
      if(Jd.le.0) return
c
c nutation-precession determined for receiving
      Kindnp = ndprec*(Ncode - 1)
      call PRCNUT(Jd,Fract)
c
c determination of geocentric coordinates of first site
c calculate diurnal polar motion
      if(mod(Jct(28)/4,2).ne.0) then
         call DIURNP(Jd,Fract,(Sidtm0+Utrec*Sidvel),Dtheta)
         Sidtm0 = Sidtm0 - Dtheta
      end if
      Sidtm = Sidtm0 + Utrec*Sidvel + Dgst
      call SITCOR(Sidtm,1,nvel,norm)
c
c determination of geocentric coordinates of second site for
c differential n-count observable
      if(nintrf.ge.0 .and. nddiff.ge.0) then
         Sidtm2 = Sidtm0 + Utrec2*Sidvel + Dgst
         if(mod(Jct(28)/4,2).gt.0) Sidtm2 = Sidtm2 - Dtheta
         call SITCOR(Sidtm2,2,nvel,norm)
      end if
c
c determine ct-at including diurnal and monthly & annual terms
      if(prmter(81).ne.0._10) then
         nterm = 3
         if(Klan.eq.17 .and. Ict(27).ne.0) nterm = 2
c fermn does not yet include diurnal "ct-at" terms in calculations
c so we need to incorporate them here
         Ctat = CTATF(jda,Frect,nterm,1)
         if(nintrf.ge.0 .and. nddiff.ge.0)
     .       ctat2 = CTATF(jda,Frect,nterm,2)
 
c ctat saved in /obscrd/, ctat2 in /coord/
         call TIMINC(jda,Frect,Jd,Fract,Ctat/Secday)
      end if
      ctrecf = Fract
      Ctrec  = ctrecf*Secday
      if(nintrf.ge.0 .and. nddiff.ge.0) Ctrec2 = Ctrec +
     .    (Utrec2 - Utrec) + (ctat2 - Ctat)
c     jd will be wrong for ctrec2 only if utrec is within 20 ms
c     of 32.15 + atuts seconds after midnight
c
c
c-----------calculate differential delays for observed object(s)
c
c           for double difference observable, do second object first
      npl  = Nplnt0
      kspt = 1
      klb  = Klanb
      npc  = Ncp0
      if(nddiff.gt.0) then
         kspt = 2
         npl  = Nplnt2
         klb  = Klans1
         npc  = Ncs1
      end if
c
c determine if observed object is a star
  100 if(npl.gt.0) then
c
c determine if object is moon or probe in cislunar space
         if(Klan.eq.17) then
            call FERMN(nvel,norm,npath,kspt)
         else if(klb.le.0) then
            if(Nspot.le.0) call SUICID(
     .'UNABLE TO CALCULATE RADAR INTERFEROMETRY OBSERVABLE, STOP IN FERC
     .TL ', 17)
c
c observed object is probe not in cislunar space
            call FERSB(nvel,norm,npath,kspt)
         else if(npc.ne.3) then
            call FERSB(nvel,norm,npath,kspt)
         else
            call FERMN(nvel,norm,npath,kspt)
         end if
 
         if(Jd.le.0) return
      else
         call FERSTR(nvel,norm,kspt)
      end if
c
c is this a double difference observable?
c if so, call ferstr, fermn, or fersb a second time
      if(nddiff.gt.0) then
         if(kspt.ne.1) then
            kspt = 1
            npl  = Nplnt0
            klb  = Klanb
            npc  = Ncp0
            goto 100
         end if
      end if
c
c
c-----------form theoretical observable------------------------------
c
c              ( dfdly and difcnt equivalenced to deriv(2,1) )
c              ( ddr              equivalenced to deriv(2,2) )
c
c           form conventional vlbi observables
      if(nintrf.lt.0) then
         if(Nice.le.0) then
            dfdly = Difdly(1,1)
 
c nddiff.eq.1:  differential vlbi observable
            if(nddiff.eq.1) dfdly = dfdly - Difdly(1,2)
         end if
         if(Nice.ge.0) then
c npath always=1 for conventional vlbi observable, so
c rates stored in difdly(2,kspt) slot
            ddr = Difdly(2,1)
            if(nddiff.eq.1) ddr = ddr - Difdly(2,2)
         end if
c
c form counted-cycle vlbi observables
      else if(npath.eq.2) then
c
c if end of count interval, store propagation media
c corrections, calculate transmitter frequencies, and
c form theoretical observable
         call FERCNT(nrvfrq)
         call PROPCO(-1,4)
      else
c
c if beginning of count interval, store propagation media
c corrections, obtain transmit frequencies, update times
c and return to fermtr with idopob=1
         call PROPCO(-2,4)
         call FERCNT(nrvfrq)
         Utrec = Utrsav + Cnttim
         npath = 2
         Jds   = Jdsav
         if(Utrec.ge.Secday) then
            Utrec = Utrec - Secday
            Jds   = Jds + 1
         end if
         idopob = 1
         Jdsav  = Jds
         Utrsav = Utrec
         return
      end if
c
c
c
c-----------calculate and apply propagation media corrections-----------
c
c           ipct= 0 apply corrections to accumulated cycle-count
c                   observable (corrections to delays calculated above
c                   with ipct= -2,-1)
c           ipct= 1 calculate corrections to differential delay
c           ipct= 2 calculate corrections to differential delay rate
c
      if(Nice.le.0) then
         ipct = 0
         if(nintrf.lt.0) ipct = 1
         call PROPCO(ipct,4)
         Deriv(2,1) = Deriv(2,1) + Sumcor(1)
 
c deriv(2,1) equivalenced to dfdly or difnct
         if(Nice.lt.0) goto 200
      end if
      ipct = 2
      call PROPCO(ipct,4)
      ddr = ddr + Sumcor(2)
c
c
c-----------apply clock corrections to observable----------------
c
  200 if(Nice.le.0) then
         if(Neqnox.gt.0) then
            xffact = 1._10
            if(nintrf.eq.0) clock = clock - eqnx(1)
 
c constant term vanishes if initial value of phase subtracted
            if(nintrf.ge.0 .and. Jct(59).eq.1) xffact = -
     .          (freqtr(1) - freqtr(2))*xfactr
c if jct(59)=0  eqnx(1-3) are corrections to observable
c (delay in seconds, phase in degrees or cycles)
c if jct(59)=1  eqnx(1-3) are clock terms (seconds), also
c applied to utrec above
            dfdly = dfdly + clock*xffact
         endif
 
c difnct= dfdly by equivalence
         if(Nrbias.gt.0) difnct = difnct + Rbsx(1)
 
c allow rbias to be used also for initial sample error for ncount
         if(Nice.lt.0) return
      endif
 
c nice always .lt.0 for counted-cycle vlbi observations
      if(Neqnox.gt.0) ddr = ddr + eqnx(2) + eqnx(3)*dutrec
c
c
      return
      end
