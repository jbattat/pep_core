      subroutine CMPAR2(mocpar)
 
      implicit none

c amuchastegui/ash - april 1970 - subroutine cmpar2

c arguments
      integer*4 mocpar

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'bdctrl.inc'
      include 'comdat.inc'
      include 'empcnd.inc'
      include 'eqenox.inc'
      include 'eqnphs.inc'
      real*4    eqnx(3)
      equivalence (eqnx, Pnox)
      integer*4 zeqnph/984/   !r8=528,r10=984
      include 'fcntrl.inc'
      include 'inodta.inc'
      integer*4 i2bod
      equivalence (Jpert, i2bod)
      include 'ltrapx.inc'
      include 'mnsprt.inc'
      integer*4 zmnspr/4080/   !r8=2072,r10=4080
      include 'mtrapx.inc'
      include 'namtim.inc'
      include 'number.inc'
      integer*2 nsite(2)
      equivalence (nsite, Nsite1)
      include 'obscrd.inc'
      include 'obstap.inc'
      include 'phase.inc'
      include 'plndta.inc'
      include 'psrstf.inc'
      include 'rdbias.inc'
      include 'rotdta.inc'
      include 'sbdta.inc'
      include 'scdta.inc'
      include 'sitcrd.inc'
      character*4 sitf4(2,2)
      equivalence (Sitf, sitf4(1,1))
      include 'skymap.inc'
      include 'skystf.inc'
      integer*4 maxnsk/80/
      include 'sptcrd.inc'
      include 'stats.inc'
      character*4 site14(2,2,2)
      equivalence (site1(1,1),site14(1,1,1))
      include 'stcord.inc'
      include 'tapdta.inc'
      include 'zeroes.inc'

c external functions
      integer*4 ITYPOB,LEG,NSCAN
c
c local
      character*8 blank8/'        '/
      character*4 blank, ap15/'AP15'/,tspot
      equivalence (blank8,blank)
 
      character*48 gmess/
     . '        NO DATA SET FOR BODY***, STOP IN CMPAR2 '/
      character*4 amper9/'&&&&'/, pound9/'####'/
      integer i,i4,j,j1,j2,jtypob,jut1,jwob,k,l,n,n1,ntapf
c
c           nk1=-1 for first observation of series, used for saved
c           partials logic.  nk1 set = 0 at end of first call to partl
c           for a series and further incremented in fermtr for proper
c           handling of quasi-vlbi partials
      Nk1 = -1
c     nk1 used to have something to do with velocity calulations in
c     radar link, so there might be incorrect comment cards in other
c     parts of pep that are no longer valid
c
c*  start=1000
c
      jtypob = ITYPOB(Ncodf)
c
c initialize bias quantities
      call ZFILL(Xpp, zeqnph)
c
c determine klan and klanb for observed body
      Mnplnt = Nplnt0 - 10
 
c npshp1=0  (done after planet tape is rewound)
      if(Nplnt0.ge.0) then
         if(Nplnt0.eq.3) goto 200
         if(Nplnt0.eq.10) then
            Klan = u_mxpl+1
            Ncp0 = 3
            Klap = u_mxpl+1
            goto 200
         else if(Nplnt0.gt.0) then
            do i = 1, Numpln
               if(Nplnt(i).eq.Nplnt0) then
                  Klap = i
                  Ncp0 = Npcent(i)
                  if(Nplnt0.le.30 .and. Ncp0.le.0 .and. 
     .             Inttyp(i).ne.0) then
                     Klan = i
                     goto 100
                  endif
                  Klanb = i
                  Klan  = 0
                  if(Ncp0.gt.0) then
                     if(Ncp0.eq.10) then
                        Klan = u_mxpl+1
                     else if(Ncp0.eq.3) then
                        Klan = u_mxpl+2
                     else
                        do j = 1, Numpln
                           if(Ncp0.eq.Nplnt(j)) then
                              Klan = j
                              goto 100
                           endif
                        end do
                        call MVC(' NO CENTRAL BODY',1,16, gmess,1)
                        i4 = Ncp0
                        call EBCDIX(i4, gmess, 17, 3)
                        goto 20
                     endif
                  endif
                  goto 100
               endif
            end do
   20       i4 = Nplnt0
            call EBCDIX(i4, gmess, 29, 3)
            call SUICID(gmess, 13)
         else
            Klan = 18
            Klap = 18
            goto 200
         endif
      endif
c
c           determine klans1 for a second observed body (ncodf.gt.20)
c                klans1 is largely a flag to indicate tape reading
c                for a second object - the object itself is nps1
c             note:  there are two other ways to set up klans1.
c                    (1) ncodf=7,8,9,13,14, spot names second object.
c                        see below at label=1460
c                    (2) ksite.gt.0, second object is on planet ksite
c                        (this logic works only for ncodf=1,2,3,4,5,6,
c                         19, and 20.)
c                          see below at label=1400
  100 if(Ncodf.gt.20 .and. Nplnt2.ne.0 .and. Nplnt2.ne.3) then
         do i = 1, Numpln
            if(Nplnt(i).eq.Nplnt2) then
               Klans1 = i
               Nps1   = Nplnt2
               Ncs1   = Npcent(i)
               goto 200
            endif
         end do
         write(Iout,110) Nplnt2
  110    format(' ***ERROR*** SECOND OBSERVED OBJECT',I3,' NOT INPUT')
         call SUICID(
     .    'SECOND OBSERVED OBJECT NOT INPUT, STOP IN CMPAR2',12)
      endif
c
c determine klanr for planet rotation
  200 if(Klan.gt.0 .and. Klan.le.u_mxpl) then
         do i = 1, Numpln
            if(Nplnt(Klan).eq.-Nplnt(i)) then
               Klanr = i
               goto 300
            endif
         end do
      endif
c*  start=1200
c
c           determine site numbers, partial derivative controls
c           nsite(1) is receive site (radar or optical)
c           nsite(2) is send site (radar only)
  300 do i = 1, 2
         Sitf(i) = blank8
      end do
      if(jtypob.eq.3) then
c
c transfer 'sending site' (indicates type of occultation)
         sitf4(1, 2) = Sita2(Ntape)
      else
         do i = 1, 2
            ntapf = Ntape
            if(i.gt.1) then
               if(Sita2(Ntape).eq.blank .or. jtypob.eq.2) goto 400
               if(Ncodf.eq.19) goto 400
               ntapf = Ntape + 7
            endif
            do k = 1, Numsit
               if(Sita1(ntapf).eq.Site(1,k)) then
                  nsite(i) = k
                  Ksite(i) = Kscrd(k)
                  Sitf(i)  = sitd(k)
                  T0sit(i) = T0site(k)
                  call LVTBDY(Lscrdx(1,i), Lscrd(1,k), Mscrdx(1,i), -6)

c lsite controls the computation of partials in sitcrd
                  do j=1,6
                     Coords(j, i) = Scord(j, k)
                     if(Iabs1.le.0 .or. Mscrdx(j,i).le.0)
     .                   Lsite(j, i) = Lscrdx(j, i)
                  end do
                  goto 350
               endif
            end do
            if(Ntape.lt.3) then
               if(Nrewna(Ntape).ge.0 .and. Nrewna(Ntape).le.1) then
                  write(Iout, 310) Ntape, Nrewna(Ntape), i, Ntape,
     .                             Sita1(ntapf)
  310             format('0NREWNA(', i1, ')=', i3,
     .' DOES NOT INDICATE THE EXISTENCE OF COORDINATES FOR NON INPUT SIT
     .E SITA', i1, '(', i1, ')=', 1A4)
                  call SUICID(' NO SITE COORDINATES, STOP IN CMPAR2',9)
               endif
               sitf4(2,i) = site14(2,i,Ntape)
            else
               sitf4(2,i) = Sita1(ntapf+1)
            endif
            sitf4(1,i) = Sita1(ntapf)
            Ksite(i)   = Kscrd1(i,Ntape)
            T0sit(i)   = T0St1(i,Ntape)
            do j=1,6
               Coords(j,i) = Scord1(j,i,Ntape)
            end do
  350    end do
      endif
c
c*  start=1400
c determine if observing site is not on earth
  400 if(Ncodf.le.20 .and. Ksite(1).gt.0 .and. Ksite(1).ne.3) then
         Nps1 = Ksite(1)
         if(Nps1.eq.10) then
            Klans1 = u_mxpl+1
            Ncs1   = 3
         else
            do i = 1, Numpln
               if(Nplnt(i).eq.Nps1) then
                  Klans1 = i
                  Ncs1   = Npcent(i)
                  goto 500
               endif
            end do
            write(Iout,410) Nps1
  410       format(' ***ERROR*** OBSERVING SITE ON PLANET',I3)
            call SUICID(
     .' OBSERVING SITE ON NON-INPUT PLANET, STOP IN CMPAR2 ',13)
         endif
      endif
c
c determine spot number, partial derivative controls
  500 call ZFILL(Rot, zmnspr)
      call LVTPRM(Lpsrx, izero2, Mpsrx, u_nmpsr)
      if(Spotf.eq.blank .or. Spotf.eq.amper9 .or. Spotf.eq.pound9)
     .    goto 1000
      do k = 1, Numpsr
         if(Spotf.eq.Sptpsr(k)) then
            Nplsr = k
            do j=1,u_nmpsr
               Psrprm(j) = Psrcn(j, k)
            end do
            call LVTPRM(Lpsrx, Lpsrcn(1,k), Mpsrx,u_nmpsr)
            Jdps0  = Jdpsr0(k)
            Plsper = Plspr(k)
            Ntyps  = Ntypsr(k)
            goto 600
         endif
      end do
      if(Ncodf.eq.18) then
         write(Iout,510) Spotf
  510    format(' ***ERROR*** OBSERVED PULSAR ',A4,' NOT INPUT')
         call SUICID('OBSERVED PULSAR NOT INPUT, STOP IN CMPAR2   ',11)
      endif
 
  600 tspot=Spotf
      if(jtypob.ne.3 .and. jtypob.ne.5) then
         if(Numspt.le.0) goto 900
 
c spot name indicates observed spot or star
         l = 1
      else
c spot name gives 2nd body for undifferenced observables
c compare from 1st non-blank character in each
         j1 = NSCAN(Spotf, 4, blank)
         n1 = 4 - j1
         do i = 1, Numpln
            j2 = NSCAN(Aplnt(i), 8, blank)
            if(j2.ge.0) then
               n = min0(n1, 8-j2)
               if(LEG(n,j1+1,Spotf,j2+1,Aplnt(i)).eq.0) then
                  Klans1 = i
                  Nps1   = Nplnt(i)
                  Ncs1   = Npcent(i)
                  goto 1200
               endif
            endif
         end do
         goto 900
      endif
  700 do k = 1, Numspt
         if(tspot.ne.Spot(k)) goto 800
         if(l.eq.2) then
            Nspot2 = k
         else
            if(Nsplnt(k).ne.Nplnt0) then
               write(Iout,710) tspot,Nsplnt(k),Nplnt0
  710          format(' ***ERROR*** SPOT ',A4,' ON BODY',I3,' NOT',I3)
               call SUICID(
     .'OBSERVED SPOT IS NOT ON CORRECT PLANET, STOP IN CMPAR2  ',14)
            endif
            Nspot = k
         endif
         call LVTBDY(Lspcdx(1,l), Lspcrd(1,k), Mspcdx(1,l), -6)
         do j = 1,6
            Spcdx(j, l) = Spcord(j, k)
            if(Iabs1.le.0 .or. Mspcdx(j,l).le.0) Lspot(j,l)= Lspcdx(j,l)
         end do
         T0spt(l)=T0spot(k)
 
         if(l.eq.2) goto 1100
         goto 1000
  800 end do
  900 write(Iout,910) tspot
  910 format(' ***ERROR*** OBSERVED SPOT ',A4,' NOT INPUT')
      call SUICID(' OBSERVED SPOT IS NOT INPUT, STOP IN CMPAR2 ', 11)
 1000 if(Spotf2.ne.blank) then
         l = 2
         tspot=Spotf2
         goto 700
      endif
 
c offset apollo 15 coordinates if alsep rather than lrrr observed
 1100 if(Ncodf.gt.3) then
         if(Spotf.eq.ap15) then
            Spcdx(2, 1) = Spcdx(2, 1) + 1.4330E-3_10
            Spcdx(3, 1) = Spcdx(3, 1) + 0.7361E-3_10
         else if(Spotf2.eq.ap15) then
            Spcdx(2, 2) = Spcdx(2, 2) + 1.4330E-3_10
            Spcdx(3, 2) = Spcdx(3, 2) + 0.7361E-3_10
         endif
      endif

c set up for analytic dissipation terms in lunar libration
c note that no calculation is needed if the coefficients are all zero
c and saved partials are available to be copied
      if(Nplnt0.eq.10) then
         if(Mrcond(14).ne.0._10 .or. Mrcond(15).ne.0._10 .or.
     .    Mrcond(16).ne.0._10) Dodiss=.true.
         do i=7,u_nmbod
            if(Lmrx(i).eq.8 .or. Lmrx(i).eq.9 .or. Lmrx(i).eq.10) then
               if(Iabs1.eq.0) then
                  Dodiss = .true.
               else
                  do j=7,u_nmbod
                     if(Mmrx(j).eq.Lmrx(i)) goto 1140
                  end do
                  Dodiss = .true.
 1140             continue
               endif
            endif
         end do
      endif
 
c check consistency of klans1
 1200 if(Klans1.gt.0) then
 
c test klans1 to see if actual 'sc' tape
         if(Klans1.eq.Klanb .or. Klans1.eq.Klan) Klans1 = 0
 
c make sure no more tapes than program can handle
         if(Ncs1.gt.0 .and. Ncs1.ne.Ncp0 .and. Ncs1.ne.Nplnt0 .and.
     .    Nps1.ne.Ncp0) then
            write(Iout,1210) Nps1,Ncs1,Nplnt0,Ncp0
 1210       format(' ***ERROR*** SECOND OBJECT',I3,' HAS CENTRAL BODY',
     .          I3,' WHICH IS NEITHER',I3,' NOR',I3)
               call SUICID(
     .           'SECOND OBJECT NOT CONCENTRIC WITH OBSERVED BODY ',12)
         endif
      endif
c
c*  start=1700
c
c determine radar bias number, partial derivative controls
      do i = 1, Numrbs
        if(Nplrbs(i).eq.Nplnt0 .and. Rdbsit(1,i).eq.sitf4(1,1) .and.
     .        Rdbsit(2,i).eq.sitf4(1,2) .and. Rdbser(i).eq.Series) then
           Nrbias = i
           do j = 1, 2
              Rbsx(j) = Rbias(j, i)
           end do
           call LVTBDY(Lrbsx, Lrbs(1,i), Mrbsx, -2)
           goto 1300
        endif
      end do
      call LVTBDY(Lrbsx, izero2, Mrbsx, -2)
c
c*  start=1800
c determine optical equinox-equator number, partial
c derivative controls
 1300 do i = 1, Numeqn
         if(Eqnsit(i).eq.sitf4(1,1) .and. Eqnser(i).eq.Series) then
            Neqnox = i
            do j = 1, 3
               eqnx(j) = deqnx(i, j)
            end do
            call LVTBDY(Leqnx, Leqn(1,i), Meqnx, -3)
            goto 1400
            endif
         end do
      call LVTBDY(Leqnx, izero2, Meqnx, -3)
c
c determine optical phase number, partial derivative
c controls
 1400 if(Iabs1.gt.0) Ncph = Ncphu
      do i = 1, Numphs
         if(Nplphs(i).eq.Nplnt0 .and. Phsit(i).eq.sitf4(1,1) .and.
     .            Phser(i).eq.Series) then
            Nphase = i
            if(Iabs1.le.0 .or. Ncphu.le.Ncphs(i)) Ncph = Ncphs(i)
            do j = 1, 9
               Aphs(j) = Aphase(j, i)
            end do
            call LVTBDY(Lphsx, Lphs(1,i), Mphsx, -9)
            goto 1500
         endif
      end do
      call LVTBDY(Lphsx, izero2, Mphsx, -9)
c
c determine sky error map for optical observations
 1500 Nskyc = 0
      if(Iabs1.gt.0) Nskyc = Mskyc
      call ZFILL(Sky, 16*maxnsk)
      do i = 1, Numstr
         if(Ctlg.eq.Ctlgnm(i)) then
            Nstar = i
            if(Nskycf(i).gt.Nskyc) Nskyc = Nskycf(i)
            do j = 1, maxnsk
               Sky(j) = Skycf(j, Nstar)
            end do
            call LVTBDY(Lsky, Lskycf(1,Nstar), Msky, -maxnsk)
            goto 1600
         endif
      end do
      call LVTBDY(Lsky, izero2, Msky, -maxnsk)
c*  start=2000
c
c error analysis for observation series of a certain type of
c a given body
 1600 if(Ncodg.le.0) then
         Ncodg = Ncodf
      else if(Ncodf.ne.Ncodg .or. Klap.ne.Klaps) then
         call ERRTOT(Ncodg)
         Ncodg = Ncodf
      endif
c
c           rewind tapes as requested
c
c           fixed tapes (nbody, moon, etc.) are assumed always needed
      if(Nrewnd.gt.0) then
 
c rewind earth rotation tape if necessary
         if(Inut.gt.0 .and. Itrwnd(Inut).gt.0) then
            rewind Inut
            Itrwnd(Inut) = 0
         endif
c
c rewind n-body tape if necessary
         if(Nbody.gt.0 .and. Libdy.ne.0 .and. Nrewnd.gt.0) then
            rewind Ibody
            Itrwnd(Ibody) = 0
            if(Kpert.gt.0) then
               rewind Kpert
               Itrwnd(Kpert) = 0
            endif
         endif
c
c rewind s-body tape if necessary
         if(i2bod.gt.0 .and. Lib2y.ne.0 .and. Nrewnd.gt.0) then
            rewind i2bod
            Itrwnd(i2bod) = 0
         endif
c
c rewind moon tape if necessary
         if(Imn.gt.0 .and. Itrwnd(Imn).gt.0) then
            rewind Imn
            Itrwnd(Imn) = 0
         endif
c
c rewind moon rotation tape if necessary
         if(Ilib.gt.0 .and. Itrwnd(Ilib).gt.0) then
            rewind Ilib
            Itrwnd(Ilib) = 0
         endif
c
c rewind earth-moon barycenter tape if necessary
         if(Iem.gt.0 .and. Itrwnd(Iem).gt.0) then
            rewind Iem
            Itrwnd(Iem) = 0
         endif
c
c rewind ut1 and wobble data sets if necessary
         if(Jct(33).gt.0) then
            jut1 = Jct(33)
            if(Itrwnd(jut1).gt.0) then
               rewind jut1
               Itrwnd(jut1) = 0
            endif
            jwob = Jct(33) + 1
            if(Itrwnd(jwob).gt.0) then
               rewind jwob
               Itrwnd(jwob) = 0
            endif
         endif
      endif
c
c           the remaining tapes may not be needed
c
c
c           rewind planet tape if necessary
      if(Klam.gt.0 .and. Klam.le.u_mxpl) then
         if(Klan.ne.Klam .or. Nrewnd.gt.0) Npshp1 = 0
         if(Jplnt.gt.0 .and. Itrwnd(Jplnt).gt.0) then
            if(Klan.ne.Klam .or. Nrewnd.gt.0) then
               rewind Jplnt
               Itrwnd(Jplnt) = 0
            endif
         endif
c
c rewind planet rotation tape if necessary
         if(Klamr.gt.0) then
            if(Jpr.gt.0 .and. Itrwnd(Jpr).gt.0) then
               if(Klanr.ne.Klamr .or. Nrewnd.gt.0) then
                  rewind Jpr
                  Itrwnd(Jpr) = 0
               endif
            endif
         endif
      endif
c
c rewind satellite-probe tape if necessary
      if(Klamb.gt.0) then
         if(Jsb.gt.0 .and. Itrwnd(Jsb).gt.0) then
            if(Klanb.ne.Klamb .or. Nrewnd.gt.0) then
               rewind Jsb
               Itrwnd(Jsb) = 0
            endif
         endif
      endif
c
c rewind second satellite or observing site on satellite tape
      if(Klams1.gt.0) then
         if(Jsc.gt.0 .and. Itrwnd(Jsc).gt.0) then
            if(Klans1.ne.Klams1 .or. Nrewnd.gt.0) then
               rewind Jsc
               Itrwnd(Jsc) = 0
            endif
         endif
      endif
c
c rewind planet tapes used for solar-system barycenter offset
c (used only for pulsar observations, so there's no overlap with klam)
      if(Klssm.gt.0) then
         do i=1,9
            if(Ssbkl(i).gt.0) then
               j=Iplss(i)
               if(j.gt.0 .and. Itrwnd(j).gt.0 .and.
     .          (Klssb.ne.Klssm .or. Nrewnd.gt.0)) then
                  rewind j
                  Itrwnd(j)=0
               endif
            endif
         end do
      endif            
c
c
      return
      end
