      subroutine PRDSET(nrmprd)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, j1, j2, k, ll, n, n1
 
c*** end of declarations inserted by spag
 
 
c        subroutine prdset - j.f.chandler - 1982 march
c        based on code from subroutine nrmict - m.e.ash, oct 1969
c     determine pointers for observation series and fill bernum common
c     called from nrmict or prdict
c
c        parameters
      integer*4 nrmprd
c  nrmprd indicates caller: 1-nrmict, 2-prdict
c
c array dimensions
      include 'globdefs.inc'

c        common
      include 'bernum.inc'
      include 'eqenox.inc'
      include 'inodta.inc'
      include 'ktrap.inc'
      include 'namtim.inc'
      include 'phase.inc'
      include 'psrstf.inc'
      include 'rdbias.inc'
      include 'skystf.inc'
      include 'sptcrd.inc'
      include 'stcord.inc'
      include 'wrkcompr.inc'

c external functions
      integer*4 LEG,NSCAN

c local
      character*4 blank/'    '/,amper9/'&&&&'/,pound9/'####'/
      character*4 ap15/'AP15'/,al15/'AL15'/
c
c*  start=1000
c determine klan and klanb for observed body
      Mnplnt = Nplnt0 - 10
      if(Nplnt0.ge.0) then
         if(Nplnt0.ne.3) then
            if(Nplnt0.eq.10) then
               Klan = u_mxpl+1
               Ncp0 = 3
               Klap = u_mxpl+1
            else if(Nplnt0.gt.0) then
               do i = 1, Numpln
                  if(Nplnt(i).eq.Nplnt0) then
                     Klap = i
                     Ncp0 = Npcent(i)
                     if(Nplnt0.le.30 .and. Ncp0.le.0 .and. 
     .                Inttyp(i).ne.0) then
                        Klan = i
                        goto 100
                     endif
                     Klanb = i
                     if(Ncp0.le.0) goto 100
                     if(Ncp0.ne.10) then
                        if(Ncp0.eq.3) goto 100
                        do j = 1, Numpln
                           if(Ncp0.eq.Nplnt(j)) then
                              Klan = j
                              goto 100
                           endif
                        end do
                        write(Iout,10) Ncp0,Nplnt0
   10                   format(' NO CENTRAL BODY', i3, ' FOR BODY', i3,
     .                         ', STOP IN PRDSET LABEL=1050 ')
                        goto 30
                     else
                        Klan = u_mxpl+1
                        goto 100
                     endif
                  endif
               end do
               write(Iout,20) Nplnt0
   20          format('  NO INPUT PLANET', i3,
     .                ', STOP IN PRDSET LABEL=1080 ')
   30          call SUICID('ERRORS IN PRDSET', 4)
            else
               Klan = u_mxpl+2
               Klap = u_mxpl+2
            endif
         endif
      endif
c
c           determine klans1 for a second observed body (ncodf.gt.20)
c             (see also below statements label=1400 and label=1500)
c                klans1 is largely a flag to indicate tape reading
c                for a second object - the object itself is nps1
  100 if(Ncodf.gt.20) then
         if(Nplnt2.gt.0 .and. Nplnt2.ne.3) then
            do i = 1, Numpln
               if(Nplnt(i).eq.Nplnt2) then
                  Klans1 = i
                  Nps1   = Nplnt2
                  Ncs1   = Npcent(i)
                  goto 200
               endif
            end do
            if(Nplnt2.ne.10) then
               call SUICID(
     .    'SECOND OBSERVED OBJECT NOT INPUT, STOP IN PRDSET LABEL=1110 '
     .    , 15)
            else
               Klans1 = u_mxpl+1
               Nps1   = 10
               Ncs1   = 3
            endif
         endif
      endif
c
c determine klanr for planet rotation
  200 if(Klan.gt.0) then
         if(Klan.le.u_mxpl) then
            do i = 1, Numpln
               if(Nplnt(Klan) + Nplnt(i).eq.0) then
                  Klanr = i
                  goto 300
               endif
            end do
         endif
      endif
c
c*  start=1200
c determine site number
  300 if(Jtypob.ne.3) then
         do i = 1, 2
            if(i.gt.1) then
               if(Jtypob.eq.2) goto 400
               if(Sitf(1,2).eq.blank) goto 400
               if(Sitf(1,2).eq.Sitf(1,1)) Nsite2=-1
            endif
            do k = 1, Numsit
               if(Sitf(1,i).eq.Site(1,k)) then
                  nsite(i)   = k
                  Sitf(2,i) = Site(2,k)
               endif
            end do
         end do
      endif
c
c*  start=1400
c determine if observing site is not on earth
  400 if(Ncodf.le.20) then
         if(Ksite(1).gt.0 .and. Ksite(1).ne.3) then
            Nps1 = Ksite(1)
            if(Nps1.ne.10) then
               if(Numpln.gt.0) then
                  do i = 1, Numpln
                     if(Nplnt(i).eq.Nps1) then
                        Klans1 = i
                        Ncs1   = Npcent(i)
                        goto 500
                     endif
                  end do
               endif
               call SUICID(
     .                'OBSERVING SITE ON NON-INPUT PLANET, STOP PRDSET '
     .                , 12)
            else
               Klans1 = u_mxpl+1
               Ncs1   = 3
            endif
         endif
      endif
c
c*  start=1500
c determine spot number or klans1 from spot name
  500 if(Spota.ne.blank) then
         if(Spota.ne.amper9 .and. Spota.ne.pound9) then
            if(Spota.eq.al15) Spota   = ap15
            if(Spota2.eq.al15) Spota2 = ap15
 
c determine pulsar number from spot name
            if(Numpsr.gt.0) then
               do i = 1, Numpsr
                  if(Spota.eq.Sptpsr(i)) then
                     Nplsr = i
                     goto 520
                  endif
               end do
            endif
  520       if(Jtypob.eq.3 .or. Jtypob.eq.5) then
c spot name gives 2nd body for undifferenced observables
c compare from 1st non-blank character in each
               j1 = NSCAN(Spota,4,blank)
               n1 = 4 - j1
               do i = 1, Numpln
                  j2 = NSCAN(Aplnt(i),8,blank)
                  if(j2.ge.0) then
                     n = min0(n1,8 - j2)
                     if(LEG(n,j1+1,Spota,j2+1,Aplnt(i)).eq.0) then
                        Klans1 = i
                        Nps1   = Nplnt(i)
                        Ncs1   = Npcent(i)
                        goto 700
                     endif
                  endif
               end do
 
c spot name indicates observed spot or star
            else if(Numspt.gt.0) then
               do k = 1, Numspt
                  if(Spota.eq.Spot(k)) then
                     if(Nsplnt(k).ne.Nplnt0) call SUICID(
     .' OBSERVED SPOT IS NOT ON CORRECT PLANET, STOP IN PRDSET ',14)
                     Nspot = k
                     goto 600
                  endif
               end do
            endif
            call SUICID(' OBSERVED SPOT IS NOT INPUT, STOP IN PRDSET ',
     .                  11)
         endif
      endif
  600 if(Spota2.ne.blank) then
         do k = 1, Numspt
            if(Spota2.eq.Spot(k)) then
               Nspot2 = k
               goto 700
            endif
         end do
         call SUICID(' OBSERVED SPOT IS NOT INPUT, STOP IN PRDSET ', 11)
      endif
c
c*  start=1600
c check consistency of klans1
  700 if(Klans1.gt.0) then
 
c test klans1 to see if actual 'sc' tape
         if(Klans1.eq.Klanb .or. Klans1.eq.Klan) Klans1 = 0
 
c make sure no more tapes than program can handle
         if(Ncs1.gt.0) then
            if(Ncs1.ne.Ncp0 .and. Ncs1.ne.Nplnt0) then
               if(Nps1.ne.Ncp0) call SUICID(
     .             'SECOND OBJECT NOT CONCENTRIC WITH OBSERVED BODY ',
     .             12)
            endif
         endif
      endif
c
c*  start=1700
c determine radar bias number
      if(Numrbs.gt.0) then
         do i = 1, Numrbs
            if(Nplrbs(i).eq.Nplnt0) then
               if(Rdbsit(1,i).eq.Sitf(1,1)) then
                  if(Rdbsit(2,i).eq.Sitf(1,2)) then
                     if(Rdbser(i).eq.Series) then
                        Nrbias = i
                        goto 800
                     endif
                  endif
               endif
            endif
         end do
      endif
c
c determine optical equinox-equator number
  800 if(Numeqn.gt.0) then
         do i = 1, Numeqn
            if(Eqnsit(i).eq.Sitf(1,1)) then
               if(Eqnser(i).eq.Series) Neqnox = i
            endif
         end do
      endif
c
c determine optical phase number
      if(Numphs.gt.0) then
         do i = 1, Numphs
            if(Nplphs(i).eq.Nplnt0) then
               if(Phsit(i).eq.Sitf(1,1)) then
                  if(Phser(i).eq.Series) then
                     Nphase = i
                     goto 900
                  endif
               endif
            endif
         end do
      endif
c
c*  start=1800
c determine star catalog
  900 if(Numstr.gt.0) then
         do i = 1, Numstr
            if(Ctlgnm(i).eq.Ctlgm) then
               Nstar = i
               goto 1000
            endif
         end do
      endif
c
c*  start=1900
c search for target body
 1000 if(Mumtar.gt.0) then
         do ll = 1, Mumtar
            do i = 1, Numpln
               if(Mtrg(ll).eq.Nplnt(i)) then
                  Klant(ll) = i
                  goto 1050
               endif
            end do
            call SUICID(' NO INPUT TARGET PLANET, STOP IN PRDSET ', 10)
 1050    end do
      endif
c
c exo-planet parameters
      do ll = 1, Mmpex
         do i = 1, Numpln
            if(Mplex(ll).eq.Nplnt(i)) then
               Klanex(ll) = i
               goto 1150
            endif
         end do
         Klanex(ll)=0
 1150 end do
c
c*  start=2000
      return
      end
