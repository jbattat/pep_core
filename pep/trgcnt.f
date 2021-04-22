      subroutine TRGCNT
 
      implicit none
 
c
c m.e.ash    june 1969     subroutine trgcnt
c calculate partial derivative control vectors for central body
c harmonics and target body init.cond.,parameters,harmonics
c
c array dimensions
      include 'globdefs.inc'
c common
      include 'ethhar.inc'
      include 'fcntrl.inc'
      include 'ltrapx.inc'
      include 'monhar.inc'
      include 'mtrapx.inc'
      include 'number.inc'
      include 'pemctl.inc'
      include 'plnhar.inc'
      include 'scoef4.inc'
      include 'sbdta.inc'
c
c local
      integer*4 i, j, klnhr, ll, n2
      integer*2 izr2/0/
c
c determine central body gravitational potential harmonic
c partial derivative controls
      if(Ict(1).gt.0 .and. Klanb.gt.0 .and. Ncp0.gt.0) then

c earth is central body
         if(Ncp0.eq.3) then
            if(Nezone.gt.25 .or. Netess.gt.20) call SUICID(
     .       ' NEZONE.GT.25.OR.NETESS.GT.20.  STOP IN TRGCNT. ',12)
            n2 = Nezone - 1
            call LVTHAR(Lczhar,Mczhar,Lezhar,1,1,1,n2,
     .       Nczone, Mczone, n2)
            n2 = (Netess*(Netess+1))/2 - 1
            call LVTHAR(Lcchar,Mcchar,Lechar,1,1,1,n2,Nctess,Mctess,n2)
            call LVTHAR(Lcshar,Mcshar,Leshar,1,1,1,n2,Nctess,Mctess,n2)

c moon is central body
         else if(Ncp0.eq.10) then
            if(Nmzone.gt.25 .or. Nmtess.gt.20) call SUICID(
     .       ' NMZONE.GT.25.OR.NMTESS.GT.20.  STOP IN TRGCNT. ',12)
            n2 = Nmzone - 1
            call LVTHAR(Lczhar,Mczhar,Lmzhar,1,1,1,n2,Nczone,Mczone,n2)
            n2 = (Nmtess*(Nmtess+1))/2 - 1
            call LVTHAR(Lcchar,Mcchar,Lmchar,1,1,1,n2,Nctess,Mctess,n2)
            call LVTHAR(Lcshar,Mcshar,Lmshar,1,1,1,n2,Nctess,Mctess,n2)

         else
c planet is central body
            do i = 1, Nmphar
               if(Nplhar(i).eq.Ncp0) then
                  klnhr = i
                  if(Npzone(klnhr).gt.25 .or. Nptess(klnhr)
     .             .gt.20) call SUICID(
     .             ' NPZONE.GT.25.OR.NPTESS.GT.20.  STOP IN TRGCNT. '
     .             ,12)
                  n2 = Npzone(klnhr) - 1
                  call LVTHAR(Lczhar,Mczhar,Lpzhar,4,klnhr,1,
     .             n2,Nczone,Mczone,n2)
                  n2 = (Nptess(klnhr)*(Nptess(klnhr)+1))/2 - 1
                  call LVTHAR(Lcchar,Mcchar,Lpchar,4,klnhr,1,
     .             n2,Nctess,Mctess,n2)
                  call LVTHAR(Lcshar,Mcshar,Lpshar,4,klnhr,1,
     .             n2,Nctess,Mctess,n2)
                  goto 100
               endif
            end do
            n2 = 0
            call LVTHAR(Lczhar,Mczhar,izr2,1,1,1,n2,Nczone,Mczone,n2)
            call LVTHAR(Lcchar,Mcchar,izr2,1,1,1,n2,Nctess,Mctess,n2)
            call LVTHAR(Lcshar,Mcshar,izr2,1,1,1,n2,Nctess,Mctess,n2)
         endif
      endif
c
c determine target body initial condition, parameter and
c gravitational potential harmonic partial derivative controls
  100 Numtar = 0
      do ll = 1, i_mxtrg
         Ntrg(ll)  = 0
         Klant(ll) = 0
         do j = 1, 30
            Ltbod(j,ll) = 0
         end do
         Ntzone(ll) = 0
         do j = 1, 4
            Ltzhar(j,ll) = 0
         end do
         Nttess(ll) = 0
         do j = 1, 5
            Ltchar(j,ll) = 0
            Ltshar(j,ll) = 0
         end do
      end do
c
      if(Ict(1).gt.0) then
         if(Klanb.gt.0) then
c
c target body partial control vectors for probe  observations
            call TRGLIC(Nkisb,Kisb,.false.)
         else if(Nplnt0.gt.0 .and. Nplnt0.ne.10) then
c target body partial control vectors for planet observations
            call TRGLIC(Nkipl,Kipl,.true.)
         else if(Nplnt0.eq.10) then
c target body partial control vectors for moon observations
            call TRGLIC(Nkimn,Kimn,.true.)
         else if(Nplnt0.eq.-4) then
c target body partial control vectors for star observations
            call TRGLIC(Nkiem,Kiem,.true.)
         endif
      endif
 
      return
      end
