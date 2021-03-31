      subroutine ANSET
 
      implicit none
c
c m.e.ash  sept 1969   subroutine anset
c check of consistency of input data and calculation of nparam
c calculate integers in /anctrl/ controlling formation of normal eqs
c
      character*4 blank/'    '/
      character*4 words(5)/' ERR','ORS ','DETE','CTED',' IN '/
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'anctrl.inc'
      integer*4 zanctr/4530/   !r8=4530,r10=4530
      include 'dtparm.inc'
      include 'eqenox.inc'
      include 'ethhar.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'monhar.inc'
      include 'namtim.inc'
      include 'phase.inc'
      include 'plnhar.inc'
      integer*4 ngdpts(4)
      equivalence (ngdpts,Scontl(1,9))
      include 'psrstf.inc'
      include 'rdbias.inc'
      include 'scoef4.inc'
      include 'skystf.inc'
      include 'sptcrd.inc'
      include 'stcord.inc'

c local
      integer   i,i1,j,jm1,n1,n11,ngd,nice,nstop
c
c zero out anctrl labelled common
      call ZFILL(Lscrd1,zanctr)
c
c calculation of nparam and check of consistency of least-
c squares adjustment control constants
      nstop  = 0
      Nparam = 0
c
c for observing sites with adjustable coordinates
      do j = 1,u_mxsit
         Lscrd1(j) = Nparam + 1
         call BCHECK(Lscrd(1,j),nice,-6)
         Lscrd2(j) = Nparam
      end do
c
c for equator-equinox corrections
      do j = 1,u_mxeqx
         Leqn1(j) = Nparam + 1
         call BCHECK(Leqn(1,j),nice,-3)
         Leqn2(j) = Nparam
      end do
c
c for sky corrections (star catalog errors)
      do j = 1,u_mxsky
         Lsky1(j) = Nparam + 1
         call BCHECK(Lskycf(1,j),nice,-80)
         Lsky2(j) = Nparam
      end do
c
c for solar system parameters
      Lprm0 = Nparam
      call ACHECK(Lprm,nice,u_nmprm)
      if(nice.gt.0) then
         nstop = nstop + nice
         write(Iout,50) nice,words,u_nmprm
   50    format(i4,5A4,'LPRM(1-',i3')')
      endif
c
c for earth-moon barycenter initial conditions and earth
c parameters
      call BCHECK(Lem,nice,u_nmbod)
      if(nice.gt.0) then
         nstop = nstop + nice
         write(Iout,100) nice,words,u_nmbod
  100    format(i4,5A4,'LEM(7-',i2,')')
      endif
c
c for earth rotation initial conditions and parameters
      call BCHECK(Ler,nice,u_nmbod)
      if(nice.gt.0) then
         nstop = nstop + nice
         write(Iout,150) nice,words,u_nmbod
  150    format(i4,5A4,'LER(7-',i2,')')
      endif
c
c for et-ut2, a1-ut1, or wobble parameters
      if(Numdt.gt.0) then
         n11 = -Numdt
         call BCHECK(Ldt,nice,n11)
         if(Jddt0.le.0) then
            call BCHECK(Ldt(201),nice,n11)
            call BCHECK(Ldt(401),nice,n11)
         endif
      endif
c
c for earth gravitational potential harmonic coefficients
      Lehar1 = Nparam + 1
      Lechr0 = Nparam
      Leshr0 = Nparam
      if(Nezone.gt.1) then
         n11 = -(Nezone - 1)
         call BCHECK(Lezhar,nice,n11)
      endif
      Lechr0 = Nparam
      Leshr0 = Nparam
      if(Netess.gt.1) then
         n11 = -((Netess*(Netess+1))/2 - 1)
         call BCHECK(Lechar,nice,n11)
         Leshr0 = Nparam
         call BCHECK(Leshar,nice,n11)
      endif
      Lehar2 = Nparam
c
c for moon initial conditions and parameters
      call BCHECK(Lmn,nice,u_nmbod)
      if(nice.gt.0) then
         nstop = nstop + nice
         write(Iout,200) nice,words,u_nmbod
  200    format(i4,5A4,'LMN(7-',i2,')')
      endif
c
c for moon rotation initial conditions and parameters
      call BCHECK(Lmr,nice,u_nmbod)
      if(nice.gt.0) then
         nstop = nstop + nice
         write(Iout,250) nice,words,u_nmbod
  250    format(i4,5A4,'LMR(7-',i2,')')
      endif
c
c for moon gravitational potential harmonic coefficients
      Lmhar1 = Nparam + 1
      Lmchr0 = Nparam
      Lmshr0 = Nparam
      if(Nmzone.gt.1) then
         n11 = -(Nmzone - 1)
         call BCHECK(Lmzhar,nice,n11)
      endif
      Lmchr0 = Nparam
      Lmshr0 = Nparam
      if(Nmtess.gt.1) then
         n11 = -((Nmtess*(Nmtess+1))/2 - 1)
         call BCHECK(Lmchar,nice,n11)
         Lmshr0 = Nparam
         call BCHECK(Lmshar,nice,n11)
      endif
      Lmhar2 = Nparam
c
c for moon spot coordinates
      Nsptm1 = 1
      Nsptm2 = 0
      do j = 1,u_mxspt
         if(Nsplnt(j).ne.10) goto 300
         Nsptm2   = Nsptm2 + 1
         Lspt1(j) = Nparam + 1
         call BCHECK(Lspcrd(1,j),nice,-3)
         Lspt2(j) = Nparam
      end do
c
c for moon radar observation biases
  300 Nrbsm1 = 1
      Nrbsm2 = 0
      do j = 1,u_mxrbs
         if(Nplrbs(j).ne.10) goto 400
         Nrbsm2   = Nrbsm2 + 1
         Lrbs1(j) = Nparam + 1
         call BCHECK(Lrbs(1,j),nice,-2)
         Lrbs2(j) = Nparam
      end do
c
c for moon optical observation phase corrections
  400 Nphsm1 = 1
      Nphsm2 = 0
      do j = 1,u_mxphs
         if(Nplphs(j).ne.10) goto 500
         Nphsm2   = Nphsm2 + 1
         Lphs1(j) = Nparam + 1
         call BCHECK(Lphs(1,j),nice,-9)
         Lphs2(j) = Nparam
      end do
c
c for sun spot coordinates
  500 Nspts1 = Nsptm2 + 1
      Nspts2 = Nsptm2
      if(Nspts1.le.u_mxspt) then
         do j = Nspts1,u_mxspt
            if(Nsplnt(j).ne.0) goto 600
            if(Spot(j).eq.blank) goto 600
            Nspts2   = Nspts2 + 1
            Lspt1(j) = Nparam + 1
            call BCHECK(Lspcrd(1,j),nice,-3)
            Lspt2(j) = Nparam
         end do
      endif
c
c for sun radar observation biases
  600 Nrbss1 = Nrbsm2 + 1
      Nrbss2 = Nrbsm2
      if(Nrbss1.le.u_mxrbs) then
         do j = Nrbss1,u_mxrbs
            if(Nplrbs(j).ne.0) goto 700
            if(Rdbsit(1,j).eq.blank) goto 700
            Nrbss2   = Nrbss2 + 1
            Lrbs1(j) = Nparam + 1
            call BCHECK(Lrbs(1,j),nice,-2)
            Lrbs2(j) = Nparam
         end do
      endif
c
c for sun optical observation phase corrections
  700 Nphss1 = Nphsm2 + 1
      Nphss2 = Nphsm2
      if(Nphss1.le.u_mxphs) then
         do j = Nphss1,u_mxphs
            if(Nplphs(j).ne.0) goto 800
            if(Phsit(j).eq.blank) goto 800
            Nphss2   = Nphss2 + 1
            Lphs1(j) = Nparam + 1
            call BCHECK(Lphs(1,j),nice,-9)
            Lphs2(j) = Nparam
         end do
      endif
c
c for planet initial conditions and parameters
  800 jm1 = 18
      do j = 1,u_mxpl
         Lpl1(j) = Nparam + 1
         call BCHECK(Lpl(1,j),nice,u_nmbod)
         Lpl2(j) = Nparam
         if(nice.gt.0) then
            nstop = nstop + nice
            write(Iout,820) nice,words,u_nmbod,j
  820       format(i4,5A4,'LPL(7-',i2,',',i2,')')
         endif
c
c for planet gravitational potential harmonic coefficients
         Klnhrx(j) = 0
         Lphar1(j) = Nparam + 1
         Lpchr0(j) = Nparam
         Lpshr0(j) = Nparam
         do i = 1,4
            if(Nplhar(i).ne.0) then
               if(Nplhar(i).eq.Nplnt(j)) then
                  Klnhrx(j) = i
                  if((Nplnt(j).ge.0) .or. (Nshape(i).le.0)) then
 
c spherical harmonics - gravitational pot. or shape
                     if(Npzone(i).gt.1) then
                        n11 = Npzone(i) - 1
                        do i1 = 1,n11
                           if(Lpzhar(i,i1).gt.0) Nparam = Nparam +
     .                         1
                        end do
                     endif
                     Lpchr0(j) = Nparam
                     Lpshr0(j) = Nparam
                     if(Nptess(i).gt.1) then
                        n11 = (Nptess(i)*(Nptess(i)+1))/2 - 1
                        do i1 = 1,n11
                           if(Lpchar(i,i1).gt.0) Nparam = Nparam +
     .                         1
                        end do
                        Lpshr0(j) = Nparam
                        do i1 = 1,n11
                           if(Lpshar(i,i1).gt.0) Nparam = Nparam +
     .                         1
                        end do
                     endif
 
c fourier series shape  nshape=1  hardwired to 122 coef.
                  else if(Nshape(i).le.1) then
                     do i1 = 1,122
                        if(Lpzhar(i,i1).gt.0) Nparam = Nparam + 1
                     end do
 
c grid shape model  nshape=2
                  else if(Nshape(i).le.2) then
                     ngd = ngdpts(i)
                     do i1 = 1,ngd
                        if(Lpzhar(i,i1).gt.0) Nparam = Nparam + 1
                     end do
                  endif
               endif
            endif
c
c
c
         end do
         Lphar2(j) = Nparam
         Lpl2(j)   = Nparam
c
c for planet spot coordinates
         Nsptp2(j) = Nsptp2(jm1)
         Nsptp1(j) = Nsptp2(j) + 1
         if(Nsptp1(j).le.u_mxspt) then
            if(Nplnt(j).gt.0) then
               n1 = Nsptp1(j)
               do i = n1,u_mxspt
                  if(Nsplnt(i).ne.Nplnt(j)) goto 850
                  Nsptp2(j) = Nsptp2(j) + 1
                  Lspt1(i)  = Nparam + 1
                  call BCHECK(Lspcrd(1,i),nice,-3)
                  Lspt2(i) = Nparam
                  Lpl2(j)  = Nparam
               end do
            endif
         endif
c
c for planet radar observation biases
  850    Nrbsp2(j) = Nrbsp2(jm1)
         Nrbsp1(j) = Nrbsp2(j) + 1
         if(Nrbsp1(j).le.u_mxrbs) then
            if(Nplnt(j).gt.0) then
               n1 = Nrbsp1(j)
               do i = n1,u_mxrbs
                  if(Nplrbs(i).ne.Nplnt(j)) goto 900
                  Nrbsp2(j) = Nrbsp2(j) + 1
                  Lrbs1(i)  = Nparam + 1
                  call BCHECK(Lrbs(1,i),nice,-2)
                  Lrbs2(i) = Nparam
                  Lpl2(j)  = Nparam
               end do
            endif
         endif
c
c for planet optical observation phase corrections
  900    Nphsp2(j) = Nphsp2(jm1)
         Nphsp1(j) = Nphsp2(j) + 1
         if(Nphsp1(j).le.u_mxphs) then
            if(Nplnt(j).gt.0) then
               n1 = Nphsp1(j)
               do i = n1,u_mxphs
                  if(Nplphs(i).ne.Nplnt(j)) goto 950
                  Nphsp2(j) = Nphsp2(j) + 1
                  Lphs1(i)  = Nparam + 1
                  call BCHECK(Lphs(1,i),nice,-9)
                  Lphs2(i) = Nparam
                  Lpl2(j)  = Nparam
               end do
            endif
         endif
  950    jm1 = j
      end do
c
c for star coordinates
      Nstar2 = Nsptp2(u_mxpl)
      Nstar1 = Nsptp2(u_mxpl) + 1
      if(Nstar1.le.u_mxspt) then
         do i = Nstar1,u_mxspt
            Lspt1(i) = Nparam + 1
            if(Nsplnt(i).lt.0) then
               Nstar2 = i
               call BCHECK(Lspcrd(1,i),nice,-3)
            endif
            Lspt2(i) = Nparam
         end do
      endif
c
c for star quantities in radar bias common
      Nrbst2 = Nrbsp2(u_mxpl)
      Nrbst1 = Nrbsp2(u_mxpl) + 1
      if(Nrbst1.le.u_mxrbs) then
         do i = Nrbst1,u_mxrbs
            Lrbs1(i) = Nparam + 1
            if(Nplrbs(i).lt.0) then
               Nrbst2 = i
               call BCHECK(Lrbs(1,i),nice,-2)
            endif
            Lrbs2(i) = Nparam
         end do
      endif
c
c for star quantities in optical phase corr.common
      Nphst2 = Nphsp2(u_mxpl)
      Nphst1 = Nphsp2(u_mxpl) + 1
      if(Nphst1.le.u_mxphs) then
         do i = Nphst1,u_mxphs
            Lphs1(i) = Nparam + 1
            if(Nplphs(i).lt.0) then
               Nphst2 = i
               call BCHECK(Lphs(1,i),nice,-9)
            endif
            Lphs2(i) = Nparam
         end do
      endif
c
c for pulsar parameters
      do i = 1,u_mxpsr
         Lpsr0(i) = Nparam
         call ACHECK(Lpsrcn(1,i),nice,16)
         if(nice.gt.0) then
            nstop = nstop + nice
            write(Iout,960) nice,words,i
  960       format(i4,5A4,'LPSRCN(.,',i2,')')
         endif
      end do
c
c stop program if any consistency errors detected
c (they all should have been caught in subroutine check of
c the input link and the program stopped then)
      if(nstop.gt.0) call SUICID(
     .' ERRORS IN CONTROL CONSTANTS FOR FORMING NORMAL EQUATIONS, STOP I
     .N ANSET',18)
 
      return
      end
