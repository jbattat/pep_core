      subroutine PRDEQS(nfirst)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ii, imdt11, imdt22, imdt33, imdt44, imdt55,
     .          imdt66, ip, j, klaps, klntrg, ll, nmpart, nshp
 
c*** end of declarations inserted by spag
 
 
c     ash/friedman   june 1969    subroutine prdeqs
c     get partial derivative pointers for prediction of new
c           observed minus theory or formation of normal eqns.

c arguments
      integer*4 nfirst

c array dimensions
      include 'globdefs.inc'

c commons
      include 'anctrl.inc'
      include 'bernum.inc'
c     nsite1 = 1,...  receiving observing site is site(1-2,nsite1)
c                     (if 0, observing site noninput)
c     nsite2 = 1,...  sending observing site is site(1-2,nsite2)
c                     (if 0, no sending site, or noninput and diff.)
c                     (if -1, sending site noninput and same)
c     nplnt0 = planet number of observed body
c           0= sun
c           1= mercury
c           2= venus
c           4= mars
c           5= jupiter
c           6= saturn
c           7= uranus
c           8= neptune
c           9= pluto
c           10=moon
c           11,...,30  other natural planets, asteroids or satellites
c           31,....    artificial space probes integrated in planet link
c     klan   = 1...u_mxpl observed or central body is nplnt(klan)
c     klan   = u_mxpl+1   observed or central body is moon
c     klan   = u_mxpl+2   observed or central body is sun
c     klan   = 0          none of these
c     klanb  = 1...u_mxpl observed body is nplnt(klanb) read in
c                         subroutine sbreed.  if npcent(klanb) not 0,
c                         then nplnt(klan)= npcent(klanb), otherwise klan=0
c                          if klanb is not 0.
c     nspot  = 1,...u_mxspt observed spot is spot(nspot)  (if 0, no spot)
c     nrbias = 1,...u_mxrbs radar biases are rbias(1-2,nrbias)
c     nphase = 1,...u_mxphs optical phase correction coefficients are
c                           aphase(1-9,nphase)
c     neqnox = 1,...u_mxeqn optical equinox-equator, etc. corrections are
c                           eqnx(1-3,neqnox)
c
      include 'dtparm.inc'
      include 'eqenox.inc'
      include 'ethhar.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'ktrap.inc'
      include 'lcntrl.inc'
      include 'monhar.inc'
      include 'phase.inc'
      include 'plnhar.inc'
      include 'psrstf.inc'
      include 'rdbias.inc'
      include 'scoef4.inc'
      include 'skystf.inc'
      include 'sptcrd.inc'
      include 'stcord.inc'
      include 'wrkcompm.inc'
c
c local
      integer*2 ione1/1/, izr2/0/, n9
c
c calculate iptr vector for first point in obs.series
c recalculate if dt's change
      if(Imdt1.eq.imdt11 .and. Imdt2.eq.imdt22 .and.
     .   Imdt3.eq.imdt33 .and. Imdt4.eq.imdt44 .and.
     .   Imdt5.eq.imdt55 .and. Imdt6.eq.imdt66 .and.
     .   nfirst.gt.0) goto 200
c
c clear ptr vector
      Lbj = 2
      call ZFILL(Ipts, 2*Nparam)
c
c partial derivatives w.r.t. site coordinates
      if(Nsite1.gt.0) then
         Jfk = Lscrd1(Nsite1) - 1
         call PRDBDY(Lscrd(1,Nsite1), Mscrd(1,1), -6)
      endif
      if(Nsite1.ne.Nsite2 .and. Nsite2.gt.0) then
         Jfk = Lscrd1(Nsite2) - 1
         call PRDBDY(Lscrd(1,Nsite2), Mscrd(1,2), -6)
      endif
c
c increment partial counter for sites which are not input
      do j=1,2
         if(nsite(j).eq.0) then
            do i=1,6
               if(Mscrd(i,j).gt.0) Lbj = Lbj + 1
            end do
         endif
      end do
c
c partial derivatives w.r.t. equator-equinox corrections for
c optical observation series
      n9 = Numeqn
      call PRDVCT(Neqnox, ione1, n9, Leqn(1,Neqnox), Meqn, Leqn1, Leqn2,
     .            -3)
c
c partial derivatives w.r.t. sky corrections (star catalogs)
      n9 = Numstr
      call PRDVCT(Nstar, ione1, n9, Lskycf(1,Nstar), Msky, Lsky1, Lsky2,
     .            -80)
c
c partial derivatives w.r.t. solar system parameters
      Jfk = Lprm0
      call PRDPRM(Lprm, Mprm, 100)
c
c partial derivatives w.r.t. earth-moon barycenter  initial
c conditions and earth parameters
      call PRDBDY(Lem, Mem, 30)
c
c partial derivatives w.r.t. earth rotation initial conditions
c and parameters
      call PRDBDY(Ler, Mer, 30)
c
c partial derivative w.r.t. et-ut2,a1-ut1, or wobble parameter
      do i = 1, Numdt
         if(Ldt(i).gt.0) then
            Jfk = Jfk + 1
            if(i.eq.Imdt1 .or. i.eq.Imdt2) then
               Lbj = Lbj + 1
               Ipts(Jfk) = Lbj
            endif
         else if(i.eq.Imdt1 .or. i.eq.Imdt2) then
            Lbj = Lbj + 1
         endif
      end do
      if(Jddt0.le.0) then
         do i = 1, Numdt
            ii = i + 200
            if(Ldt(ii).gt.0) then
               Jfk = Jfk + 1
               if(ii.eq.Imdt3 .or. ii.eq.Imdt4) then
                  Lbj = Lbj + 1
                  Ipts(Jfk) = Lbj
               endif
            else if(ii.eq.Imdt3 .or. ii.eq.Imdt4) then
               Lbj = Lbj + 1
            endif
         end do
         do i = 1, Numdt
            ii = i + 400
            if(Ldt(ii).gt.0) then
               Jfk = Jfk + 1
               if(ii.eq.Imdt5 .or. ii.eq.Imdt6) then
                  Lbj = Lbj + 1
                  Ipts(Jfk) = Lbj
               endif
            else if(ii.eq.Imdt5 .or. ii.eq.Imdt6) then
               Lbj = Lbj + 1
            endif
         end do
      endif
      if(Numdt.le.0 .and. Mumdt.gt.0) then
         if(Imdt1.gt.0) Lbj = Lbj + 1
         if(Imdt2.gt.0) Lbj = Lbj + 1
         if(Imdt3.gt.0) Lbj = Lbj + 1
         if(Imdt4.gt.0) Lbj = Lbj + 1
         if(Imdt5.gt.0) Lbj = Lbj + 1
         if(Imdt6.gt.0) Lbj = Lbj + 1
      endif
c
c partial derivatives w.r.t. earth gravitational potential
c harmonic coefficients with central body the earth
      if(Klanb.gt.0 .and. Ncp0.eq.3) then
 
c zonal harmonics
         call PRDHAR(Lezhar, 1, ione1, Mczhar, Mczone)
 
c tesseral cosine harmonics
         Jfk = Lechr0
         call PRDHAR(Lechar, 1, ione1, Mcchar, Mctess)
 
c tesseral sine harmonics
         Jfk = Leshr0
         call PRDHAR(Leshar, 1, ione1, Mcshar, Mctess)
      endif
      Jfk = Lehar2
c
c partial derivatives w.r.t. moon initial conditions and
c parameters
      call PRDBDY(Lmn, Mmn, 30)
c
c partial derivatives w.r.t. moon rotation initial conditions
c and parameters
      call PRDBDY(Lmr, Mmr, 30)
c
c partial derivatives w.r.t. moon gravitational potential
c harmonic coefficients
      if(Nplnt0.eq.10 .or. Nplnt2.eq.10) then
c moon is observed body
c zonal harmonics
         call PRDHAR(Lmzhar, 1, ione1, Mszhar, Mszone)
 
c tesseral cosine harmonics
         Jfk = Lmchr0
         call PRDHAR(Lmchar, 1, ione1, Mschar, Mstess)
 
c tesseral sine harmonics
         Jfk = Lmshr0
         call PRDHAR(Lmshar, 1, ione1, Msshar, Mstess)
      else if(Klanb.gt.0 .and. Ncp0.eq.10) then
c moon is central body for probe
c zonal harmonics
         call PRDHAR(Lmzhar, 1, ione1, Mczhar, Mczone)
 
c tesseral cosine harmonics
         Jfk = Lmchr0
         call PRDHAR(Lmchar, 1, ione1, Mcchar, Mctess)
 
c tesseral sine harmonics
         Jfk = Lmshr0
         call PRDHAR(Lmshar, 1, ione1, Mcshar, Mctess)
      endif
      Jfk = Lmhar2
c
c partial derivatives w.r.t. initial conditions and
c parameters for observing body not the earth
      if(Klans1.gt.0) then
         Jfk = Lpl1(Klans1) - 1
         call PRDBDY(Lpl(1,Klans1), Msc, 30)
      endif
c
c partial derivatives w.r.t. planet initial conditions
c and parameters
      if(Klan.gt.0 .and. Klan.le.u_mxpl) then
         Jfk = Lpl1(Klan) - 1
         call PRDBDY(Lpl(1,Klan), Mpl, 30)
c
c partial derivatives w.r.t. planet rotation initial
c conditions and parameters to be inserted here
         if(Klanr.gt.0) then
c jfk=lpl1(klanr)-1
c call prdbdy(lpl(1,klanr),mpr,30)
c
c partial derivatives wrt planet shape harmonic coefficients
            Jfk  = Lphar1(Klanr) - 1
            nshp = Nshape(Klnhrx(Klanr))
            if(nshp.eq.1) then
 
c fourier model - nshp=1
               call PRDHAR(Lpzhar, 4, Klnhrx(Klanr), Mszhar, Mnfour)
            else if(nshp.eq.2) then
 
c grid model - nshp=2
               call PRDGRD(Lpzhar, 4, Klnhrx(Klanr), Mszhar, Mngd)
            else if(nshp.eq.3) then
c external model - nshp=3 do nothing
            else
c spherical harmonics - nshp=0
c zonal harmonics
               call PRDHAR(Lpzhar, 4, Klnhrx(Klanr), Mszhar, Mszone)
 
c tesseral harmonics
               Jfk = Lpchr0(Klanr)
               call PRDHAR(Lpchar, 4, Klnhrx(Klanr), Mschar, Mstess)
               Jfk = Lpshr0(Klanr)
               call PRDHAR(Lpshar, 4, Klnhrx(Klanr), Msshar, Mstess)
            endif
c
c
c shape harmonics not input, might be on tape
         else if(Mngd.gt.0) then
 
            call PRDGRD(izr2, 1, izr2, Mszhar, Mngd)
         else
            call PRDHAR(izr2, 1, izr2, Mszhar, Mszone)
            call PRDHAR(izr2, 1, izr2, Mschar, Mstess)
            call PRDHAR(izr2, 1, izr2, Msshar, Mstess)
         endif
c
c partial derivatives w.r.t. planet gravitational potential
c harmonic coefficients with central body the planet
         if(Klanb.le.0) goto 100
         Jfk = Lphar1(Klan) - 1
 
c zonal harmonics
         call PRDHAR(Lpzhar, 4, Klnhrx(Klan), Mczhar, Mczone)
 
c tesseral cosine harmonics
         Jfk = Lpchr0(Klan)
         call PRDHAR(Lpchar, 4, Klnhrx(Klan), Mcchar, Mctess)
 
c tesseral sine harmonics
         Jfk = Lpshr0(Klan)
         call PRDHAR(Lpshar, 4, Klnhrx(Klan), Mcshar, Mctess)
      endif
c
c partial derivatives w.r.t. satellite-probe initial
c conditions and parameters
      if(Klanb.gt.0) then
         Jfk = Lpl1(Klanb) - 1
         call PRDBDY(Lpl(1,Klanb), Msb, 30)
      endif
 
  100 klaps = Klap
 
c maybe nplnt0=-4 should already mean klap=19???
      if(Nplnt0.lt.0) Klap = 19
      if(Klap.gt.0) then
c
c spot coordinates
         call PRDVCT(Nspot, Nsptp1(Klap), Nsptp2(Klap), Lspcrd(1,Nspot),
     .               Mspcrd, Lspt1, Lspt2, -3)
c
c second spot coordinates
         if(Nspot2.gt.0 .and. Nplnt2.eq.Nplnt0) then
            Jfk = Lspt1(Nspot2) - 1
            call PRDBDY(Lspcrd(1,Nspot2), Mspcrd(1,2), -3)
         endif
c
c radar biases
         call PRDVCT(Nrbias, Nrbsp1(Klap), Nrbsp2(Klap), Lrbs(1,Nrbias),
     .               Mrbs, Lrbs1, Lrbs2, -2)
c
c optical phase corrections
         call PRDVCT(Nphase, Nphsp1(Klap), Nphsp2(Klap), Lphs(1,Nphase),
     .               Mphs, Lphs1, Lphs2, -9)
 
c maybe should leave klap???
         Klap = klaps
      endif
c
c partial derivatives w.r.t. target body initial conditions,
c parameters, and gravitional potential harmonic coefficients
      do ll = 1, Mumtar
         if(Mtrg(ll).gt.0) then
c
c search for target planet
            klntrg = Klant(ll)
c
c partial derivatives w.r.t. target planet init.cond.
            Jfk = Lpl1(klntrg) - 1
            call PRDBDY(Lpl(1,klntrg), Mtbod(1,ll), 30)
c
c search for target planet harmonics
            if(Mtzone(ll).gt.0 .or. Mttess(ll).gt.0) then
               if(Klnhrx(klntrg).le.0) call SUICID(
     .   'NO INPUT TARGET PLANET HARMONICS, STOP IN PRDEQS', 12)
c
c partial derivatives w.r.t. target planet gravitational
c potential harmonic coefficients
c zonal harmonics
               Jfk = Lphar1(klntrg) - 1
               call PRDHAR(Lpzhar, 4, Klnhrx(klntrg), Mtzhar(1,ll),
     .                     Mtzone(ll))

c tesseral cosine harmonics
               Jfk = Lpchr0(klntrg)
               call PRDHAR(Lpchar, 4, Klnhrx(klntrg), Mtchar(1,ll),
     .                     Mttess(ll))

c tesseral sine harmonics
               Jfk = Lpshr0(klntrg)
               call PRDHAR(Lpshar, 4, Klnhrx(klntrg), Mtshar(1,ll),
     .                     Mttess(ll))
            endif
         endif

      end do
c
c partial derivatives w.r.t. exo-planet parameters
      do ll=1,Mmpex
         if(Mplex(ll).gt.0) then
            klntrg = Klanex(ll)
            Jfk = Lpl1(klntrg) - 1
            call PRDBDY(Lpl(1,klntrg), Mpex(1,ll), 30)
         endif
      end do
c
c partial derivatives w.r.t. pulsar parameters
      if(Nplsr.gt.0) then
         Jfk = Lpsr0(Nplsr)
         call PRDPRM(Lpsrcn(1,Nplsr), Mpsrx, 16)
      endif
c
c test for error
      if(Lbj.ne.Numpar) then
         write(Iout, 150) Lbj, Numpar
  150    format('0PRDEQS    LBJ=', i4, ' IS NOT EQUAL TO NUMPAR=', i4)
         call SUICID(
     .' NUMBER OF PARTIAL DERIVATIVES NOT CORRECT, STOP IN PRDEQS  ',15)
      endif
c
c           setup iptr vector from ipts vector
c ipts(i)= derivative number on observation library tape for row i of
c          normal equations, i=1,nparam
c iptr(i)= row number in normal equations of derivative i on observation
c          library tape, i=1,numpar
c
      call ZFILL(Iptr, 2*Numpar)
      do i = 1, Nparam
         if(Ipts(i).gt.0) then
            ip = Ipts(i)
            Iptr(ip) = i
         endif
      end do
c
c save dt ptrs to avoid repetition
      nfirst = 1
      nmpart = Numpar
      imdt11 = Imdt1
      imdt22 = Imdt2
      imdt33 = Imdt3
      imdt44 = Imdt4
      imdt55 = Imdt5
      imdt66 = Imdt6
c
c see if number of partials is same
  200 if(Numpar.ne.nmpart) call SUICID(
     .'NUMBER OF PARTIAL DERIVATIVES NOT SAME THROUGHOUT OBS.SERIES, STO
     .P IN PRDEQS', 19)
c
c if topo grid model is used, need 4 new ptrs for each obs.
      if(Mngd.gt.0) call GRDSID
      return
      end
