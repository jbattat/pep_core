      subroutine FRMMVE(rst,ntop)
 
      implicit none
c
c m.e.ash   feb 1970    subroutine frmmve
c move saved vector into restored vector (ntop positive)
c calling program is subroutine nrmfrm
c
c arguments
      integer ntop
      real*10 rst(1000)
c     nrst = index for restored vector (jumps around in routine)
c     nsav = index for saved vector (increases monotonically in routine)
c     rst  = output restored right side or row of normal equations
c     sav  = input saved right side or row of normal equations
c     ntop = positive, maximum length of right side or row to be
c            restored
c
c array dimensions
      include 'globdefs.inc'

c        commons
      include 'anctrl.inc'
      include 'dtparm.inc'
      include 'ethhar.inc'
      include 'lcntrl.inc'
      include 'mcnfrm.inc'
      include 'mhrfrm.inc'
      include 'monhar.inc'
      include 'namtim.inc'
      include 'plnhar.inc'
      integer*4 ngdpts(4)
      equivalence (ngdpts,Scontl(1,9))
      include 'psrstf.inc'
      include 'psrstm.inc'
      include 'restor.inc'
      include 'scoef4.inc'
      include 'zeroes.inc'
 
c local
      integer   i,k,klan,klmhar,klnhar,m1,mphase,mrbias,mspot,
     .          n1,nshp
      integer*2 mpl10/10/
      integer*2 mpl0/0/
      integer*2 onemns/-1/
c
c*  start=100
      Nsav=0
c
c observing site coordinates
      call FRMSIT(rst,ntop)
c
c equinox-equator-declination corrections
      call FRMEQN(rst,ntop)
c
c sky corrections (star catalog error models)
      call FRMSTR(rst,ntop)
c
c solar system parameters
      Nrst=Lprm0
      call FRMPRM(rst,Lprm,Mprm,u_nmprm,ntop)
      if(Nrst.ge.ntop) goto 750
c
c earth-moon barycenter init.cond, and earth parameters
      call FRMBDY(rst,Lem,Mem,u_nmbod,ntop)
      if(Nrst.ge.ntop) goto 750
c
c earth rotation initial conditions and parameters
      call FRMBDY(rst,Ler,Mer,u_nmbod,ntop)
      if(Nrst.ge.ntop) goto 750
c
c et-ut2 or a1-ut1 and wobble parameters
      if(Numdt.eq.0 .or. Mumdt.eq.0) then
         if(Numdt.ne.Mumdt) call SUICID(
     .   'NUMDT.NE.MUMDT FOR SAVED NORMAL EQUATIONS, WARNING IN FRMMVE',
     .    -15)
      else if(Mumdt.ne.Numdt) then
         call SUICID(
     .'MUMDT.NE.NUMDT FOR ET-UT2 PARAMETERS IN SAVED NORMAL EQUATIONS, S
     .TOP IN FRMMVE  ',20)
      endif
      call FRMBDY(rst,Ldt,Mdt,-600,ntop)
      if(Nrst.ge.ntop) goto 750
c
c earth gravitational potential harmonic coefficients
c zonal
      n1=Nezone-1
      m1=Mezone-1
      call FRMHAR(rst,Lezhar,Mezhar,1,1,1,n1,m1,ntop)
      if(Nrst.ge.ntop) goto 750
 
c tesseral cosine
      n1=(Netess*(Netess+1))/2-1
      m1=(Metess*(Metess+1))/2-1
      call FRMHAR(rst,Lechar,Mechar,1,1,1,n1,m1,ntop)
      if(Nrst.ge.ntop) goto 750
 
c tesseral sine
      call FRMHAR(rst,Leshar,Meshar,1,1,1,n1,m1,ntop)
      if(Nrst.ge.ntop) goto 750
c
c moon initial conditions and parameters
      call FRMBDY(rst,Lmn,Mmn,u_nmbod,ntop)
      if(Nrst.ge.ntop) goto 750
c
c moon rotation initial conditions and parameters
      call FRMBDY(rst,Lmr,Mmr,u_nmbod,ntop)
      if(Nrst.ge.ntop) goto 750
c
c moon gravitational potential harmonic coefficients
c zonal
      n1=Nmzone-1
      m1=Mmzone-1
      call FRMHAR(rst,Lmzhar,Mmzhar,1,1,1,n1,m1,ntop)
      if(Nrst.ge.ntop) goto 750
 
c tesseral cosine
      n1=(Nmtess*(Nmtess+1))/2-1
      m1=(Mmtess*(Mmtess+1))/2-1
      call FRMHAR(rst,Lmchar,Mmchar,1,1,1,n1,m1,ntop)
      if(Nrst.ge.ntop) goto 750
 
c tesseral sine
      call FRMHAR(rst,Lmshar,Mmshar,1,1,1,n1,m1,ntop)
      if(Nrst.ge.ntop) goto 750
c
c moon spot coordinates
      mspot=0
      call FRMSPT(rst,mspot,mpl10,ntop)
      if(Nrst.ge.ntop) goto 750
c
c moon radar observation biases
      mrbias=0
      call FRMRBS(rst,mrbias,mpl10,ntop)
      if(Nrst.ge.ntop) goto 750
c
c moon optical observation phase corrections
      mphase=0
      call FRMPHS(rst,mphase,mpl10,ntop)
      if(Nrst.ge.ntop) goto 750
c
c sun spot coordinates
      call FRMSPT(rst,mspot,mpl0,ntop)
      if(Nrst.ge.ntop) goto 750
c
c sun radar observation biases
      call FRMRBS(rst,mrbias,mpl0,ntop)
      if(Nrst.ge.ntop) goto 750
c
c sun optical observations phase corrections
      call FRMPHS(rst,mphase,mpl0,ntop)
      if(Nrst.ge.ntop) goto 750
c
c*  start=500
c start of planet loop
      do k=1,Mumpln
c
c find planet for restored normal equations
         klan=0
         do i=1,Numpln
            if(Nplnt(i).eq.Mplnt(k)) then
               klan=i
               goto 10
            endif
         end do
c
c planet initial conditions and parameters
   10    if(klan.gt.0) then
            Nrst=Lpl1(klan)-1
            call FRMBDY(rst,Lpl(1,klan),Mpl(1,k),u_nmbod,ntop)
         else
            call FRMBDY(rst,izero2,Mpl(1,k),u_nmbod,ntop)
         endif
c
c see if there are gravitational potential harmonics
c for the planet
         klmhar=0
         do i=1,Mumphr
            if(Mplhar(i).eq.Mplnt(k)) then
               klmhar=i
               goto 20
            endif
         end do
   20    if(klmhar.gt.0) then
            klnhar=0
            if(klan.gt.0) then
               do i=1,Nmphar
                  if(Nplhar(i).eq.Nplnt(klan)) then
                     klnhar=i
                     goto 30
                  endif
               end do
   30          if(klnhar.gt.0 .and.
     .          Nshape(klnhar).ne.mnshap(klmhar)) call SUICID(
     .          ' INPUT SHAPE MODEL DOESNT MATCH SNE. FRMMVE ',11)
            endif
            nshp=mnshap(klmhar)+1
            if(nshp.eq.2) then
c
c fourier
               m1=122
               if(klnhar.gt.0) then
                  n1=122
                  call FRMHAR(rst,Lpzhar,Mpzhar,klnhar,klmhar,4,
     .             n1,m1,ntop)
               else
                  call FRMHAR(rst,izero2,Mpzhar,1,klmhar,4,0,m1,ntop)
               endif
            else if(nshp.eq.3) then
c
c grid
               m1=Mngd(klmhar)
               if(klnhar.gt.0) then
                  n1=ngdpts(klnhar)
                  call FRMHAR(rst,Lpzhar,
     .             Mpzhar,klnhar,klmhar,4,n1,m1,ntop)
               else
                  call FRMHAR(rst,izero2,Mpzhar,1,klmhar,4,0,m1,ntop)
               endif
            else
c
c planet gravitational potential harmonic coefficients
c zonal
               m1=Mpzone(klmhar)-1
               if(klnhar.gt.0) then
                  n1=Npzone(klnhar)-1
                  call FRMHAR(rst,Lpzhar,Mpzhar,klnhar,klmhar,4,
     .             n1,m1,ntop)
               else
                  call FRMHAR(rst,izero2,Mpzhar,1,klmhar,4,0,m1,ntop)
               endif
 
c tesseral cosine
               m1=(Mptess(klmhar)*(Mptess(klmhar)+1))/2-1
               if(klnhar.gt.0) then
                  n1=(Nptess(klnhar)*(Nptess(klnhar)+1))/2-1
                  call FRMHAR(rst,Lpchar,Mpchar,klnhar,klmhar,4,
     .             n1,m1,ntop)
               else
                  call FRMHAR(rst,izero2,Mpchar,1,klmhar,4,0,m1,ntop)
               endif
 
c tesseral sine
               if(klnhar.gt.0) then
                  call FRMHAR(rst,Lpshar,Mpshar,klnhar,klmhar,4,
     .             n1,m1,ntop)
               else
                  call FRMHAR(rst,izero2,Mpshar,1,klmhar,4,0,m1,ntop)
               endif
            endif
         endif
c
c*  start=570
c do not call spot or bias routines for planet rotation
         if(Mplnt(k).gt.0) then
c
c planet spot coordinates
            call FRMSPT(rst,mspot,Mplnt(k),ntop)
c
c planet radar observation biases
            call FRMRBS(rst,mrbias,Mplnt(k),ntop)
c
c planet optical observation phase corrections
            call FRMPHS(rst,mphase,Mplnt(k),ntop)
         endif
c
c end of planet loop
      end do
c
c*  start=600
c star coordinates
      call FRMSPT(rst,mspot,onemns,ntop)
c
c star quantities in radar bias common
      call FRMRBS(rst,mrbias,onemns,ntop)
c
c star quantities in optical phase corr.common
      call FRMPHS(rst,mphase,onemns,ntop)
c
c pulsar parameters
      do k=1,Mumpsr
c
c find corresponding input pulsar
         do i=1,Numpsr
            if(Sptpsr(i).eq.Sptps1(k)) then
               Nrst=Lpsr0(i)
               call FRMPRM(rst,Lpsrcn(1,i),Mpsrcn(1,k),u_nmpsr,ntop)
               goto 40
            endif
         end do
c
c no match
         call FRMPRM(rst,izero2,Mpsrcn(1,k),u_nmpsr,ntop)
   40 end do

  750 return
      end
