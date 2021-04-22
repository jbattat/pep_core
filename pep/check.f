      subroutine CHECK(nstop)
 
      implicit none
c
c m.e.ash  sept 1969   subroutine check
c check of consistency of input data and calculation of nparam

c arguments
      integer*4 nstop

c array dimensions
      include 'globdefs.inc'

c commons
      include 'dtparm.inc'
      include 'empcnd.inc'
      include 'eqenox.inc'
      include 'ethhar.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'maxmt0dt.inc'
      include 'lcntrl.inc'
      include 'monhar.inc'
      include 'namtim.inc'
      include 'obsdta.inc'
      include 'phase.inc'
      include 'plnhar.inc'
      integer*4 ngdpts(4)
      equivalence (Scontl(1,9),ngdpts)
      include 'psrstf.inc'
      include 'rdbias.inc'
      include 'scoef4.inc'
      include 'skystf.inc'
      include 'sptcrd.inc'
      include 'stcord.inc'
 
c temporary storage
      common/WRKCOM/ Buff(1000)
      real*10 Buff
 
c local
      integer   i,i1,j,n11,ngd,nice,nphs01,nphs02,nrbs01,
     .          nrbs02,nspt01,nspt02,numdt1
      character*4 blank/'    '/
      character*4 words(5)/' ERR','ORS ','DETE','CTED',' IN '/
c
c calculation of nparam and check of consistency of least-
c squares adjustment control constants
      nstop  = 0
      Nparam = 0
      if(Ict(1).gt.0) then
c
c for duplicate planet numbers
         do i = 1,Numpln
            do j = 1,Numpln
               if(i.eq.j) goto 50
               if(Nplnt(i).ne.0 .and. Nplnt(j).ne.0) then
                  if(Nplnt(i).eq.Nplnt(j)) then
                     write(Iout,10) i,Nplnt(i),j,Nplnt(j)
                     if(Mout.gt.0) write(Mout,10) i,Nplnt(i),j,Nplnt(j)
   10                format(' NPLNT(',i2,')=',i3,' AND NPLNT(',i2,
     .                      ')=',i3,
     .                      ' NOT ALLOWED IN LEAST SQUARES ANALYSIS')
                     nstop = nstop + 1
                  endif
               endif
            end do
   50    end do
c
c for observing sites with adjustable coordinates
         do j = 1,u_mxsit
            call BCHECK(Lscrd(1,j),nice,-6)
         end do
c
c for equator-equinox corrections
         do j = 1,u_mxeqx
            call BCHECK(Leqn(1,j),nice,-3)
         end do
c
c for sky corrections (star catalog errors)
         do j = 1,u_mxsky
            call BCHECK(Lskycf(1,j),nice,-80)
         end do
c
c for solar system parameters
         call ACHECK(Lprm,nice,u_nmprm)
         if(nice.gt.0) then
            nstop = nstop + nice
            write(Iout,60) nice,words,'LPRM(1-',u_nmprm
            if(Mout.gt.0) write(Mout,60) nice,words,'LPRM(1-',u_nmprm
   60       format(i4,5A4,a,i3')')
         endif
c
c for earth-moon barycenter initial conditions and earth
c parameters
         call BCHECK(Lem,nice,u_nmbod)
         if(nice.gt.0) then
            nstop = nstop + nice
            write(Iout,60) nice,words,'LEM(7-',u_nmbod
            if(Mout.gt.0) write(Mout,60) nice,words,'LEM(7-',u_nmbod
         endif
c
c for earth rotation initial conditions and parameters
         call BCHECK(Ler,nice,u_nmbod)
         if(nice.gt.0) then
            nstop = nstop + nice
            write(Iout,60) nice,words,'LER(7-',u_nmbod
            if(Mout.gt.0) write(Mout,60) nice,words,'LER(7-',u_nmbod
         endif
c
c for et-ut2, a1-ut1 or wobble parameters
         if(Numdt.gt.0) then
            n11 = -Numdt
            call BCHECK(Ldt,nice,n11)
            numdt1 = Numdt - 1
            if(Jddt0.le.1) then
               if(Jddt0.le.0) then
                  call BCHECK(Ldt(201),nice,n11)
                  call BCHECK(Ldt(401),nice,n11)
               endif
               do i = 1,numdt1
                  i1 = i + 1
                  if(Jddt(i1).le.Jddt(i)) then
                     write(Iout,110) i1,Jddt(i1)
                     if(Mout.gt.0) write(Mout,110) i1,Jddt(i1)
  110                format(' JDDT(',i3,')=',i8,
     .                      ' NOT INCREASING FOR A1-UT1 OR WOBBLE')
                     nstop = nstop + 1
                  endif
               end do
            else
               do i = 1,numdt1
                  i1 = i + 1
                  if(Jddt(i1).ge.Jddt(i)) then
                     write(Iout,120) i1,Jddt(i1)
                     if(Mout.gt.0) write(Mout,120) i1,Jddt(i1)
  120                format(' JDDT(',i3,')=',i8,
     .                      ' NOT DECREASING FOR ET-UT2')
                     nstop = nstop + 1
                  endif
               end do
            endif
         endif
c
c for earth gravitational potential harmonic coefficients
         if(Nezone.gt.1) then
            n11 = -(Nezone - 1)
            call BCHECK(Lezhar,nice,n11)
         endif
         if(Netess.gt.1) then
            n11 = -((Netess*(Netess+1))/2 - 1)
            call BCHECK(Lechar,nice,n11)
            call BCHECK(Leshar,nice,n11)
         endif
c
c for moon initial conditions and parameters
         call BCHECK(Lmn,nice,u_nmbod)
         if(nice.gt.0) then
            nstop = nstop + nice
            write(Iout,60) nice,words,'LMN(7-',u_nmbod
            if(Mout.gt.0) write(Mout,60) nice,words,'LMN(7-',u_nmbod
         endif
c
c for moon rotation initial conditions and parameters
         call BCHECK(Lmr,nice,u_nmbod)
         if(nice.gt.0) then
            nstop = nstop + nice
            write(Iout,60) nice,words,'LMR(7-',u_nmbod
            if(Mout.gt.0) write(Mout,60) nice,words,'LMR(7-',u_nmbod
         endif
c
c for moon gravitational potential harmonic coefficients
         if(Nmzone.gt.1) then
            n11 = -(Nmzone - 1)
            call BCHECK(Lmzhar,nice,n11)
         endif
         if(Nmtess.gt.1) then
            n11 = -((Nmtess*(Nmtess+1))/2 - 1)
            call BCHECK(Lmchar,nice,n11)
            call BCHECK(Lmshar,nice,n11)
         endif
c
c for moon spot coordinates
         nspt02 = 0
         do j = 1,u_mxspt
            if(Nsplnt(j).ne.10) goto 200
            nspt02 = nspt02 + 1
            call BCHECK(Lspcrd(1,j),nice,-6)
         end do
c
c for moon radar observation biases
  200    nrbs02 = 0
         do j = 1,u_mxrbs
            if(Nplrbs(j).ne.10) goto 250
            nrbs02 = nrbs02 + 1
            call BCHECK(Lrbs(1,j),nice,-2)
         end do
c
c for moon optical observation phase corrections
  250    nphs02 = 0
         do j = 1,u_mxphs
            if(Nplphs(j).ne.10) goto 300
            nphs02 = nphs02 + 1
            call BCHECK(Lphs(1,j),nice,-9)
         end do
c
c for sun spot coordinates
  300    nspt01 = nspt02 + 1
         if(nspt01.le.u_mxspt) then
            do j = nspt01,u_mxspt
               if(Nsplnt(j).ne.0) goto 350
               if(Spot(j).eq.blank) goto 350
               nspt02 = nspt02 + 1
               call BCHECK(Lspcrd(1,j),nice,-6)
            end do
         endif
c
c for sun radar observation biases
  350    nrbs01 = nrbs02 + 1
         if(nrbs01.le.u_mxrbs) then
            do j = nrbs01,u_mxrbs
               if(Nplrbs(j).ne.0) goto 400
               if(Rdbsit(1,j).eq.blank) goto 400
               nrbs02 = nrbs02 + 1
               call BCHECK(Lrbs(1,j),nice,-2)
            end do
         endif
c
c for sun optical observation phase corrections
  400    nphs01 = nphs02 + 1
         if(nphs01.le.u_mxphs) then
            do j = nphs01,u_mxphs
               if(Nplphs(j).ne.0) goto 450
               if(Phsit(j).eq.blank) goto 450
               nphs02 = nphs02 + 1
               call BCHECK(Lphs(1,j),nice,-9)
            end do
         endif
c
c for planet initial conditions and parameters
  450    do j = 1,u_mxpl
            call BCHECK(Lpl(1,j),nice,u_nmbod)
            if(nice.gt.0) then
               nstop = nstop + nice
               write(Iout,460) nice,words,u_nmbod,j
               if(Mout.gt.0) write(Mout,460) nice,words,u_nmbod,j
  460          format(i4,5A4,'LPL(7-',i2,',',i2,')')
            endif
c
c for planet gravitational potential harmonic coefficients
c or planet shape coefficients
            do i = 1,4
               if(Nplhar(i).ne.0) then
                  if(Nplhar(i).eq.Nplnt(j)) then
                     if((Nplnt(j).ge.0) .or. (Nshape(i).le.0))
     .                  then
 
c spherical harmonic coefs-grav. potential or shape
                        if(Npzone(i).gt.1) then
                           n11 = Npzone(i) - 1
                           do i1 = 1,n11
                              if(Lpzhar(i,i1).gt.0)
     .                            Nparam = Nparam + 1
                           end do
                        endif
                        if(Nptess(i).gt.1) then
                           n11 = (Nptess(i)*(Nptess(i)+1))/2 - 1
                           do i1 = 1,n11
                              if(Lpchar(i,i1).gt.0)
     .                            Nparam = Nparam + 1
                           end do
                           do i1 = 1,n11
                              if(Lpshar(i,i1).gt.0)
     .                            Nparam = Nparam + 1
                           end do
                        endif
 
c non-spherical harmonic shape models
                     else if(Nshape(i).le.1) then
c fourier series shape model-nshape=1
c current model has 122 coefficients
                        do i1 = 1,122
                           if(Lpzhar(i,i1).gt.0) Nparam = Nparam +
     .                         1
                        end do
                     else if(Nshape(i).le.2) then
c altitude grid local shape model-nshape=2
c ngdpts=number of grid points
                        ngd = ngdpts(i)
                        do i1 = 1,ngd
                           if(Lpzhar(i,i1).gt.0) Nparam = Nparam +
     .                         1
                        end do
                     endif
                  endif
               endif
 
c no shape model yet for nshape gt2
            end do
c
c for planet spot coordinates
            nspt01 = nspt02 + 1
            if(nspt01.le.u_mxspt) then
               if(Nplnt(j).gt.0) then
                  do i = nspt01,u_mxspt
                     if(Nsplnt(i).ne.Nplnt(j)) goto 480
                     nspt02 = nspt02 + 1
                     call BCHECK(Lspcrd(1,i),nice,-6)
                  end do
               endif
            endif
c
c for planet radar observation biases
  480       nrbs01 = nrbs02 + 1
            if(nrbs01.le.u_mxrbs) then
               if(Nplnt(j).gt.0) then
                  do i = nrbs01,u_mxrbs
                     if(Nplrbs(i).ne.Nplnt(j)) goto 500
                     nrbs02 = nrbs02 + 1
                     call BCHECK(Lrbs(1,i),nice,-2)
                  end do
               endif
            endif
c
c for planet optical observation phase corrections
  500       nphs01 = nphs02 + 1
            if(nphs01.le.u_mxphs) then
               if(Nplnt(j).gt.0) then
                  do i = nphs01,u_mxphs
                     if(Nplphs(i).ne.Nplnt(j)) goto 550
                     nphs02 = nphs02 + 1
                     call BCHECK(Lphs(1,i),nice,-9)
                  end do
               endif
            endif
  550    end do
c
c for star coordinates
         nspt01 = nspt02 + 1
         if(nspt01.le.u_mxspt) then
            do i = nspt01,u_mxspt
               if(Nsplnt(i).lt.0) then
                  nspt02 = i
                  call BCHECK(Lspcrd(1,i),nice,-6)
               endif
            end do
         endif
c
c for star quantities in radar bias common
         nrbs01 = nrbs02 + 1
         if(nrbs01.le.u_mxrbs) then
            do i = nrbs01,u_mxrbs
               if(Nplrbs(i).lt.0) then
                  nrbs02 = i
                  call BCHECK(Lrbs(1,i),nice,-2)
               endif
            end do
         endif
c
c for star quantities in optical phase corr.common
         nphs01 = nphs02 + 1
         if(nphs01.le.u_mxphs) then
            do i = nphs01,u_mxphs
               if(Nplphs(i).lt.0) then
                  nphs02 = i
                  call BCHECK(Lphs(1,i),nice,-9)
               endif
            end do
         endif
c
c for pulsar parameters
         do i = 1,u_mxpsr
            call ACHECK(Lpsrcn(1,i),nice,u_nmpsr)
            if(nice.gt.0) then
               nstop = nstop + nice
               write(Iout,60) nice,words,'LPSRCN(.,',i
               if(Mout.gt.0) write(Mout,60) nice,words,'LPSRCN(.,',i
            endif
         end do
      endif
c
c check of consistency of planet control constants
      call PCHECK(nstop)
c
c check if maximum number of observation library o-c tapes
c or saved normal equations tapes is exceeded
      if((Numobt.gt.10) .or. (Nummt0.gt.maxmt0)) then
         write(Iout,600) Numobt,Nummt0
         if(Mout.gt.0) write(Mout,600) Numobt,Nummt0
  600    format(' NUMOBT=',i3,'  NUMMT0=',i3,
     .          '  MAXIMUM ALLOWED IS 10, ERROR')
         nstop = nstop + 1
      endif
c
c determine precession model if not specified
      if(Jct(21).eq.0) then
         if(Ercond(28).eq.0._10 .or. Ercond(28).eq.5026.75_10) then
c old precession model
            Jct(21)=-1
         endif
      endif
      if(Jct(21).lt.0 .and. Jct(13).gt.0) then
         write(Iout,610)
         if(Mout.gt.0) write(Mout,610)
  610    format(' PRE-1976 PRECESSION MODEL NOT IMPLEMENTED FOR J2000',
     .    ' EPOCH')
         nstop=nstop+1
      endif

c further checks to be inserted here
c initialize filter direct access data sets
      if(Ict(42).gt.0) call PREPDA(Buff)
 
      return
      end
