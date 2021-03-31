      subroutine DIFNOM(ic, irc, dif)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   klan, klmhar, klnhar, m1, mphase, mrbias,
     .          mspot, n1, nshp
      real*10 pper1s, ppers
 
c*** end of declarations inserted by spag
 
 
c       subroutine difnom - j.f.chandler - 1980 may
c           based on subr. frmmve - m.e.ash, 1970 feb
c       compare input nominal values with those just read from sne
c       data set.  optionally, copy sne values into current
c       storage.  form and return difference vector (saved-input) of
c       all adjustable parameters.
c       return count of non-zero differences
c arguments
      integer ic,irc
      real*10 dif(1)

c array dimensions
      include 'globdefs.inc'
c commons 
      include 'anctrl.inc'
      include 'dtparm.inc'
      include 'empcnd.inc'
      include 'ethhar.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'mcnfrm.inc'
      include 'mhrfrm.inc'
      include 'monhar.inc'
      include 'namtim.inc'
      include 'param.inc'
      include 'plnhar.inc'
 
      integer*4 ngdpts(4)
      equivalence (ngdpts, Scontl(1,9))
 
      include 'psrstf.inc'
      include 'psrstm.inc'
      include 'restor.inc'
 
c nsav= index for difference vector (jumps around in routine)
c (if ntop positive)
c sav = difference vector for adjustable parameters
c
      include 'scoef4.inc'
      include 'wrkcomrs.inc'
 
      integer*2 mpl10/10/, mpl0/0/, mplstr/-4/
      character*8 blank/' '/, emnam(5)/' EMBARY ',' EROTAT ','  MOON  ',
     .          ' MROTAT ', ' EARTH  '/
      character*8 plsr/'PLSR....'/
      character*6 qcond/'COND  '/, harn(3)/'ZHAR  ', 'CHAR  ', 'SHAR  '/
      integer*4 ndif(2)
      logical*4 cdf(5)
      integer   i, j, k
 
c
c initialize vector pointer and counters
      Nsav    = 0
      ndif(1) = 0
      ndif(2) = 0
      cdf(1)  = mod(ic, 2).eq.1
      cdf(2)  = mod(ic/2, 2).eq.1
      cdf(3)  = mod(ic/4, 2).eq.1
      call ZFILL(Sav, 16*Nparam)
      call DIFPTS
c
c observing site coordinates
      call DIFSIT(ndif, cdf)
c
c equinox-equator-declination corrections
      call DIFEQN(ndif, cdf)
c
c sky corrections (star catalog error models)
      call DIFSTR(ndif, cdf)
c
c solar system parameters
      Nsav = Lprm0
      call DIFPRM(prmter, Prmtr1, Lprm, 100, ndif, cdf, blank, 'PRMTER')
c
c earth-moon barycenter init.cond, and earth parameters
      call DIFBDY(Econd, Econd1, Lem, 30, ndif, cdf, emnam(1), qcond)
c
c earth rotation initial conditions and parameters
      call DIFBDY(Ercond, Ercnd1, Ler, 30, ndif, cdf, emnam(2), qcond)
c
c et-ut2 or a1-ut1 and wobble parameters
      if(Numdt.ne.0 .and. Mumdt.ne.0)
     .    call DIFPM4(Dt, Dt1, Ldt, 600, ndif, cdf, blank, 'DT    ')
c
c earth gravitational potential harmonic coefficients
c zonal
      Nsav = Lehar1 - 1
      n1   = Nezone - 1
      m1   = Mezone - 1
      call DIFHAR(Ezhar, Ezham, Lezhar, n1, m1, 1, ndif, cdf, emnam(5),
     .            harn(1))
 
c tesseral cosine
      n1 = (Netess*(Netess+1))/2 - 1
      m1 = (Metess*(Metess+1))/2 - 1
      call DIFHAR(Echar, Echam, Lechar, n1, m1, 1, ndif, cdf, emnam(5),
     .            harn(2))
 
c tesseral sine
      call DIFHAR(Eshar, Esham, Leshar, n1, m1, 1, ndif, cdf, emnam(5),
     .            harn(3))
c
c moon initial conditions and parameters
      call DIFBDY(Mcond, Mcond1, Lmn, 30, ndif, cdf, emnam(3), qcond)
c
c moon rotation initial conditions and parameters
      call DIFBDY(Mrcond, Mrcnd1, Lmr, 30, ndif, cdf, emnam(4), qcond)
c
c moon gravitational potential harmonic coefficients
c zonal
      n1 = Nmzone - 1
      m1 = Mmzone - 1
      call DIFHAR(Mzhar, Mzham, Lmzhar, n1, m1, 1, ndif, cdf, emnam(3),
     .            harn(1))
 
c tesseral cosine
      n1 = (Nmtess*(Nmtess+1))/2 - 1
      m1 = (Mmtess*(Mmtess+1))/2 - 1
      call DIFHAR(Mchar, Mcham, Lmchar, n1, m1, 1, ndif, cdf, emnam(3),
     .            harn(2))
 
c tesseral sine
      call DIFHAR(Mshar, Msham, Lmshar, n1, m1, 1, ndif, cdf, emnam(3),
     .            harn(3))
c
c moon spot coordinates
      mspot = 0
      call DIFSPT(mspot, mpl10, ndif, cdf)
c
c moon radar observation biases
      mrbias = 0
      call DIFRBS(mrbias, mpl10, ndif, cdf)
c
c moon optical observation phase corrections
      mphase = 0
      call DIFPHS(mphase, mpl10, ndif, cdf)
c
c sun spot coordinates
      call DIFSPT(mspot, mpl0, ndif, cdf)
c
c sun radar observation biases
      call DIFRBS(mrbias, mpl0, ndif, cdf)
c
c sun optical observations phase corrections
      call DIFPHS(mphase, mpl0, ndif, cdf)
c
c start of planet loop
      do k = 1, Mumpln
c
c find planet for restored normal equations
         klan = 0
         if(Numpln.le.0) go to 100
         do i = 1, Numpln
            if(Nplnt(i).eq.Mplnt(k)) then
               klan = i
c
c planet initial conditions and parameters
               Nsav = Lpl1(klan) - 1
               call DIFBDY(Pcond(1,klan), Pcond1(1,k), Lpl(1,klan),
     .                     30, ndif, cdf, Aplnt(klan), qcond)
c see if there are gravitational potential harmonics
c for the planet
               klmhar = 0
               do j = 1, Mumphr
                  if(Mplhar(j).eq.Mplnt(k)) then
                     klmhar = j
                     go to 10
                  endif
                  end do
   10          if(klmhar.gt.0) then
                  klnhar = 0
                  do j = 1, Nmphar
                     if(Nplhar(j).eq.Nplnt(klan)) then
                        klnhar = j
                        go to 20
                     endif
                     end do
 
   20             nshp = Nshape(klnhar) + 1
                  if(nshp.eq.1) then
c
c planet gravitational potential harmonic coefficients
c zonal
                     m1 = Mpzone(klmhar) - 1
                     n1 = Npzone(klnhar) - 1
                     call DIFHAR(Pzhar(klnhar,1), Pzham(klmhar,1),
     .                           Lpzhar(klnhar,1), n1, m1, 4,
     .                           ndif, cdf, Aplnt(klan), harn(1))
 
c tesseral cosine
                     m1 = (Mptess(klmhar)*(Mptess(klmhar)+1))/2 - 1
                     n1 = (Nptess(klnhar)*(Nptess(klnhar)+1))/2 - 1
                     call DIFHAR(Pchar(klnhar,1), Pcham(klmhar,1),
     .                           Lpchar(klnhar,1), n1, m1, 4,
     .                           ndif, cdf, Aplnt(klan), harn(2))
 
c tesseral sine
                     call DIFHAR(Pshar(klnhar,1), Psham(klmhar,1),
     .                           Lpshar(klnhar,1), n1, m1, 4,
     .                           ndif, cdf, Aplnt(klan), harn(3))
                  else if(nshp.eq.2) then
c
c fourier
                     m1 = 122
                     n1 = 122
                     call DIFHAR(Pzhar(klnhar,1), Pzham(klmhar,1),
     .                           Lpzhar(klnhar,1), n1, m1, 4,
     .                           ndif, cdf, Aplnt(klan), 'FOUR. ')
                  else if(nshp.eq.3) then
c
c grid
                     m1 = Mngd(klmhar)
                     n1 = ngdpts(klnhar)
                  endif
               endif
c call difgrd(pzhar(klnhar,1),pzham(klmhar,1),lshp4(klnhar,1),
c 1 n1,m1,4,ndif,cdf,aplnt(klan),'grid  ')
c
c do not call spot or bias routines for planet rotation
               if(Mplnt(k).gt.0) then
c
c planet spot coordinates
                  call DIFSPT(mspot, Mplnt(k), ndif, cdf)
 
c planet radar observation biases
                  call DIFRBS(mrbias, Mplnt(k), ndif, cdf)
 
c planet optical observation phase corrections
                  call DIFPHS(mphase, Mplnt(k), ndif, cdf)
               endif
               go to 50
            endif
 
c no match
            end do
 
c end of planet loop
   50    end do
c
c star coordinates
  100 call DIFSPT(mspot, mplstr, ndif, cdf)
 
c star quantities in radar bias common
      call DIFRBS(mrbias, mplstr, ndif, cdf)
 
c star quantities in optical phase corr.common
      call DIFPHS(mphase, mplstr, ndif, cdf)
c
c pulsar parameters
      do k = 1, Mumpsr
c
c find corresponding input pulsar
         if(Numpsr.le.0) go to 200
         do i = 1, Numpsr
            if(Sptpsr(i).eq.Sptps1(k)) then
               plsr(5:8)  = Sptpsr(i)
               Nsav       = Lpsr0(i)
               ppers      = Psrcn(6, i)
               Psrcn(6,i) = ppers + Plspr(i)
               pper1s     = Psrcn1(6, k)
               Psrcn1(6,k)= pper1s + Plspr1(k)
               call DIFPRM(Psrcn(1,i), Psrcn1(1,k), Lpsrcn(1,i), 16,
     .                     ndif, cdf, plsr, 'CON   ')
               Psrcn(6,i)  = ppers
               Psrcn1(6,k) = pper1s
               go to 150
            endif
            end do
  150    end do
c
c return differences
  200 irc = ndif(1)
      if(Line.gt.56) call NEWPG
      write(Iout, 300) ic, ndif
  300 format('-SAVED NORMAL EQUATION NOMINAL PARAMETER VALUE SUMMARY'/
     .       5x, 'IC=', z8, ', DISCREPANCIES:', i4, ' FREE AND', i4,
     .       ' FIXED PARAMETERS')
      Line = Line + 4
      if(irc.gt.0) then
         do i = 1, Nparam
            dif(i) = Sav(i)
            end do
      endif
 
      return
      end
