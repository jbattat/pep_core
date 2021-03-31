      subroutine PRNFIL(noprnt, nstop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   j, j0, k, kj, linc, lm, n, noprnt, nstop
 
c*** end of declarations inserted by spag
 
 
c
c           j.f.chandler - 1980 aug - subroutine prnfil
c        based on z. goldberg's subr. finprt
c
c        suroutine to print out input parameters for kalman filtering.
c
c        common
      include 'fcntrl.inc'
      include 'filtda.inc'
      include 'filtds.inc'
      include 'filnit.inc'
      include 'filtim.inc'
      include 'inodta.inc'
 
      integer*2 j0s(6), j1s(6)
      character*1 blank/' '/
 
      if(Ict(42).gt.0) then
 
c check for correct number of input epochs
         if(Nepcht.le.Nepoch) then
            call PAGCHK(60, 2, 0)
            write(Iout, 20) Nepcht
   20       format('0***', i5, ': TOO FEW FILTER EPOCHS INPUT, ERROR',
     .             ' IN PRNFIL ***')
            nstop = nstop + 1
         endif
 
c see if printout desired
         if(noprnt.le.0) then
 
            call PAGCHK(60, 5, 0)
            write(Iout, 40) Ict(5), Ict(42), Jct(56)
   40       format('-KALMAN FILTER INPUT & CONTROL PARAMETERS       ',
     .          '   ICT(5)=', i3, '   ICT(42)=', i3, '   JCT(56)=', i3)
 
            write(Iout, 60) Insne, Filter, Outsne, Smat, Iconof, Wzero
   60       format('0DATA SETS:', 7x, 'INSNE=', i3, '   FILTER=', i3,
     .             '   OUTSNE=', i3, '   SMAT=', i3, '   ICONOF=', i3,
     .             '   WZERO=', i3)
c
c print solution controls, lfilt
            call PAGSET('FILTER SPAN SMOOTHING CONTROLS  ', -8)
            if(Line.gt.55) call NEWPG
            call PAGHED(0)
            j0 = 1
            kj = 0
            do j = 2, 251
               if(Lfilt(j).ne.Lfilt(j0) .or. j.ge.251) then
                  kj = kj + 1
                  j0s(kj) = j0
                  j1s(kj) = j - 1
                  j0 = j
                  if(kj.ge.6 .or. j.ge.251) then
                     call PAGCHK(59, 1, 1)
                     write(Iout, 70) (blank, j0s(k), j1s(k), Lfilt(j0s(
     .                               k)), k = 1, kj)
   70                format(6(a1,'  LFILT(',i3,'-',i3,')=',i3))
                     kj = 0
                  endif
               endif
            end do
c
c print control constants, filprm
            do j = 1, 46
               lm = 51 - j
               if(Filprm(lm).ne.0._10) goto 80
            end do
   80       linc = 3 + (lm - 1)/4
            call PAGCHK(60, linc, 0)
 
            write(Iout, 100) (blank, j, Filprm(j), j = 1, lm)
  100       format('0FILTER CONTROL CONSTANTS'/a1, '  TIMEZ', i2, ' =',
     .             1pd21.14, a1, '  DELTA', i2, ' =', d21.14, a1,
     .             ' DELPRD', i2, ' =', d21.14, a1, 'FILPRM(', i2, ')=',
     .             d21.14, a1/(1x,4('FILPRM(',i2,')=',d21.14,a1)))
c
c print control integers, fict
            lm   = 50
            linc = 3 + (lm - 1)/10
            call PAGCHK(60, linc, 0)
 
            write(Iout, 120) (blank, j, Fict(j), j = 1, lm)
  120       format('0FILTER CONTROL INTEGERS'/10(a1,'FICT(',i2,')=',i3))
c
c print test quantities, feps
            lm   = 20
            linc = 3 + (lm - 1)/5
            call PAGCHK(60, linc, 0)
 
            write(Iout, 140) (blank, j, Feps(j), j = 1, lm)
  140       format('0FILTER TESTING QUANTITIES'/5(a1,' FEPS(',i2,')=',
     .             1pe15.8))
c
c print prediction controls, lprdct
            call PAGSET('FILTER POST-FIT PREDICTION CONTROLS ', -9)
            if(Line.gt.55) call NEWPG
            call PAGHED(0)
            j0 = 1
            kj = 0
            do j = 2, 251
               if(Lprdct(j).ne.Lprdct(j0) .or. j.ge.251) then
                  kj = kj + 1
                  j0s(kj) = j0
                  j1s(kj) = j - 1
                  j0 = j
                  if(kj.ge.6 .or. j.ge.251) then
                     call PAGCHK(59, 1, 1)
                     write(Iout, 150) (blank, j0s(k), j1s(k), Lprdct(
     .                                j0s(k)), k = 1, kj)
  150                format(6(a1,'  LPRDCT(',i3,'-',i3,')=',i3))
                     kj = 0
                  endif
               endif
            end do
c
c print smear matrix weights, sweigh
            call PAGSET('FILTER SMEAR MATRIX WEIGHTS ', -7)
            if(Line.gt.55) call NEWPG
            call PAGHED(0)
            j0 = 1
            kj = 0
            do j = 2, 251
               if(Sweigh(j).ne.Sweigh(j0) .or. j.ge.251) then
                  kj = kj + 1
                  j0s(kj) = j0
                  j1s(kj) = j - 1
                  j0 = j
                  if(kj.ge.3 .or. j.ge.251) then
                     call PAGCHK(59, 1, 1)
                     write(Iout, 160) (blank, j0s(k), j1s(k), Sweigh(
     .                                j0s(k)), k = 1, kj)
  160                format(3(a1,'  SWEIGH(',i3,'-',i3,')=',1pd22.15))
                     kj = 0
                  endif
               endif
            end do
c
c print parameter names, pnames
            linc = 3 + (Npnp - 1)/5
            call PAGCHK(60, linc, 0)
            write(Iout, 180) (j, (Pnames(n,j),n=1,2), j = 1, Npnp)
  180       format('0PROCESS NOISE PARAMETER NAMES'/5(i5,'. ',a8,1x,a8)
     .             )
c
c print parameter values, pnprms
            linc = 3 + (Npnp - 1)/4
            call PAGCHK(60, linc, 0)
            write(Iout, 200) (blank, j, Pnprms(j), j = 1, Npnp)
  200       format('0PROCESS NOISE PARAMETER VALUES'/4(a1,'PNPRMS(',i2,
     .             ')=',1pd20.13))
c
c print filter epochs, fep
            lm   = Nepoch + 1
            linc = 3 + (lm - 1)/4
            call PAGCHK(60, linc, 0)
            write(Iout, 220) (blank, j, Fep(j), j = 1, lm)
  220       format('0FILTER EPOCHS'/4(a1,'  FEP(',i2,')= ', f17.9,4x))
 
            call PAGCHK(60, 3, 0)
            write(Iout, 240) Npnp, Nepoch, Maxe, Maxp
  240       format('0OTHER FILTER INPUT QUANTITIES'/' NPNP=', i3,
     .             '  NEPOCH=', i3, '  MAXE=', i3, '  MAXP=', i3)
         endif
      else if(noprnt.le.0) then
         call PAGCHK(60, 2, 0)
         write(Iout, 250)
  250    format('0THERE IS NO KALMAN FILTERING')
      endif
 
      return
 
      end
