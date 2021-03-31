      subroutine SOLPPR(scale)
 
      implicit none
c
c           r.w.babcock  june 30, 1983
c
c     solve for parameters which have been partially prereduced away
c     using formalism of memo 75-3, distribution/r.d.reasenberg,
c     29 july 1975.

c           parameters
      real*10 scale(1000)
c
c array dimensions
      include 'globdefs.inc'

c commons
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'iptrst.inc'
      include 'mcnfrm.inc'
      include 'nrmmat.inc'
      include 'restor.inc'
      real*10 fba(1000)
      equivalence (Sav,fba)
      include 'rtsidesl.inc'
      common/WRKCOM/ Iptr(1000),Oldsol(1000),Zbar(1000),Dvars(1000)
      real*10 Oldsol,Zbar,Dvars
      integer*2 Iptr
 
c local
      real*10 ai, DOTN, dum, sigd, term
      integer   i, ia, iaprio, ic, id, idum, imats, ippr, ivectk, izn, 
     .          j, k, kc, krow, l, lc, nflgd, nflgz, nparm
      equivalence (ai,ia)
      character*4 subscr/'(  )'/
      character*8 zerneg(2)/'  ZERO', 'NEGATIVE'/
      integer*4 ierzn(2)
 
      ierzn(1) = 0
      ierzn(2) = 0
      do i = 1, Nummt0
         imats = Imat0(i)
         call EBCDIX(i,subscr,2,2)
         call FRMHED(imats,'IMAT0', subscr, 2, ippr, 1)
         if(ippr.ne.1) call SUICID(
     .    'NORMAL EQUATIONS NOT PREREDUCED, STOP IN SOLPPR ',12)
         call ZFILL(Zbar,16*Nparam)
         call ZFILL(Iptr,2*Mparam)
         do j = 1, Mparam
            ia     = j
            Sav(j) = ai
         end do
         call FRMMVE(Zbar,Nparam)
         do j = 1, Nparam
            ai = Zbar(j)
            if(ia.ne.0) Iptr(ia) = j
         end do
c
c read first record of sne
         read(imats) (idum,j = 1,5),(dum,j = 1,5),iaprio,ivectk
c
c read pointer group for pre-reduced sne
         read(imats) nparm,Ncparm,Ndparm,Znsqpp,Nserpp,Nauxpp,
     .               (Icrest(j),j = 1,Ncparm),
     .               (Idrest(j),j = 1,Ndparm)
c
c skip over reduced normal equations
         call BSKIP(imats,Ncparm)
c
c skip over reduced a priori information
         if(iaprio.gt.0) call BSKIP(imats,Ncparm)
c
c copy (subset of) solution vector in order of c matrix
c
         call ZFILL(Oldsol,16*Ncparm)
         do j = 1, Ncparm
            ic = Iptr(Icrest(j))
            Icrest(j) = ic
            if(ic.ne.0) Oldsol(j) = Solut(ic)
         end do
c
c restore d variances
         read(imats) nflgd,(Dvars(j),j = 1,Ndparm)
         if(nflgd.eq.-1) then
c
c restore zbar
            read(imats) nflgz,(Zbar(j),j = 1,Ndparm)
         else
            nflgz = nflgd
            do j = 1, Ndparm
               Zbar(j) = Dvars(j)
            end do
         endif
         if(nflgd.ne.-1) call SUICID(
     .       'MISSING VARIANCES FOR UNINTERESTING PARAMETERS IN SOLPPR'
     .       , -14)
         if(nflgz.ne.Ndparm) call SUICID(
     .   'INVALID SOLN FLAG FOR UNINTERESTING PARAMETERS, STOP SOLPPR '
     .   , 15)
c
c restore fbar-adjoint row-by-row
         do j = 1, Ndparm
            read(imats) idum,(fba(k),k = 1,Ncparm)
            id = Iptr(Idrest(j))
            if(id.ne.0) then
               if(Solut(id).ne.0._10) then
                  call PAGCHK(60,2,0)
                  write(Iout,10) id
   10             format('0WARNING: PARAMETER', i4,
     .' (SEE SOLUTION) APPEARS IN BOTH SAVED SOLN AND PPR NEQ''S, OR TWO
     .PPR NEQ''S')
                  if(Mout.ne.0) write(Mout,10) id
               endif
               Solut(id) = Zbar(j) - DOTN(Oldsol,fba,Ncparm)
               if(nflgd.eq.-1) then
                  sigd = Dvars(j)
c note: fbar-adjoint is properly scaled, but b is not
c it's less work to apply the scale factors to fba, though.
                  do k = 1, Ncparm
                     kc = Icrest(k)
                     if(kc.gt.0) fba(k) = fba(k)*scale(kc)
                  end do
c multiply b by row of fba, left and right, to get scalar
c use symmetry of b to halve the work
                  do k = 1, Ncparm
                     kc = Icrest(k)
                     if(kc.gt.0) then
                        krow = (kc*(kc-1))/2
                        do l = 1, Ncparm
                           lc = Icrest(l)
                           if(lc.gt.0 .and. lc.le.kc) then
                              term = B(krow + lc)*fba(k)*fba(l)
                              if(l.ne.k) term = term + term
                              sigd = sigd + term
                           endif
                        end do
                     endif
                  end do
                  Sigma(id) = SQRT(ABS(sigd))
                  if(sigd.le.0._10) then
                     call PAGCHK(60,1,0)
                     izn = 1
                     if(sigd.lt.0._10) izn = 2
                     write(Iout,20) id,zerneg(izn),j
   20                format(' VARIANCE', i4, ' IS ', a8, ' (', i4,
     .                      'TH OF UNINTERESTING PARAMETERS)')
                     ierzn(izn) = ierzn(izn) + 1
                  endif
               endif
            endif
         end do
 
         rewind imats
         Itrwnd(imats) = 0
      end do
 
      do i = 1, 2
         if(ierzn(i).ne.0) then
            call PAGCHK(60,1,0)
            write(Iout,40) ierzn(i),zerneg(i)
   40       format(i5,' UNINTERESTING VARIANCES WERE ', a8)
         endif
      end do
 
      return
      end
