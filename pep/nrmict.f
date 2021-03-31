      subroutine NRMICT
 
      implicit none

c
c m.e.ash    oct 1969     subroutine nrmict
c main program for reading observed minus theory and partials tapes
c and incrementing normal equations
 
c array dimensions
      include 'globdefs.inc'

c commons
      include 'aprtbf.inc'
      include 'bernum.inc'
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'filmis.inc'
      include 'filtim.inc'
      include 'inodta.inc'
      include 'ktrap.inc'
      integer*2 numobs, mspcx1(3),mspcx2(3)
      equivalence (Num2,numobs),(Mspcrd(1,1),mspcx1(1)),
     1 (Mspcrd(1,2),mspcx2(1))
      include 'lcntrl.inc'
      include 'namobs.inc'
      include 'namtim.inc'
      include 'nrmmat.inc'
      include 'nrmwgt.inc'
      include 'obsdta.inc'
      include 'prdmov.inc'
      include 'rtside.inc'
      include 'wrkcompr.inc'
      include 'zeroes.inc'

c external functions
      integer*4 LEG
c
c local
      real*10 ersav(3)
      character*2 ovrid(2)/'  ',' +'/
      integer*4 icnt(4)
      integer*4 nfvctk/-1/
      real*4    a(3)
      logical   lhsflg, rhsflg
      real*10 weight
      real*10 dist, ersav1, fnobs
      integer*4 i, iaprio, incore, ipoch, ivectk, keasmt, l1, l2,
     .          measav, nbytes, ncodsv, nfirst, nplnsv
      integer*4 j
c
c filter quantities
      real*10 jdf
 
c
c initialize filter quantities
      Twrite = 0
      Nor    = 0
      incore = -1
      do i = 1, 400
         Nwrite(i) = 0
      end do
c
c see if normal equations are to be saved for each
c observation series
      if(Imat1.gt.0) call SAVHED(Imat1,'NRMEQSER')
c
c iabs1=input obs-theory and parital derivitive tapes
c imat1=output saved normal equations tape
c
c initialization
      Mesmt1 = 0
      Neasmt = 0
      ncodsv = 0
      ivectk = 1
      if(Jct(60).gt.0) ivectk = 0
      iaprio = 0
      nplnsv = 0
      Line   = 66
      call PAGSET('FORM NORMAL EQS FROM OBS-TH ', 7)
      call PAGSET(
     .'NTAP NSEQ TYPOBS PLANET NPLNT  SITE1 SERIES  SITE2  SPOT     ERRO
     .R WEIGHTS     IOBCON INCLUDED   DELETED   AVERAGE   RT.MEAN SQUARE
     . ', -33)
      do i = 1, 3
         Ermes1(i) = 0._10
      end do
 
c
c niobc =-1 do not read iobcon
c niobc = 0 read first record of series but not interior records
c niobc = 1 read interior records of series but not first record
c
      Niobc = 0
      if(Iobcon.le.0) Niobc = -1
c
c increment tape counter
      Jtape = 0
      do while( .true. )
         call PRDOBS(1,1)
         if(Iabs1.le.0) then
c
c error analysis for observations that went into forming
c normal equations
            call PAGCHK(60,9,0)
            fnobs = Mesmt1
            if(Mesmt1.le.0) fnobs = 1._10
            a(1)  = Ermes1(1)/fnobs
            a(2)  = Ermes1(2)/fnobs
            fnobs = Ermes1(3)/fnobs
            a(3)  = SQRT(fnobs)
            write(Iout,20) Mesmt1,a,fnobs,Ermes1(3)
   20       format('0ERROR ANALYSIS FOR THE', i8, ' MEASUREMENTS USED TO
     . INCREMENT THE NORMAL EQUATIONS FROM OBS-TH TAPES'/
     .      '          AVERAGE (OBS-TH)/ERROR', 1pe14.5/
     .      '       AVERAGE ABS(OBS-TH)/ERROR', 1pe14.5/
     .      ' ROOT MEAN SQUARE (OBS-TH)/ERROR', 1pe14.5/
     .      '     AVERAGE ((OBS-TH)/ERROR)**2', 1pe14.5/
     .      '         SUM ((OBS-TH)/ERROR)**2', 1pe14.5)
            keasmt = Neasmt - Mesmt1
            write(Iout,40) keasmt,Neasmt
   40       format(/i8,' OF THE TOTAL', i8,
     .             ' MEASUREMENTS WERE DELETED *')
            if(Niobc.ge.0) then
               rewind Iobcon
               Itrwnd(Iobcon) = 0
            endif
c
c terminate saved normal equations tape
            if(Imat1.gt.0) then
               write(Imat1) izr4,(Zero(1),i = 1,18)
               rewind Imat1
            endif
c
c first write last epoch to direct access, then
c read direct access and rewrite into filter input sne
            if(Ict(42).gt.0) then
               call WFILDA(B,Side,incore,Nparam,Sigma,lhsflg,
     .                     rhsflg)
               call WINSNE(B,B,Side,Sigma,Nparam,lhsflg,rhsflg,
     .                     Mesmt1, Ermes1)
            endif
c
c print out timer information
            call TIMRIT(
     .'FORMING NORMAL EQNS FROM OBSERVED MINUS THEORY AND PARTIAL DERIVA
     .TIVE TAPES ', 19)
 
            return
         endif
 
         call PAGHED(0)
c
c if ict(80)=1, move in correct initial conditions of
c embary,moon and earth,moon rotation
         if((Ict(80).gt.0) .and. (Iterat.le.1)) then
            call LIBCND(Jdem9,Econd9,Lem,Jdem0,Econd)
            call LIBCND(Jdmn9,Mcond9,Lmn,Jdmn0,Mcond)
            call LIBCND(Jder9,Ercnd9,Ler,Jder0,Ercond)
            call LIBCND(Jdmr9,Mrcnd9,Lmr,Jdmr0,Mrcond)
         endif
 
         do while( .true. )
c
c read first record of observation series
            call PRDOBS(1,3)
            if(Ncodf.le.0) goto 200
            nfirst = 0
            if((Nplnt0.ne.nplnsv) .or. (Ncodf.ne.ncodsv)) then
               call PAGCHK(59,1,1)
               write(Iout,50)
   50          format('    ')
               nplnsv = Nplnt0
               ncodsv = Ncodf
            endif
            do i = 1, 4
               icnt(i) = 0
            end do
c
c if ict(80)=1, move in correct initial conditions of
c planet and/or probe
            if((Ict(80).gt.0) .and. (Iterat.le.1)) then
               if(Klan.gt.0 .and. Klan.le.u_mxpl) then
                  call LIBCND(Jdpp0,Ppcond,Lpl(1,Klan),
     .             Jdpl0(Klan),Pcond(1,Klan))
               endif
               if(Klanb.gt.0)
     .             call LIBCND(Jdbb0,Bbcond,Lpl(1,Klanb),
     .             Jdpl0(Klanb),Pcond(1,Klanb))
               if(Klans1.gt.0)
     .             call LIBCND(Jevv0(1),Vvcone(1,1),Lpl(1,Klans1),
     .             Jdpl0(Klans1),Pcond(1,Klans1))
            endif
c
c save statistics
c get correct observation type name
            do i = 1, 3
               ersav(i) = Ermeas(i)
               measav   = Measmt
            end do
            if(Ncodg.eq.4) then
               if(Klanb.gt.0 .and.
     .             (Ncp0.eq.3 .or. Ncp0.eq.10)) Ncodg = 20
               if(Ksite(1).gt.0 .and. Ksite(1).ne.3)
     .             Ncodg = 21
            endif
 
            do while( .true. )
c
c read observation record
               call PRDOBS(1,4)
               if(Ncode.le.0) then
c
c printout statistics for observation series
                  call PAGCHK(58,1,1)
                  measav = Measmt - measav
                  Mesmt1 = Mesmt1 + measav
                  do i = 1, 3
                     ersav(i)  = Ermeas(i) - ersav(i)
                     Ermes1(i) = Ermes1(i) + ersav(i)
                  end do
                  fnobs  = measav
                  ersav1 = 0._10
                  if(fnobs.gt.0._10) then
                     ersav1 = ersav(1)/fnobs
                     fnobs  = SQRT(ersav(3)/fnobs)
                  endif
                  write(Iout,60) Ntape,Nseq,Typobs(Ncodg),Ncodf,
     .                            Plnnam, Nplnt0,
     .                            (Sitf(i,1),i = 1,2),Series,
     .                            (Sitf(i,2),i = 1,2),Spota,
     .                            (Erwgt2(i),i = 1,2),
     .                            ovrid(Irid1),ovrid(Irid2),icnt,
     .                            ersav1, fnobs
   60             format(i4,i5,a6,i2,1x,2A4,i3,1x,2A4,1x,
     .                   1A4, 1x, 2A4, 1x, 1A4, 1p, 2E12.5, 2A2,
     .                   2I5, 1x, 2I5, 1p, 2D13.5)
                  if(Ncodf.gt.20) then
                     write(Iout,70) Nplnt2,Spota2
   70                format(26x,i3,24x,1A4)
                     Line = Line + 1
                  endif
                  Neasmt = icnt(1) + icnt(2) + icnt(3) + icnt(4)
     .                     + Neasmt
c
c save normal equations for this series
                  if(Imat1.gt.0 .and. measav.gt.0) then
                     write(Imat1) measav,Ntape,Nseq,Erwgt2,
     .                     ersav, Zero(1),Zero(1),iaprio,
     .                     ivectk, (Zero(1),i = 1,10)
                     if(ivectk.gt.0) write(Imat1) nfvctk,
     .                  (Vectk(i),i = 1,Nparam)
                     write(Imat1) Nparam,
     .                     (Side(i),i = 1,Nparam)
                     call FWSIG(Imat1,Nparam,B,Sigma)
                     call NRMSET(0)
                  endif
                  goto 120
               else
                  dist = Save(3)
                  if(dist.lt.0._10 .or. dist.gt.1.E5_10) dist = 0._10
c
c delete bad observations
                  if(Nice.le.0) then
                     if(ABS(Deriv(2,1)).lt.Eps(Jacc+1)
     .                   *Deriv(1,1)) then
                        icnt(1) = icnt(1) + 1
                     else
                        Num1    = 2
                        icnt(3) = icnt(3) + 1
                     endif
                  endif
 
                  if(Nice.ge.0) then
                     if(ABS(Deriv(2,Num2)).lt.Eps(Jacc+Num2)
     .                   *Deriv(1,Num2)) then
                        icnt(2) = icnt(2) + 1
                     else
                        Num2 = 1
                        if(Nice.gt.0) Num2 = 0
                        icnt(4) = icnt(4) + 1
                     endif
                  endif
c
c get pointers
                  call PRDEQS(nfirst)
 
c if there are correlated observations to be used
c and extra partials,  iptr must be increased
                  if((Ict(19).ne.0) .and. (Next.ne.0)) then
 
c read in parameter names from ibuf1
                     rewind Ibuf1
                     read(Ibuf1) ((Nm(i,j),i=1,2),j = 1,Nparam)
                     rewind Ibuf1
 
c compare exnams to these names and add to iptr
                     nbytes = 16
                     do i = 1, Next
                        Iptr(Numpar + i) = 0
                        l1 = nbytes*(i - 1) + 1
                        do j = 1, Nparam
                           l2 = nbytes*(j - 1) + 1
                           if(LEG(nbytes,l1,Exnams,l2,Nm).eq.0) then
                              Iptr(Numpar + i) = j
                              goto 80
                           endif
                        end do
   80                end do
                  endif
c
c
c ict(17)= 0 form both the left hand side (lhs) and the right hand
c         side (rhs) of the normal equations
c ict(17)= 1 form the rhs only
c ict(17)=-1 form the lhs only
c
                  lhsflg = Ict(17).le.0
                  rhsflg = Ict(17).ge.0
 
c
c in old pep, call to nrmcfz used to be here
c
 
c if obs deleted, get another
                  if(Num1.le.Num2) then
c
c if filter set jdf = jd.fract of obs
                     if(Ict(42).gt.0) then
                        jdf = Jd + Ctrecf - dist/86400._10
                        if(jdf.ge.Fep(1)) then
c
c calculate index to epoch for this obs
                           do i = 1, Nepoch
                              if(jdf.lt.Fep(i+1)) then
                                 ipoch = i - 1
                                 if((ipoch+1).gt.Nepoch)
     .                            goto 90
                                 if(incore.eq.-1) incore = ipoch
c
c this obs in range so bump counters
                                 Nwrite(ipoch+1)= Nwrite(ipoch+1) + 1
                                 Twrite = Twrite + 1
c
c see if this epoch in core
                                 if(ipoch.ne.incore) then
c
c write out normal eqn for old incore epoch
                                    call WFILDA(B,Side,incore,
     .                                 Nparam, Sigma, lhsflg, rhsflg)
                                    incore = ipoch
c
c see if this epoch not written out yet
                                    if(Nwrite(ipoch+1).gt.1) then
c
c not in core and already written so read back in
                                      call RFILDA(B,Side,ipoch,
     .                                   Nparam, Sigma, lhsflg,
     .                                   rhsflg)
                                    else
                                      call NRMSET(1)
                                    endif
                                 endif
                                 goto 100
 
                              endif
 
                           end do
 
                        endif
c
c adjust counters to reflect obs out of range
   90                   if(Nice.le.0) then
                           icnt(1) = icnt(1) - 1
                           icnt(3) = icnt(3) + 1
                           Nor     = Nor + 1
                        endif
                        if(Nice.ge.0) then
                           icnt(2) = icnt(2) - 1
                           icnt(4) = icnt(4) + 1
                           Nor     = Nor + 1
                        endif
 
                        goto 110
 
                     endif
c
c apply tape-wide weights
  100                do i = Num1, Num2
                        Deriv(1,i) = Deriv(1,i)/Wgtobs(Jtape)
                     end do
c
c increment normal equations
c note: the following loop used to be the routine nrmsid
c
                     do j = Num1, Num2
                        Measmt    = Measmt + 1
                        weight    = Deriv(2,j)/Deriv(1,j)
                        Ermeas(1) = Ermeas(1) + weight
                        Ermeas(2) = Ermeas(2) + ABS(weight)
                        Ermeas(3) = Ermeas(3) + weight**2
                        Wobsth    = 1._10/Deriv(1,j)**2
                        Wobst1    = 1._10/Deriv(1,j)
 
c     increment r. side and coeff. matrix of normal eqs.
                        if(Nmp2.le.2 .or. Ict(19).le.0) then
                           call NRMCOF(Deriv(1,j),Deriv(1,j),
     .                        Numpar, Nmp2, ict(19),lhsflg,rhsflg,B)
                        else
                           call NRMCOF(Deriv(1,j),Dvx(1,j),Numpar,
     .                        Nmp2, ict(19),lhsflg,rhsflg,B)
                        endif
                     end do
 
                  endif
               endif
  110       end do
  120    end do
  200 end do
      end
