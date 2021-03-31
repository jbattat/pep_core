      subroutine BDREED(jd, mplnt, kcall)
 
      implicit none
 
c     subroutine BDREED - M.E.Ash  1969 Aug
c     n-body tape is read either forward or backward in time
c     observations in error records are skipped if ICT(31) so indicates
c     copy data to individual body arrays as indicated

c arguments
      integer*4 jd,kcall
      integer*2 mplnt

c        jd   =julian date to be in middle of second record
c        mplnt=planet number whose reed routine is calling bdreed
c        kcall=caller's flag

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'b2ydta.inc'
      include 'bdctrl.inc'
      include 'bddta.inc'
      include 'bdydta.inc'
      include 'comdat.inc'
      include 'comnut.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'number.inc'
      include 'sbdta.inc'
      include 'scdta.inc'
      include 'tapdta.inc'
      include 'tapdte.inc'
      include 'tapdtm.inc'
      include 'tapdtp.inc'
      include 'trpcom.inc'
 
c local
      integer*4 i,ibdbad(3),ivl,ivlcm,j,j1,j2,j3,jdbdst,jdmt,jdx(3),
     . k,kbd,kc,l,l1,l2,l3,lock,m,m1,m2,m3,mm,mvl,n,nbdrec
c
c          kcall indicates caller:
c        1 - em, 2 - mn, 3 - pl, 4 - sb, 5 - sc, 6 - so
c        needed mainly just to distinguish 3 from 5
c
c           test to see if jd is on tape
      if(jd.gt.Jdbd1 .and. jd.lt.Jdbd2) then
         if(Kpert.le.0) goto 200
         if(jd.gt.Jdcom1 .and. jd.lt.Jdcom2) goto 200
         write(Iout, 20) jd, Jdcom1, Jdcom2
   20    format(i17, ' NOT ON CENTER-OF-MASS DATA SET (', i7, '-',
     .          i7, ')')
         goto 110
         endif
      write(Iout, 100) jd, Jdbd1, Jdbd2
  100 format(i17, ' NOT ON N-BODY DATA SET (', i7, '-', i7, ')')
  110 jd = 0
      return
c
c get correct records of n-body tape into storage
  200 n = Ibdsgn*(jd - Jdbd(2))
      if(n.lt.0) then
c
c correct records are behind on tape
         n = 4 - (n + 1)/Intbd5
         do i = 1, n
            nbdrec = nbdrec - 1
            backspace Ibody
            if(Kpert.gt.0) backspace Kpert
         end do
         goto 500
      else if(n.eq.0) then
         goto 1000
      else
c
c correct records are ahead on tape
         n = n/Intbd5
         if(n.eq.0) goto 1000
         mm = n - 3
         if(mm.lt.0) then
c
c correct middle record is no more than two ahead on tape
c storage must be shifted
            l2 = 5
            mm = 2
            l3 = n*5
            j3 = n*10
            m3 = n*40
            Jdbd(1)   = Jdbd(n + 1)
            Frct(1)   = Frct(n + 1)
            Ivel(1)   = Ivel(n + 1)
            Mvel(1)   = Mvel(n + 1)
            Ivcm(1)   = Ivcm(n + 1)
            Jdcom(1)  = Jdcom(n + 1)
            Frcom(1)  = Frcom(n + 1)
            ivl       = Ivel(1)
            mvl       = Mvel(1)
            ivlcm     = Ivcm(1)
            ibdbad(1) = ibdbad(n + 1)
            if(ibdbad(1).gt.0) goto 400
         else if(mm.eq.0) then
            goto 500
         else
            do i = 1, mm
               nbdrec = nbdrec + 1
               read(Ibody, err=210)
               if(Kpert.gt.0) read(Kpert)
               goto 240
  210          read(Ibody)
 
c second read of error record might not be needed if system changes
               write(Iout, 220) nbdrec, Ibody
  220          format(' **** ERROR RECORD', i6,
     .                ' SKIPPED ON N-BODY DATA SET', i3,
     .                ' IN BDREED LABEL=42 ****')
  240       end do
            goto 500
         endif
      endif
  300 l1 = l2 - 4
      j2 = l2*2
      j1 = j2 - 9
      m2 = l2*8
      m1 = m2 - 39
      do k = l1, l2
         l = k + l3
         do j = 1, Nbdy2
            do i = 1, ivl
               Body(i, k, j) = Body(i, l, j)
            end do
         end do
         if(Kpert.gt.0) then
            do i = 1, ivlcm
               Comcrd(i, k) = Comcrd(i, l)
            end do
         endif
      end do
      do k = j1, j2
         l = k + j3
         do i = 1, ivl
            Merc(i, k) = Merc(i, l)
         end do
      end do
      do k = m1, m2
         l = k + m3
         do i = 1, mvl
            Mon(i, k) = Mon(i, l)
         end do
         Psid(k) = Psid(l)
         Epsd(k) = Epsd(l)
         do i = 1, 3
            Librt(k, i) = Librt(l, i)
         end do
      end do
  400 if(n.eq.2) goto 600
      Jdbd(2)  = Jdbd(3)
      Frct(2)  = Frct(3)
      Ivel(2)  = Ivel(3)
      Mvel(2)  = Mvel(3)
      Ivcm(2)  = Ivcm(3)
      Jdcom(2) = Jdcom(3)
      Frcom(2) = Frcom(3)
      ivl      = Ivel(2)
      mvl      = Mvel(2)
      ivlcm    = Ivcm(2)
      n  = 2
      l2 = 10
      mm = 3
      ibdbad(2) = ibdbad(3)
      if(ibdbad(2).gt.0) goto 600
      goto 300
c
c read n-body data into storage
      entry BDRED1(jd)
      nbdrec = 0
c           on first entry, center-of-mass tape is not yet synched to
c           the n-body tape.  the two must have matching tabular
c           intervals and record epochs, but need not have the same
c           time spans -- we must line them up and mark them locked
c           together.  use jdcom(1) as flag for 1st c-of-m tape read.
c           (we also have jd=0.)
      Jdcom(1) = -1
      lock     = 0
 
  500 mm    = 1
  600 do l  = mm, 3
         l2 = l*5
         l1 = l2 - 4
         j2 = l*10
         j1 = j2 - 9
         m2 = l*40
         m1 = m2 - 39
  650    ibdbad(l) = 0
         Jdbd(l)   = 0
         nbdrec    = nbdrec + 1
         read(Ibody, err=800) Jdbd(l), Frct(l), ivl, mvl,
     .                          ((Merc(i,j),i=1,ivl), j = j1, j2),
     .                          (((Body(i,j,k),i=1,ivl),j=l1,l2),
     .                          k = 1, Nbdy2),
     .                          ((Mon(i,j),i=1,mvl), j = m1, m2),
     .                          (Psid(j), Epsd(j), j = m1, m2),
     .                          ((Librt(j,i),i=1,3), j = m1, m2)
         Ivel(l) = ivl
         Mvel(l) = mvl
         if(Kpert.le.0) goto 850
 
c must read from c-of-m tape 1st time so we can compare
         if(lock.eq.0 .and. Jdcom(1).lt.0 .or. lock.eq.1) goto 750
 
c not locked yet, see if times match
  700    i = (Jdbd(l) - Jdcom(l))*Ibdsgn
         if((i/Intbd5)*Intbd5.ne.i) call SUICID(
     .'CENTER-OF-MASS TABULAR POINTS OUT OF PHASE WITH N-BODY, STOP IN B
     .DREED  ', 18)
 
c if not locked, go read from the appropriate tape
         if(i.lt.0) goto 650
         if(i.eq.0) goto 850
  750    read(Kpert) Jdcom(l), Frcom(l), ivlcm,
     .               ((Comcrd(i,j),i=1,ivlcm), j = l1, l2)
         Ivcm(l) = ivlcm
         if(lock.eq.0) goto 700
         if(Jdbd(l).ne.Jdcom(l)) call SUICID(
     .       'CENTER-OF-MASS TAPE LOST SYNCH, STOP BDREED ', 11)
         goto 850
  800    ibdbad(l) = 1
         read(Ibody)
  850    lock = 1
      end do
 
      Nvels(8) = 9999
      Nbtrp(8) = 0
      Nvels(9) = 9999
      Nbtrp(9) = 0
      Nmmm=0
c
c reconstruct dates of start of bad records
      if(ibdbad(1).gt.0 .and. Jdbd(1).le.0) then
         if(ibdbad(2).le.0 .or. Jdbd(2).gt.0) then
            Jdbd(1) = Jdbd(2) - ISIGN(Intbd5, Ibdsgn)
            Jdbd(3) = Jdbd(2) + ISIGN(Intbd5, Ibdsgn)
            goto 900
         else if(ibdbad(3).le.0 .or. Jdbd(3).gt.0) then
            Jdbd(2) = Jdbd(3) - ISIGN(Intbd5, Ibdsgn)
            Jdbd(1) = Jdbd(2) - ISIGN(Intbd5, Ibdsgn)
            goto 900
         endif
         if(jd.le.0) call SUICID(
     .   ' 3 BAD RECORDS IN A ROW ON N-BODY TAPE, STOP IN BDREED  ', 14)
         Jdbd(1) = jdbdst + ISIGN(Intbd5*(nbdrec-3), Ibdsgn)
      endif
      Jdbd(2) = Jdbd(1) + ISIGN(Intbd5, Ibdsgn)
      Jdbd(3) = Jdbd(2) + ISIGN(Intbd5, Ibdsgn)
  900 if(jd.gt.0) goto 200
      jdbdst = Jdbd(1)
      if(Ibdsgn.gt.0) Jdbd1 = Jdbd(2)
      if(Ibdsgn.lt.0) Jdbd2 = Jdbd(2)
      return
c
c see if jd of observation is in bad records
 1000 do i = 1, 3
         if(ibdbad(i).gt.0) then
            write(Iout, 1020) jd, Ibody, ibdbad
 1020       format(i17,
     .             ' DELETED BECAUSE OF BAD RECORD ON N-BODY DATA SET',
     .             i3, ' WITH IBDBAD=', 3I2, '  ********')
            jd = -1
            if(Ict(35).eq.0) goto 1300
            if(Ict(35).gt.0 .and. (Nplnt0.le.30.or.Ict(35).gt.1)) return
            call SUICID(
     .         ' ERROR RECORD ON N-BODY DATA SET, STOP IN BDREED', 12)
         endif
      end do
c
c search for body of interest
      if(mplnt.eq.0) return
      do i = 1, Nbdy
         if(mplnt.eq.Npl(i)) then
            kbd = i - 1
            if(i.eq.Nbdy) goto 1600
            goto 1200
         endif
      end do
      write(Iout, 1100) jd, mplnt, (Npl(i), i = 1, Nbdy)
 1100 format('0INCONSISTENCY IN BDREED WITH JD=', i7, ', MPLNT=', i3,
     .       ', NPL=', 10I3)
      call SUICID(' PLANET NOT ON N-BODY DATA SET, STOP IN BDREED  '
     .            , 12)
c
c check if satellite or probe is body of interest
 1200 if(mplnt.gt.30 .or. Ncp(kbd+1).gt.0) call SUICID(
     .' SATELLITE OR PROBE NOT ALLOWED ON N-BODY DATA SET, STOP IN BDREE
     .D  ', 17)
      
      kc = ((jd-Jdbd(1))*Ibdsgn)/Intb(kbd+1)
      j1 = 0
c effectively equivalence (merc(1,31),body(1,1,1))
      if(kbd.gt.0) j1 = 15*(kbd+1)

      if(kcall.eq.6) then
         Ttb2(2) = Jdbd(1) + Ibdsgn*kc*Intb(kbd+1)
         Idxb2 = j1 + kc - 4
         return

      else if(kcall.eq.1) then
c
c earth-moon barycenter is body of interest
         if(Jdbd(1).eq.Jdem(1)) return
         do i = 1, 3
            Fem(i)  = 0._10
            Jdem(i) = Jdbd(i)
         end do
         do j = 1, 15
            do i = 1, 6
               Embary(i, 1, j) = Body(i, j, kbd)
            end do
         end do
         goto 1800
      else
c
c copy proper dates into jdx
         if(kcall.eq.5) then
            do i = 1, 3
               jdx(i) = Jdsc(i)
            end do
         else if(kcall.eq.4) then
            do i = 1, 3
               jdx(i) = Jdsb(i)
            end do
         else
            do i = 1, 3
               jdx(i) = Jdp(i)
            end do
         endif
 
         if(kbd.gt.0) then
c
c planet is body of interest
            if(Jdbd(1).eq.jdx(1)) return
            do i = 1, 3
               jdx(i) = Jdbd(i)
            end do
            goto 1500
         else
c
c mercury is body of interest
            n = (Jdbd(2) - Jdbd(1))/2
            if(Ibdsgn.gt.0) then
               if(jd.ge.jdx(2) .and. jd.lt.jdx(3)) return
               if(jd.lt.(Jdbd(2)+n)) goto 1400
            else
               if(jd.le.jdx(2) .and. jd.gt.jdx(3)) return
               if(jd.gt.(Jdbd(2)+n)) goto 1400
            endif
            jdx(1) = Jdbd(2)
            jdx(2) = Jdbd(2) + n
            jdx(3) = Jdbd(3)
            j1     = 10
            goto 1500
         endif
      endif
 1300 if(Ncodf.le.3) then
         call SUICID(
     .        ' ERROR RECORD ON N-BODY DATA SET, STOP IN BDREED', 12)
      endif
      return
 
 1400 jdx(1) = Jdbd(1) + n
      jdx(2) = Jdbd(2)
      jdx(3) = Jdbd(2) + n
      j1     = 5
 
c determine where to put info
 1500 if(kcall.eq.5) then
 
c caller was screed
         do j = 1, 15
            l = j1 + j
            do i = 1, 6
               Satprc(i, 1, j) = Merc(i, l)
            end do
         end do
         do i = 1, 3
            Fsc(i)  = 0._10
            Jdsc(i) = jdx(i)
         end do
      else if(kcall.eq.4) then
 
c caller was sbreed
         do j = 1, 15
            l = j1 + j
            do i = 1, 6
               Satprb(i, 1, j) = Merc(i, l)
            end do
         end do
         do i = 1, 3
            Fsb(i)  = 0._10
            Jdsb(i) = jdx(i)
         end do
      else
 
c caller was plreed
         do j = 1, 15
            l = j1 + j
            do i = 1, 6
               Planet(i, 1, j) = Merc(i, l)
            end do
         end do
 
c copy out record dates
         do i = 1, 3
            Fp(i)  = 0._10
            Jdp(i) = jdx(i)
         end do
      endif
      goto 1800
c
c moon is body of interest
 1600 n = (Jdbd(2) - Jdbd(1))/5
      if(Ibdsgn.gt.0) then
         if(jd.ge.Jdm(2) .and. jd.lt.Jdm(3)) return
      else if(jd.le.Jdm(2) .and. jd.gt.Jdm(3)) then
         return
      endif
      l = iabs(jd - Jdbd(2))
      m = 0
      do i = 1, 5
         m = m + n
         if(l.lt.iabs(m)) goto 1700
      end do
 1700 jdmt = Jdbd(2) + m - 3*n
      m1   = 24 + 2*iabs(m)
      do i = 1, 3
         jdmt    = jdmt + n
         Jdm(i)  = jdmt
         Fm(i)   = 0._10
         Jder(i) = Jdm(i)
      end do
      do j = 1, 24
         l = m1 + j
         do i = 1, 6
            Moon(i, 1, j) = Mon(i, l)
         end do
         Psidx(j) = Psid(l)
         Epsdx(j) = Epsd(l)
         do i = 1, 3
            Librat(j, i) = Librt(l, i)
         end do
      end do
c
c set interpolation flags to signal tape motion
 1800 Nvels(kcall) = 9999
      Nbtrp(kcall) = 0
 
      return
      end
