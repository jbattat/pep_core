      subroutine SBRAAD(jd,fract,*)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 fract, repoch, xx, yy, zz, zz1
      integer   i, ientr, iord, j, jd, jord, k, l, m, m1, m2, mm, n, 
     .          nrec, numpt, numpt1
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash, l.friedman- sept. 1969 --subroutine sbraad
c satellite-probe tape is read either forward or backward in time
c in the variable tabular interval mode
c modified for kalman filter by paul macneil dec., 1977
c
c           jd   = julian day no.
c           fract= fraction of day from midnight
c          *    = error exit for end of file on jsb
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'fcntrl.inc'
      include 'flcomp.inc'
      include 'inodta.inc'
      include 'number.inc'
      include 'sbdta.inc'
      include 'sbdtavtp.inc'
      include 'spqind.inc'
      include 'trpcom.inc'
 
      logical*4 frstrd
      integer*4 iend/0/
 
      ientr  = 0
      frstrd = .false.
c
c*  start=100
c see if the correct records are in storage
  100 if(Ifiltr.gt.0) then
         if(jd + fract.gt.Tnext) goto 500
      endif
      xx = jd - Jdsb(2)
      xx = xx + (fract - Fsb(2))
 
c xx= given time - time of record 2
      yy = jd - Jdsb(5)
      yy = yy + (fract - Fsb(5))
 
c yy= given time - time of record 5
      if(xx*yy.le.0._10) goto 1000
c
c see if correct records are ahead on tape
      if(ABS(xx).lt.ABS(yy)) then
c
c*  start=500
c correct records are behind on tape
         mm = 1
         zz = xx + hc1(2)
 
c zz= given time - time of record 1
         n = 8
         if(zz*xx.gt.0.0_10) then
            zz1 = zz + hc1(1)
 
c zz1= given time - time of record 0
            n = 9
            if(zz*zz1.gt.0.0_10) then
 
c correct records are back an unknown number of records
               n = -zz1/hc1(1)
               n = n/2 + 12
               if(n.gt.nrec/3) then
                  if(Ifiltr.gt.0) then
                     if(n.gt.nrec .and. n.lt.nrec + 5) n = nrec
                     if(n.le.nrec) goto 150
                     call SUICID('NO FILTER REWIND, STOP IN SBRAAD', 8)
                  endif
                  rewind Jsb
                  write(Iout,110) jd,fract,Jsb,Nplnt0
  110             format(i17,f14.13,' TOO FAR BEHIND ON TAPE', i3,
     .                   ' FOR SATELLITE-PROBE', i4, ', REWIND TAPE')
                  do i = 1, 2
                     read(Jsb,err=120)
                     goto 130
  120                read(Jsb)
  130             end do
                  goto 200
               endif
            endif
         endif
c
c*  start=800
c backspace records
  150    if(nrec.lt.n) then
c
c backspacing out of span
            write(Iout,160) jd,fract
  160       format(' JD.F=', i8, f9.8,
     .             ', BACKSPACE OUT OF DATA SPAN IN SBRAAD')
            jd = 0
            goto 1000
         else
            do i = 1, n
               backspace Jsb
            end do
            nrec = nrec - n
 
c must be within proper filter span, reset end counter
            iend = 0
            goto 300
         endif
      else
c
c see if we move 5 tabular points in storage
         zz = yy - hc1(6)
 
c zz= given time - time of record 6
         if(zz*yy.gt.0._10) then
c
c move 4 points in storage and then read more
            m1 = 4
            m2 = 2
            mm = 5
         else
            m1 = 5
            m2 = 1
            mm = 6
         endif
c
c*  start=200
c move tabular points in storage
         do l = 1, m1
            m = l + m2
            Jdsb(l)   = Jdsb(m)
            Fsb(l)    = Fsb(m)
            iorder(l) = iorder(m)
            jorder(l) = jorder(m)
            numprt(l) = numprt(m)
            hc1(l)    = hc1(m)
            hc2(l)    = hc2(m)
            iord      = iorder(l)
            do j = 1, iord
               do i = 1, 3
                  sprb(l,j,i) = sprb(m,j,i)
               end do
            end do
            if(numprt(l).gt.0) then
               numpt = numprt(l)
               jord  = jorder(l)
               do k = 1, numpt
                  do j = 1, jord
                     do i = 1, 3
                        dsprb(l,j,i,k) = dsprb(m,j,i,k)
                     end do
                  end do
               end do
            endif
         end do
         goto 400
      endif
c
c read variable tabular points into storage
      entry SBRAD1
      frstrd = .true.
      ientr  = 1
  200 nrec   = 0
      iend   = 0
  300 mm     = 1
  400 if(Ifiltr.le.0) goto 900
      if(frstrd) goto 700
      if((jd+fract).le.Tnext) goto 900
  500 if(Nlflag.eq.1) then
         write(Iout,550)
  550    format(
     .        ' WARNING: SBRAAD READING BEYOND END OF LAST FILTER SPAN'
     .  )
         Line = Line + 1
         goto 900
      endif
 
c loop through remainder of this integration segment
  600 do while( iend.lt.3 )
         read(Jsb) Jdsb(6),Fsb(6)
         if((Jdsb(6)+Fsb(6)).gt.Tnext) iend = iend + 1
      end do
 
c end of remainder reading loop
      Recalc = .true.
c
c*  start=1000
c start next filter span
  700 iend = 0
 
c read next header
      read(Jsb) Tprev,Tthis,Tnext,Iepoch,Nlflag
 
c skip reading nrcfil for now
      Nrcfil = -1
      if(Line.gt.56) call OBSPAG
      write(Iout,800) Tprev,Tthis,Tnext,Iepoch,Nlflag,Nrcfil
  800 format(' FILTER HEADER/TRAILER: TPREV=', f15.6, ' TTHIS=', f15.6,
     .       ' TNEXT=', f15.6, ' IEPOCH=', i4/' NLFLAG=', i2,
     .       ' NRCFIL=', i3)
      Line = Line + 2
      nrec = 0
 
c read first six records of this span
      mm = 1
c
c*  start=2000
  900 do l = mm, 6
         read(Jsb,end=1100) Jdsb(l),Fsb(l),iord,jord,numprt(l),numpt1,
     .                         ((sprb(l,j,i),i=1,3),j = 1,iord),
     .                         (((dsprb(l,j,i,k),i=1,3),j=1,jord),
     .                         k = 1, numpt1), hc1(l),hc2(l)
         nrec = nrec + 1
         iorder(l) = iord
         jorder(l) = jord
         if(Ifiltr.gt.0) then
            if((Jdsb(l)+Fsb(l)).gt.Tnext) iend = iend + 1
 
c end of segment, start reading new segment
            if(iend.ge.3) then
               Recalc = .true.
               goto 700
            endif
         endif
      end do
      if(ientr.ne.2) then
         if(.not. (frstrd)) goto 100
c
c*  start=9000
c get a sharper criterion for saying point is before start
         if(Idirsb.le.0) then
            Frsb2 = Fsb(2)
            Jdsb2 = Jdsb(2)
         else
            Frsb1 = Fsb(2)
            Jdsb1 = Jdsb(2)
         endif
      endif
 
 1000 return
 1100 return 1
c
c entry for checking epoch against filter span limit
c if beyond, then skip remaining records immediately
      entry SBRAD2(repoch)
      if(Ifiltr.le.0) return
      if(repoch.lt.Tnext) return
      ientr = 2
      goto 600
 
      end
