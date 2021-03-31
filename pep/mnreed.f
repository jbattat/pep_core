      subroutine MNREED(jd)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ivl, j, jdmnst, k, l, l1, l2, l3, mm, n,
     .          nmnrec
 
c*** end of declarations inserted by spag
 
 
c m.e.ash   may 1966   subroutine mnreed
c moon tape is read either forward or backward in time

c arguments
      integer*4 jd

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdat.inc'
      include 'comnut.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'number.inc'
      include 'plndta.inc'
      include 'tapdta.inc'
      include 'tapdtm.inc'
      include 'trpcom.inc'

c external functions
      integer*4 JBDTST

c local
      integer*4 intmn5/4/
      integer*2 np10/10/, imnbad(3)
c
c test to see if moon is on n-body tape
      if(JBDTST(np10).le.0) then
c
c test to see if jd is on tape
         if(jd.gt.Jdm1) then
            if(jd.lt.Jdm2) goto 100
         endif
         write(Iout,50) jd,Jdm1,Jdm2
   50    format(i17,' NOT ON MOON DATA SET (', i7, '-', i7, ')')
         jd = 0
      else
         call BDREED(jd,np10,2)
      endif
      return
c
c get correct records of moon tape into storage
  100 n = Idirmn*(jd - Jdm(2))
      if(n.lt.0) then
c
c correct records are behind on tape
         n = 4 - (n + 1)/4
         do i = 1, n
            nmnrec = nmnrec - 1
            backspace Imn
         end do
         goto 400
      else if(n.eq.0) then
         goto 800
      else
c
c correct records are ahead on tape
         n = n/4
         if(n.eq.0) goto 800
         mm = n - 3
         if(mm.lt.0) then
c
c correct middle record is no more than two ahead on tape
c storage must be shifted
            l1     = 1
            l2     = 8
            l3     = 8*n
            mm     = 2
            Jdm(1) = Jdm(n + 1)
            Jder(1)   = Jder(n + 1)
            Imvel(1)  = Imvel(n + 1)
            Fm(1)     = Fm(n + 1)
            ivl       = Imvel(n + 1)
            imnbad(1) = imnbad(n + 1)
            if(imnbad(1).gt.0) goto 300
         else if(mm.eq.0) then
            goto 400
         else
            do i = 1, mm
               nmnrec = nmnrec + 1
               read(Imn,err=110)
               goto 140
  110          read(Imn)
 
c second read of error record might not be needed if system changes
               write(Iout,120) Imn
  120          format(' **** ERROR RECORD SKIPPED ON MOON DATA SET',
     .                i3, ' IN MNREED 42 ****')
  140       end do
            goto 400
         endif
      endif
  200 do k = l1, l2
         l = k + l3
         Psidx(k) = Psidx(l)
         Epsdx(k) = Epsdx(l)
         do i = 1, 3
            Librat(k,i) = Librat(l,i)
         end do
         do j = 1, Iparm
            do i = 1, ivl
               Moon(i,j,k) = Moon(i,j,l)
            end do
         end do
      end do
  300 if(n.eq.2) goto 500
      Jdm(2)   = Jdm(3)
      Jder(2)  = Jder(3)
      Imvel(2) = Imvel(3)
      Fm(2)    = Fm(3)
      ivl      = Imvel(3)
      n  = 2
      l1 = 9
      l2 = 16
      mm = 3
      imnbad(2) = imnbad(3)
      if(imnbad(2).gt.0) goto 500
      goto 200
c
c read moon data into storage
      entry MNRED1(jd)
      nmnrec = 0
  400 mm     = 1
  500 do l   = mm, 3
         l2  = l*8
         l1  = l2 - 7
         imnbad(l) = 0
         Jdm(l)    = 0
         nmnrec    = nmnrec + 1
         read(Imn,err=550) Jdm(l),Fm(l),ivl,
     .                        (((Moon(i,j,k),i=1,ivl),j=1,Iparm),
     .                        k = l1, l2),
     .                        (Psidx(k),Epsdx(k),k = l1,l2),
     .                        ((Librat(k,i),i=1,3),k = l1,l2)
         Imvel(l) = ivl
         Jder(l)  = Jdm(l)
         goto 600
  550    imnbad(l) = 1
         read(Imn)
  600 end do
      Nvels(2) = 9999
      Nbtrp(2) = 0
c
c
c reconstructs dates of start of bad records
      if(imnbad(1).gt.0 .and. Jdm(1).le.0) then
         if(imnbad(2).le.0 .or. Jdm(2).gt.0) then
            Jdm(1) = Jdm(2) - ISIGN(intmn5,Idirmn)
            Jdm(3) = Jdm(2) + ISIGN(intmn5,Idirmn)
            goto 700
         else if(imnbad(3).le.0 .or. Jdm(3).gt.0) then
            Jdm(2) = Jdm(3) - ISIGN(intmn5,Idirmn)
            Jdm(1) = Jdm(2) - ISIGN(intmn5,Idirmn)
            goto 700
         else if(jd.le.0) then
            call SUICID(
     .            ' 3 BAD RECORDS ON MOON DATA SET, STOP IN MNREED 908 '
     .            , 13)
            goto 700
         else
            Jdm(1) = jdmnst + ISIGN(intmn5*(nmnrec-3),Idirmn)
         endif
      endif
      Jdm(2) = Jdm(1) + ISIGN(intmn5,Idirmn)
      Jdm(3) = Jdm(2) + ISIGN(intmn5,Idirmn)
  700 if(jd.gt.0) goto 100
      jdmnst = Jdm(1)
      return
c
c see if jd of observations is in bad records
  800 do i = 1, 3
         if(imnbad(i).gt.0) then
            write(Iout,820) jd,Imn,imnbad
  820       format(i17,
     .             '  DELETED BECAUSE OF BAD RECORD ON MOON DATA SET',
     .             i3, '  WITH IMNBAD=', 3I2, '  ********')
            jd = -1
            if(Ict(35).lt.0) then
               goto 900
            else if(Ict(35).eq.0) then
               if(Ncodf.le.3) goto 900
            else
               if(Nplnt0.le.30) return
               if(Ict(35).gt.1) return
               goto 900
            endif
         endif
      end do
 
      return
 
900   call SUICID('ERROR RECORD ON MOON DATA SET, STOP IN MNREED   ',12)
 
      end
