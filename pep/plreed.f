      subroutine PLREED(jd)
 
      implicit none

c
c m.e.ash, a.d.rasinski   may 1966   subroutine plreed
c reference to s-body tape added 1977 jul - j.f.chandler
c planet tape is read either forward or backward in time

c parameters
      integer*4 jd
 
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdat.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'tapdta.inc'
      include 'tapdtp.inc'
      include 'trpcom.inc'

c local variables
      integer   i, ivl, j, jdpst, k, l, l1, l2, mm, n, nplrec
      integer*2 iplbad(3)
c
c test to see if planet is on n-body tape
      if(Jplnt.lt.0) then
         call B2REED(jd,0._10,Nplnt(Klan),4)
         return
      else if(Jplnt.eq.0) then
         call BDREED(jd,Nplnt(Klan),3)
         return
      else
 
         if(jd.le.Jdp1) goto 900
         if(Jdp2.le.jd) goto 900
      endif
  100 n = Jdp(2) - jd
      if(Idirpl.lt.0) then
 
         if(n.lt.0) then
c
c backspace logic
            n = -n
            goto 200
         else if(n.eq.0) then
            goto 700
         else
            n = Jdp(3) - jd
         endif
 
      else if(n.lt.0) then
         n = jd - Jdp(3)
      else if(n.eq.0) then
         goto 700
      else
         goto 200
      endif
      if(n.lt.0) goto 700
      if(n.ne.0) then
 
         n = n/Intp5 - 1
         if(n.lt.0) then
         else if(n.eq.0) then
 
            Jdp(1)    = Jdp(3)
            Ipvel(1)  = Ipvel(3)
            Fp(1)     = Fp(3)
            iplbad(1) = iplbad(3)
            if(iplbad(1).le.0) then
               ivl = Ipvel(1)
               do k = 1, 5
                  n = k + 10
                  do j = 1, Iparp
                     do i = 1, ivl
                        Planet(i,j,k) = Planet(i,j,n)
                     end do
                  end do
               end do
            endif
            mm = 2
            goto 400
         else
            n = n - 1
            if(n.gt.0) then
               do i = 1, n
                  read(Jplnt,err=110)
                  nplrec = nplrec + 1
                  goto 130
  110             read(Jplnt)
 
c second read of error record might not be needed if system changes
                  write(Iout,120) Jplnt
  120             format(
     .                  ' **** ERROR RECORD SKIPPED ON PLANET DATA SET'
     .                  , i3, ' IN PLREED 31 ****')
  130          end do
            endif
            goto 300
         endif
      endif
 
      Jdp(1)    = Jdp(2)
      Ipvel(1)  = Ipvel(2)
      Fp(1)     = Fp(2)
      iplbad(1) = iplbad(2)
      Jdp(2)    = Jdp(3)
      Ipvel(2)  = Ipvel(3)
      Fp(2)     = Fp(3)
      iplbad(2) = iplbad(3)
      if(iplbad(1).le.0 .or. iplbad(2).le.0) then
         ivl = max0(Ipvel(1),Ipvel(2))
         do k = 1, 10
            n = k + 5
            do j = 1, Iparp
               do i = 1, ivl
                  Planet(i,j,k) = Planet(i,j,n)
               end do
            end do
         end do
      endif
      mm = 3
      goto 400
  200 n = (n - 1)/Intp5 + 4
      do i = 1, n
         nplrec = nplrec - 1
         backspace Jplnt
      end do
      goto 300
c
c
c read planet data into storage
      entry PLRED1(jd)
      nplrec = 0
  300 mm     = 1
  400 do l   = mm, 3
         l2  = l*5
         l1  = l2 - 4
         iplbad(l) = 0
         Ipvel(l)  = 0
         Jdp(l)    = 0
         nplrec    = nplrec + 1
         read(Jplnt,err=450) Jdp(l),Fp(l),ivl,
     .                          (((Planet(i,j,k),i=1,ivl),j=1,Iparp),
     .                          k = l1, l2)
         Ipvel(l) = ivl
         goto 500
  450    iplbad(l) = 1
         read(Jplnt)
  500 end do
c
c reconstructs dates of start of bad records
      if(iplbad(1).gt.0 .and. Jdp(1).le.0) then
         if(iplbad(2).le.0 .or. Jdp(2).gt.0) then
            Jdp(1) = Jdp(2) - ISIGN(Intp5,Idirpl)
            Jdp(3) = Jdp(2) + ISIGN(Intp5,Idirpl)
            goto 600
         else if(iplbad(3).le.0 .or. Jdp(3).gt.0) then
            Jdp(2) = Jdp(3) - ISIGN(Intp5,Idirpl)
            Jdp(1) = Jdp(2) - ISIGN(Intp5,Idirpl)
            goto 600
         else if(jd.le.0) then
            call SUICID(
     .        ' 3 BAD RECORDS ON PLANET DATA SET, STOP IN PLREED 908   '
     .        , 14)
            goto 600
         else
            Jdp(1) = jdpst + ISIGN(Intp5*(nplrec-3),Idirpl)
         endif
      endif
      Jdp(2) = Jdp(1) + ISIGN(Intp5,Idirpl)
      Jdp(3) = Jdp(2) + ISIGN(Intp5,Idirpl)
  600 if(jd.gt.0) goto 100
      jdpst = Jdp(1)
      return
c
c see if jd of observations is in bad record
  700 do i = 1, 3
         if(iplbad(i).gt.0) then
            write(Iout,720) jd,Jplnt,iplbad
  720       format(i17,
     .             ' DELETED BECAUSE OF BAD RECORD ON PLANET DATA SET',
     .             i3, '  WITH IPLBAD=', 3I2, '  ********')
            jd = -1
            if(Ict(35).lt.0) then
               goto 800
            else if(Ict(35).eq.0) then
               if(Ncodf.le.3)  goto 800
            else
               if(Nplnt0.le.30) return
               if(Ict(35).gt.1) return
               goto 800
            endif
         endif
      end do
 
      return
 
  800 call SUICID(' ERROR RECORD ON PLANET DATA SET, STOP IN PLREED',12)
 
  900 write(Iout,1000) jd,Jdp1,Jdp2
 1000 format(i17,' IS NOT ON PLANET DATA SET (', i7, '-', i7, ')')
      jd = 0
 
      end
