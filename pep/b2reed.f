      subroutine B2REED(jd,fract,mplnt,kcall)
 
      implicit none

c j.f.chandler  1977 jun   subroutine b2reed
c read s-body tape either forwards or backwards in time
c
c arguments
      integer*4 jd,kcall
      real*10 fract
      integer*2 mplnt
c jd,fract  epoch to be in middle of second record
c mplnt     planet number of body of interest
c supply data to calling routine (indicated by kcall)
c kcall = 1 sbreed
c kcall = 2 screed
c kcall = 3 szreed
c kcall = 4 plreed
 

c array dimensions
      include 'globdefs.inc'
c common 
      include 'b2dta.inc'
      include 'b2ydta.inc'
      include 'comdat.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      integer*4 i2bod
      equivalence (i2bod,Jpert)
      include 'sbdta.inc'
      include 'scdta.inc'
      include 'tapdta.inc'
      include 'tapdtp.inc'
 
c local variables
      real*10 dim,dim5,dt,fk,tb,tb21
      integer   i,iv1,iv4,iv8,ivl,j,j1,jj,k,kb,kc,
     .          l,l1,l2,l3,l4,l5,l6,m,m1,m4,m8,mm,n,nb2rec
      integer*2 ib2bad(3)
c
c test to see if jd is on tape
      if(jd.ge.Jdb21 .and. jd.le.Jdb22) then
         tb = jd + fract
         goto 300
      endif
  100 write(Iout,200) jd,Jdb21,Jdb22
  200 format(i17,' NOT ON S-BODY DATA SET (',i7,'-',i7,')')
      jd = 0
      return
 
c get correct records of s-body tape into storage
  300 fk = ((jd-Jdb2(2)) + fract - Frb2(2))/B2int5
      if(fk.lt.0._10) then
c
c correct record is behind on tape
         n = fk
         n = 4 - n
         if(n.gt.nb2rec) goto 100
         do i = 1,n
            nb2rec = nb2rec - 1
            backspace i2bod
         end do
         goto 600
      endif
 
      n = fk
      if(n.eq.0) goto 550
      m = n - 3
      if(m.gt.0) then
 
c correct records are ahead on tape
         do i = 1,m
            nb2rec = nb2rec + 1
            read(i2bod,err=350)
            goto 360
  350       read(i2bod)
  360    end do
         goto 600
      endif
 
      if(m.eq.0) goto 600
c correct record is no more than two ahead of present middle
c shift storage and then read the balance
      mm = 4 - n
      m1 = 40*n
      m4 = 10*n
      m8 = 5*n
      Jdb2(1)   = Jdb2(n+1)
      Frb2(1)   = Frb2(n+1)
      Ivl1(1)   = Ivl1(n+1)
      Ivl4(1)   = Ivl4(n+1)
      Ivl8(1)   = Ivl8(n+1)
      ib2bad(1) = ib2bad(n+1)
      if(ib2bad(1).gt.0) goto 500
      l2 = 40
      l4 = 10
      l6 = 5
 
  400 l1 = l2 - 39
      l3 = l4 - 9
      l5 = l6 - 4
      if(Na1.gt.0) then
         ivl = Ivl1(n)
         do i = 1,ivl
            do j = l1,l2
               jj = j + m1
               do k = 1,Na1
                  Bod1(i,j,k) = Bod1(i,jj,k)
               end do
            end do
         end do
      endif
      if(Na4.gt.0) then
         ivl = Ivl4(n)
         do i = 1,ivl
            do j = l3,l4
               jj = j + m4
               do k = 1,Na4
                  Bod4(i,j,k) = Bod4(i,jj,k)
               end do
            end do
         end do
      endif
      if(Na8.gt.0) then
         ivl = Ivl8(n)
         do i = 1,ivl
            do j = l5,l6
               jj = j + m8
               do k = 1,Na8
                  Bod8(i,j,k) = Bod8(i,jj,k)
               end do
            end do
         end do
      endif
  500 if(n.eq.2) goto 700
 
c shift other record, too
      Jdb2(2)   = Jdb2(3)
      Frb2(2)   = Frb2(3)
      Ivl1(2)   = Ivl1(3)
      Ivl4(2)   = Ivl4(3)
      Ivl8(2)   = Ivl8(3)
      ib2bad(2) = ib2bad(3)
      if(ib2bad(2).gt.0) goto 700
      l2 = 80
      l4 = 20
      l6 = 10
      n  = 2
      goto 400
 
c
c test if any bad records now
  550 do i = 1,3
         if(ib2bad(i).gt.0) then
            write(Iout,560) jd,i2bod,ib2bad
  560       format(i17,' DELETED - BAD RECORD ON S-BODY DATA SET',
     .              i3, ' WITH BAD=',3I2,' * * * *')
            jd = -1
 
c test to see if run should be aborted
            if(Ict(36).lt.0) then
            else if(Ict(36).eq.0) then
               if(Ncodf.gt.3) return
            else
               if(Ict(36).gt.1 .or. mplnt.le.30) return
            endif
            call SUICID('ERRORS ON S-BODY TAPE, STOP IN B2REED   ',10)
 
         endif
      end do
 
c find planet number on tape
      if(mplnt.eq.0) return
      do i = 1,Nast
         if(mplnt.eq.Np2(i)) then
            kb = i - 1
            goto 330
         endif
      end do
      write(Iout,320) mplnt,i2bod
  320 format('0* * * PLANET NUMBER',i3,
     .       ' NOT FOUND ON S-BODY DATA SET',i3,' * * * *')
      call SUICID('BODY MISSING FROM S-BODY TAPE   ',8)
 
c determine which class the body of interest belongs to
  330 if(kb.ge.Na1) then
         kc = kb - Na1
         if(kc.ge.Na4) then
            kc  = 480 + 120 + (kc - Na4)*15
            dim = D(3)
            l   = 7
         else
            kc  = 480 + kc*30
            dim = D(2)
            l   = 4
         endif
      else
         kc  = kb*120
         dim = D(1)
         l   = 1
      endif
 
c test if date is in storage already
      dt = (tb - Tgo(kcall))/dim
      if(dt.ge.0._10 .and. dt.lt.5._10) return
      ivl  = max0(Ivl1(l),Ivl1(l+1),Ivl1(l+2))
      dim5 = dim*5._10
      tb21 = Jdb2(1) + Frb2(1)
      j1   = (tb - tb21)/dim
      if(kcall.ne.5) j1 = j1/5*5
      tb = tb21 + dim*j1
 
c now tb is epoch of desired middle record
      Ttb2(1)    = tb - dim5
      Ttb2(2)    = tb
      Ttb2(3)    = tb + dim5
      Tgo(kcall) = tb
      j1    = j1 - 4 + kc
      Idxb2 = j1
 
c
c now move data to proper storage
 
      if(kcall.eq.1) then
 
c satellite observed body
         do i = 1,3
            Isbvel(i) = ivl
            Jdsb(i)   = Ttb2(i)
            Fsb(i)    = Ttb2(i) - Jdsb(i)
         end do
         do j = 1,15
            do i = 1,ivl
               Satprb(i,1,j) = Bod1(i,j1,1)
            end do
            j1 = j1 + 1
         end do
 
      else if(kcall.eq.2) then
 
c satellite observing body
         do i = 1,3
            Iscvel(i) = ivl
            Jdsc(i)   = Ttb2(i)
            Fsc(i)    = Ttb2(i) - Jdsc(i)
         end do
         do j = 1,15
            do i = 1,ivl
               Satprc(i,1,j) = Bod1(i,j1,1)
               end do
            j1 = j1 + 1
         end do
 
      else if(kcall.eq.4) then
 
c minor planet
         do i = 1,3
            Ipvel(i) = ivl
            Jdp(i)   = Ttb2(i)
            Fp(i)    = Ttb2(i) - Jdp(i)
         end do
         do j = 1,15
            do i = 1,ivl
               Planet(i,1,j) = Bod1(i,j1,1)
               end do
            j1 = j1 + 1
         end do
 
      else if(kcall.eq.5) then
         return
      else
         call SUICID('ILLEGAL CALL TO B2REED  ',6)
      endif
 
      return
 
c
c enter here to read first three data records
c
      entry B2RED1(jd)
 
      nb2rec = 0
  600 mm     = 1
  700 do l = mm,3
         l2 = l*40
         l1 = l2 - 39
         l4 = l*10
         l3 = l4 - 9
         l6 = l*5
         l5 = l6 - 4
         if(Na1.le.0) l2 = l1
         if(Na4.le.0) l4 = l3
         if(Na8.le.0) l6 = l5
         ib2bad(l) = 0
         Jdb2(l)   = 0
         Frb2(l)   = 0._10
         nb2rec    = nb2rec + 1
         read(i2bod,err=750) Jdb2(l),Frb2(l),iv1,iv4,iv8,
     .                          (((Bod1(i,j,k),i=1,iv1),j=l1,l2),
     .                          k = 1,max0(1,0+Na1)),
     .                          (((Bod4(i,j,k),i=1,iv4),j=l3,l4),
     .                          k = 1,max0(1,0+Na4)),
     .                          (((Bod8(i,j,k),i=1,iv8),j=l5,l6),
     .                          k = 1,max0(1,0+Na8))
         Ivl1(l) = iv1
         Ivl4(l) = iv4
         Ivl8(l) = iv8
         goto 800
  750    ib2bad(l) = 1
         read(i2bod)
  800 end do
 
      if(jd.gt.0) goto 300
 
      return
      end
