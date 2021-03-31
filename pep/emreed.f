      subroutine EMREED(jd,fract)
 
      implicit none
 
c
c m.e.ash, a.d.rasinski   may 1966    subroutine emreed
c earth-moon barycenter tape is read either forward or backward in
c time
c
c arguments
      integer*4 jd
      real*10 fract

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdat.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'number.inc'
      include 'plndta.inc'
      include 'tapdta.inc'
      include 'tapdte.inc'

c local variables 
      real*10 xdist,femst
      integer*4 i,ivl,j,jdemst,k,l,l1,l2,mm,n,nemrec
      integer*2 np3/3/, iembad(3)
c external functions
      integer*4 JBDTST

c
c test to see if embary is on n-body tape
      if(JBDTST(np3).le.0) then
c
c test for day,time on tape
         if(jd.gt.Jdem1 .or. (jd.eq.Jdem1.and.fract.gt.Frem1)) then
            if(jd.lt.Jdem2 .or. (jd.eq.Jdem2.and.fract.lt.Frem2))
     .       goto 100
         endif
         write(Iout,50) jd,fract,Jdem1,Frem1,Jdem2,Frem2
   50    format(i17,f5.4,' NOT ON EARTH-MOON BARYCENTER DATA SET (',i7,
     .    f5.4,'-',i7,f5.4,')')
         jd = 0
         return
      else
         call BDREED(jd,np3,1)
         return
      endif
  100 xdist = Jdem(2)-jd
      xdist = xdist + (Fem(2)-fract)
      xdist = xdist/Emint5
c
c determine if data must be read in
      if(xdist.gt.0._10) goto 200
      if(xdist.gt.-1._10) goto 700
      n=-xdist
      if(n.le.2) then
c correct records are immediately ahead. shift storage and replenish
         l2=3-n
         l1=5*n
         do l=1,l2
            Jdem(l)   = Jdem(l+n)
            Ievel(l)  = Ievel(l+n)
            Fem(l)    = Fem(l+n)
            iembad(l) = iembad(l+n)
            if(iembad(l).le.0) then
               ivl = Ievel(l)
               mm  = 5*l
               do k=mm-4,mm
                  do j = 1, Iparem
                     do i = 1, ivl
                        Embary(i,j,k) = Embary(i,j,k+l1)
                     end do
                  end do
               end do
            endif
         end do
         mm=l2+1
         goto 400
      endif
c desired records are n ahead on tape. skip if necessary
      if(n.eq.3) goto 300
      n=n-3
      do i = 1, n
         nemrec = nemrec + 1
         read(Iem,err=110)
         goto 130
  110    read(Iem)
 
c second read of error record might not be needed if system changes
         write(Iout,120) Iem
  120    format(' **** ERROR RECORD SKIPPED ON EMBARY DATA SET',i3,
     .    ' IN EMREED ****')
  130 end do
      goto 300

c correct records are behind on tape. backspace and fill again
  200 n=xdist
      if(n.eq.xdist) n=n-1
      n=n+4
      do i = 1, n
         nemrec = nemrec - 1
         backspace Iem
      end do
      goto 300
c
c read embary data into storage
      entry EMRED1(jd)
      nemrec = 0
  300 mm     = 1
  400 do l = mm, 3
         l2 = l*5
         l1 = l2 - 4
         iembad(l) = 0
         Ievel(l)  = 0
         Jdem(l)   = 0
         nemrec    = nemrec + 1
         read(Iem,err=450) Jdem(l),Fem(l),ivl,
     .                        (((Embary(i,j,k),i=1,ivl),j=1,Iparem),
     .                        k = l1, l2)
         Ievel(l) = ivl
         goto 500
  450    iembad(l) = 1
         read(Iem)
  500 end do
c
c reconstruct dates of start of bad records
      if(iembad(1).gt.0 .and. Jdem(1).le.0) then
         if(iembad(2).le.0 .or. Jdem(2).gt.0) then
            Fem(1) =Jdem(2)+Fem(2)-Emint5
            Jdem(1)=Fem(1)
            Fem(1) =Fem(1)-Jdem(1)
            Fem(3) =Jdem(2)+Fem(2)+Emint5
            Jdem(3)=Fem(3)
            Fem(3) =Fem(3)-Jdem(3)
            goto 600
         else if(iembad(3).le.0 .or. Jdem(3).gt.0) then
            Fem(2) =Jdem(3)+Fem(3)-Emint5
            Jdem(2)=Fem(2)
            Fem(2) =Fem(2)-Jdem(2)
            Fem(1) =Jdem(2)+Fem(2)-Emint5
            Jdem(1)=Fem(1)
            Fem(1) =Fem(1)-Jdem(1)
            goto 600
         else if(jd.le.0) then
            call SUICID(
     .        ' 3 BAD RECORDS ON EMBARY DATA SET, STOP IN EMREED   '
     .        , 13)
            goto 600
         else
            Fem(1) = jdemst+femst + Emint5*(nemrec-3)
            Jdem(1)=Fem(1)
            Fem(1) =Fem(1)-Jdem(1)
         endif
      endif
      Fem(2) =Jdem(1)+Fem(1)+Emint5
      Jdem(2)=Fem(2)
      Fem(2) =Fem(2)-Jdem(2)
      Fem(3) =Jdem(2)+Fem(2)+Emint5
      Jdem(3)=Fem(3)
      Fem(3) =Fem(3)-Jdem(3)
  600 if(jd.gt.0) goto 100
      jdemst = Jdem(1)
      femst  = Fem(1)
      return
c
c see if jd of observations is in bad record
  700 do i = 1, 3
         if(iembad(i).gt.0) then
            write(Iout,720) jd,Iem,iembad
  720       format(i17,
     .             ' DELETED BECAUSE OF BAD RECORD ON EMBARY DATA SET',
     .             i3, '  WITH IEMBAD=', 3I2, '  ********')
            jd = -1
            if(Ict(35).lt.0) then
               call SUICID(
     .            ' ERROR RECORD ON EMBARY DATA SET, STOP IN EMREED'
     .            , 12)
            else if(Ict(35).eq.0) then
               goto 800
            else
               if(Nplnt0.le.30) return
               if(Ict(35).gt.1) return
               call SUICID(
     .            ' ERROR RECORD ON EMBARY DATA SET, STOP IN EMREED'
     .            , 12)
            endif
         endif
      end do
      return
  800 if(Ncodf.le.3) then
         call SUICID(
     .            ' ERROR RECORD ON EMBARY DATA SET, STOP IN EMREED'
     .            , 12)
      endif
 
      return
      end
