      subroutine RTREED(jd,fract,nrss,nplntr)
 
      implicit none
c
c m.e.ash, r.w.king   september 1972   subroutine rtreed
c earth,moon or planet rotation tape is read either forward or
c backward in time
c constant tabular interval, five tabular points per record
c
c jd   = julian day no.
c fract= fraction of day from midnight
      integer*4 jd,nrss,nplntr
      real*10 fract

c array dimensions
      include 'globdefs.inc'

c common
      include 'inodta.inc'
      include 'rotdta.inc'
      include 'trpcom.inc'
 
c printx,erintx= tabular interval for moon/planet, earth
      real*10 xdist,xx
      integer*4 i,ivl,j,k,l,l1,l2,l3,mm,n
c
c determine rotation body tape to be read
      if(nplntr.eq.3) then
c
c-----------------------------------------------------------------------
c
c           tape to be read is jer for earth rotation
c
c           test to see if jd is on tape
         if(jd.lt.Jder1) goto 600
         if(jd.ne.Jder1) goto 800
         if(fract.gt.Frer1) goto 800
         goto 600
c
c-----------------------------------------------------------------------
c
c           tape to be read is jpr for moon or planet rotation
c
c           test to see if jd is on tape
      else if(jd.lt.Jdpr1) then
      else if(jd.eq.Jdpr1) then
         if(fract.gt.Frpr1) goto 300
      else
         goto 300
      endif
  100 write(Iout,200) jd,fract
  200 format(i17,f17.13,' NOT ON PLANET OR MOON ROTATION DATA SET')
      jd = 0
      return
  300 if(jd.lt.Jdpr2) then
      else if(jd.eq.Jdpr2) then
         if(fract.ge.Frpr2) goto 100
      else
         goto 100
      endif
c
c get correct records of tape into storage
c determine number of tabular points between
c jd, fract and jdpr(2), fpr(2)
  400 xx    = jd - Jdpr(2)
      xx    = xx + (fract - Fpr(2))
      xdist = xx/Print
      if(xdist.lt.0) then
      else if(xdist.eq.0) then
c
c correct records are behind on tape
         if(xx*Print.ge.0) return
      else
c
c correct records are ahead on tape
         n = xdist/5.0_10
         if(n.eq.0) return
         mm = n - 3
         if(mm.lt.0) then
c
c correct middle record is no more than two ahead on tape
c storage must be shifted
            l1 = 1
            l2 = 5
            l3 = 5*n
            mm = 2
            Jdpr(1)   = Jdpr(n + 1)
            Iprvel(1) = Iprvel(n + 1)
            ivl       = Iprvel(1)
            Fpr(1)    = Fpr(n + 1)
            do while( .true. )
               do k = l1,l2
                  l = k + l3
                  do j = 1,Iparpr
                     do i = 1,ivl
                        Plnmon(i,j,k) = Plnmon(i,j,l)
                     end do
                  end do
               end do
               if(n.eq.2) goto 500
               Jdpr(2)   = Jdpr(3)
               Fpr(2)    = Fpr(3)
               Iprvel(2) = Iprvel(3)
               ivl       = Iprvel(2)
               n  = 2
               l1 = 6
               l2 = 10
               mm = 3
            end do
         else if(mm.ne.0) then
            do i = 1,mm
               read(Jpr,err=410)
               goto 440
  410          read(Jpr)
 
c second read of error record might not be needed if system changes
               write(Iout,420) Jpr
  420          format(
     . ' **** ERROR RECORD SKIPPED ON PLANET OR MOON ROTATION DATA SET'
     . ,i3,' IN RTREED 42 ****')
  440       end do
         endif
         mm = 1
         goto 500
      endif
      n = xdist/5._10
      n = 4 - n
      do i = 1,n
         backspace Jpr
      end do
c
c read moon or planet rotation data into storage
      entry RTRED1(jd,nplntr)
      if(nplntr.eq.3) goto 1000
      mm    = 1
  500 do l  = mm,3
         l2 = l*5
         l1 = l2 - 4
         read(Jpr) Jdpr(l),Fpr(l),ivl,
     .             (((Plnmon(i,j,k),i=1,ivl),j=1,Iparpr),k = l1,l2)
         Iprvel(l) = ivl
      end do
      if(jd.le.0) return
      goto 400
  600 write(Iout,700) jd,fract
  700 format(i17,f17.13,' NOT ON EARTH ROTATION DATA SET')
      jd = 0
      return
  800 if(jd.lt.Jder2) then
      else if(jd.eq.Jder2) then
         if(fract.ge.Frer2) goto 600
      else
         goto 600
      endif
c
c get correct records of earth rotation tape into storage
c determine number of tabular points between jd, fract and
c jder(2),fer(2)
  900 xx    = jd - Jder(2)
      xx    = xx + (fract - Fer(2))
      xdist = xx/Erint
      if(xdist.lt.0) then
      else if(xdist.eq.0) then
c
c correct records are behind on tape
         if(xx*Erint.gt.0) return
      else
c
c correct records are ahead on tape
         n = xdist/5.0_10
         if(n.eq.0) return
         mm = n - 3
         if(mm.lt.0) then
c
c correct middle record is no more than two ahead on tape
c storage must be shifted
            l1 = 1
            l2 = 5
            l3 = 5*n
            mm = 2
            Jder(1)   = Jder(n + 1)
            Iervel(1) = Iervel(n + 1)
            ivl       = Iervel(1)
            Fer(1)    = Fer(n + 1)
            do while( .true. )
               do k = l1,l2
                  l = k + l3
                  do j = 1,Iparer
                     do i = 1,ivl
                        Earth(i,j,k) = Earth(i,j,l)
                     end do
                  end do
               end do
               if(n.eq.2) goto 1100
               Jder(2)   = Jder(3)
               Fer(2)    = Fer(3)
               Iervel(2) = Iervel(3)
               ivl       = Iervel(2)
               n  = 2
               l1 = 6
               l2 = 10
               mm = 3
            end do
         else if(mm.ne.0) then
            do i = 1,mm
               read(Jer,err=910)
               goto 940
  910          read(Jer)
 
c second read of error record might not be needed if system changes
               write(Iout,920) Jer
  920          format(
     .          ' **** ERROR RECORD SKIPPED ON EARTH ROTATION DATA SET'
     .          ,i3,' IN RTREED 142 ****')
  940       end do
         endif
         goto 1000
      endif
      n = xdist/5._10
      n = 4 - n
      do i = 1,n
         backspace Jer
      end do
c
c read earth rotation data into storage
 1000 mm    = 1
 1100 do l  = mm,3
         l2 = l*5
         l1 = l2 - 4
         read(Jer) Jder(l),Fer(l),ivl,
     .             (((Earth(i,j,k),i=1,ivl),j=1,Iparer),k = l1,l2)
         Iervel(l) = ivl
      end do
      if(jd.gt.0) goto 900
 
      return
      end
