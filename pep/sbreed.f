      subroutine SBREED(jd,fract)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 epoch, fract, xdist, xx
      integer   i, ivl, j, jd, k, l, l1, l2, l3, mm, n
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   november 1967    subroutine sbreed
c reference to s-body tape added 1977 jul - j.f.chandler
c satellite-probe tape is read either forward or backward in time
c constant tabular interval, five tabular points per record
c
 
c jd   = julian day no.
c fract= fraction of day from midnight

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'number.inc'
      include 'sbdta.inc'
      include 'trpcom.inc'
c
c test to see if satellite or probe is on n-body tape
      if(Jsb.lt.0) then
         call B2REED(jd,fract,Nplnt0,1)
         goto 900
      else if(Jsb.eq.0) then
         call BDREED(jd,Nplnt0,4)
         goto 900
c
c test to see if jd is on tape
      else if(jd.lt.Jdsb1) then
      else if(jd.eq.Jdsb1) then
         if(fract.gt.Frsb1) goto 300
      else
         goto 300
      endif
  100 write(Iout,200) jd,fract
  200 format(i17,f14.13,' NOT ON SATELLITE-PROBE DATA SET')
      jd = 0
      goto 900
  300 if(jd.lt.Jdsb2) then
      else if(jd.eq.Jdsb2) then
         if(fract.ge.Frsb2) goto 100
      else
         goto 100
      endif
c
c see if this is nordsieck variable tabular interval type
      if(Ksb(88).le.0) then
         call SBRAAD(jd,fract,*700)
         goto 900
      endif
c
c get correct records of satellite-probe tape into storage
c determine number of tabular points between jd, fract and
c jdsb(2), fsb(2)
  400 xx    = jd - Jdsb(2)
      xx    = xx + (fract - Fsb(2))
      xdist = xx/Sbint
      if(xdist.lt.0) then
      else if(xdist.eq.0) then
c
c correct records are behind on tape
         if(xx*Sbint.ge.0) goto 900
      else
c
c correct records are ahead on tape
         n = xdist/5.0_10
         if(n.eq.0) goto 900
         mm = n - 3
         if(mm.lt.0) then
c
c correct middle record is no more than two ahead on tape
c storage must be shifted
            l1 = 1
            l2 = 5
            l3 = 5*n
            mm = 2
            Jdsb(1)   = Jdsb(n + 1)
            Isbvel(1) = Isbvel(n + 1)
            ivl       = Isbvel(1)
            Fsb(1)    = Fsb(n + 1)
            do while( .true. )
               do k = l1, l2
                  l = k + l3
                  do j = 1, Iparsb
                     do i = 1, ivl
                        Satprb(i,j,k) = Satprb(i,j,l)
                     end do
                  end do
               end do
               if(n.eq.2) goto 600
               Jdsb(2)   = Jdsb(3)
               Fsb(2)    = Fsb(3)
               Isbvel(2) = Isbvel(3)
               ivl       = Isbvel(2)
               n  = 2
               l1 = 6
               l2 = 10
               mm = 3
            end do
         else if(mm.ne.0) then
            do i = 1, mm
               read(Jsb,err=410)
               goto 440
  410          read(Jsb)
 
c second read of error record might not be needed if system changes
               write(Iout,420) Jsb
  420          format(
     .         ' **** ERROR RECORD SKIPPED ON SATELLITE-PROBE DATA SET'
     .         , i3, ' IN SBREED 42 ****')
  440       end do
         endif
         goto 500
      endif
      n = xdist/5.0_10
      n = 4 - n
      do i = 1, n
         backspace Jsb
      end do
c
c read satellite-probe data into storage
      entry SBRED1(jd)
  500 mm    = 1
  600 do l  = mm, 3
         l2 = l*5
         l1 = l2 - 4
         read(Jsb) Jdsb(l),Fsb(l),ivl,
     .             (((Satprb(i,j,k),i=1,ivl),j=1,Iparsb),k = l1,l2)
         Isbvel(l) = ivl
         if(mod(Jct(6)/1024,2).eq.1) then
            if(Line.gt.56) call OBSPAG
            write(Iout,610) l,Jdsb(l),Fsb(l)
  610       format(' SBREED: replacing storage record',i2,' with JD.F=',
     .       i10,f15.14)
            Line=Line+1
         endif
      end do
      if(jd.le.0) goto 900
      goto 400
 
c error exit for end-of-file on jsb in sbraad
  700 write(Iout,800) jd,fract
  800 format(i17,f17.13,
     .       ' TOO NEAR THE END OF SATELLITE-PROBE DATA SET')
      jd = 0
 
  900 return
      entry SBRED2(epoch)
      if(Ksb(88).le.0) call SBRAD2(epoch)
      return
      end
