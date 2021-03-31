      subroutine SCREED(jd,fract)
 
      implicit none
 
 
c*** start of declarations inserted by spag
 
c*** end of declarations inserted by spag
 
 
c     m.e.ash    may 1970        subroutine screed
c        reference to s-body tape added 1977 jul - j.f.chandler
c     observing body  tape is read either forward or backward in time
c     constant tabular interval, five tabular points per record

c arguments
      integer*4 jd
      real*10 fract
c jd   = julian day no.
c fract= fraction of day from midnight
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'inodta.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'scdta.inc'
      include 'trpcom.inc'

c local variables
      real*10 xdist,xx
      integer*4 i,ivl,j,k,l,l1,l2,l3,mm,n

c test to see if observing body is on n-body tape
      if(Jsc.lt.0) then
         call B2REED(jd,fract,Nplnt(Klans1),2)
         return
      else if(Jsc.eq.0) then
         call BDREED(jd,Nplnt(Klans1),5)
         return
c
c test to see if jd is on tape
      else if(jd.lt.Jdsc1) then
      else if(jd.eq.Jdsc1) then
         if(fract.gt.Frsc1) goto 300
      else
         goto 300
      endif
  100 write(Iout,200) jd,fract
  200 format(i17,f14.13,' NOT ON OBSERVING BODY DATA SET')
      jd = 0
      return
  300 if(jd.lt.Jdsc2) then
      else if(jd.eq.Jdsc2) then
         if(fract.ge.Frsc2) goto 100
      else
         goto 100
      endif
c
c see if this is nordsieck variable tabular interval type
      if(Ksc(88).le.0) then
         call SCRAAD(jd,fract)
         return
      endif
c
c get correct records of satellite-probe tape into storage
c determine number of tabular points between jd, fract and
c jdsc(2), fsc(2)
  400 xx    = jd - Jdsc(2)
      xx    = xx + (fract - Fsc(2))
      xdist = xx/Scint
      if(xdist.lt.0) then
      else if(xdist.eq.0) then
c
c correct records are behind on tape
         if(xx*Scint.ge.0) return
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
            Jdsc(1)   = Jdsc(n + 1)
            Iscvel(1) = Iscvel(n + 1)
            ivl       = Iscvel(1)
            Fsc(1)    = Fsc(n + 1)
            do while( .true. )
               do k = l1, l2
                  l = k + l3
                  do j = 1, Iparsc
                     do i = 1, ivl
                        Satprc(i,j,k) = Satprc(i,j,l)
                     end do
                  end do
               end do
               if(n.eq.2) goto 600
               Jdsc(2)   = Jdsc(3)
               Fsc(2)    = Fsc(3)
               Iscvel(2) = Iscvel(3)
               ivl       = Iscvel(2)
               n  = 2
               l1 = 6
               l2 = 10
               mm = 3
            end do
         else if(mm.ne.0) then
            do i = 1, mm
               read(Jsc,err=410)
               goto 440
  410          read(Jsc)
 
c second read of error record might not be needed if system changes
               write(Iout,420) Jsc
  420          format(
     .         ' **** ERROR RECORD SKIPPED ON OBSERVING BODY DATA SET '
     .         , i3, ' IN SCREED 42 ****')
  440       end do
         endif
         goto 500
      endif
      n = xdist/5.0_10
      n = 4 - n
      do i = 1, n
         backspace Jsc
      end do
c
c read observing body data into storage
      entry SCRED1(jd)
  500 mm    = 1
  600 do l  = mm, 3
         l2 = l*5
         l1 = l2 - 4
         read(Jsc) Jdsc(l),Fsc(l),ivl,
     .             (((Satprc(i,j,k),i=1,ivl),j=1,Iparsc),k = l1,l2)
         Iscvel(l) = ivl
      end do
      if(jd.gt.0) goto 400
 
      return
      end
