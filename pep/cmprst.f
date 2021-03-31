      subroutine CMPRST
 
      implicit none
c
c m.e.ash    march 1969    subroutine cmprst
c checkpoint restart for compar link (position observation
c tapes)
c

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdat.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obsdta.inc'
      include 'zeroes.inc'

c external functions
      integer*4 LEG

c local
      integer   i,ict33,mtabs1,mtabs2,mtobs,nrwb,nseqb,
     .          nseqb1,nseqb2,ntabs1,ntabs2,ntapb,ntobs
      integer*2 ncod
      character*1 htype,colon/':'/

c     ict(31) = 0            no checkpoint restart in compar link
c     ict(31) = negative,    checkpoint restart in compar at start of
c                            observation library tape
c     ict(31) = positive,    checkpoint restart in compar in middle of
c                            observation library tape
c     in either of the above two cases, jtape counter in compar is set
c     to iabs(ict(31))
c     ict(32) = value of ntape for compar checkpoint restart
c     ict(33) = value of nseq for compar checkpoint restart
c               if ict(33).gt.0
c     ict(33) = value of nseq for series before compar checkpoint
c               restart if ict(33).lt.0
c     iobs = observation card data set (sysin if 5, in which case it is
c            assumed that the person restarting compar has the obser-
c            vation cards beginning at the restart point)
c     iabs1= input observation library tape
c     iabs2= output observation library tape
c
      ict33 = Ict(33)
      if(Ict(33).lt.0) ict33 = iabs(ict33) + 1
      ntabs1 = 0
      if(Iabs1.le.0) ntabs1   = 1
      if(Ict(31).lt.0) ntabs1 = 1
      ntabs2 = 0
      if(Iabs2.le.0) ntabs2   = 1
      if(Ict(31).lt.0) ntabs2 = 1
      ntobs = 0
      if(Iobs.le.0) ntobs = 1
      if(Iobs.eq.5) ntobs = 1
c
c skip first 2 records of output tape
      if(ntabs2.gt.0) goto 300
      call INOUTV(Iabs2)
      read(Iabs2,err=100)
      goto 200
  100 read(Iabs2)
  200 read(Iabs2)
c
c skip first record of observation series
  300 if(ntabs1.le.0) then
         read(Iabs1) nseqb1
         if(nseqb1.ge.ict33) then
            ntabs1 = 1
            backspace Iabs1
         endif
      endif
      if(ntabs2.le.0) then
         read(Iabs2) nseqb2
         if(nseqb2.ge.ict33) then
            ntabs2 = 1
            backspace Iabs2
         endif
      endif
      if(ntobs.le.0) then
         if(Jct(69).gt.0) then
            read(Iobs,320) ncod,htype,nrwb,ntapb,nseqb
  320       format(i2,3x,a1,57x,i2,i3,i4)
         else
            read(Iobs,340) ncod,nrwb,ntapb,nseqb
  340       format(i2,68x,i2,i3,i5)
         endif
         if(ntapb.ge.Ict(32) .and. nseqb.ge.ict33) then
            ntobs = 1
            backspace Iobs
         else
            if(Jct(69).gt.0 .and. LEG(1,1,htype,1,colon).eq.0)
     .          read (Iobs,340)
            if(ncod.gt.20) read(Iobs,340)
            if(nrwb.lt.0 .or. nrwb.gt.1) read(Iobs,340)
         endif
      endif
      if(ntabs1.gt.0 .and. ntabs2.gt.0 .and. ntobs.gt.0) then
c
c write message
         call PAGCHK(60,3,0)
         write(Iout,350) (Ict(i),i = 31,33),Iobs,Iabs1,Iabs2
  350    format(/
     .' TAPES POSTIONED FOR CHECKPOINT RESTART IN COMPAR LINK   JTAPE=IC
     .T(31)=',i3,'   NTAPE=ICT(32)=',i4,'   NSEQ=ICT(33)=',
     .i6/' IOBS=',i3,'  IABS1=',i3,'  IABS2=',i3)
 
         return
      else
c
c skip records of observation series
         mtabs1 = ntabs1
         mtabs2 = ntabs2
         mtobs  = ntobs
      endif
  400 if(mtabs1.le.0) then
         read(Iabs1,err=500) ncod
         if(ncod.le.0) mtabs1 = 1
      endif
      goto 600
  500 read(Iabs1)
  600 if(mtabs2.le.0) then
         read(Iabs2,err=900,end=700) ncod
         if(ncod.le.0) mtabs2 = 1
      endif
      goto 1000
  700 ntabs2 = 1
      backspace Iabs2
      write(Iabs2) izero,izero
      call PAGCHK(60,2,0)
      write(Iout,800) nseqb2,Iabs2
  800 format(/' END ACCEPTED IN SERIES NSEQ=',i6,' ON IABS2=',i3,
     .       ' ZERO RECORD WRITTEN FOR RESTART ON NEW SERIES')
      mtabs2 = 1
      goto 1000
  900 read(Iabs2)
 1000 if(mtobs.le.0) then
         read(Iobs,1050,err=1100) ncod
 1050    format(i1,79x)
         if(ncod.le.0) mtobs = 1
      endif
      goto 1200
 1100 read(Iobs,1050)
c
c test for end of skipping records of observation series
 1200 if(mtabs1.le.0) goto 400
      if(mtabs2.le.0) goto 600
      if(mtobs.le.0) goto 1000
c
c have we arrived at series just before restart
      if(Ict(33).le.0) then
         if(ntabs1.le.0) then
            if(nseqb1 + Ict(33).eq.0) ntabs1 = 1
         endif
         if(ntabs2.le.0) then
            if(nseqb2 + Ict(33).eq.0) ntabs2 = 1
         endif
         if(ntobs.le.0) then
            if(nseqb + Ict(33).eq.0) ntobs = 1
         endif
      endif
      goto 300
      end
