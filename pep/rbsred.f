      subroutine RBSRED(in0,nstop,init)
 
      implicit none

c
c k.m.becker     september 1967    subroutine rbsred
c radar observation biases are initialized and read
c
c parameters
      integer*4 in0,nstop
      logical*4 init
c
c        modified for *command july 1978  r.b. goldstein
c
c
c         setup and read quantities into temporary storage, then
c     rearrange into labelled common,such that corrections to a given
c     planet are grouped together.  the order of the planets is as
c     follows
c
c              10 = moon
c               0 = sun
c        nplnt(i) = planet(aplnt(i)), i=1,Numpln
c
c     with other (positive) planet numbers following in increasing
c     numerical order, and finally,
c
c              -4 = star
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'inodta.inc'
      include 'namtim.inc'
      include 'rdbias.inc'
c
c for j=1,...,numrbs we have
c rdbsit(1,j)=first four characters of 8 character receiving site name
c rdbsit(2,j)=first four characters of 8 character sending site name
c rdbser(j)  =four character observation series name
c nplrbs(j)  =planet number of observed body
c rbias(1,j) =constant bias in time delay (sec)
c rbias(2,j) =constant bias in doppler shift(cyc/sec)
c lrbs(i,j)=0 rbias(i,j) not adjusted, i=1,2
c lrbs(i,j)=1 rbias(i,j) adjusted, i=1,2
c
c                 vi. radar biases               (rbsred)
c  for each radar site and observation series
c  columns
c    1- 4  receiving site name first 4 characters            (1a4)
c    6- 9  sending   site name first 4 characters            (1a4)
c   11-14  observation series name                           (1a4)
c   15-17  planet number                                     (i3)
c   22-23  lrbs=1,0 time delay bias adjusted or not          (i2)
c   24-25  lrbs=1,0 doppler    bias adjusted or not          (i2)
c   26-27  flag, if '##' then use new format, otherwise not
c   26-37  radar bias in time delay                          (e12.5)
c            or 28-41   (e14.7) ##
c   38-49  radar bias in doppler                             (e12.5)
c            or 42-53   (e14.7) ##
c  limitation  no more than u_mxrbs radar sites with biases allowed
c
c           shared work space in input link
      common /WRKCOM/ Rdbstt, Rdbsrt, Rbiast, Lrbst, Nplrbt
      character*4 Rdbstt(2,u_mxrbs),Rdbsrt(u_mxrbs)
      real*4 Rbiast(2,u_mxrbs),rbs(2)
      integer*2 Lrbst(2,u_mxrbs),Nplrbt(u_mxrbs),lrb(2),npl
      integer*2 k, m, dumi2
      equivalence (Numrbs,k)
      character*4 arcb/'AREC'/
      character*4 srarec/'6764'/
      character*4 blank/'    '/, rec, sen, rser
      integer*4 nnrbs,nplcop,i,j,l,n,ns
      character*80 rcard
c
c initialize site data
      do j = 1, u_mxrbs
         Rdbser(j) = blank
         Rdbsrt(j) = blank
         Nplrbs(j) = 0
         Nplrbt(j) = 0
         do i = 1, 2
            Rdbsit(i,j) = blank
            Rdbstt(i,j) = blank
            Rbias(i,j)  = 0.0
            Rbiast(i,j) = 0.0
            Lrbs(i,j)   = 0
            Lrbst(i,j)  = 0
         end do
      end do
      Numrbs = 0
c
c set up standard receiving site names, sending site names and
c observation series names
      m = 0
c
c set up standard names
      do l = 1, 2
         m = m + 1
         do i = 1, 2
            Rdbstt(i,m) = arcb
         end do
         Rdbsrt(m) = srarec
         Nplrbt(m) = l
      end do
 
      Numrbs = m
      if(.not. (init)) then
c
c spool radar bias cards from in to in0 with a-format printout
         call PEPTIC(In,Iout,in0,4,'RADAR BIAS CARDS', nstop, 1)
c
c read radar bias cards
         k = m
         do while( .true. )
            read(in0,15) rcard
   15       format(a)
            if(rcard(26:27).eq.'##') then
               read(rcard,18) rec,sen,rser,npl,lrb,rbs
   18          format(a4,1x,a4,1x,a4,i3,4x,2I2,2x,2E14.7)
            else
               read(rcard,20) rec,sen,rser,npl,lrb,rbs
   20          format(a4,1x,a4,1x,a4,i3,4x,2I2,2E12.5)
            endif
            if(rec.eq.blank) goto 100
c
c search to see if standard site and series
            if(k.gt.0) then
 
              do i = 1, k
                if(rec.eq.Rdbstt(1,i) .and. sen.eq.Rdbstt(2,i) .and.
     .              rser.eq.Rdbsrt(i) .and. npl.eq.Nplrbt(i)) then
                      ns = i
                      goto 40
                   endif
                end do
              endif
            k = k + 1
            if(k.gt.u_mxrbs) then
               write(Iout,30) rec,sen,rser,k
   30          format('0ERROR IN RBSRED', 10x, 3A6, '  K=', i4)
               call SUICID('TOO MANY RADAR BIASES INPUT ',7)
            endif
            Rdbstt(1,k) = rec
            Rdbstt(2,k) = sen
            Rdbsrt(k)    = rser
            Nplrbt(k)    = npl
            ns   = k
   40       do i = 1, 2
               Lrbst(i,ns)  = lrb(i)
               Rbiast(i,ns) = rbs(i)
            end do
         end do
      endif
c
c rearrange data to labelled common so that corrections to a given
c planet are together
  100 if(Numrbs.gt.0) then
         nnrbs = Numrbs
         l = 0
 
c find moon radar bias corrections
         call REDCPY(dumi2,Nplrbt,Rdbsrt,Rdbstt,Rbiast,Lrbst,
     .               2,2,nnrbs,l,10,
     .               dumi2,Nplrbs,Rdbser,Rdbsit,Rbias,Lrbs)
c
c find sun  radar bias corrections
         call REDCPY(dumi2,Nplrbt,Rdbsrt,Rdbstt,Rbiast,Lrbst,
     .               2,2,nnrbs,l,0,
     .               dumi2,Nplrbs,Rdbser,Rdbsit,Rbias,Lrbs)
c
c find planet radar bias corrections
         do n = 1,Numpln
            nplcop = Nplnt(n)
            if(nplcop.gt.0)
     .        call REDCPY(dumi2,Nplrbt,Rdbsrt,Rdbstt,Rbiast,Lrbst,
     .                    2,2,nnrbs,l,nplcop,
     .                    dumi2,Nplrbs,Rdbser,Rdbsit,Rbias,Lrbs)
         end do
c
c arrange other planets in special sequence as described at beginning
         do j = 1,30
            call REDCPY(dumi2,Nplrbt,Rdbsrt,Rdbstt,Rbiast,Lrbst,
     .       2,2,nnrbs,l,j,dumi2,Nplrbs,Rdbser,Rdbsit,Rbias,Lrbs)
         end do
 
         call REDCPY(dumi2,Nplrbt,Rdbsrt,Rdbstt,Rbiast,Lrbst,
     .               2,2,nnrbs,l,-4,
     .               dumi2,Nplrbs,Rdbser,Rdbsit,Rbias,Lrbs)
 
         if(l.ne.Numrbs) call SUICID(
     .     ' ERROR IN NUMBER OF RBIAS CORRECTIONS,STOP IN RBSRED',13)
      endif
 
      rewind in0
      return
      end
