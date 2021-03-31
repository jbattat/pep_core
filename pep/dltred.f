      subroutine DLTRED(in0, nstop, init)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, in0, jct69, LEG, nstop
      real      temp
 
c*** end of declarations inserted by spag
 
 
c m.e. ash  aug 1967  subroutine dltred
c read in data controlling observation series changes and error
c weighting and also data controlling dummy observations, then
c write into disk data set iobcon.
      logical   init
c
c        modified for *command july 1978  r.b. goldstein
c        modified to allow blank fields in data editing mode
c           1980 march - z.m.goldberg and j.f. chandler
c        * * * warning * * * there must be something non-blank in
c           cols 1-16, or peptic will treat it as all blank * * *
c
c        * * * warning * * * site, spot, and series names may not be
c           changed to blanks, even by intent
c
c     1983 feb - alternate header card format added
c
c        common
      include 'dltflg.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obsdta.inc'
 
c temporary work space
      common/WRKCOM/ Freq,Freq2,Card(10),Card2(9),Cardb(9),Blank(9),
     .        Ctlg,Erwgt
      character*8 card,card2,cardb,blank,ctlg
      character*4 site1,series,site2,spot,spot2
      real*10 freq,freq2
      real*4    Erwgt(2), acctim, fdev
      integer*4 ntape, nseq
      integer*2 ncodf, nplnt, itime, nrewnd, nplnt2
c
c     ncodf = code indicating type of observation series
c           ncodf=1  radar observation series of time delay and
c                    doppler shift
c           ncodf=4  meridian circle observation series of right
c                    ascension and declination referred to the true
c                    equinox and equator of date
c           ncodf=5  photographic observation series of right
c                    ascension and declination referred to the mean
c                    equinox and equator of reference epoch
c           ncodf=11 interferometry observation series of differential
c                    delay and differential delay rate
c           ncodf=12 quasi-vlbi observation series of differential
c                    n-count or differential delay increment
c           ncodf.gt.20 observation series of type ncodf - 20 of two
c                    observed objects
c***note. description here of first card of series incomplete
c***see comment cards in subroutine prmred for first card of series
c***for observation card data set iobs
c     nplnt= planet number of observed body
c           0=sun
c           1=mercury
c           2=venus
c           4=mars
c           5=jupiter
c           6=saturn
c           7=uranus
c           8=neptune
c           9=pluto
c          10=moon
c          11-30= other possible natural planets,satellites or asteroids
c           31,...= artificial space probes
c     site1  = first four characters of eight character receiving
c              observing site name for radar and optical observations
c     series = four characters giving observation series name
c     site2  = first four characters of eight characters sending
c              observing site name for radar observations
c     spot   = four characters giving name of spot actually observed
c              on given body nplnt, if not the center of the body or
c              sub-radar point.
c     erwgt(i)=factor by which quoted error of measurement i is
c              multiplied
c           i=1    radar measurement of time delay or optical
c                  measurement of right ascension
c           i=2    radar measurement of doppler shift or optical
c                  measurement of declination
c           if erwgt(1).le.-1.e3 and erwgt(2).le.-1.e3, then whole
c              observation series is skipped
c     acctim = time accuracy constant for observation series
c              for radar observations: accuracy constant for time
c              delay iteration usually 10**-6 (seconds of delay)
c              for meridian circle observations: accuracy constant for
c              meridian crossing iteration, usually 10**-3 (radians of
c              Earth rotation).
c              for transit/occultation observations: accuracy constant
c              for alignment iteration, usually 10**-3 (seconds of time)
c     itime  = time code for observation series
c            = 10*c + f  (for c .ge. 0)
c            = c         (for c .lt. 0)
c         For radar observations, c is forced to be non-negative, so
c         f may be -1 in this case (means same as f=2).
c         Only last two digits of year are given for observations, so
c         an indication must be given at start of series what first
c         two digits are.  A series must not cross the start of a
c         century year, except for dummy observations.
c                 c is the century code, i.e., YEAR/100 - 19
c              1 21st century observations
c              0 20th century observations
c             -1 19th century observations
c             -2 18th century observations,etc.
c                 f is a flag for the type of observable
c                 (treated as 0 for c .lt. 0) 
c              0 optical: before atomic time
c              0 radar: observation time tag is UT2
c              1 atomic time exists, offset to be computed
c                observation time tag is UTC
c              2 AT - UT offset suppied with data
c              3 radar: observation time tag is UTC send time
c              5 normal point pseudo-obs at coordinate time
c                  these values of itime are necessary because each
c                  observation card only has the last two digits of the
c                  year and an indication must be made at the start of
c                  the series what the first two digits are
c     ctlg  =  name of reference star catalog used in reducing
c              photographic observations to r.a. and dec.
c     fdev  =  fractional offset in units of 10**-10 from a.1
c              time of unit of time used in time delay measurement
c              for radar observation series
c     freq  =  frequency of sending radar site for radar
c              observation series
c     nrewnd=  0 at start of observation series only those tapes are
c              rewound which are not used in the series
c     nrewnd= -1 same as nrewnd=0 with site card(s) read after 1st card
c     nrewnd=  1   at start of observation series all tapes are rewound
c     nrewnd=  2 same as nrewnd=1 with site card(s) read after 1st card
c     ntape =  observation library tape number
c              if 0, signal to program to generate theoretical value of
c              observations and partial derivatives for given
c              observation series.
c             must have ntape greater or equal to previous ntape read in
c     nseq  = sequence number of observation series on given observation
c             library tape.
c             must have nseq greater than previous nseq read in
c             for same value of ntape.
c
c           variables describing dummy observation or alterations of
c           observation errors
      integer*4 jd(2), intday, intsec
      integer*2 ihr(2), imin(2)
      real*4    sec(2), erobs(2)
c
c     if ntape is greater than zero, the observations of the given
c           observation series lying between julian date jd(1) at
c           hour ihr(1), minute imin(1), second sec(1) and julian date
c           jd(2) at hour ihr(2), minute imin(2), second sec(2) have the
c           errors of their first and second measurements multiplied by
c           erobs(1) and erobs(2), respectively.
c     intday and intsec have no meaning if ntape is greater than zero.
c
c     if ntape is less than or equal to zero, then the theoretical
c     values and partial derivatives of the theoretical values of the
c     observations of the given observation series lying between julian
c     date jd(1) at hour ihr(1), minute imin(1), second sec(1) and
c     julian date jd(2) at hour ihr(2), minute imin(2), second sec(2)
c     are calculated at intervals of intday days & intsec seconds.  for
c     use in forming the coefficient matrix of the normal equations, the
c     errors of these dummy observations are assumed to be erobs(1) and
c     erobs(2) for measurements 1 and 2 respectively.
c
c     for both the options ntape.le.0 or ntape.gt.0 above, if erobs(i)
c     is zero,  theoretical value of measurement i is calculated but
c     not used in normal equations, whereas if erobs(i) is .le. -1.e03,
c     the theoretical value of measurement i is not calculated (i=1,2).
c
c
c                 ix    observation library series  (dltred)
c
c  changes to observation series
c  for each series
c card 1
c  columns
c    1- 2  ncodf (type of series)                                i2
c    3- 5  nplnt (observed planet no.)                           i3
c    6- 6  htype (flag for extra header cards)                   a1
c    7-10  receive site name -first 4 characters                 a4
c   12-15  series name       -first 4 characters                 a4
c   17-20  sending site name -first 4 characters                 a4
c   22-25  spot name         -first 4 characters                 a4
c   26-31  errwgt(1) error weight for 1st observable             e6.0
c   32-37  errwgt(2) error weight for 2nd observable             e6.0
c   38-43  acctim (time accuracy constant)                       e6.0
c   44-45  itime  (time code for series)                         i2
c   46-52  fdev   (frac.offset from a.1 time in units of 1e-10)  f7.2
c   53-70  freq   (frequency of sending site)                    d18.11
c   71-72  nrewnd                                                i2
c   73-75  ntape  (observation library tape)                     i3
c   76-80  nseq   (sequence number of observation series)        i5
c  notes-must put in order of increasing nseq
c       -this format corresponds to record 3 of obs.lib.
c       -if errwgt(i).le.-1.e3 then whole obs.series is deleted!
c       -if jct(69).gt.0 the 1st card has a different format from
c        column 46 onward, and there may be 1 or more extra cards.
c       -the format may also be selected via '#80' or '#72' cards
c        inserted in front of any series (for old and new formats,
c        respectively).  in any case, the format remains selected
c        until changed by another '#' card.  the new format has:
c  columns
c   46-63  freq   (frequency of sending site)                    d18.11
c   64-65  nrewnd                                                i2
c   66-68  ntape  (observation library tape)                     i3
c   69-72  nseq   (sequence number of observation series)        i4
c       -if jct(69).gt.0 and htype=':' there is an extra header
c        card ('colon card') with the following format:
c  columns
c    1- 1  htype (further card indicator)                        a1
c    2- 8  fdev  (frac.offset in time)                           f7.2
c    9-16  ctlg  (name of star catalog)                          a8
c
c       -if ncodf.gt.20 an extra card is read with format:
c  columns
c    1- 2  nplnt2 (planet no. of second observed object)         i2
c             (or, if jct(69).gt.0)
c    1- 3  nplnt2 (planet no. of second observed object)         i3
c    4- 7  second spot name -first 4 characters                  a4
c    8-25  freq2  (reference frequency for second object)        d18.11
c card 2     (may be more than one of these cards)
c  columns
c    1- 7  jd(1)      if ntape.le.0 then obs.from time jd(1)     i7
c    8-10  ihr(1)                                     ihr(1)     i3
c   11-13  imin(1)                                   imin(1)     i3
c   14-21  sec(1)                                     sec(1)     f8.4
c   25-31  jd(2)      to time jd(2)                              i7
c   32-34  ihr(2)            ihr(2)                              i3
c   35-37  imin(2)          imin(2)                              i3
c   38-45  sec(2)            sec(2) are calculated with assumed  f8.4
c   49-56  errobs(1)  errors errobs(1)                           e8.1
c   57-64  errobs(2)         errobs(2)                           e8.1
c   65-66  intday     at intervals of intday days                i2
c   67-72  intsec                     intsec seconds             i6
c        if ntape.gt.0 then obs. are calculated with values of errobs
c                     if errobs(i).le.0 then meas(i) is deleted
c                     in the normal equations.
c
c        a blank card separates series
c        2 blank cards signifies the end
c
      character*4 ablnk/'    '/
      real*10 zro8/0._10/
      integer*4 zro4
      integer*2 zro2
      equivalence(zro8, zro4, zro2)
      character*1 htype, colon/':'/
      character*4 num80/'#80'/
 
      if( init ) then
c
c*  start=200
         write(Iobcon) (zro8, i = 1, 4), zro2, zro4, zro8, zro2, zro8,
     .                 zro2, zro4, zro8, zro8, (zro8, i = 1, 4)
         go to 400
      else
         call PEPTIC(In, Iout, in0, 9,
     .               'OVERRIDE ERRORS AND DUMMY OBS. CARDS', nstop, 2)
 
c create long string of blanks
         call MOVEBL(ablnk, 4, Blank, 72)
         jct69 = Jct(69)
      endif
  100 do while( .true. )
c
c read variables describing observation series
c and write onto disk
c initialize in-core copy
         call MOVEBL(ablnk, 4, Card2, 2*72)
         read(in0, 150) Card
  150    format(10A8)
         if( LEG(1,1,Card,1,num80) .ne. 0 ) then
            backspace in0
            fdev   = 0.
            Ctlg   = Blank(1)
            nplnt2 = 0
            spot2  = ablnk
            Freq2  = 0._10
            if( jct69 .gt. 0 ) then
 
c new format header card a
               read(in0, 160) ncodf, nplnt, htype, site1, series,
     .                        site2, spot, Erwgt, acctim, itime, Freq,
     .                        nrewnd, ntape, nseq
  160          format(i2, i3, a1, 3(a4,1x), a4, 3E6.0, i2, d18.11, i2,
     .                i3, i4)
               if( LEG(1,1,htype,1,colon) .eq. 0 ) then
 
c colon card
                  read(in0, 170) htype, fdev, Ctlg
  170             format(a1, f7.2, a8)
                  backspace in0
                  read(in0, 150) Cardb
               endif
 
               if( ncodf .le. 20 ) go to 300
               read(in0, 180) nplnt2, spot2, Freq2, htype
  180          format(i3, a4, d18.11, a1)
            else
 
c old format header card a
               read(in0, 190) ncodf, nplnt, site1, series, site2, spot,
     .                        Erwgt, acctim, itime, fdev, Freq, nrewnd,
     .                        ntape, nseq
  190          format(i2, i3, 4(1x,a4), 3E6.0, i2, f7.2, d18.11, i2,
     .                i3, i5)
               call MVC(Card, 46, 7, Cardb, 2)
               call MVC(Card, 53, 23, Card, 46)
               call MVC(Card, 77, 4, Card, 69)
               if( ncodf .le. 20 ) go to 300
               read(in0, 200) nplnt2, spot2, Freq2
  200          format(i2, 1x, a4, d18.11)
            endif
            backspace in0
            read(in0, 150) Card2
            go to 300
         else
            jct69 = 1
            if( LEG(3,1,Card,1,num80) .eq. 0 ) jct69 = 0
         endif
      end do
 
c set modification flags
  300 Gncode    = ncodf .ne. 0
      Gnplt1    = nplnt .ne. 0
      Gsite1    = site1 .ne. ablnk
      Gseres    = series .ne. ablnk
      Gsite2    = site2 .ne. ablnk
      Gspot     = spot .ne. ablnk
      Gerwgt(1) = LEG(6, 26, Card, 1, Blank) .ne. 0
      Gerwgt(2) = LEG(6, 32, Card, 1, Blank) .ne. 0
      Gacctm    = LEG(6, 38, Card, 1, Blank) .ne. 0
      Gitime    = LEG(2, 44, Card, 1, Blank) .ne. 0
      Gfreq     = LEG(18, 46, Card, 1, Blank) .ne. 0
      Gnrwnd    = LEG(2, 64, Card, 1, Blank) .ne. 0
      Gfdev     = LEG(7, 2, Cardb, 1, Blank) .ne. 0
      Gctlg     = LEG(8, 9, Cardb, 1, Blank) .ne. 0
      Gnplt2    = nplnt2 .ne. 0
      Gspot2    = spot2 .ne. ablnk
      Gfreq2    = LEG(18, 8, Card2, 1, Blank) .ne. 0
      write(Iobcon) ncodf, nplnt, site1, series, site2, spot, Erwgt,
     .              acctim, itime, fdev, Freq, nrewnd, ntape, nseq,
     .              nplnt2, spot2, Freq2, Ctlg, Gncode, Gnplt1, Gsite1,
     .              Gseres, Gsite2, Gspot, Gerwgt, Gacctm, Gitime,
     .              Gfdev, Gfreq, Gnrwnd, Gnplt2, Gspot2, Gfreq2,
     .              Gctlg, zro8
      if( ncodf .gt. 0 .or. nseq .ne. 0 ) then
c
c read site variables
         if(nrewnd.lt.0 .or. nrewnd.gt.1) then
            call SUICID('NREWND NOT 0 OR 1, STOP IN DLTRED   ',9)
         endif
         do while( .true. )
c
c read variables describing dummy observations or
c alterations of observation errors and write onto disk.
            read(in0, 320) (jd(i), ihr(i), imin(i), sec(i), i = 1, 2),
     .                     (erobs(i), i = 1, 2), intday, intsec
  320       format(2(i7,2I3,f8.4,3x), 2E8.1, i2, i6)
c adjust to non-negative fraction of day, less than 1
c (assuming input is not too strange)
            do i = 1, 2
               temp = 3600*ihr(i) + 60*imin(i) + sec(i)
               if( temp .lt. 0. ) then
                  if( sec(i) .lt. 0. ) then
                     sec(i)  = sec(i) + 60.
                     imin(i) = imin(i) - 1
                  endif
                  if( imin(i) .lt. 0 ) then
                     imin(i) = imin(i) + 60
                     ihr(i)  = ihr(i) - 1
                  endif
                  if( ihr(i) .lt. 0 ) then
                     ihr(i) = ihr(i) + 24
                     jd(i)  = jd(i) - 1
                     go to 340
                  endif
               endif
               if( temp .ge. 86400. ) then
                  if( sec(i) .ge. 60. ) then
                     sec(i)  = sec(i) - 60.
                     imin(i) = imin(i) + 1
                  endif
                  if( imin(i) .ge. 60 ) then
                     imin(i) = imin(i) - 60
                     ihr(i)  = ihr(i) + 1
                  endif
                  if( ihr(i) .ge. 24 ) then
                     ihr(i) = ihr(i) - 24
                     jd(i)  = jd(i) + 1
                  endif
               endif
  340       end do
            write(Iobcon) (jd(i), ihr(i), imin(i), sec(i), i = 1, 2),
     .                    (erobs(i), i = 1, 2), intday, intsec
            if( jd(1) .le. 0 ) go to 100
         end do
      endif
c
c no more data cards
  400 endfile Iobcon
      rewind Iobcon
      rewind in0
      return
      end
