      subroutine PEPTIN(input, print, error)

      implicit none


c*** start of declarations inserted by spag
      integer   i, input, irc, ISCAN, j, k, nesf, NSCAN

c*** end of declarations inserted by spag


      integer*4 print, error
c        subroutine peptin  --  j.f.chandler  august 1978
c        cms interface added 1980 dec
c        this routine reads cards from the input stream until it
c        encounters a *include card.  then it reads from the requested
c        external source file either until the end or an optional
c        maximum number of card images.  anyway, it then switches back
c        to the usual input stream.  note: a *include found within an
c        external source file has the same effect as one in the
c        primary input stream except that it causes an immediate end of
c        that file, that is, there is no push-down facility.  any
c        errors in syntax or missing files or i/o errors on external
c        source files are noted and tallied in 'error'.  the syntax:
c        *include ddname(member)
c        there is another format available for cms files:
c        *include 'name type mode'
c        the ddname or fileid may begin in any column up to 16.
c
c        any cards with a '$' in column 1 will be printed and ignored.
c
c        input - fortran unit of the primary input stream
c        print - fortran unit for messages
c        error - cumulative error count from caller
c
      include 'crdbuf.inc'
      character*1 cardc(80)
      equivalence (cardc, Card80)

      character*8 erc(8)/'NOT FND', 'NO DD  ', 'RECFM  ', 'NO OPN ',
     .          'I/0 ERR', 'I/O ERR', 'FILEID ', ' ????? '/, qqfm,
     .          fname, mname
      character*8 qcms/' FILEID '/, qfile/' FILE DD'/, qmem/' MEMBER '/
      character*8 sincld/'*INCLUDE'/
c
c*  start=100
      if(Ieof.eq.1) return
      if(Useit) goto 900
  100 if(Lcnt.eq.0) goto 700

c currently reading from a file, see if reached limit
  200 if(nesf.gt.1) then
         call CMSGET(Card80, 0, irc)
      else
         call PDSGET(Card80, irc)
      endif
      if(irc.eq.0) Lcnt = Lcnt + 1
      if(irc.lt.0) goto 500
      if(irc.eq.0) goto 900

c error on read attempt
  300 write(print, 400) erc(irc)
  400 format(81x, ' **PDS/CMS ERROR IN PEPTIN (', a7, ')**')
      error = error + 1
c*  start=200
c signal end of reading from latest file
  500 write(print, 600) Lcnt
  600 format(82x, 'END OF INCLUSION AFTER', i4, ' CARDS')
      Lcnt = 0

c get a card from input
  700 read(input, 800, end=1200) Card80
  800 format(a80)
  900 if(cardc(1).eq.'$') then
c ignore a comment card
         write(print, 920) Card80
  920    format(1x, a80)
         Useit = .false.
         goto 100
      endif

      if(Card8(1).ne.sincld .and. .not.Useit) call PLOG(Card80)
      Useit = .false.
      if(Card8(1).ne.sincld) return
c*  start=300
c set up source file
      if(Lcnt.gt.0) write(print, 950) mname
  950 format(82x, 'WARNING. FILE ', a8, ' ENDS AT *INCLUDE')
      write(print, 1000) Card80
 1000 format(1x, a80, ' START OF INCLUSION')

c find file ddname
      qqfm = qfile

c find next non-blank
      k = NSCAN(cardc(9), 8, ' ')
      if(k.ge.0) then
         k = 9 + k
         if(cardc(k).eq.'''') then
c*  start=400
c cms file
            qqfm = qcms
            k    = k + 1
            i    = ISCAN(cardc(k), 80 - k, '''')
            if(i.gt.0) then
               call CMSOPN(cardc(k), i, irc)
               nesf = 2

c set up mname in case of error
               j = NSCAN(cardc(k), i, ' ')
               if(j.ge.0) then
                  i = ISCAN(cardc(k+j), 8, ' ')
                  if(i.lt.0 .or. i.gt.7) i = 7
                  call MOVEBL(cardc(k+j), i + 1, mname, 8)
               endif
               goto 1100
            endif
         else
            i = ISCAN(cardc(k), 9, '(')
            if(i.gt.0) then
               call MOVEBL(cardc(k), i, fname, 8)
               i = i + k + 1

c find member name
               qqfm = qmem
               j    = ISCAN(cardc(i), 9, ')')
               if(j.gt.0) then
                  call MOVEBL(cardc(i), j, mname, 8)

c open pds member
                  call PDSOPN(fname, mname, irc)
                  nesf = 1
                  goto 1100
               endif
            endif
         endif
      endif
      write(print, 1050) qqfm
 1050 format(81x, ' BAD', a8, 'NAME.  ERROR IN PEPTIN')
      error = error + 1
      goto 100
 1100 if(irc.gt.0) goto 300
      Lcnt = 0
      goto 200
c*  start=500
c reached end of input stream
 1200 Ieof = 1

c*  start=900
      return
      end
