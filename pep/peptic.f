      subroutine PEPTIC(input,print,output,nword,title,nstop,lend)

      implicit none

c
c     *end or eof may define the end of any input sub-stream
c    n.brenner/r.reasenberg  5/74  subroutine peptic
c        pep translator of input controls
c     program to replace symbolic names delimited by / with given
c     translation strings.
c
c     card images are read from an input unit (like 5), converted and
c     both printed (on 6, say) and sent to an output unit (1,say).
c
c     the input arguments are--
c                            -- input, print, and output, which are the
c     numbers of the input unit, the print and the output units.
c                            -- the title and its length, nword
c tf revised march 1975      -- nstop which is pep's input error count
c     which causes a stop at the end of the input link for nstop.gt.0
c tf revised march 1975      -- error is set to 0 in this routine each
c     time this routine is entered and counts errors detected by this
c     routine.  nstop is updated.
c
c    input commands that start in column 1 with a *
c
c        *diction  defines the start of a dictionary building
c        sequence of cards.  each such card has a name (column 1
c        thru 8) and a length (column 10,11 ),  it may have an
c        alpha string (of length = len bytes) starting in column 13.
c        if len is zero, the alpha string is the name of a previously
c        defined dictionary entry.  if len is negative, the aplha string
c        is on the next card and has length -len.
c        format(a8,i3,1x,8a8,a4)
c
c        *dictend  defines the end of the dictionary building string
c
c        *dictcle  clears the dictionary
c
c     *end may define the end of a non-namelist input stream for all
c     current calls to this routine.  two blank cards are output on
c     data set 'output'.
c
c     &end defines the end of a namelist input stream.  the &end may
c     occur anywhere in columns 2 thru 80, except after '$'.
c
c     one blank card ends the non-namelist input streams except for
c     dltred which requires two blank cards in a row, except that the
c     very first card being blank ends dltred stream.  a card is
c     considered to be blank if the first 16 columns are blank.
c     at the end of each sub-stream either two blank cards or '&end'
c     is written on 'output' for good measure.
c
c           t.forni  march 1975  modified this routine so that
c           old input stream will always work with no extra cards
c        modified 1978 aug. - interfaced with * commands and external
c        source files for input
c
c           lend indicates what is considered the end of a particular
c           input stream
c     for  lend=-1    an end of file  or  *end card (not used as yet)
c     for  lend= 0    an  &end  anywhere in cols 2 thru 80
c     for  lend= 1    1 blank card
c     for  lend= 2    2 blank cards in a row, except that the very first
c                     card being blank ends the stream
c           note:  a card is considered to be blank if the first 16
c                  columns are blank.  this is only tested when lend=1,2
c
c

c arguments
      integer*4 input,print,output,nword,nstop,lend
      character*4 title(nword)

c common
      include 'crdbuf.inc'
      character*1 cardc(80)
      equivalence (cardc, Card8)

c local variables
      integer   i,iam,icol1,idict,iend,jcard,jstore,lcard,len,nrd
      character*1 delim/'/'/
      character*8 scmd(5)/'*DICTION','*DICTEND','*DICTCLE','*END    ',
     .          '**END   '/
      character*1 bstar
      equivalence (scmd(1),bstar)
      character*80 cards
      integer   pmod, pold/0/
      logical   goon

      character*8 name,store(1000),dname(100)
c set up maximum length of dictionary arrays and initial pointers
c user must include explicit *dictclear card to clear the dictionary
      integer   dlen(100), dloc(100)
      integer*4 ndict/1/, maxdic/99/, istore/2/, maxsto/990/

c set up all blank default dictionary entry for unknown names
      data dlen(1)/1/, dloc(1)/1/, store(1)/' '/, dname(1)/' '/

      character*8  blnk/'        '/
      integer*4 iblnk, itf, itst, error
c define namelist recognition character
      include 'nmlchrdt.inc'
c external functions
      integer*4 ENDTST,ISCAN

      error = 0
      iblnk = 0
      itf   = 0
      itst  = 0
      jcard = 0

      goon = .true.

      if(.not. Useit) then
         call PEPTIN(input, print, error)
      else
         if(Card8(1).eq.scmd(5)) goto 1600

c should reread card8, is it ok?
         Useit    = .false.
         if(cardc(1).eq.bstar) then
            do i = 1, 4
               if(Card8(1).eq.scmd(i)) goto 100
            end do

c error - unknown command - skip it
            write(print, 20) Card80
   20       format(1x, a80, ' UNKNOWN * COMMAND. ERROR PEPTIC')
            error = error + 1
            call PEPTIN(input, print, error)
         endif
      endif
  100 do while(Ieof.ne.1)
         jcard = jcard + 1

c convert tabs to blanks
         do i=1,80
            if(cardc(i).eq.char(9)) cardc(i)=' '
         end do

c check the first column of the card for a *
         if(cardc(1).ne.bstar) then
c
            if(itf.ne.1) then
               itf = 1
               write(print, 110) Card80, title
  110          format(1x, a80, 1x, 12A4, 1A2)
            endif

c save card as it came
            call MVC(Card8, 1, 80, cards, 1)
            pmod = 0

c remove $ sign  and columns to the right of $ sign if any
            lcard = ISCAN(cardc(1), 80, '$') + 1
            if(lcard.lt.1) then
               lcard = 80
            else if(lcard.eq.1) then

c print card with $ in column 1
               if(jcard.gt.1) write(print, 110) Card80
               call PEPTIN(input, print, error)
               goto 200
            else
               call MOVEBL(' ', 1, cardc(lcard), 81-lcard)
            endif
            if(lend.eq.0) then

c check for '&end'
               if(ENDTST(Card8).gt.0) goon = .false.
            else if(lend.ne.-1) then
c
c test on number of blank cards = lend for end of particular
c input stream
c blank card if first 16 columns are blank
               if(Card8(1).ne.blnk .or. Card8(2).ne.blnk) iblnk=0
               if(Card8(1).eq.blnk .and. Card8(2).eq.blnk) iblnk=iblnk+1
               if(iblnk.eq.lend) goon = .false.

c detect single blank card as end
               if(jcard.eq.1 .and. iblnk.eq.1) goon = .false.
            endif
            if(.not.goon .or. itst.eq.0) goto 1400
c
c     scan the card for delimiters flanking a symbol of up to 8 long
c     if such is found, blank it out and cover it up with the translate
c     string from the dictionary, which is placed starting at the first
c     delim and going as long as it is, covering up stuff on the card
            icol1 = 1
            do while( .true. )
               iam = ISCAN(cardc(icol1), lcard-icol1, delim) + 1
               if(iam.le.0) goto 1400
               pmod = pmod + 1
               iend = ISCAN(cardc(iam+icol1), lcard+1-icol1-iam, delim)
               if(iend.lt.0) iend = 8
               call MOVEBL(cardc(iam+icol1), iend, name, 8)
               do idict = 1, ndict
                  if(name.eq.dname(idict)) goto 120
               end do
               idict = 1
               write(print, 320) name
               error = error + 1
  120          call MOVEBL(store(dloc(idict)), dlen(idict),
     .                     cardc(iam+icol1-1), iend + 2)
               icol1 = icol1 + iam + iend + 1
               if(icol1.ge.lcard) goto 1400
            end do
         else

c found a * in column 1
            itst = 1
            if(itf.ne.1) then
               itf = 1
               write(print, 130) title
  130          format(82x, 12A4, 1A2)
            endif

c check list of * commands
c *DICTION
            if(Card8(1).eq.scmd(1)) then
               write(print, 140) Card80
  140          format(1x, a80, ' START OF  DICTIONARY')
               goto 300
c *DICTEND
            else if(Card8(1).eq.scmd(2)) then
               goto 400
c *DICTCLE
            else if(Card8(1).eq.scmd(3)) then
               ndict  = 1
               istore = 2
               write(print, 150) Card80
  150          format(1x, a80, ' CLEAR THE DICTIONARY')
               call PEPTIN(input, print, error)
               goto 200
            endif
c
c have a **END card
            if(Card8(1).eq.scmd(5)) Useit = .true.

c have a *END or **END
            if(Card8(1).eq.scmd(4) .or. Card8(1).eq.scmd(5)) then
               write(print, 110) Card80
               goto 1600
            endif

c send the card back to input for testing
            Useit = .true.
            goto 1600
         endif
  200 end do
      goto 1600
c
c read a dictionary entry card
  300 call PEPTIN(input, print, error)
      if(Ieof.eq.1) goto 1200
      name = Card8(1)
      if(name.eq.scmd(2)) goto 400

c the length of the translate string is in bytes
      call MVC(Card8, 13, 68, store(istore), 1)
      i = 9
      call DECODI(Card8, i, 12, len)
      nrd = (7 + iabs(len))/8

      if(len.lt.0) then
c the length of the translate string is negative if the string is
c the next card
         len = -len
         call PEPTIN(input, print, error)
         if(Ieof.eq.1) goto 1200
         call MVC(Card8, 1, len, store(istore), 1)
         goto 600
      else if(len.eq.0) then
c The length is zero when the translate string is another symbol.
c Search the dictionary for a given name symbol.  If not found,
c the default all-blank entry is returned.
         do idict = 1, ndict
            if(store(istore).eq.dname(idict)) goto 340
         end do
         error = error + 1
         write(print, 320) store(istore)
  320    format(81x, ' UNKNOWN NAME ', a8, ' ERROR PEPTIC')
         idict = 1
  340    len   = dlen(idict)
         nrd   = (len+7)/8
         ndict = ndict + 1
         if(ndict.gt.maxdic) goto 800
         dloc(ndict) = dloc(idict)
         goto 700
      else
         goto 600
      endif

  400 write(print, 500) Card80
  500 format(1x, a80, ' END OF DICTIONARY')
      call PEPTIN(input, print, error)
      goto 100
  600 ndict = ndict + 1
      if(ndict.gt.maxdic) goto 800
      dloc(ndict) = istore
      istore = istore + nrd

c check for overrun
      if(istore.ge.maxsto) goto 800
  700 dname(ndict) = name
      dlen(ndict)  = len
      goto 1000
  800 write(print, 900) ndict, istore
  900 format(' DICTIONARY OVERFLOW.  NDICT=', i5, ', ISTORE=', i6, 56x,
     .       ' ERROR PEPTIC')
      ndict  = maxdic
      istore = maxsto
      nrd    = 0
      error  = error + 1
 1000 jstore = dloc(ndict)
      write(print, 1100) name, len, (store(jstore+i-1), i = 1, nrd)
 1100 format(' DICTIONARY: ', a8, i3, 1x, 12A8)
      goto 300
 1200 write(print, 1300)
 1300 format(81x, ' EOF WHILE READING DICTIONARY.  ERROR PEPTIC')
      error = error + 1
      goto 1900
c
c print the translated card
 1400 if(pmod.gt.pold) then
         write(print, 1450) cards, Card80
 1450    format(1x, a80, ' CHANGED TO'/1x, a80)
      else
         if(jcard.gt.1) write(print, 110) cards
      endif
      if(lend.eq.0) then
         if(goon .and. lcard.gt.1 .and.
     .        (cardc(2).ne.'&' .or. ISCAN(cardc(2),79,'=').gt.0)) then
            do i = lcard, 2, -1
               if(cardc(i).eq.',') goto 1460
               if(cardc(i).ne.' ') then
                  if(i.lt.80 .and. cardc(i).ne.'=') cardc(i+1) = ','
                  goto 1460
               endif
            end do
         endif
 1460    call NMLFIX(Card8)
      endif
      write(output, 1500) Card80
 1500 format(a80)

      if(goon) then
         call PEPTIN(input, print, error)
         goto 100
      endif
 1600 if(lend.gt.0) write(output, 1700)
 1700 format(/)
      if(lend.eq.0) write(output, 1800) Nmlchr
 1800 format(' ', a1, 'END')
c
      if(error.eq.0) goto 2100
 1900 write(print, 2000) error
 2000 format(81x, ' ERROR COUNT =', i5, '  PEPTIC')
      nstop = nstop + error

 2100 endfile output
      rewind output
      return
      end
