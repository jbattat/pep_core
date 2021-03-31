      subroutine PRMTRN
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   ieqm, indref, io, iprm, iprm1, j, jeq, jeqn, jeqz, jn,
     .          jqf, jsv, jsv1, k, kn, ktyp, ll, lpgr, mode, mpteq
      integer   n, neq, next, ni, nie, nnm, npmo, npteq, nstop
 
c*** end of declarations inserted by spag
 
 
c       subr. prmtrn - j.f.chandler - 1979 nov 2
c       process 'multpar' control cards saved on ibuf5.
c       convert parameter names to indices into the total parameter set
c       write output arrays back out onto ibuf5.
c        derived from:   p. macneil  june 1976  subr. analix
c                           and  april 1977  subr. pprctl
c
c       the purpose of the input controls is, in essence, to assemble
c       a set of parameter lists for the following tasks:
c       partial prereduction, parameter constraints, and multi-set runs.
c       for ppr and constraints the lists start empty, but multi-set
c       lists may be cumulative, i.e., start from a previous list.
c    1. the ppr list indicates which parameters are "uninteresting" and
c       should be reduced out of the normal equations.  obviously, this
c       list must not be empty or all-inclusive.
c    2. the constraint lists indicate which parameters to tie
c       together.  for any particular solution there may be zero or
c       more such lists.
c    3. each multi-set list contains the parameter set to be solved
c       for in an extra solve/adjust cycle at the end of processing.
c
c       ppr and multi-set lists are specified by successive additions
c       and deletions from the given starting point.  each new
c       addition or deletion is applied in the input order and may
c       duplicate or counteract a previous addition or deletion.
c       the parameter constraints, however, are specified completely
c       in one list of parameters for each constraint.
c
c           input stream syntax
c       the following annotated skeleton stream shows all of the
c       possible items that may follow the '*multpar' card.  the
c       meanings are context-dependent, so the order is important,
c       but any of the items may be omitted if so desired.
c       parameter name cards contain one name per card (one or two
c       words each).  refer to subroutine namprm for a description
c       of the names.  words with imbedded blanks should be enclosed
c       in single quotes.  the names need not be aligned to any
c       particular column, but #-commands should be left-justified.
c                  (first, ppr requests)
c     #include     (or #delete)
c     (parameter)
c     (name)
c     (cards)
c     . . .
c                  (next, multi-sets and constraints, if any.
c                   note: 1st #-card other than #include or #delete
c                   signals the end of ppr requests.)
c     #equate      (optional - global parameter constraints)
c     (parameter)
c     (name)
c     (cards)
c     . . .
c     #sne         (any #equate after this is applied after
c                   saving the normal equations.)
c     . . .
c     #corrow      (begin list of highly interesting parameters --
c                   print entire row of correlation matrix for each)
c     (parameter)
c     (name)
c     (cards)
c     . . .
c     #ref         (define first parameter set)
c     #name xxxx   (establish name 'xxxx' for this run)
c     #insert      (or #delete - introduce parameter names to be
c                   inserted or deleted for this run)
c     (parameter)
c     (name)
c     (cards)
c     . . .
c           (any number of #insert or #delete groups may be used)
c     #equate      (optional - constraints for just this run)
c     (parameter)
c     (name)
c     (cards)
c     . . .
c     #ref yyyyy   ('yyyyy' may be the name of a previous run or
c                   blank - the net result of #inserts and #deletes
c                   for that run is used as the basis for this one,
c                   or, if blank, then start with all parameters.)
c     . . .
c
c       note: the commands #insert and #include are equivalent
c
c       after the usual parameter adjustment, extra solutions are
c       performed (one for each #ref) with various subsets of the
c       overall parameter set.  in addition, each solution may be
c       constrained by tying groups of parameters rigidly together.
c       the equations are actually reduced by adding together rows
c       and columns, and the adjustment is copied out from the 1st
c       ('master') parameter of each group.  these operations are
c       performed in the input order, and so, groups may be tied
c       together by including their masters in a later group.
c
c       note: the total set of normal equations may be constrained
c       either before or after saving them (or both).  any
c       constraints applied before saving will also hold for all
c       subsequent multi-set solutions.
c
c       special pseudo-names may be used in the lists for the
c       #equate, #insert, and #delete commands.  these names may
c       be used in place of either of the usual two words of
c       a full parameter name.
c
c           all
c     represents any valid name.  its use is logically equivalent to
c     a list of individual names as follows:
c     'aaaaa all'   - all names of the form 'aaaaa zzzzz'
c     'all   aaaaa' - all names of the form 'zzzzz aaaaa'
c
c           each
c     used only with #equate.  if the master name of a list has
c     the form 'aaaaa each', the list is processed once for
c     each parameter 'aaaaa bbbbb' with 'bbbbb' substituted for the
c     2nd word of each item in the list.
c
c       the special command '#delete all' (or '#insert all') may be
c       used to delete (or insert) all parameters.
c       similarly, '#equate each' may be used to introduce a list
c       of parameter first names.
c
c        commons
      include 'aprtbf.inc'
      include 'correl.inc'
      include 'maxcrrdt.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
c matrix used as work space
c - some arrays could be moved if necessary
      common /NRMMAT/ Names(2,1000),Name(100),Ipartn(1000),
     .        Icorr(1000),Iii(1000,100),Ieqp(2000),
     .        Nm(2,2000),Kptrn(1000),Kptro(1000)
      character*8 names,name,nm
      integer*2 Iii, Ieqp, jptr(1000),Ipartn,Icorr,Kptrn,Kptro
      equivalence (Iii(1,1),jptr(1))
c note: no provision is made for catching an overflow of ieqp
c because 'nm' array follows and can be overwritten
c note: ipartn(.)=iii(.,-1), icorr(.)=iii(.,)
c
      logical*4 cando, pprstf
      character*8 key(10) /' ','#REF','#NAME','#INSERT','#DELETE',
     1 '#EQUATE', '#SNE', '#INCLUDE', '#CORROW', ' ' /
      character*8 ablnk
      equivalence (key,ablnk)
      character*8 nall/'ALL'/, neach/'EACH'/
      integer*2 ione/1/, itwo/2/, ithree/3/, ifour/4/, izero/0/, npm
      character*8 qtype(2)/'MULTISET','  PPR   '/
      character*8 an1,an2
      integer   ii, i, jj
c
c
c initialize
      Line = 60
c        types of input cards processed for each value of next (other
c        types printed out, then ignored):
c        next=1   #ref cards only or #equate for overall soln.
c                   or #corrow
c        next=2   #corrow, #name or #ref
c        next=3   #insert, #delete, #equate, #corrow or #ref
c        next=4   parameter name (for insertion), #insert,
c                   #delete, #equate, #corrow or #ref
c        next=5   parameter name (for deletion), #insert,
c                   #delete, #equate, #corrow or #ref
c        next=6   parameter name (for equation), #insert,
c                   #delete, #equate, #corrow or #ref
c        next=7   parameter name (for correlation print) or #ref
c
c        read input requests (if any exist)
      nstop = 0
      if(Ibuf5.le.0 .or. Itrwnd(Ibuf5).eq.0) goto 3600
      read(Ibuf5) nnm,((Nm(j,i),j=1,2),i = 1,nnm)
      rewind Ibuf5
      Itrwnd(Ibuf5) = -1
c
c read parameter names generated by namprm
      read(Ibuf1) (Names(1,i),Names(2,i),i = 1,Nparam)
      rewind Ibuf1
      Itrwnd(Ibuf1) = -1
c
c initialize loop to fill iii and ieqp
      Ncorp  = -1
      j      = 0
      jeqn   = 1
      pprstf = .false.
      mode   = 1
      n      = -1
      next   = 1
      if(Jct(53).gt.0) then
 
c start by looking for ppr requests
         mode = 2
         next = 3
         goto 700
      endif
  100 do while( .true. )
c
c*  start=100
c start loop to fill iii         --------------------
c increment j
         j    = j + 1
         ktyp = 1
         if(j.gt.nnm) goto 2300
 
         do kn = 1, 10
            if(Nm(1,j).eq.key(kn)) then
               ktyp = kn
               if(ktyp.eq.8) ktyp = 4
               if(mode.eq.2 .and. ktyp.ne.4 .and. ktyp.ne.5)
     .             then
 
c end of ppr requests
                  mode = 1
                  next = 1
               endif
               if(ktyp.eq.1) goto 2300
               if(ktyp.eq.2) then
c
c*  start=200
c process #ref card
                  if(next.eq.6) goto 900
                  goto 300
               else if(ktyp.eq.3) then
c end of #ref processing
c
c*  start=300
c process #name card
                  if(next.ne.2) goto 1800
                  n = n + 1
                  Name(n) = Nm(2,j)
                  next    = 3
                  if(indref.eq.0) goto 500
                  do i = 1, Nparam
                     Iii(i,n) = Iii(i,indref)
                  end do
                  goto 200
               else if(ktyp.eq.4 .or. ktyp.eq.5 .or. ktyp.eq.8)
     .                  then
c end of #name processing
c
c*  start=400
c process #insert (#delete) card
                  if(next.lt.3 .or. next.gt.6) goto 1800
                  if(next.eq.6) goto 900
                  goto 600
               else if(ktyp.eq.6) then
c end of #sne processing
c
c*  start=600
c process #equate card
                  if(next.eq.2 .or. next.gt.6) goto 1800
 
c check for legality of #equate before #ref
                  if(next.eq.1 .and. n.gt.0) goto 1800
                  if(next.ne.6) goto 1500
                  goto 900
               else if(ktyp.eq.7) then
c end of #insert (#delete) processing
c
c*  start=500
c process #sne card
                  if(n.ne.-1) goto 1800
                  if(next.eq.6) goto 900
                  goto 800
               else if(ktyp.eq.9) then
c end of bad card processing
c
c*  start=1100
c process #corrow card
                  if(next.eq.6) goto 900
                  goto 2200
               endif
            endif
         end do
c card not blank and not a control card, should then be a
c parameter name card
         ktyp = 0
c end of #equate processing
c
c*  start=800
c process parameter names
         if(next.ne.6) then
            if(next.lt.4 .or. next.gt.7) goto 1800
            goto 1600
c*  start=900
c           process parameter for #equate
c      jqf = -1: using 'each' option
c             0: first name card
c             1: ordinary names
         else if(jqf.lt.0) then
         else if(jqf.eq.0) then
 
c first of the set, check for 'each'
            jqf = -1
            if(Nm(2,j).ne.neach) then
               jqf = 1
               goto 1600
            endif
         else
            goto 1600
         endif
  200 end do
  300 if(n.eq.-1) n = 0
      if(n.eq.0) then
 
c no reference set name, use all as reference set
         indref = 0
      else if(Nm(2,j).eq.ablnk) then
         indref = 0
      else
 
c identify reference parameter set
         do i = 1, n
            if(Nm(2,j).eq.Name(i)) then
               indref = i
               goto 400
            endif
         end do
         next = 1
         goto 1800
      endif
  400 if(n.gt.99) then
         next = 1
         goto 1800
      else
         next = 2
         goto 100
      endif
  500 do i = 1, Nparam
         Iii(i,n) = 1
      end do
      goto 100
  600 next = ktyp
      if(Nm(2,j).ne.nall) goto 100
      next = 3
      if(ktyp.eq.4) goto 500
 
c '#delete all' command
  700 do i = 1, Nparam
         Iii(i,n) = 0
      end do
      goto 100
  800 next = 1
      n    = 0
      goto 100
c
c process a batch of #equate parameter names for 'each' option
c this code is executed as soon as the next command is found
  900 iprm  = Nparam
      cando = .false.
      if(jqf.ge.0) goto 1100
 
c process 'each' option
      jn = j - 1
      if(jn.le.jsv) goto 1200
      jsv1 = jsv + 1
 
c initialize scan for occurences of first name
      iprm = 0
 1000 if(iprm.lt.Nparam) then
         iprm1 = iprm + 1
         do i = iprm1, Nparam
            if(Nm(1,jsv).eq.Names(1,i)) then
 
c found match of first name
               iprm = i
               jeq  = jeq + 1
               Ieqp(jeq) = iprm
               an2 = Names(2,iprm)
 
c scan for matches to second name
               do jj = jsv1, jn
                  an1 = Nm(1,jj)
                  do ii = 1, Nparam
                     if(an1.eq.Names(1,ii) .and.
     .                   an2.eq.Names(2,ii)) then
                        jeq = jeq + 1
                        Ieqp(jeq) = ii
                        goto 1010
                     endif
                  end do
 1010          end do
               goto 1100
            endif
 
c no more matches of first name
         end do
      endif
c finished all processing for last set of names
c return to process next command
      if(cando) goto 1400
      ll = (j - jsv)/3 + 1
      assign 1300 to lpgr
      goto 1900
 
c finished a series of equates, update pointers
 1100 if(jeq.gt.jeqn + 2) then
         Ieqp(jeqn + 1) = jeq - jeqn - 1
         jeqn  = jeq + 1
         cando = .true.
      endif
 
c set up for next series of equates
 1200 Ieqp(jeqn) = n
      jeq = jeqn + 1
      goto 1000
 1300 write(Iout,2100) qtype(mode),
     .                  (Nm(1,jj-1),Nm(2,jj-1),jj = jsv,j)
      if(Mout.gt.0) write(Mout,2100) qtype(mode),
     .                        (Nm(1,jj-1),Nm(2,jj-1),jj = jsv,j)
      nstop = nstop + 1
 1400 if(ktyp.eq.1) goto 2400
      if(ktyp.eq.2) goto 300
      if(ktyp.eq.3) goto 1800
      if(ktyp.eq.4 .or. ktyp.eq.5 .or. ktyp.eq.8) goto 600
      if(ktyp.eq.7) goto 800
      if(ktyp.eq.9) goto 2200
c
c initialize for reading in #equate parameter names
 1500 Ieqp(jeqn) = n
      jeq  = jeqn + 1
      next = 6
      jqf  = 0
      if(Nm(2,j).eq.neach) jqf = -1
      jsv = j + 1
      goto 100
 1600 cando = .false.
      do i = 1, Nparam
         if(Names(1,i).eq.Nm(1,j) .or. Nm(1,j).eq.nall) then
            if(Names(2,i).eq.Nm(2,j) .or. Nm(2,j).eq.nall) then
               if(next.eq.6) then
 
c weed out duplicates
                  if(jeq.gt.jeqn + 1 .and. i.eq.Ieqp(jeqn+2))
     .                goto 1700
                  jeq = jeq + 1
                  Ieqp(jeq) = i
               else
                  Iii(i,n) = 5 - next
               endif
               if(Nm(1,j).ne.nall .and. Nm(2,j).ne.nall) then
                  if(mode.eq.2) pprstf = .true.
                  goto 100
               else
                  cando = .true.
               endif
            endif
         endif
 1700 end do
      if(cando) then
         if(mode.eq.2) pprstf = .true.
         goto 100
      endif
c end of parameter name processing
c
c*  start=1000
c process bad input cards
 1800 ll = 1
      assign 2000 to lpgr
 1900 Line = Line + ll
      if(Line.gt.58) then
         write(Iout,1950) Heding,Date,Npage
 1950    format('1PPR,MULTI-SET & CONSTRAINT CONTROLS', t40, 18A4, 1x,
     .          2A4, ' PAGE', i5/)
         Npage = Npage + 1
         Line  = 2 + ll
      endif
      go to lpgr
 2000 write(Iout,2100) qtype(mode),(Nm(i,j),i = 1,2)
      if(Mout.gt.0) write(Mout,2100) qtype(mode),
     .                        (Nm(i,j),i = 1,2)
 2100 format(' THE FOLLOWING ', a8, ' REQUEST CAN NOT BE HONORED:',
     .       (t55,3(a8,1x,a8,7x)))
      nstop = nstop + 1
      goto 100
 2200 if(next.eq.7 .or. n.gt.0) then
         next = 1
         goto 1800
      else
         n     = 0
         next  = 7
         Ncorp = 0
         goto 700
      endif
c end of #corrow processing
c*  start=2000
c processing for end of request data set
 2300 if(next.eq.6) goto 900
 2400 jeqz = jeqn - 1
      Ieqp(jeqn) = 9999
      jeqn = 1
 
c process global #equates first
      ni   = -1
      nie  = ni
      npmo = 0
      if(n.lt.0) n = 0
c
c write ppr controls, if any, before all others
      if(.not. pprstf) goto 2700
      npm = Nparam
      write(Ibuf5) ithree,npm,(Ipartn(i),i = 1,npm)
      ll = (npm - 1)/50 + 2
      assign 2500 to lpgr
      goto 1900
 2500 write(Iout,2600) (Ipartn(i),i = 1,npm)
 2600 format('0PRE-REDUCTION FLAGS:', (t25,5(i3,9I2)))
c
c*  start=2100
 2700 if(ni.lt.1) goto 3000
      npm = 0
      do i = 1, Nparam
         if(Iii(i,ni).ne.0) then
            npm = npm + 1
            jptr(npm) = i
         endif
      end do
      if(npm.eq.0) goto 3500
      write(Ibuf5) ione,npm,(jptr(i),i = 1,npm)
      if(jeqz.le.0) goto 3000
      ll = 2 + (npm - 1)/25
      assign 2800 to lpgr
      goto 1900
 2800 write(Iout,2900) nie,(jptr(i),i = 1,npm)
 2900 format('0POINTERS FOR RUN', i4, ':', (t25,25I4))
c
c*  start=2200
c write out corresponding #equate arrays
 3000 npteq = 0
 3100 if(Ieqp(jeqn).lt.ni) goto 3400
      if(Ieqp(jeqn).eq.ni) then
         neq = Ieqp(jeqn + 1) + 1
         write(Ibuf5) itwo,(Ieqp(jeqn+i),i=1,neq)
         ll = 3 + (neq - 2)/5
         assign 3200 to lpgr
         goto 1900
      else
c
c*  start=2300
c end of related #equates
         write(Ibuf5) itwo,izero,izero
c
c form inverse array
c
c set up null operation list  (npm=0)
         npm = 0
         do i = 1, Nparam
            Kptrn(i) = i
         end do
 
c apply #equates in reverse order  (flag npm non-zero)
         do while( npteq.ne.0 )
            npm   = Nparam
            mpteq = npteq + 1
            npteq = Ieqp(npteq)
            neq   = Ieqp(mpteq)
            ieqm  = Kptrn(Ieqp(mpteq+1))
 
c place ptrs to master parameter in each slave position
            do i = 2, neq
               Kptrn(Ieqp(mpteq+i)) = ieqm
            end do
         end do
         if(ni.gt.-1) then
 
c merge ptr arrays from pre- and post-sne
            if(npmo.ne.0) then
               do i = 1, Nparam
                  io = Kptro(i)
                  Kptrn(i) = Kptrn(io)
               end do
               npm = Nparam
            endif
c write out final list of inverse #equate for
c redidstribution of adjustments
            write(Ibuf5) ifour,npm,(Kptrn(i),i=1,npm)
         else
 
c initial set, save ptr array for later merging
            do i = 1, npm
               Kptro(i) = Kptrn(i)
            end do
            npmo = npm
         endif
c
c increment parameter set number (external)
         nie = nie + 1
         goto 3500
      endif
 3200 write(Iout,3300) nie,
     .                  (Ieqp(jeqn+i),(Names(k,Ieqp(jeqn+i)),k=1,2),
     .                  i=2,neq)
 3300 format('0 EQUATE PARAMETERS FOR RUN', i4, ':'/5(i5,'.',2(1x,a8)))
 
c relink list backwards for later forming inverse array
      Ieqp(jeqn) = npteq
      npteq = jeqn
 
c advance pointer to next array
 3400 jeqn = jeqn + Ieqp(jeqn+1) + 2
      goto 3100
 
 3500 ni = ni + 1
      if(ni.le.n) goto 2700
 
c end of multiparameter output
      write(Ibuf5) ione,izero,izero
      endfile Ibuf5
      rewind Ibuf5
 
c signal ibuf5 contains translated controls
      Itrwnd(Ibuf5) = -2
c end  of loop to fill iii       --------------------
c
c*  start=3000
c collect correlation print requests
      if(Ncorp.ge.0) then
         do i = 1, Nparam
            if(Icorr(i).ne.0) then
               Ncorp = Ncorp + 1
               Iicorr(Ncorp) = i
               if(Ncorp.ge.maxcrr) goto 3600
            endif
         end do
      endif
c
c*  start=9000
 3600 if(nstop.eq.0) return
      if(mod(Jct(20)/2,2).eq.0) return
      call SUICID('STOP DUE TO ERRORS IN PRMTRN', 7)
      return
      end
