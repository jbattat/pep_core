      subroutine PERDTA(in0)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, in0
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   june 1966    subroutine perdta
c assign standard values to peripheral data set numbers
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'aprtbf.inc'
      include 'inodta.inc'
      include 'maxmt0dt.inc'
      include 'obsdta.inc'
      include 'plndta.inc'
c
c           setup /inodta/
c     in to be determined from typewriter or card reader input (later)
c     in =5    setup in subroutine prmred
c     in0=98   setup in subroutine prmred or read from title card
c     iout =6  setup in subroutine prmred
      Jout   = 0
      Kout   = 0
      Lout   = 0
      Mout   = 0
      Nout   = 0
      Ipunch = 7
      Jpunch = 0
      Kpunch = 0
      Igraph = 97
      Intern = 99
      Ibuf   = 4
      Imat   = 0
      Jmat   = 0
      Ipert  = 90
      Jpert  = 0
      Kpert  = 0
      Nummt0 = 0
      do i = 1, maxmt0
         Imat0(i) = 0
      end do
      Imat1 = 0
      Imat2 = 0
      Imat3 = 0
      Imat4 = 0
      Ictat = 0
      Ieng  = 0
      Jeng  = 0
      Keng  = 0
c
c setup /plndta/
      Iem  = 13
      Imn  = 20
      Inut = 0
      Ilib = 0
      do i = 1, u_mxpl
         Iplnt(i) = 0
      end do
      Iplcon = 1
      Ivcnt  = 0
c
c setup /obsdta/
      Iobs = 5
      do i = 1, 10
         Iobs0(i) = 0
         Iobs1(i) = 0
         Iobs2(i) = 0
      end do
      Iobcon = 2
      Numobt = 1
c
c setup /aprtbf/
      Ipert1 = 0
      Ipert2 = 0
      Jpert1 = 0
      Jpert2 = 0
      Kpert1 = 0
      Kpert2 = 0
      Ibuf1  = 0
      Ibuf2  = 0
      Ibuf3  = 0
      Ibuf4  = 0
      Ibuf5  = 0
      Ibuf6  = 0
      Ibuf7  = 0
      Ibuf8  = 0
      Ibuf9  = 0
      Ibuf10 = 0
      do i = 1, 30
         Epsa(i) = 0.0E0
      end do
      Epsa(30) = 5.0E-5
c
c change tables for error monitor
      call ERRSET(207, 300, 15, 0, 0, 209)
      call ERRSET(211, 1, 1, 0, 0, 214)
      call ERRSET(215, 1, 1, 0, 0, 215)
 
c call errset (215,300,300,0,0,215)
      call ERRSET(217, 1, 1, 0, 0, 220)
      call ERRSET(221, 1, 1, 0, 0, 225)
c     call errset (221,300,300,0,0,225)
c
c     ihc207 - ihc209 set to unlimited errors, up to 15 messages.
c
c     ihc211 - ihc214 , ihc217 - ihc220  program stopped for any
c     occurance, except err option in read statement takes precidence
c     if it applies.
c
c     ihc215, ihc221 - ihc225  set to unlimited allowable errors.
c     in subroutine prntot the program is stopped if any errors occured
c     (unlimited allowed to this point to catch all input errors). if no
c     such errors, errset called to set number of allowed errors to none
c     and program continues.
c
      return
      end
