      subroutine ACMIN(in0,nstop,init)
 
      implicit none

c     d. white  april 1973  subroutine acmin
c
c        modified for * commands rbg/jfc
c
c The first card in the a priori input stream to PEP is a title
c card.  Then come one or more a priori input groups.  An input
c group begins with a namelist (NMLST4) which has as members:
c    N - I*4 - the number of parameters for which there are a priori
c            values or sigmas.
c    DIAGON - L*4 - .true. if the input B matrix is in diagonal form,
c            otherwise matrix is assumed to be lower diagonal half.
c    COVAR - L*4 - .true. if the input B matrix is a covariance matrix,
c            otherwise assumed that matrix is coefficient (or
c            information) matrix.
c    APEST - L*4 - .true. if an a priori estimate vector is included in
c            this group, otherwise assumed none such.
c    SQRTB - L*4 - .true. (.and.diagon) B is the sqrt of the required
c            value, i.e., either sigma or 1/sigma
c    OFFSET - L*4 - .true. if a priori estimates are of offset from
c            first parameter in packet.  In this case, the first
c            parameter name is a reference parameter, others are
c            offsets.  The X matrix may contain actual values, or X(1)
c            may be zero with offsets in X(2-N); only the differences
c            are significant.  The B matrix may contain sigma's (if
c            COVAR is .true.) or 1/sigma's (if COVAR is false); SQRTB
c            is used, APEST and DIAGON are ignored.
c            Note: offset a priori information generates a singular
c            a priori matrix.  For this reason, a priori nominals are
c            not computed, and iteration is prohibited.
c    W - R*8 - weight factor by which the input matrix (coeficient or
c            covariance) is multiplied.  If zero, the matrix is not
c            used -- NAMRED and MATCH are not called.
c    B - R*8 - input B matrix whose form, content, and dimension are
c            described above.
c    X - R*8 - input a priori parameter values (if APEST .true.) for
c            parameters whose sigmas are in the B matrix.
c
c  The parameters for which a covariance or coefficient matrix and a
c  priori estimates are input are indicated by name cards immediately
c  following the namelist.  One parameter name goes on each card and
c  the names are in the same order as the values in the X and B
c  matrices.  Most parameters have two part names; these are described
c  in subroutine NAMPRM.  The names can be placed on the name cards in
c  free field format, at least one blank separating the first part
c  from the second.  The name cards have to be ended by a blank card.
c  Any number of groups of namelists and names can be input and the
c  program stops looking when a NMLST4 with N = 0 is found.
c
c  Data sets needed are IBUF1 for parameter names and IBUF2 for
c  the a priori B and U matrices.  Optional printout of input
c  values is made on LOUT (if LOUT .ne. 0).
c
c parameters
      integer*4 in0, nstop
      logical*4 init
c
c common
      include 'aprtbf.inc'
      include 'crdbuf.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'matzro.inc'
      include 'numnum.inc'
c
c work space for input link
      common/WRKCOM/ Names(2,1000),Apcoef(5050),Xap(100)
      character*8 Names
      real*10 Apcoef,Xap
c
c local and arrays equivalenced with work space
      integer i,ieofs,ij,im,j,jlim,n,n1,nprd
c maximum parameters in one a priori coeff matrix = 100
c 100 * (100 + 1) / 2 = 5050
      character*8 apnams(2,100)
      real*10 b(5050)
      integer*4 numpar,nappar
      equivalence (Apcoef,b),(nappar,n)
      real*10 x(100),w,sumsq
      equivalence (x,Xap)
      real*10 iskale(1000),xnom(1000),rhs(1000)
      logical*4 diagon,covar,sqrtb,offset
      logical*4 apest,anyest,dg,alldg,anyoff
      character*80 aptitl
c
c namelists
      namelist /NMLST4/nappar, n, diagon, covar, Apcoef, b, x, Xap,
     .         apest, w, sqrtb, offset
c
c initialization
      anyest = .false.
      anyoff = .false.
      alldg  = .true.
      do i = 1, 1000
         iskale(i) = 1._10
         rhs(i)    = 0._10
      end do
      call ZFILL(Coeff, 16*Nrmsiz)
      sumsq = 0._10
c
c rewind data sets for repeated input streams
      if(Lout.gt.0) rewind Intern
      if(Ibuf1.gt.0) rewind Ibuf1
      if(Ibuf2.gt.0) rewind Ibuf2
c
c get names of parameters being adjusted
      if(Ict(44).ne.0 .or. Ibuf1.ge.1)
     .   call NAMPRM(numpar, Names, iskale, xnom)
c
c see if any apriori input
      if(Ict(44).eq.0 .or. init) return
      ieofs = Ieof
c
c determine input data set for apriori
      im = In
      if(Ict(44).gt.1) im = Ict(44)
 
c if necessary, turn off any pending external source file
      if(im.ne.In) Lcnt = 0
c
c read and print title
      call PEPTIN(im, Iout, nstop)
      aptitl = Card80
      write(Iout, 100) aptitl
  100 format(1x, A80, ' A PRIORI INPUT TITLE CARD')
      do while(.true.)
c
c loop over sets of nmlst4 and names until nappar = 0
c first initialize nmlst4 parameters
         nappar = 0
         w     = 1._10
         sqrtb = .false.
         diagon = .false.
         covar  = .false.
         apest  = .false.
         offset = .false.
         call ZFILL(Xap, 16*100)
         call ZFILL(Apcoef, 16*5050)
c
c spool nmlst4
         call PEPTIC(im, Iout, in0, 9,
     .               'A PRIORI B AND X MATRICES  &NMLST4  ', nstop, 0)
         read(in0, NMLST4)
         rewind in0
         if(nappar.eq.0) goto 200
c
c note if any parm estimates are input
         if(w.ne.0._10) anyest = anyest .or. (apest .and..not.offset)
         if(w.ne.0._10) anyoff = anyoff .or. offset
         dg = diagon
         if(nappar.eq.1 .or. (nappar.eq.2.and.Apcoef(2).eq.0._10) .or.
     .       w.eq.0._10) dg = .true.
         alldg = alldg .and. dg
c
c spool a priori parameter name cards until 1 blank card
         call PEPTIC(im, Iout, in0, 6, 'A PRIORI PARAMETER NAMES',
     .               nstop, 1)
         if(w.ne.0._10) then
c
c now read names ignoring blank records
            call NAMRED(in0, apnams, nappar, nprd)
c
c check for mismatched number of parameters
            if(nprd.ne.nappar) goto 110
            if(nappar.lt.100) then
               ij   = nappar
               jlim = 1
               if(.not.diagon) ij = nappar*(nappar + 1)/2
               n1 = nappar + 1
               do i = n1, 100
                  if(apest .and. Xap(i).ne.0._10) goto 110
                  if(.not.diagon) jlim = i
                  do j = 1, jlim
                     ij = ij + 1
                     if(Apcoef(ij).ne.0._10) goto 110
                  end do
               end do
            endif
            goto 120
 
  110       call SUICID('EXTRANEOUS A PRIORI INFORMATION ', -8)
c
c match apriori names with main names and insert coefficients
  120       call MATCH(Coeff,Names,numpar,xnom,Apcoef,
     .                 apnams,w,nappar,Xap,diagon,covar,
     .                 apest,offset,nstop,iskale,rhs,sqrtb,sumsq)
         endif
      end do
 
c
c write out compressed coeff matrix in format of saved norm eqn
c note: acwrap clobbers 'names'
  200 call ACWRAP(Coeff,rhs,xnom,numpar,aptitl,anyest,
     .            alldg,nstop,anyoff,sumsq)
      if(im.ne.In) then
 
c return to ordinary input
         Ieof = ieofs
         Lcnt = 0
      endif
      if(anyest .and. Ict(44).eq.-1) call SUICID(
     . 'APRIORI PARAMETER VALUES SHOULD NOT BE USED WITH ICT(44)=-1 ',
     . -15)
      return
      end
