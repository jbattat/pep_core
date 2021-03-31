      subroutine PCOPS(m1,ivec,iabs1i)
 
      implicit none

c          j.f.chandler  -  1977 sep
c          subroutine to oversee copying of partials from old obslib
c          tape.  the decision whether to copy or compute is based on
c          the value of ict(4).
c  ict(4)=-3  do quick copy of all partials.
c         -2  do quick copy, provided that all partials exist on old
c             tape.  if some do not, act as if ict(4)=0.
c     * * * note * * in quick copy this routine is never called.
c         -1  copy only if new l vector is not set. (not implemented)
c          0  copy if old partial exists (except partials w.r.t. probe
c             parameters and certain other cases subject to change).
c          1  always compute partials.  iabs1 is set to zero in obsred
c             to indicate this option (or if there is no old tape)
c
c arguments to pcops
      character*4 ivec
      integer*4 m1,iabs1,iabs1i
c   m1 - index of 1st element to use in mvec (list mode only)
c   ivec - 4-character name of vector (e.g. 'prm ' for lprm,mprm)
c   iabs1 - copying flag, non-zero iff copying is allowed.
c          pcops must be called once for each new l vector.  note - even
c          if the l vector contains both 'slots' and a 'list' , only
c          one setup is needed.  for example, lpl(1-6) are slots for
c          initial condtions, and lpl(7-30) are a list of other
c          parameters.   the setup is done as follows:
c     l=0
c     call pcops(7,'pl  ',iabs1)
c          then process partials in two batches, e.g.:
c     iflag=0
c     call pcopy(l,6,iflag,1,lpl,mpl)
c          and then:
c     iflag=1
c     call pcopy(l,30,iflag,1,lpl,mpl)

c arguments to pcopy
      integer*4 l,lim,iflag,iswtch
      integer*2 lvec(100),mvec(100)

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'lfix.inc'
      include 'ltrapobs.inc'
      integer*2 kind
      equivalence (Numpar,kind)
      include 'mtrap.inc'
      include 'partcm.inc'
 
c local
      integer*4 m
      character*4 mesg(11)/'   L','****',' LAC','KS I','TEM ','OF M',
     .          '****',', ST','OP I','N PC','OPY '/
 
      iabs1 = iabs1i
      mesg(2) = ivec
      mesg(7) = ivec
      m = m1
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c          actual processing is done with calls to 'pcopy'
c          the control vector lvec is processed, starting at l.
c
c     input to pcopy
c   l - where to start in lvec (begin at lvec(l+1))
c   lim - highest allowed index for lvec (i.e. size of lvec)
c   iflag - mode indicator:
c        0  'slot' mode.  lvec consists of ones and zeroes.
c        1  'list' mode.  lvec contains a list of parameter numbers,
c           ended by first zero.
c   iswtch - copy switch:
c        1  Copy if possible.
c        0  Do not copy partials, but return a value indicating whether
c           the old partial could be copied.
c        -1 Copy after all.  Setup was done with the previous call, when
c           iswtch was 0 (and the copying was deferred).
c   lvec - l vector for output partials - deriv(296,2)
c   mvec - m vector for old obslib tape partials - vired(296,2)
c
c     output:
c          indices l, m (internal to pcopy), kind, and mind are updated
c          and partials are copied until either (a) the l vector is
c          finished or (b) a partial can't be copied for some reason.
c          (that includes the case where input iswtch=0)
c   iflag - completion indicator:
c        0  lvec not finished.  must call pcopy again after computing
c           latest partial (i.e. lvec(l) )
c        1  lvec is finished.  either l.ge.lim, or end of list.
c   iswtch - not changed unless input value was zero or -1.
c        0  Partial cannot be copied, must compute.
c        1  Partial can be copied if necessary (i.e., if it can't be
c           computed).  In order to copy this partial after all and
c           rejoin the logical flow, just reset iflag to indicate mode,
c           set iswtch to -1, and go back to the loop calling PCOPY.
c
      entry PCOPY(l,lim,iflag,iswtch,lvec,mvec)
 
      if(iswtch.lt.0) then
c partial can't be computed, come here to copy after all and resume
         iswtch = 0
         goto 300
      endif
 
  100 if(l.ge.lim) goto 900
      l = l + 1
      if(lvec(l).le.0) goto 150
      kind = kind + 1
 
c check for override of copying
      if(iabs1.le.0) goto 920
      if(iflag.ne.0) goto 400
 
c check if old partial exists
      if(mvec(l).le.0) goto 920
  110 Mind = Mind + 1
 
      if(iswtch.ne.0) goto 300
 
c option for probes - always try to compute partials fresh
      iswtch = 1
      goto 920
 
  150 if(iflag.gt.0) goto 900
      if(iabs1.le.0 .or. mvec(l).le.0) goto 100
 
c something missing from lvec
  160 call SUICID(mesg,11)
 
c old partial exists, copy directly
  300 Deriv(kind,1) = Vired(Mind,1)
      Deriv(kind,2) = Vired(Mind,Mun2)
      Lold(kind)    = 1
 
c go back for next partial
      goto 100
 
c process here when a zero means end-of-list (e.g. lprm)
c check for end of old partials
  400 if(mvec(m).le.0 .or. lvec(l).lt.mvec(m)) goto 920
      if(lvec(l).gt.mvec(m)) goto 160
 
c proper old partial found
      m = m + 1
      goto 110
c                return with proper flag
c finished lvec
  900 iflag = 1
      return
c not finished
  920 iflag = 0
      return
      end
