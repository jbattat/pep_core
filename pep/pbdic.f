      subroutine PBDIC(kp, lplnt, kplnt, l, iq2, klak, iswtch)
 
      implicit none

c        j.f.chandler , 1977 aug
c        logic control for finding partials w.r.t. initial conditions
c        from integration tape for body nplnt(klak).

c arguments
      character*2 iq2
      integer*4 lplnt,kplnt,l,iswtch
      integer*2 klak,kp(30)
c        the integration control integer array 'kp' is examined to see
c        if i.c. partial number 'l' exists.  if so, tape ptr 'kplnt'
c        and kp ptr 'lplnt' are updated, and control returns.
c        'iswtch' indicates whether the corresponding observation
c        partial can be copied from saved values (1 if so, 0 if not).
c        if the partial exists neither place then the offending body and
c        kp are printed and processing halted.  if the partial was not
c        integrated but can be copied, then iswtch is set equal to -1
c        and control returns.
c               note: 'iq2' is a 2-letter abbreviation for the body.
c        special values for klak: -3 earth-moon barycenter,
c                                 -2 moon
c                                 -1 earth rotation
c                                  0 moon rotation

c array dimensions
      include 'globdefs.inc'

c commons
      include 'inodta.inc'
      include 'namtim.inc'
 
      character*8 aemmn(4)/' EMBARY ','  MOON  ','  EARTH ','  MOON  '/,
     1 rotn/'ROTATION'/,blank/'        '/, rname,aname
      integer   i, lt
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(kp(1).ne.0) then
         lt = l + 1
         do i = lplnt, lt
            if(kp(i).ge.0) kplnt = kplnt + 1
            end do
         lplnt = lt + 1
         if(kp(lt).ge.0) return
      endif
      if(iswtch.eq.0) then
         rname = blank
         if(klak.gt.0) then
 
c planet other than earth or moon
            aname = Aplnt(klak)
            if(Nplnt(klak).lt.0) rname = rotn
         else
 
c earth or moon (or earth or moon rotation)
            aname = aemmn(klak + 4)
            if(klak.ge.-1) rname = rotn
         endif
         write(Iout, 50) (iq2, i, kp(i), i = 1, 7), iq2, l, aname,
     .                   rname
   50    format(7('   K',a2,'(',i1,')=',i3)/' BUT L', a2, '(', i1,
     .')>0. INCONSISTENT CONTROL CONSTANT FOR PARTIAL W.R.T. INITIAL CON
     .DITION FOR ', a8, 1x, a8)
         if(klak.gt.0) write(Iout, 100) klak, Nplnt(klak)
  100    format(10x, '(  NPLNT(', i2, ')=', i3, '  )')
         call SUICID('CONTROL CONSTANT INCONSISTENCY, STOP IN PBDIC   ',
     .               12)
      else
 
c partial can be copied from old obslib vector
         iswtch = -1
      endif
      return
      end
