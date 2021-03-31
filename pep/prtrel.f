      subroutine PRTREL(ltest,nkb,kb,lparb,itype,kick,iswtch)
 
      implicit none

c  j.f. chandler  -  1977 sep
c  compute partial with respect to relativity motion factor, including
c  scale factor from 'hippo'.  if control 'ltest' does not indicate
c  a partial w.r.t. body rel. factor, then either a) suicid with usual
c  message if iswtch=0 or b) set iswtch to -1 and return if iswtch>0.

c arguments
      integer*4 lparb,itype,kick,iswtch
      integer*2 ltest,nkb,kb(99)

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'ltrapobs.inc'
      include 'obscrd.inc'
      include 'partcm.inc'
      integer*4 kbods(8)
      equivalence (kbods,Kembry)

c local
      integer*4 i,jtype,ksave

c for messages
      character*8 empp(5)/'EARTH','MOON','PLANET','PROBE','PROBEC'/,
     .          mesg(9)/'CAN''T CO','MPUTE PA','RTIALS W','. R. T. ',
     .          '********','PARAMETE','RS, STOP',' IN PRTR','EL'/
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if(ltest.eq.24) goto 200
 
c not relativity motion factor (or can't find partial on tape
  100 if(iswtch.gt.0) then
c
c not found, but further action is possible
         iswtch = -1
         return
      else
         mesg(5) = empp(itype)
         call SUICID(mesg,17)
      endif
 
c look for partial on body tape
  200 do i = 8,nkb
         if(kb(i).eq.0) goto 100
         if(kb(i).lt.31) then
         else if(kb(i).eq.31) then
            jtype = itype
 
c get correct index into array of partials ptrs.
            if(itype.eq.5) jtype = 8
            ksave = kbods(jtype)
            kbods(jtype) = i - (7 - lparb)
            call CPARTL(itype,2,kick)
 
c restore partials ptr
            kbods(jtype) = ksave
 
c scale output partials
            if(Nice.le.0) Deriv(Numpar,1) = Deriv(Numpar,1)
     .          /Hippo(1)
            if(Nice.ge.0) Deriv(Numpar,2) = Deriv(Numpar,2)
     .          /Hippo(1)
            return
         else
            goto 100
         endif
      end do
      goto 100
      end
