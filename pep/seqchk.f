      subroutine SEQCHK(nntp,nochk,ierr)
 
      implicit none

c       subroutine seqchk - j.f. chandler - 1979 sep
c           check for input planet numbers and series sequencing
c           for compar link.

c arguments
      integer*4 nntp,ierr
      logical*4 nochk
c     nntp= index into arrays (1-iobcon,2-iobs,3-iabs1)
c     nochk=logical flag for suppressing sequence check
c     ierr= returned error flag for non-input planets
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'comdat.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'namtim.inc'
      include 'obstap.inc'
      include 'stats.inc'
 
c local variables
      character*4 tapnm(4)/'OBCN','OBS ','ABS1','ABS2'/
      integer*4 i,nntpi,nntps,npct,npltt
 
      ierr  = 0
      npct  = 0
      npltt = Nplnta(nntp)
 
c check for existence of body
  100 if(npltt.gt.0 .and. npltt.ne.3 .and. npltt.ne.10) then
         do i = 1, Numpln
            if(npltt.eq.Nplnt(i)) goto 200
         end do
 
c body not found, signal error
         write(Iout,150) npct,npltt,tapnm(nntp)
  150    format('-+++ NPLNT', i1, '=', i3, ' FROM II', a4,
     .          ' NOT AN INPUT NPLNT +++')
         ierr = 1
      endif
 
c found that body, now find nplnt2
  200 if(npct.ne.2) then
         npct  = 2
         npltt = Nplta2(nntp)
         goto 100
 
c see if any missing bodies
      else if(ierr.eq.0) then
         nntps = nntp
         if(nochk) goto 400
         goto 300
      else
         if(nntp.ne.3) call SUICID('PLANET NOT INPUT, STOP IN SEQCHK',8)
         return
      endif
 
      entry SEQCK1(nntp,nntpi)
 
c enter here for just sequence check
      nntps = nntpi
  300 if(Ict(8).le.0) then
         if(Ntapa(nntp).lt.Ntaps(nntps)) then
            write(Iout,320) Ntapa(nntp),Ntaps(nntps),tapnm(nntps)
  320       format('0TAPE SEQUENCE NUMBER', i4, ' .LT.', i4, ' FOR II',
     .             a4)
            call SUICID('IMPROPER SERIES SEQUENCING, STOP IN SEQCHK  ',
     .                  11)
         else if(Ntapa(nntp).eq.Ntaps(nntps)) then
            if(Nseqa(nntp).le.Nseqs(nntps)) then
               write(Iout,330) Nseqa(nntp),Nseqs(nntps),tapnm(nntps)
  330          format('0SERIES SEQUENCE NUMBER', i6, ' .LE.', i6,
     .                ' FOR II', a4)
               call SUICID(
     .                    'IMPROPER SERIES SEQUENCING, STOP IN SEQCHK  '
     .                    , 11)
            endif
         endif
      endif
 
c passed all tests, save new numbers
  400 Ntaps(nntps) = Ntapa(nntp)
      Nseqs(nntps) = Nseqa(nntp)
      return
      end
