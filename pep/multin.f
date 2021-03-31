      subroutine MULTIN(in0, nstop, init)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   in0, j, nprd, nstop
 
c*** end of declarations inserted by spag
 
 
c        program to spool and read multiple parameter set requests
c        p. macneil june 1976
c        modified by p. macneil april, 1977 to also handle
c        partial prereduction of normal equations control (request)
c        cards
c
      logical*4 init
 
c commons
      include 'aprtbf.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
 
c shared work space
      common /WRKCOM/ Nm(2, 2000)
      character*8 Nm
 
      integer*4 maxreq/2000/
 
      if(Ibuf5 .gt. 0 .and. Itrwnd(Ibuf5) .ne. 0) then
         if(Itrwnd(Ibuf5) .gt. 0) rewind Ibuf5
         Itrwnd(Ibuf5) = 0
      endif
      if( .not. (init)) then
         if(Jct(51) .le. 0 .and. Jct(53) .le. 0) then
 
c * command input stream implies useful info here
            if(Jct(27) .eq. 0) return
         endif
         call PEPTIC(In, Iout, in0, 9,
     .               'PPR, MULTI-SET, CONSTRAINT REQUESTS ', nstop, 1)
c
c read requests in free field format
         call NAMRED(in0, Nm, maxreq, nprd)
         if(Ibuf5 .le. 0) then
 
c error - forgot to define ibuf5
            nstop = nstop + 1
            write(Iout, 20)
   20       format(82x, 'IBUF5=0, ERROR IN MULTIN ***')
c
c save input controls on ibuf5 for later translation
         else if(nprd .gt. maxreq) then
            nstop = nstop + 1
            write(Iout, 40) maxreq
   40       format(82x, 'MORE THAN', i5,
     .             ' REQUESTS, ERROR IN MULTIN ***')
 
         else if(nprd .gt. 0) then
            write(Ibuf5) nprd, (Nm(1,j), Nm(2,j), j = 1, nprd)
            endfile Ibuf5
            rewind Ibuf5
 
c signal ibuf5 contains information
            Itrwnd(Ibuf5) = -1
         endif
      endif
      return
      end
