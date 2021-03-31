      subroutine INITAL(called, omitac, initf, in0, nstop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, in0, nstop
 
c*** end of declarations inserted by spag
 
 
      logical*1 called(20)
      logical   omitac, initf
c
c r.b. goldstein  july 1978
c
c common
      include 'fcntrl.inc'
 
      if(  .not. called(2) ) call SITRED(in0, nstop, initf)
      if(  .not. called(3) ) call SPTRED(in0, nstop, initf)
      if(  .not. called(4) ) call RBSRED(in0, nstop, initf)
      if(  .not. called(5) ) call EQNRED(in0, nstop, initf)
      if(  .not. called(6) ) call PHSRED(in0, nstop, initf)
      if(  .not. called(7) ) call DLTRED(in0, nstop, initf)
 
      do i = 2, 7
         called(i) = .true.
      end do
 
      if(  .not. (omitac) ) then
         if(  .not. called(9) ) call ACMIN(in0, nstop, initf)
         if(  .not. called(10) ) call MULTIN(in0, nstop, initf)
 
         do i = 9, 10
            called(i) = .true.
         end do
         if( Ict(42) .gt. 0 ) then
            if(  .not. called(8) ) call FILTIN(in0, nstop, initf)
            if(  .not. called(13) ) call FILTPN(in0, nstop, initf)
            if(  .not. called(14) ) call FILTEP(in0, nstop, initf)
 
            called(8)  = .true.
            called(13) = .true.
            called(14) = .true.
         endif
      endif
 
      return
      end
