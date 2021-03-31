      subroutine ROTLOG(ilpl, lrotdp, mpl, iabs1)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, iset, j
 
c*** end of declarations inserted by spag
 
 
      integer*4 iabs1
      integer*2 ilpl(1), mpl(1)
      logical*1 lrotdp(12)
 
      logical*1 anyp, jmp
      integer*2 lpl(30)
      integer*2 map(12)/1, 2, 6, 7, 3, 4, 5, 8, 9, 10, 11, 12/
c
c        logic to set up lrotdp from lpl.
c        lrotdp is .true. or .false. coresponding to calculate or
c        not to calculate a partial. the order of lrotdp (and also
c        of partial calculations) is as follows:
c             phi0, p, i0, psi0, mu, delta0, alpha0
c             pdot,a,b,c,d (terms added oct. 1978)
c        logic is as follows:
c             set all lrotdp that are specified by lpl
c             if any partial at all is desired, set lrotdp(1)
c             if partial wrt alpha0 or delta0 is desired, set lrotdp(3)
c                  and lrotdp(4)
c             if partial wrt mu is desired, set lrotdp(3) and lrotdp(4)
c
c
c        set up vector lpl, internal to rotlog, which indicates
c        the parameters for which partials calculations are
c        needed. construction of lpl depends on ict(4), iabs1,
c        and mpl vector.
c
      do i = 1, 20
         lpl(i) = ilpl(i)
      end do
      if( iabs1 .gt. 0 ) then
c
c        here ict(4)=0 or -1. cannot yet implement the -1
c        option since do not have original l vector. (only the
c        merged vector is available)
c
         do i = 1, 30
            if( mpl(i) .ne. 0 ) then
               do j = 1, 30
                  if( lpl(j) .eq. mpl(i) ) then
                     lpl(j) = 0
                     go to 50
                  endif
               end do
            endif
   50    end do
      endif
c
c
      do i = 1, 7
         lrotdp(i) = .false.
      end do
      anyp = .false.
      do i = 1, 30
         if((lpl(i) .ge. 6) .and. (lpl(i) .le. 17) ) then
            anyp = .true.
            iset = lpl(i) - 5
            iset = map(iset)
            lrotdp(iset) = .true.
         endif
      end do
      if( anyp ) then
         lrotdp(1) = .true.
         jmp =  .not. (lrotdp(5) .or. lrotdp(6) .or. lrotdp(7))
         if(  .not. (jmp) ) then
            lrotdp(3) = .true.
            lrotdp(4) = .true.
         endif
      endif
 
      return
      end
