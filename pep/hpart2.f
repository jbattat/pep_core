      subroutine HPART2(kick)
 
      implicit none
c
c r.w.king    june 1980    subroutine hpart2
c increment moon rotation partials w.r.t. j2 of the moon by lunar
c orbit partials.
c

c array dimensions
      include 'globdefs.inc'

c common
      include 'number.inc'
      include 'partcm.inc'
      include 'pemctl.inc'
      include 'tapdta.inc'
c
c local
      integer   i,iflag,j,k,kick,nspt,nstart
      integer*2 ltest
c
c     at present this routine is written very specifically to handle
c     the compbination of partials w.r.t. j2 of the moon from a
c     rotation integration and an orbit integration.  hence, the calling
c     arguments are minimal, and the tests for whether to call hpart2 in
c     hpartl includes three specific checks.  this was done in order
c     to minimize the number of conditions under which the subroutine
c     is called.  in the future, the routine could be generalized by
c     adding more arguments to the call list and executing the call
c     under more circumstantces.
c
c
c           initialization
      ltest  = 1031
      nstart = 8
      Kmon   = Lparm
      iflag  = -1
      nspt   = 1
      if(Nspot2.gt.0) nspt = 2
c
c determine nstart in the km array and kmon, the index to
c the j2 partials on the moon tape
      call PBDPRM(Nkimn,Kimn,nstart,Kmon,ltest,iflag)
      if(iflag.le.0) then
c
c is the partial w.r.t. j2 on the moon tape?  if so, it
c must be the first zonal
         if(Kimn(nstart).eq.2) then
c
c increment derpr by the orbit partials
c if the j2 partial is present, kmon must point to it
            call CPARTL(2,1,kick)
            do k = 1,nspt
               do j = 1,Mouse
                  do i = 1,Index
                     Derpr(i,j,k) = Derpr(i,j,k) - Dermn(i,1)
                  end do
               end do
            end do
            Ivze(1) = 1
         endif
      endif
 
      return
      end
