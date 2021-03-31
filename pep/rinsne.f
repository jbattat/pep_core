      subroutine RINSNE(b, side, buff, nparam, kool, ipoch)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, irow, j, k, kode
 
c*** end of declarations inserted by spag
 
 
c
c d. white  april 1974  subroutine rinsne
c
c read one epoch filter input saved norm eqn
c kool=0 => store rhs, kool=1 => store lhs
c
c parameters
      integer*4 nparam, ipoch
      integer*4 kool
      real*10   b(1)
      real*10 side(nparam), buff(nparam)
c
c common
      include 'filtds.inc'
c
c local
      integer*2 nobs
c
c id
      read(Insne) kode, Jd1fil, Jd2fil, i, nobs
c
c epoch check
      if( i .ne. ipoch ) call SUICID('EPOCH CHECK IN RINSNE   ', 6)
c
c if no obs, then no data
      if( nobs .le. 0 ) return
c
c rhs
      if( mod(kode,2) .eq. 1 ) then
         read(Insne) buff
         if( kool .lt. 1 ) then
            do i = 1, nparam
               side(i) = side(i) + buff(i)
            end do
         endif
      endif
c
c lhs
      if( kode .gt. 1 ) then
         do i = 1, nparam
            read(Insne) buff
            if( kool .gt. 0 ) then
               irow = i*(i - 1)/2
               do j = 1, i
                  k    = irow + j
                  b(k) = b(k) + buff(j)
               end do
            endif
         end do
      endif
      return
      end
