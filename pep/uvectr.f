      subroutine UVECTR(ictl, x, r, u, du)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 DOT, s
      integer   i, j, nuv
 
c*** end of declarations inserted by spag
 
 
      integer*4 ictl
      real*10 x(6), r, u(3), du(3)
c
c        r.b. goldstein  r.king  march 1978
c        routine to calculate any or all of the following:
c             r : magnitude of input vector
c             u : unit vector
c             du: derivative of unit vector (wrt time)
c        the input vector, x, contains position and velocity
c
c        ictl: binary coded control integer
c             1 bit=1:  calculate r
c             2 bit=1:  calculate u
c             4 bit=1:  calculate du
c
c        if any of the bits are not set, and the associated quantity
c        is needed further on in the routine (i.e. r is needed
c        to calculate u) it is assumed that
c        the correct quantities were sent throuth the calling
c        statement.
c
c        common
      include 'fcntrl.inc'
      include 'inodta.inc'
 
      logical*4 prt
 
      prt = mod(Jct(6)/2, 2) .eq. 1
      nuv = 0
c
c compute r
c
      if( mod(ictl,2) .ne. 0 ) r = SQRT(DOT(x,x))
c
c compute u
c
      if( mod(ictl/2,2) .ne. 0 ) then
         do i = 1, 3
            u(i) = x(i)/r
         end do
         nuv = 3
      endif
c
c compute du
c
      if( mod(ictl/4,2) .ne. 0 ) then
         s = DOT(u, x(4))
         do i = 1, 3
            j     = i + 3
            du(i) = (x(j) - s*u(i))/r
         end do
         nuv = 6
      endif
 
      if( prt ) then
         if( Line .ge. 56 ) call OBSPAG
         write(Iout, 50) r, (x(i), i = 1, nuv)
   50    format(' UVECTR: R=', 1pd22.15, '  X=', (t40,3D23.15))
         Line = Line + 1 + nuv/4
      endif
 
      return
      end
