      subroutine SOLCTL(b, getipt)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, iaprio, ict76o, ict76t, niobc, nsize
 
c*** end of declarations inserted by spag
 
 
c r.reasenberg/d.white   jan 1974   subroutine solctl
c increment normal equations with saved normal equations and solve
c or restore saved solution (for maximum liklihood estimator
c or kalman-bucy filter)
c
c parameter is b matrix (coefficient or covariance)
c b's dimension set in block data
      real*10 b(1)
      logical*4 getipt
c
c common
      include 'aprtbf.inc'
      include 'ciptrx.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obsdta.inc'
      include 'rtsidesl.inc'
      include 'scail.inc'
c
c local
c scale is used after iptr is finished with
      integer*2 icv(12), iptr(1000)
      equivalence(Scale, iptr)
c
c if only restoring solution, by pass form & solve
      if((Iterat .le. 1) .and. (Ict(5) .gt. 1) ) return
      call PAGSET('RESTORING NORMAL EQUATIONS  ', 7)
      do i = 1, 12
         icv(i) = 0
      end do
      icv(7) = 2
c
c           set icv(6) for dipsnec control from ict(76)
c ict(76)= 0 no dipsnec operations
c          2 flag differences between adjustable nominals
c          4 flag differences between all nominals
c          6 copy nominals from imat0(1)
c       10+n same as n, but also correct rhs for differences
c
      ict76t = Ict(76)/10
      ict76o = Ict(76) - ict76t*10
      if( ict76o .ge. 2 ) icv(6) = 2
      if( ict76o .ge. 4 ) icv(6) = icv(6) + 4
      if( ict76o .ge. 6 ) icv(6) = icv(6) + 1
      if( ict76t .le. 0 .and. Ict(76) .gt. 0 ) icv(6) = icv(6) + 64
 
c if generating pprne from old sne, bypass restoration
      if(  .not. ((Jct(53) .ge. 2) .and.  .not. (Pprdon)) ) then
c
c if multiple parameter set solution from nrm eqs
c generated this run, restore these normal equations
         if(  .not. (getipt) ) then
            if( Jct(51) .le. 0 .or.  .not. Hedskp ) then
c
c initialize for series by series rhs save
               if( Ict(18) .ne. 0 ) call SETSD(Ibuf3, Ibuf4, Ibuf1)
c
c setup nrmfrm controls
               niobc = 0
               if( Iobcon .le. 0 ) niobc = -1
c
c restore saved normal equations from previous runs
               if(  .not. ((Ict(5) .le. -1) .or. ((Ict(5).ge.1) .and. (
     .             Iterat.gt.1))) ) then
                  icv(1) = 0
                  icv(2) = niobc
                  call NRMFRM(icv, iptr, b)
               endif
c
c restore normal equations saved this run
               if(  .not. (((Ict(5).ge.1) .and. (Iterat.le.1)) .or.
     .             (Imat1 .le. 0)) ) then
                  icv(1) = 1
                  icv(2) = -1
                  call NRMFRM(icv, iptr, b)
               endif
c
c add a priori normal equations
               if( Ict(44) .ne. 0 ) then
                  icv(1) = -1
                  icv(2) = -1
                  icv(4) = Ibuf2
                  call NRMFRM(icv, iptr, b)
               endif
c
c add rhs up for series by series adjust
               if( Ict(18) .ne. 0 ) call ADDSD(Solut, Sigma)
               go to 50
            endif
         endif
         niobc  = 0
         icv(1) = 2
         icv(2) = -1
         icv(5) = 1
         if( getipt ) icv(5) = 2
         call NRMFRM(icv, iptr, b)
         icv(5) = 0
c was solctl entered only to get iptr for multiple
c parameter set solutions
         if( getipt ) return
c
c save total normal eqn and change scale
   50    nsize  = Nparam
         iaprio = mod(icv(7)/4, 2)
         call PAGSET('LEAST SQUARES ANALYSIS  ', 6)
         call NRMSAV(Scale, Hedskp, nsize, iaprio, iptr)
c
c generate partially pre-reduced normal equations
         if( Jct(53) .le. 0 ) then
c
c shall there be a matrix inversion and solution
            if( Ict(15) .eq. 0 ) then
c
c invert coefficient matrix, determine solution, calculate
c correlations
               call SNORML(Scale)
c
c fix scale of solutions
               do i = 1, Nparam
                  Sigma(i) = Sigma(i)*Scale(i)
                  Solut(i) = Solut(i)*Scale(i)
               end do
               return
            else
               call PAGCHK(60, 3, 0)
               write(Iout, 60) Ict(15)
   60          format(
     .           '-MATRIX INVERSION AND SOLUTION CANCELLED BY ICT(15)='
     .           , i3)
               Ict(1) = Iterat - 1
               if( Ict(11) .lt. -2 .and. Ict(10) .gt. -2 ) Ict(1)
     .             = Iterat
               return
            endif
         endif
      endif
      call PPRCTL
      Ict(1) = Iterat
      return
      end
