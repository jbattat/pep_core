      subroutine ANALIX
 
      implicit none

c        p. macneil  june 1976  subroutine analix
c        main program for controlling solution runs with multiple
c        parameter sets (maximum liklihood estimator)
c     logic for translating input names moved to new routine prmtrn
c     j.f.chandler - 1979 november
c
c array dimensions
      include 'globdefs.inc'
c
c        commons
      include 'aprtbf.inc'
      include 'ciptrx.inc'
      include 'inodta.inc'
      include 'fcntrl.inc'
      include 'mcnfrm.inc'
      include 'nrmmat.inc'
      integer*2 jptr(1000)
      equivalence (B,jptr)
      include 'rtside.inc'
c        information flow:
c                         in nrmfrm (regular pass):
c                              iptr --> iptrx (common);
c                         then in analix:
c                              iptrx (common) --> jti (local);
c                         then in analix:
c                              jti + (input) jptr --> iptrx (common);
c                         then in nrmfrm (for multipasses):
c                              (common) iptrx --> iptr (used in
c                              calling seq. to nrmadd).
c           local variables
      integer*4 i,ipm,ncol,ntot,numpar
      integer*2 jti(1000),ityp,npm
 
      if(Iterat.gt.1 .and. Jct(51).gt.0) call SUICID(
     .  'ITERAT GREATER THAN 1 IN ANALIX, MULTIPLE PARAMETER SET RUN '
     .  , 15)
c
c see if ibuf5 contains information
      if(Ibuf5.gt.0 .and. Itrwnd(Ibuf5).ne.0) then
 
c set up iptr if not already done
         if(Imat2.ne.0) then
            call NRMSET(0)
            call SOLCTL(B,.true.)
         endif
c
c initialize for loop through iptrx
         Hedskp = .true.
         do while( .true. )
c
c start of loop through iptrx
            read(Ibuf5) ityp,npm,(jptr(i),i = 1,npm)
            if(Itrwnd(Ibuf5).le.0) Itrwnd(Ibuf5) = 2
            if(ityp.eq.1) then
               if(npm.eq.0) then
 
                  rewind Ibuf5
                  Itrwnd(Ibuf5) = -2
                  if(Jct(52).gt.0) then
                     endfile Ibuf6
                     rewind Ibuf6
                     if(Itrwnd(Ibuf6).gt.0) Itrwnd(Ibuf6)
     .                   = -Itrwnd(Ibuf6)
                  endif
c end of loop through iptrx
c
c summary printout
                  ncol   = 2*Nsoltn + 2
                  numpar = Nparam
                  ntot   = Nsoltn
                  Nsoltn = 0
                  if(Jct(52).gt.0)
     .                call PRNSUM(ntot,B,numpar,ncol)
                  return
               else
                  if(Nsoltn.le.0) then
                     do i = 1, Nparam
                        jti(i) = 0
                     end do
 
c invert iptrx array
                     do i = 1, Mparam
                        if(Iptrx(i).gt.0) jti(Iptrx(i)) = i
                     end do
                  endif
                  Nsoltn = Nsoltn + 1
 
c now invert jptr array into iptrx
                  do i = 1, Mparam
                     Iptrx(i) = 0
                  end do
                  do i = 1, npm
                     ipm = jptr(i)
                     if(ipm.gt.0 .and. jti(ipm).gt.0)
     .                   Iptrx(jti(ipm)) = ipm
                  end do
 
c prepare normal equations
                  call NRMSET(0)
 
c add other info in saved norm eqn form and solve
                  call SOLCTL(B,.false.)
c
c write data for summary printout
                  if(Jct(52).gt.0) then
                     write(Ibuf6) Nsoltn,
     .                            (Side(i),Sigma(i),i = 1,Nparam)
                     Itrwnd(Ibuf6) = 1
                  endif
 
c adjust parameters and prepare other applications of solut.
                  call SOLUSE(B,.false.)
               endif
            endif
         end do
      else
         write(Iout,50) Ibuf5
   50    format('- NO MULTI-PARAMETER SET INFORMATION ON IBUF5=', i3,
     .          4x, '* * * WARNING IN ANALIX * * *'/)
      endif
c
c
      return
      end
