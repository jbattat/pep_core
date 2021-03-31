      subroutine SOLUSE(b,keepit)
 
      implicit none
c
c r.reasenberg/d.white   jan 1974   subroutine soluse
c adjust parameters and prepare other applications of solution
c
c
c parameter is b matrix (coefficient or covariance)
c b's dimension set in block data
c if keepit true then keep this solution
      logical*4 keepit
      real*10 b(1)
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'rtsidesl.inc'
      include 'scail.inc'
c
c local
      logical*4 oldsol
      integer*4 imatx,i,ij,j,j1,l
c
c restore solution and scale it if necessary
      oldsol = .false.
      if((Iterat.le.1) .and. (Ict(5).gt.1)) then
         oldsol = .true.
         call SOLFRM(Jmat,Scale,1)
c
c solve for "reduced away" parameters
         if(Ict(5).eq.4) call SOLPPR(Scale)
      endif
c
c store adjusted process noise ic's on direct access
      if(Ict(42).gt.0) call ICNSTR(Solut)
c
c adjust parameters
      call ADJUST(keepit)
c
c while solution in core, form rhs contribution of a priori
c parameter values for next iteration.  (destroys b)
      if((Ict(44).ge.1) .and. (Iterat.lt.Ict(1))) call FORMU(b)
c
c series by series adjustment
      if(.not. (oldsol .or. Ict(15).ne.0 .or. Ict(18).eq.0))
     .    call ADJSER(Solut,Iout,b,b(Nparam*23+1),Sigma)

      if(Jct(51).lt.0) return

c Enter here to restore the solution, if possible, after other activity
c has overlaid it.  Remember whether an old solution was restored for
c the purposes of this run, and restore it again if necessary.

      entry SOLUS2(b,keepit)
c
c if no uncertainty of prediction then leave
      if(Ict(1).lt.Iterat) return
      if(Ict(10).le.-2 .or. Ict(14).le.0) return
c
c if restored solution + reduced param's then don't have covar.
      if(Iterat.ne.1 .or. Ict(5).ne.4) then
c
c if restored solution then only scale covariance matrix
         if(oldsol .and. Ict(44).le.0 .and. Jct(51).ge.0) goto 300
c
c if did an extra half-iteration, then must restore solution anyway
         imatx=Imat
         if(oldsol) imatx=Jmat
c
c if imat not on write message
         if(imatx.gt.0) then
c
c restore covariance matrix for uncertainty of prediction
            call SOLFRM(imatx,Scale,2)
            goto 300
         endif
      endif
      call PAGCHK(60,2,0)
      write(Iout,100) Imat,Ict(5),Ict(10),Ict(14)
  100 format(/' IMAT=', i2, ' ICT(5)=', i2, ' ICT(10)=', i2,
     .       ' ICT(14)=', i2,
     .' COVARIANCE MATRIX NOT SAVED FOR UNCERTAINTY IN PREDICTION, PROGR
     .AM WILL BE STOPPED AFTER ADJUSTMENTS')
      if(Mout.gt.0) write(Mout,200)
  200 format(' COVARIANCE MATRIX NOT SAVED FOR UNCERTAINTY IN PREDICT')
      Ict(1)  = Iterat - 1
      Ict(14) = 0
      return
c
c change scale of covariance matrix
  300 l = 0
      do i = 1, Nparam
 
c test for information in this row
         if(b(l+i).eq.1._10 .and. Solut(i).eq.0._10
     .    .and. Sigma(i).eq.1._10) then
            ij = l
            j1 = 1
            do j = 1, Nparam
               ij = ij + j1
               if(j.ne.i .and. b(ij).ne.0._10) goto 350
               if(j.ge.i) j1 = j
            end do
 
c no information in this row, clear diagonal element
            b(l + i) = 0._10
         endif
 
  350    do j    = 1, i
            l    = l + 1
            b(l) = b(l)*Scale(i)*Scale(j)
         end do
      end do
      return
      end
