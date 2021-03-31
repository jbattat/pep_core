      subroutine SNORML(scale)
 
      implicit none

c     ash/amuchastegui   nov 1968    subroutine snorml
c     inversion of coefficient matrix and solution of double precision
c     normal equations with extended precision intermediate steps
c     using fact that coefficient matrix is symmetric
c
c arguments
      real*10 scale(1000)
c
c common
      include 'aprtbf.inc'
      include 'correl.inc'
      include 'fcntrl.inc'
      include 'filtim.inc'
      include 'inodta.inc'
      include 'nrmmat.inc'
      real*10 qb(1001)
      equivalence (qb,B)
      include 'rtside.inc'
 
c miscellaneous quantities (dimensions set for order 1000)
      common /WRKCOM/ Pvrow(2000),Pvcol(2000),Pvrwb(2),Iuse(1000),
     . qtmp(1000)
      real*10 Pvcol, Pvrow, Pvrwb
      integer*2 Iuse
      real*10 qtmp
      real*10 f(1000)
      equivalence (f,Pvrow(1001))
 
c local
      integer*4 nfvctk/-1/, iaprio/0/, izr4/0/
      integer*4 elim,i,i0,ier,iercod,ierinv,iersol,ivectk,j,k,
     .          length,n0velim,n2,nsize
      real*4 one/1.0/
      real*10 zero/0._10/
      logical*4 print
c
c*  start=100
      ierinv = 0
      iersol = 0
      ier    = 0
      iercod = 0
      nsize  = Nparam
      print = Fict(6).eq.0
c
c save coefficient matrix and right side on disk
      call FWSIG(Ibuf,-nsize,B,Pvrwb)
      write(Ibuf) nsize,(Side(j),j = 1,nsize)
      endfile Ibuf
      rewind Ibuf
c
c printout coefficient matrix
      if(print) then
         call NRMRIT(' (SCALED) SYMMETRIC COEFFICIENT ', B, nsize,
     .               Measmt, Ict(47).gt.0)
 
c print out elapsed times
         call TIMRIT('SAVING AND WRITING OUT NORMAL EQUATIONS ', 10)
      endif
c
c        solve normal equations and determine inverse of coefficient
c        matrix with double precision arithmetic using fact that
c        matrix is symmetric
c        gauss-jordan method with pivots on diagonal
c
      if(Ict(46).le.0) then
c double precision direct inversion
         call SYMINV(B,Side,Nparam,1,Pvrow,Pvcol,Iuse,Pvrwb,
     .               ierinv)

      else if(Ict(46).eq.1) then
c direct inversion with extended precision on intermediate steps
         n2 = Nparam*2
         do i = 2, n2, 2
            Pvcol(i) = 0.0_10
            Pvrow(i) = 0.0_10
         end do
         Pvrwb(2) = 0.0_10
         call SYMINX(B,Side,Nparam,1,Pvrow,Pvcol,Iuse,Pvrwb,
     .               ierinv)
      else if(Ict(46).eq.2) then
c direct inversion with full extended precision, then convert back
c this could be pushed into a separate subroutine to keep r*16 confined
         do i=1,Nparam
            qtmp(i) = Side(i)
         end do
         length = (Nparam*(Nparam+1))/2
         do i=1,length
            qb(length+1-i) = B(length+1-i)
         end do
         call SYMINQ(qb, qtmp, Nparam, 1, Pvrow, Pvcol, Iuse, Pvrwb,
     .            ierinv)
         do i=1,length
            B(i)=qb(i)
         end do
         do i=1,Nparam
            Side(i)=qtmp(i)
         end do
      else if(Ict(46).ge.10) then
         call DERANK(B,Side,Nparam,Ict(46)-10,Pvrow,Pvcol,Iuse,ierinv)
c at this point, B contains not only the scaled covariance, but also the
c eigenvector matrix. Pvrow has the reciprocals of the eigenvalues.
c Pvcol has the adjustments to the eigenparameters
         if(Jout.gt.0) then
            elim=MIN(Nparam-1,Ict(46)-10)
            write(Jout,20) Nparam, elim
   20       format('1ADJUSTMENTS TO',i5' SCALED EIGENPARAMETERS (',
     .       i5,' ELIMINATED)'/
     .       '0 NUM      ADJUSTMENT     SIGMA       FRACT')
            do i=2,4
               qtmp(i)=0._10
            end do
            do i=1,Nparam
               if(Pvrow(i).gt.0._10.and.i.le.Nparam-elim) then
                  qtmp(9)=SQRT(Pvrow(i))
                  qtmp(1)=Pvcol(i)/qtmp(9)
                  qtmp(2)=qtmp(2)+qtmp(1)
                  qtmp(3)=qtmp(3)+ABS(qtmp(1))
                  qtmp(4)=qtmp(4)+qtmp(1)**2
               else
                  qtmp(1)=0._10
                  qtmp(9)=0._10
               endif
               write(Jout,30) i,Pvcol(i),qtmp(9),qtmp(1)
   30          format(i5,'. ',1pe16.8,1pe10.3,0pf10.3)
            end do
            do i=2,4
               qtmp(i)=qtmp(i)/(Nparam-elim)
            end do
            qtmp(10)=Pvrow(Nparam)/Pvrow(1)
            qtmp(11)=Pvrow(Nparam-elim)/Pvrow(1)
            write(Jout,40) Nparam-elim,qtmp(2),qtmp(3),SQRT(qtmp(4)),
     .       qtmp(4),qtmp(10),qtmp(11)
   40       format(/' THE',i5,' EIGENPARAMETERS USED HAD THE FOLLOWING',
     .       ' FRACTIONAL STATISTICS'/
     .       10x, 'AVERAGE FRACT', 1pe14.5/ 5x, 'AVERAGE ABS(FRACT)',
     .       1pe14.5/1x, 'ROOT MEAN SQUARE FRACT', 1pe14.5/7x,
     .       'AVERAGE FRACT**2', 1pe14.5/
     .       ' THE RATIO OF LARGEST TO SMALLEST EIGENVALUE',1pe14.5/
     .       ' THE RATIO OF LARGEST TO SMALLEST USED      ',1pe14.5)
            if(elim.gt.0) then
c compute the potential contributions of the eliminated eigenvalues
c version 1 - RSS of one-sigma offsets to all such
c version 2 - direct sum of the actual omitted adjustments
               write(Jout,44) Nparam
   44          format('1BIASES IN SIGMA UNITS TO THE',i5,
     .          ' CONVENTIONAL PARAMETERS DUE TO ELIMINATION'/
     .          '0 NUM   RSS(ONE-SIGMA)    ABS(ACTUAL)')
               n0velim=Nparam*(Nparam+Nparam-elim)
               call PRODCT(B(n0velim+1),Pvcol(Nparam-elim+1),f,
     .          Nparam,elim,1)
               do i=1,Nparam
                  qtmp(9)=SQRT(B((i*(i+1))/2))
                  qtmp(1)=0._10
                  do j=1,elim
                     qtmp(1)=qtmp(1)+Pvrow(Nparam-elim+j)*
     .                B(n0velim+j*Nparam-Nparam+i)**2
                  end do
                  qtmp(1)=SQRT(qtmp(1))/qtmp(9)
                  qtmp(2)=ABS(f(i))/qtmp(9)
                  write(Jout,45) i,qtmp(1),qtmp(2)
   45             format(i5,'.',1p2e16.8)
               end do
            endif
            write(Jout,'('' END OF EIGENVALUE REPORT'')')
         endif
      else
c other inversion methods to be implemented here
      endif
c
c*  start=1000
c redistribute 'equated' solutions, covariances
      if(Ibuf5.gt.0) call CNSTRI(Ibuf5,Nparam,B,Side,scale,Iuse)
c
c write out page heding and timer information
      if(print) then
         if(Ict(46).lt.10) then
            call TIMRIT(
     .   'GAUSS-JORDAN DIRECT INVERSION, SOLUTION OF SYMMETRIC MATRIX '
     .   , 15)
         else
            call TIMRIT(
     .   'MATRIX INVERSION, SOLUTION BY EIGENVECTOR DECOMPOSITION ',14)
         endif
      endif
      if(ierinv.le.0 .and. Ict(49).gt.0) then
c
c iterative cleanup of covariance matrix
         call SNORM2(f,Pvrow,Pvcol,Pvcol(1001),B,Nparam,Ict(49),
     .               Eps(11),iercod)
      else if(print) then
         call PAGCHK(60,2,0)
         write(Iout,50) Nparam
   50    format('0NO CLEAN-UP OF INVERSE OF ORDER', i4)
      endif
c
c*  start=2000
c printout and save covariance matrix
      if(print) call NRMRIT('COVARIANCE (INVERSE COEFFICIENT)', B,
     .                        nsize, Measmt, Ict(48).ge.0)
c
c save solution vector and covariance matrix
      if(Imat.gt.0) then
         call SAVHED(Imat,'SOLUTCOV')
         ivectk = 1
         if(Jct(60).gt.0) ivectk = 0
         write(Imat) Measmt,izr4,izr4,one,one,Ermeas,zero,zero,
     .               iaprio, ivectk, (zero,i = 1,10)
         if(ivectk.gt.0) write(Imat) nfvctk,
     .                             (Vectk(i),i = 1,Nparam)
         write(Imat) Nparam,(Side(i),i = 1,Nparam),
     .               (scale(i),i = 1,Nparam)
         call FWSIG(Imat,-nsize,B,Sigma)
 
c end file imat
         rewind Imat
         call PAGCHK(60,2,0)
         write(Iout,100) Imat
  100    format(
     .       '0SOLUTION VECTOR AND COVARIANCE MATRIX SAVED ON DATA SET'
     .       , i3)
 
c write out timer information
         call TIMRIT('  WRITING OUT COVARIENCE MATRIX ', 8)
      else if(print) then
         call PAGCHK(60,2,0)
         write(Iout,150)
  150    format('0SOLUTION VECTOR AND COVARIANCE MATRIX NOT SAVED')
         call TIMRIT('  WRITING OUT COVARIENCE MATRIX ', 8)
      endif
c
c
c*  start=3000
      if(Ict(18).ne.0) call SNORM6(Pvrow,Pvcol,B,scale)
c
c write out error matrix f
      if(Ict(49).ge.0) then
         call MULERR(B,f,Pvrow,Pvcol,Nparam)
      else if(print) then
         call PAGCHK(60,2,0)
         write(Iout,200)
  200    format('0ERROR MATRIX NOT CALCULATED')
      endif
c
c*  start=4000
c           calculate standard deviations and correlations,
c           printout correlation matrix, read normal equations
c           back into storage
      call SNORM4(Pvcol,Pvcol,B,Sigma,nsize,Ict(48).lt.-1,print,Ncorp,
     .            Iicorr,Ibuf1,ier)
 
      i0 = 0
      do i = 1, nsize
         read(Ibuf) j,(B(i0+k),k = 1,i)
         i0 = i0 + i
      end do
      read(Ibuf) j,(Pvrow(i),i = 1,nsize)
      rewind Ibuf
c
c clean up solution, calculate error in solution
      call SNORM5(Pvrow,Pvcol,iersol)
c
c error message
      if((iersol+ier+iercod+ierinv).gt.0) then
         Ict(1) = Iterat - 1
         call PAGCHK(60,3,0)
         write(Iout,250) iersol,ier,iercod,ierinv
  250    format(
     .'0SOMETHING WRONG WITH INVERSE AND/OR SOLUTION OF NORMAL EQUATIONS
     ., PROGRAM WILL BE STOPPED AFTER ADJUSTEMENTS TO PARAMETERS'/
     .    ' DIAGNOSTIC CODES IERSOL,IER,IERCOD,IERINV:',4I6)
      endif
c
c*  start=9000
      return
      end
