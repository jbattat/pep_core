      subroutine SOLFRM(imats,scale,key)
 
      implicit none

c     m.e.ash   july 1970   subroutine solfrm
c     restore solution from saved solution of normal equations
c     also restore standard deviations
c
c     special option to skip restoring covariance and to allow an input
c     parameter set larger than saved set if ict(5)=4
c
c calling parameters
      integer*4 imats,key
      real*10 scale(1000)
c IMATS - Dataset containing saved solution and covariance.
c SCALE - Vector of scale factors applied to standard deviations.
c KEY   - Control switch passed to FRMHED to control reading of saved
c         headers from IMATS.
c         1=>read and print headers and adopt IC values
c         2=>read and print but don't adopt IC's
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'mcnfrm.inc'
      include 'nrmmat.inc'
      include 'numnum.inc'
      include 'restor.inc'
      include 'rtside.inc'

c shared external work area
      common/WRKCOM/ Iptr
      integer*2 Iptr(1000)
c
c local
      integer*4 i,iaprio,icnt,ippr,irow,ivectk,j,jr,
     .          jrow,length,lparam,mtst,nsav0,nseq1,ntape1
      real*4    erwgt1(2)
      integer*2 ierr(2)
      character*8    zerneg(2)/'ZERO','NEGATIVE'/
      real*10 ai
      integer*4 ia
      equivalence (ai,ia)
      character*8    tmesg(6)/'RESTORIN','G SOLUTI','ON AND S',
     .    'IGMAS OF',
     .          ' ORDER  ','*****   '/
 
      length = (Nparam*(Nparam+1))/2
      if(length.gt.Nrmsiz .or. (Ict(46).eq.2.and.length*2.gt.Nrmsiz))
     .     call SUICID('TOO MANY PARAMETERS, STOP IN SOLFRM ',9)
      if(Ict(46).ge.10 .and. 2*Nparam*Nparam.gt.Nrmsiz) call SUICID
     . ('TOO MANY PARAMETERS FOR RANK REDUCTION, STOP IN SOLFRM  ',14)
 
      call PAGSET('RESTORING SOLN. AND COVAR.  ',7)
c read records indicating the parameters which are in the
c saved solution and covariance matrix
      call FRMHED(imats,'IMATS','    ',key,ippr,1)
      if(ippr.ne.0) call SUICID(
     . 'PPR FORMAT NOT VALID, STOP IN SOLFRM    ',10)
      if(Ict(5).ne.4) then
         if(Mparam.ne.Nparam) call SUICID(
     .       ' MPARAM NOT EQUAL TO NPARAM, STOP IN SOLFRM ',11)
      endif
      Weight = 1.0_10
      lparam = max0(Mparam,Nparam)
      if(lparam.gt.999)
     .     call SUICID('MORE THAN 999 PARAMETERS, STOP IN SOLFRM',10)
c
c read first record of total series for saved solution and
c covariance matrix
      read(imats) Measmt,ntape1,nseq1,(erwgt1(i),i = 1,2),
     .            (Ermeas(i),i = 1,3),Sumzns,Sumzsm,iaprio,ivectk
      call PAGCHK(60,4,0)
      write(Iout,100) imats,ntape1,nseq1,Measmt,Ermeas
  100 format('-RESTORING SOLUTION AND STANDARD DEVIATIONS'/' IMATS=',
     .       i3,'  NTAPE=',i4,'  NSEQ=',i5,'  MEASMT=',i8,
     .       '  ERMEAS= (',1pd12.5,',',1pd12.5,',',1pd12.5,')')
c
c setup iptr vector
c mparam found in frmhed read number 2
      if(Ict(5).eq.4) then
         call ZFILL(Iptr,2*Mparam)
         call ZFILL(Sigma,16*Nparam)
      endif
      do i = 1,Mparam
         ia     = i
         Sav(i) = ai
      end do
      call FRMMVE(Sigma,Nparam)
      icnt = 0
      do i = 1,Nparam
         ai = Sigma(i)
         Sigma(i) = 0.0_10
         if(ia.eq.0) then
            if(Ict(5).eq.4) goto 200
            call SUICID('PARAMETER SETS NOT IDENTICAL, STOP IN SOLFRM',
     .                  11)
         endif
         Iptr(ia) = i
         icnt     = icnt + 1
  200 end do
      if(icnt.ne.Mparam) call SUICID(
     .'RESTORED SOLUTION HAS PARAMETERS NOT FOUND IN INPUT, WARNING IN S
     .OLFRM  ',-18)
c iptr(ia) is the location in side (or in b) of the
c element ia of the vector sav
      if(ivectk.gt.0) then
         read(imats) mtst,(Sav(i),i = 1,Mparam)
         if(mtst.ne.-1)
     .        call SUICID('INVALID ROW COUNTER, STOP SOLFRM',8)
         if(Ict(5).eq.4) call ZFILL(Vectk,16*Nparam)
         do i = 1,Mparam
            ia = Iptr(i)
            if(ia.gt.0) Vectk(ia) = Sav(i)
         end do
      endif
c
c read solution and scale vectors
      read(imats) Mparam,(Sav(i),i = 1,Mparam),
     .            (B(i),i = 1,Mparam)
c
c must clear arrays if missing parameters allowed
      if(Ict(5).eq.4) then
         call ZFILL(Sigma,16*Nparam)
         call ZFILL(Side,16*Nparam)
         call ZFILL(scale,16*Nparam)
      endif
c
c fix order of solution and scale vectors
      do i = 1,Mparam
         ia = Iptr(i)
         if(ia.gt.0) then
            Side(ia)  = Sav(i)*B(i)
            scale(ia) = B(i)
         endif
      end do
 
      if(Ict(5).eq.4) call ZFILL(B,16*length)
 
c increment saved row counter
      Nsav = 0
      do while( .true. )
         Nsav = Nsav + 1
         if(Nsav.gt.Mparam) then
c
c end of restoring quantities
            rewind imats
            Itrwnd(imats) = 0
c
c calculate sigmas (square roots of diagonal elements of
c covariance matrix)
            ierr(1) = 0
            ierr(2) = 0
            do j = 1,Nparam
               if(Sigma(j).lt.0) then
                  Sigma(j) = -SQRT(-Sigma(j))
                  jr = 2
               else if(Sigma(j).eq.0) then
                  jr = 1
                  if(Ict(5).eq.4) goto 220
               else
                  Sigma(j) = SQRT(Sigma(j))
                  goto 220
               endif
               call PAGCHK(58,1,0)
               write(Iout,210) j,zerneg(jr)
  210          format(' DIAGONAL ELEMENT',i4,' IS ',a8)
               ierr(jr) = ierr(jr) + 1
  220          Sigma(j) = Sigma(j)*scale(j)
            end do
c
c see if any diagonal elements were zero or negative
            if(ierr(1) + ierr(2).ne.0) then
               write(Iout,230) ierr
  230          format(' THERE ARE',i4,' ZERO AND',i4,
     .               ' NEGATIVE DIAGONAL ELEMENTS OF COVARIANCE MATRIX'
     .               /
     .       ' PROGRAM WILL BE STOPPED AFTER ADJUSTMENTS TO PARAMETERS'
     .      )
               Ict(1) = Iterat - 1
               Line   = Line + 2
            endif
c
c print out timer information
            call EBCDI(Mparam,tmesg(6),5)
            call TIMRIT(tmesg,12)
            return
         else
c
c read row of saved covariance matrix
            read(imats) nsav0,(Sav(i),i = 1,Mparam)
            if(nsav0.ne.Nsav) call SUICID(
     .          ' NSAV0 NOT EQUAL TO NSAV, STOP IN SOLFRM',10)
c
c restore row of saved covariance matrix
            irow = Iptr(Nsav)
            if(irow.ne.0) then
               Sigma(irow) = Sav(Nsav)
               jrow = (irow*(irow-1))/2
               do i = 1,Mparam
                  ia = Iptr(i)
                  if(ia.le.irow) then
                     ia    = ia + jrow
                     B(ia) = Sav(i)
                  endif
               end do
            endif
         endif
      end do
      end
