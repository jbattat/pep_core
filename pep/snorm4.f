      subroutine SNORM4(pvcol,pnames,b,sigma,nparam,nocorr,print,ncorp,
     .                  iicorr,ibuf1,ier)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 adia, apvcol, hist
      integer   i, i0, ibuf1, idx, ier, ii, ij, iok, ix, j, jj,
     .          k, kk, ksml, kz, ldj, linexp, mlrg, mlrg1
      integer   nlrg, nlsv, nmov, nparam
 
c*** end of declarations inserted by spag
 
 
c m.ash   march 1972  subroutine snorm4
c calculate standard deviations and correlations
c printout correlation matrix
c read normal equations back into storage
      real*10 pvcol(2000),b(1),sigma(nparam)
      character*8 pnames(2000)
      logical*4 nocorr, print
      integer*2 ncorp, iicorr(1000)
c N O T E: pvcol and pnames arrays are expected to share storage
c  pvcol - work array
c  pnames- character work array
c  b     - input: covariance, returned as correlations
c  sigma - returned vector of std. devs.
c  nparam- number of equations
c  nocorr- if true, then skip correlation matrix
c  print - if true, then print normally, else, suppress all but errors
c  ncorp - extra number of correlation matrix rows to print
c  iicorr- array of matrix row numbers (see ncorp)
c  ibuf1 - if nonzero, dataset number to read parameter names from
c  ier   - returned count of errors
c
      include 'inodta.inc'
 
      real*4    btst
      real*4    blrg(50)
      integer*2 ilrg(50),jlrg(50)
      character*8 qis/' IS '/, qcont/' (CONT) '/
      character*72 phdr/'THE LOWER DIAGONAL HALF OF THE CORRELATION MATR
     .IX OF ORDER              '/
c
c dimension of ihist must be.gt.nhist
      integer*4 ihist(21),nhist/20/
c
c* start=200
c calculate square root of diagonal elements of covariance
c matrix (standard deviations of parameters estimates)
      ier = 0
      ldj = 0
      do i = 1, nparam
         ldj = ldj + i
         if(b(ldj).ge.0._10) then
            sigma(i) = SQRT(b(ldj))
         else
            if(ier.le.0) then
               call PAGSET('ERROR MESSAGES FOR THE COVARIANCE MATRIX',
     .                     -10)
               call PAGE(0,1)
            endif
            call PAGCHK(60,1,1)
            write(Iout,20) i
   20       format(' DIAGONAL ELEMENT ', i3, ' IS NEGATIVE')
            sigma(i) = 0.0_10
            ier = ier + 1
         endif
      end do
c
c check  if any diagonal elements negative
c
      if(ier.eq.0) then
c
c*  start=400
c calculate and print correlation matrix
c
         hist = nhist
         do i = 1, nhist
            ihist(i) = 0
         end do
         ihist(nhist + 1) = 0
         do i = 1, 50
            blrg(i) = 0.0
            ilrg(i) = 0
            jlrg(i) = 0
         end do
         nlrg = 0
         mlrg = 0
         if(print .and. .not.nocorr) then
            call EBCDIX(nparam,phdr,61,4)
            call MVC(qis,1,8,phdr,65)
            call PAGSET(phdr,-18)
            call PAGSET(-1,Jout)
            call PAGE(0,1)
            call MVC(qcont,1,8,phdr,65)
            call PAGSET(phdr,-18)
            call PAGSET(-1,Jout)
         endif
         ii = 0
         do i = 1, nparam
            adia     = sigma(i)
            i0       = ii
            ii       = ii + i
            pvcol(i) = b(ii)
            linexp   = 2 + (i - 1)/16
            jj       = 0
            do j = 1, i
               ij = i0 + j
               jj = jj + j
c special treatment for correlations of unity
c (because of sqrt truncation)
               if(j.ne.i) then
                  if(pvcol(i)*pvcol(j).ne.b(ij)*b(ij)) then
                     b(ij)  = b(ij)/(adia*sigma(j))
                     apvcol = ABS(b(ij))
                     if(apvcol.le.1._10) then
                        idx = apvcol*hist + 1
                     else
 
c it should always skip this
                        call PAGCHK(59,1,1)
                        write(Iout,30) i,j,b(ij)
   30                   format(' CORRELATION(', i3, i4, ') =', f25.17,
     .                         ' OUT OF RANGE')
                        if(Jout.gt.0) write(Jout,30) i,j,b(ij)
 
                        idx = nhist + 1
                     endif
                     goto 40
                  endif
               endif
               b(ij) = 1._10
               if(j.eq.i) goto 60
               apvcol     = 1._10
               idx        = apvcol*hist + 1
   40          ihist(idx) = ihist(idx) + 1
               if(.not. print) goto 100
               if(apvcol.gt.0.5_10) then
                  mlrg = mlrg + 1
                  if(apvcol.ne.1._10) then
                     btst = apvcol
                     nlsv = nlrg
                     if(nlrg.lt.50) nlrg = nlrg + 1
                     if(nlsv.gt.0) then
 
c find slot for inserting this correlation
                        do ksml = 1, nlsv
                           if(btst.gt.abs(blrg(ksml))) then
 
c shift lower values (if any)
                              nmov = nlrg - ksml
                              if(nmov.gt.0) then
                                 kz = nlrg
                                 do kk = 1, nmov
                                    blrg(kz) = blrg(kz - 1)
                                    ilrg(kz) = ilrg(kz - 1)
                                    jlrg(kz) = jlrg(kz - 1)
                                    kz = kz - 1
                                 end do
                              endif
                              goto 50
                           endif
                        end do
                     endif
 
c none smaller, insert at end
                     ksml = nlsv + 1
                     if(ksml.gt.nlrg) goto 60
   50                blrg(ksml) = b(ij)
                     ilrg(ksml) = i
                     jlrg(ksml) = j
                  endif
               endif
 
c*  start=600
   60       end do
 
            if(print .and. .not.nocorr) then
               call PAGCHK(60,linexp,1)
               write(Iout,70) i,(b(i0+k),k = 1,i)
   70          format('0', i3, (t5,16F8.4))
               if(Jout.gt.0) write(Jout,70) i,(b(i0+k),k = 1,i)
            endif
  100    end do
c
c*  start=700
c extra print
         if(.not. print) return
         if(ncorp.gt.0) then
            call PAGSET('SELECTED ROWS OF THE CORRELATION MATRIX ', -10)
            call PAGSET(-1,Jout)
            call PAGHED(0)
            do k = 1, ncorp
               iok = 0
               i   = iicorr(k)
               ij  = i*(i - 1)/2
               ix  = 1
               do j = 1, nparam
                  ij = ij + ix
                  pvcol(j) = b(ij)
                  if(b(ij).ne.0._10 .and. j.ne.i) iok = 1
                  if(j.ge.i) ix = j
               end do
               if(iok.ne.0) then
                  call PAGCHK(60,linexp,1)
                  write(Iout,70) i,(pvcol(j),j=1,nparam)
                  if(Jout.gt.0) write(Jout,70) i,(pvcol(j),j=1,nparam)
               endif
            end do
         endif
         if(nlrg.gt.0) then
c
c*  start=800
c print out large correlations
            call PAGCHK(60,nlrg + 8,0)
            mlrg1 = (nparam*(nparam-1))/2
 
c read parameter names from ibuf1
            call ZFILL(pnames,16*nparam)
            if(ibuf1.gt.0) then
               rewind ibuf1
               j = 2*nparam
               read(ibuf1) (pnames(i),i=1,j)
               rewind ibuf1
            endif
            write(Iout,120) mlrg,mlrg1,nlrg,
     .                       (i,jlrg(i),ilrg(i),blrg(i),
     .                       pnames(2*jlrg(i)-1),pnames(2*jlrg(i)),
     .                       pnames(2*ilrg(i)-1),pnames(2*ilrg(i)),
     .                       i = 1, nlrg)
            if(Jout.gt.0) write(Jout,120) mlrg,mlrg1,nlrg,
     .                          (i,jlrg(i),ilrg(i),blrg(i),
     .                          pnames(2*jlrg(i)-1),pnames(2*jlrg(i)),
     .                          pnames(2*ilrg(i)-1),pnames(2*ilrg(i)),
     .                          i = 1, nlrg)
  120       format(/i8,' OF THE', i8,
     .             ' CORRELATIONS ARE GREATER THAN 0.5'/'-THE LARGEST',
     .             i3, ' CORRELATIONS IN ABSOLUTE VALUE ARE'/(i4,'. (',
     .             i3,',',i3,') =',f14.10,3x,2(a8,1x),' AND ',2(1x,a8)))
            write(Iout,140) ihist(nhist + 1),(ihist(i),i = 1,nhist)
            if(Jout.gt.0) write(Jout,140) ihist(nhist + 1),
     .                          (ihist(i),i = 1,nhist)
  140       format('0THE DISTRIBUTION OF THE ABSOLUTE VALUES OF THE ',
     .             'CORRELATIONS WITH ', i3,
     .             ' OUT OF THE RANGE 0.0 TO 1.0 IS', /(4(2x,5I6)))
         else
            call PAGCHK(60,2,0)
            write(Iout,160)
  160       format('0 ALL CORRELATIONS ARE LESS THAN 0.5')
            if(Jout.gt.0) write(Jout,160)
         endif
      else
         call PAGCHK(60,3,0)
         write(Iout,200) ier
  200    format('0A TOTAL OF', i4, ' DIAGONAL ELEMENTS ARE NEGATIVE'/
     .          ' CORRELATION MATRIX CANNOT BE CALCULATED')
      endif
c
c*  start=1000
c
c print square root of diagonal elements of inverse matrix
      linexp = 3 + (nparam - 1)/8
      call PAGSET(-1,0)
      call PAGCHK(60,linexp,0)
      write(Iout,300) nparam,(sigma(i),i = 1,nparam)
  300 format('0THE', i4,
     .       ' STANDARD DEVIATIONS OF THE PARAMETER ESTIMATES ARE', /,
     .       (4x,1p,8D16.8))
 
      call TIMRIT(
     .'  CALCULATING & WRITING OUT CORRELATION MATRIX & STD DEVIATIONS '
     ., 16)
c
c*  start=2000
c
      return
      end
