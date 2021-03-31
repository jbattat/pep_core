      subroutine PRNSUM(ntotal,b,numpar,ncol)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 dum
      integer   i, ind, indsol, j, k, l, m, mind, ncol, nind, nm1, nnow,
     .          nsolut, ntodo, ntotal, numpar
 
c*** end of declarations inserted by spag
 
 
c
c subroutine prnsum to print summary for multiple parameter sets
c or multiple epochs
c p. macneil august, 1976
c maximum of 50 solutions
c
c
c common
      include 'aprtbf.inc'
      include 'fcntrl.inc'
      include 'filstf.inc'
      real*10 w(300,7)
      equivalence (w(1,1),G(1,1))
      include 'filtds.inc'
      include 'filtim.inc'
      include 'inodta.inc'
      include 'revptr.inc'
 
c external work area
      common/WRKCOM/ Names(2,1000),Dw(7),Nsol(50)
      real*10 Dw
      character*8 Names
      integer*4 Nsol
 
      real*10 b(ncol,numpar)
c        use of b(i,j) i=1,ncol (=2*ntotal+2), j=1,nparam
c        b(i,j) =
c                  xnom(j), i = 1
c                  scale(j), i = 2
c                  adjustment(j), i = 3,5,7,....
c                  sigma(j), i = 4,6,8,.....
      integer*2 n(2)
      real*10 oldval(2),sig(2),frac(2),adj(2),nval(2)
      character*4 qsol/' SOL'/, qnum/'  # '/
      character*8 name(2,2)
 
      if(Jct(52).le.0) return
      if(ntotal.le.0 .or. ntotal.gt.50)
     .     call SUICID(' NTOTAL OUT OF RANGE, STOP IN PRNSUM', 9)
c
c read in data from previous solutions
      nind = 0
      do i = 3, ncol, 2
         nind = nind + 1
         read(Ibuf6) Nsol(nind),(b(i,j),b(i+1,j),j = 1,Nparam)
      end do
      rewind Ibuf6
c
c read parameter names, nominals, scales
      read(Ibuf1) (Names(1,j),Names(2,j),j=1,Nparam)
      read(Ibuf1) (b(1,j),b(2,j),j=1,Nparam)
      rewind Ibuf1
c scale nominals, adjustments, sigmas, but not scales
c
      do j = 1, Nparam
         if(b(2,j).ne.1.0_10) then
            do i = 1, ncol
               if(i.ne.2) b(i,j) = b(i,j)*b(2,j)
            end do
         endif
      end do
c
c initialize for print loops
      nsolut = ntotal
      ntodo  = numpar
      Line   = 60
c
c accumulate w
      if(.not. (Iterat.le.1 .and. Jct(56).lt.2 .and.
     .    .not. Filflg(1) .and. Fict(8).le.0)) then
         rewind Iconof
 
c skip epochs
         read(Iconof)
         nm1 = Nepoch - 1
 
         do i = 1, Nepoch
            do j = 1, Npnp
               w(i,j) = 0.0_10
            end do
         end do
c
c read and sum forward w values
         do i = 1, nm1
            read(Iconof) (dum,j = 1,Npnp),(Dw(j),j = 1,Npnp)
            do j = 1, Npnp
               w(Nepoch,j) = w(Nepoch,j) + Dw(j)
            end do
         end do
c
c accumulate net backward w
         do i = 1, nm1
            l = nm1 + 1 - i
            read(Iconof) (dum,j = 1,Npnp),(Dw(j),j = 1,Npnp)
            do j = 1, Npnp
               w(l,j) = w(l + 1,j) + Dw(j)
            end do
         end do
         rewind Iconof
      endif
c
c loop through parameters
      do i = 1, Nparam, 2
 
c how many solutions per line (1 or 2)?
         nnow = 2
         if(ntodo.le.2) nnow = ntodo
         ntodo = ntodo - nnow
         if(Line + nsolut.ge.54) then
            call NEWPG
            if(Jct(51).gt.0) write(Iout,20)
   20       format('0MULTIPLE SOLUTION SUMMARY:', /)
            if(Ict(42).gt.0) write(Iout,40)
   40       format(
     .'0KALMAN FILTER SOLUTION SUMMARY   (NOTE: NOMINAL VALUES APPLY TOE
     .POCH 1)')
            Line = 4
         endif
c
c print parameter header
         do k = 1, nnow
            ind  = i + k - 1
            n(k) = ind
            name(1,k) = Names(1,ind)
            name(2,k) = Names(2,ind)
            oldval(k)  = b(1,ind)
         end do
         write(Iout,50) (qsol,n(k),name(1,k),name(2,k),oldval(k),
     .                   k = 1, nnow)
   50    format(2(a4,i7,2x,a8,1x,a8,' =',1pd23.15,11x))
c
c print column headings
         write(Iout,100) (qnum,k = 1,nnow)
  100    format(2(a4,2x,'ADJUSTMENT',6x,'NEW VALUE',14x,'SIGMA',8x,
     .          'FRACT',3x))
         Line = Line + 2
c
c loop through solutions
         do j = 1, nsolut
 
c find values to be printed
            indsol = Nsol(j)
            ind    = 2*(j - 1) + 3
            do k = 1, nnow
               l = i + k - 1
               adj(k)  = b(ind,l)
               nval(k) = oldval(k) + adj(k)
 
               if((Ict(42).gt.0) .and. (Jct(56).gt.1)) then
                  m = Rptr(l)
 
c process noise parameter, correct for offset of nominal
                  if(m.gt.0) nval(k) = nval(k) + w(indsol,m)
               endif
 
               mind    = ind + 1
               sig(k)  = b(mind,l)
               frac(k) = 1.0E19_10
               if(sig(k).ne.0.0_10) frac(k) = adj(k)/sig(k)
            end do
c
c print summary line
            call PAGCHK(58,1,0)
            write(Iout,120) (indsol,adj(k),nval(k),sig(k),frac(k),
     .                       k = 1, nnow)
  120       format(2(i3,'.',1pd16.8,d23.15,d11.3,0pf9.3,3x))
         end do
         write(Iout,120)
         Line = Line + 1
      end do
      return
      end
