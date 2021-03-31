      subroutine ADDINV(h,hinv,g,n,np)
 
      implicit none
c
c d. white  subroutine addinv  december 1974
c
c setup for inversion (scaling) and cleanup
c
c parameters
      integer*4 n,np
      real*10   h(n,n),hinv(n,n),g(np,n)
c
c common
      include 'fcntrl.inc'
      include 'filtim.inc'
      include 'inodta.inc'
c
c local
      integer*4 fict4,i,indic,j
      real*10   rmax,rterm,cmax,cterm,d
      real*10   seps/1.0e-55_10/
      integer*4 me/-1/
c
c see if any scaling or printing
      fict4 = Fict(4)
      if(fict4.ne.0) then
c
c print unscaled h
         if(mod(fict4,2).ne.0)
     .        call XPMPO(n,hinv,n,0,'H UNSCALED',10)
c
c make scale factors
c g(i,1) has row scale factor for row i
c g(i,2) has column scale factor for column i
         if(fict4.ge.2) then
            do i = 1,n
               rmax = abs(h(i,1))
               cmax = abs(h(1,i))
               do j = 2,n
                  rterm = abs(h(i,j))
                  cterm = abs(h(j,i))
                  if(rterm.gt.rmax) rmax = rterm
                  if(cterm.gt.cmax) cmax = cterm
               end do
               g(i,1) = sqrt(rmax)
               g(i,2) = sqrt(cmax)
            end do
c
c print out scale factors
            if(fict4.ge.3) then
               call XPMPO(n,g(1,1),0,0,'ROW SCALE FACTORS',17)
               call XPMPO(n,g(1,2),0,0,'COLUMN SCALE FACTORS',20)
            endif
c
c scale hinv (and h, holding for cleanup)
            do j = 1,n
               cterm = g(j,2)
               do i = 1,n
                  hinv(i,j) = hinv(i,j)/(g(i,1)*cterm)
                  h(i,j)    = hinv(i,j)
               end do
            end do
c
c print out scaled h
            if(fict4.ge.3) call XPMPO(n,hinv,n,0,'H SCALED',8)
         endif
      endif
c
c call to dpinv replaced by
      call ECROUT(h,hinv,d,seps,n,me,n,n,indic)
      write(Iout,100) indic
  100 format(' INDIC =',i2)
      if(d.ne.0.0) then
c
c call for cleanup (g(*,3), g(*,4), g(*,5) work)
c run with fict(3) = -1 to avoid error print
         if(Fict(3).ge.0) call ADDCLN(h,hinv,g(1,3),g(1,4),g(1,5),n)
c
c print out scaled inverse
         if(fict4.ge.3) call XPMPO(n,hinv,n,0,'HINV SCALED',11)
      else
         write(Iout,150)
  150    format('-*** BAD MATRIX INVERSION IN ADDINV ***')
         Line = 60
      endif
c
c unscale hinv (reversing sense of scale factors)
      if(fict4.ge.2) then
         do j = 1,n
            cterm = g(j,1)
            do i = 1,n
               hinv(i,j) = hinv(i,j)/(g(i,2)*cterm)
            end do
         end do
      endif
c
c print out unscaled inverse
      if(mod(fict4,2).ne.0)
     .     call XPMPO(n,hinv,n,0,'HINV UNSCALED',13)
      return
      end
