      subroutine FRMSTR(rst,ntop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, k
 
c*** end of declarations inserted by spag
 
c  subr. frmstr - j.f.chandler - 1983 feb
c  sky correction (star catalog error) logic routine for forming
c  normal equations from saved normal equations.
c
c arguments
      real*10 rst(1000)
      integer*4 ntop

c array dimensions
      include 'globdefs.inc'

c        common
      include 'anctrl.inc'
      include 'restor.inc'
      include 'skystf.inc'
      include 'skystm.inc'
 
      do k = 1, Mumstr
         do j = 1, Numstr
            if(Ctlgnm(j).eq.Ctlgn1(k)) then

c found match, copy pointers
               Nrst = Lsky1(j) - 1
               call FRMBDY(rst,Lskycf(1,j),Msky(1,k),-80,ntop)
               goto 50
            endif
         end do
c no match, update saved equation counters
         do i = 1, 80
            if(Msky(i,k).gt.0) Nsav = Nsav + 1
         end do

   50 end do
 
      return
      end
