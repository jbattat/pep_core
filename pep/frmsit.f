      subroutine FRMSIT(rst,ntop)
 
      implicit none
c
c m.e. ash february 1970 subroutine frmsit
c observing site logic routine for forming normal equations
c from saved equations
c
c
c arguments
      real*10 rst(1000)
      integer*4 ntop

c array dimensions
      include 'globdefs.inc'

c commons
      include 'anctrl.inc'
      include 'restor.inc'
      include 'stcord.inc'
      include 'stcrdm.inc'

c local
      integer   i,j,k

      do k = 1, Mumsit
         do j = 1, Numsit
            if(Site1(1,k).eq.Site(1,j)) then
               Nrst = Lscrd1(j) - 1
               call FRMBDY(rst,Lscrd(1,j),Mscrd(1,k),-6,ntop)
               goto 50
            endif
         end do

         do i = 1,6
            if(Mscrd(i,k).gt.0) Nsav = Nsav + 1
         end do
   50 end do
 
      return
      end
