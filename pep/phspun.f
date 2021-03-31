      subroutine PHSPUN(ncard)
 
      implicit none
c
c k.m.becker   june 1968   subroutine phspun
c punch adjusted optical observation phase corrections
c
c arguments
      integer*4 ncard

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'inodta.inc'
      include 'phase.inc'

c local
      integer   i,istr,j,k,nc
 
      if(Numphs.le.0) return
      istr = 0
      do j = 1,Numphs
         do k = 1,9
            if(Lphs(k,j).ne.0) then
               if(istr.eq.0) then
                  istr = 1
                  write(Ipunch,10)
   10             format('*PHASES')
                  ncard = ncard + 1
               endif
               write(Ipunch,20) Phsit(j),Phser(j),Nplphs(j),
     .          Ncphs(j),(Lphs(i,j),i = 1,9),(Aphase(i,j),i = 1,3)
   20          format(a4,1x,a4,2I3,9I2,3x,1p,3E12.5)
               ncard = ncard + 1
               if(Ncphs(j).gt.3) then
                  nc = Ncphs(j)
                  write(Ipunch,30) (Aphase(i,j),i = 4,nc)
   30             format(1p,6E12.5)
                  ncard = ncard + 1
               endif
               goto 100
            endif
         end do
  100 end do
 
      return
      end
