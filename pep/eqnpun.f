      subroutine EQNPUN(ncard)
 
      implicit none
c
c k.m.becker   june 1968   subroutine eqnpun
c punch adjusted equinox-equator-latitude corrections
c
c arguments
      integer*4 ncard

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'eqenox.inc'
      include 'inodta.inc'

c local
      integer*4 i,istr,j,k

      if(Numeqn.le.0) return
      istr = 0
      do j = 1,Numeqn
         do k = 1,3
            if(Leqn(k,j).ne.0) then
               if(istr.eq.0) then
                  istr = 1
                  write(Ipunch,10)
   10             format('*EECORR')
                  ncard = ncard + 1
               endif
               write(Ipunch,20) Eqnsit(j),Eqnser(j),
     .                           (Leqn(i,j),i = 1,3),Denox(j),
     .                           Dequat(j),Dlat(j)
   20          format(a4,1x,a4,6x,3(1x,i1),15x,1p,3E12.5)
               ncard = ncard + 1
               goto 100
            endif
         end do
  100 end do
 
      return
      end
