      subroutine RBSPUN(ncard)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, istr, j, k, ncard
 
c*** end of declarations inserted by spag
 
 
c
c k.m.becker   june 1968   subroutine rbspun
c punch adjusted radar observation biases
c

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'inodta.inc'
      include 'rdbias.inc'
 
      if(Numrbs.le.0) return
      istr = 0
      do j = 1, Numrbs
         do k = 1, 2
            if(Lrbs(k,j).ne.0) then
               if(istr.eq.0) then
                  istr = 1
                  write(Ipunch,10)
   10             format('*BIASES')
                  ncard = ncard + 1
               endif
               write(Ipunch,20) (Rdbsit(i,j),i = 1,2),Rdbser(j),
     .                           Nplrbs(j),(Lrbs(i,j),i = 1,2),
     .                           (Rbias(i,j),i = 1,2)
   20          format(a4,1x,a4,1x,a4,i3,4x,2I2,1p,2E12.5)
               ncard = ncard + 1
               goto 100
            endif
         end do
  100 end do
 
      return
      end
