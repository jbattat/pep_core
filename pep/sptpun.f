      subroutine SPTPUN(ncard)
 
      implicit none
c
c k.m.becker   june 1968   subroutine sptpun
c punch adjusted spot coordinates
c

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'inodta.inc'
      include 'sptcrd.inc'

c local
      integer   i,istr,j,jd,k,ncard,nsv
 
      if(Numspt.le.0) return
      istr = 0
      do j = 1,Numspt
         do k = 1,6
            if(Lspcrd(k,j).ne.0) then
               if(istr.eq.0) then
                  istr = 1
                  write(Ipunch,10)
   10             format('*SPOTS')
                  ncard = ncard + 1
               endif
               nsv=k
               do i=k+1,6
                  if(Lspcrd(i,j).gt.0 .or. Spcord(i,j).ne.0._10) nsv=i
               end do
c if a future upgrade allows for cylindrical coordinates, as with sites,
c then this punch statement must choose a format (f16.9 vs f16.11) for
c the third coordinate accordingly
               if(nsv.le.3) then
                  write(Ipunch,20) Spot(j),Nsplnt(j),
     .             (Spcord(i,j),i = 1,3),(Lspcrd(i,j),i = 1,3)
   20             format(1a4,i3,1x,f18.11,2f18.13,8x,3i2,2x,'##')
                  ncard = ncard + 1
               else
                  jd=T0spot(j)
                  write(Ipunch,21) Spot(j),Nsplnt(j),
     .             (Spcord(i,j),i=1,3),(Lspcrd(i,j),i=1,3),
     .             Spot(j),(Spcord(i,j),i=4,6),jd,(Lspcrd(i,j),i=4,6)
   21             format(1a4,i3,1x,f18.11,2f18.13,' 6',6x,3i2,2x,'##'/
     .                a4,4x,3f18.11,i8,3i2,2x,'##')
                  ncard = ncard + 1
               endif
               goto 100
            endif
         end do
  100 end do
 
      return
      end
