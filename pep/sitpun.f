      subroutine SITPUN(ncard)
 
      implicit none
 
c
c k.m.becker   june 1968   subroutine sitpun
c punch adjusted observing site coordinates
c
c
c argument
      integer*4 ncard

c array dimensions
      include 'globdefs.inc'

c common
      include 'inodta.inc'
      include 'stcord.inc'

c local
      integer*4 i,istr,j,jd,k,nsv
 
      if(Numsit.le.0) return
      istr = 0
      do j=1,Numsit
         do k=1,6
            if(Lscrd(k,j).ne.0) then
               if(istr.eq.0) then
                  istr = 1
                  write(Ipunch,10)
   10             format('*SITES')
                  ncard = ncard + 1
               endif
               nsv=k
               do i=k+1,6
                  if(Lscrd(i,j).gt.0 .or. Scord(i,j).ne.0._10) nsv=i
               end do
c longitude can have 3 integer digits plus a sign, and so can squeeze
c 11 decimal digits into the 16-byte field.  same for latitude.
c radii on planets other than earth might not fit and so the distance
c coordinates get only 9 decimal digits.
               if(nsv.le.3) then
                  if(Kscrd(j).eq.0) then
                     write(Ipunch,20) (Site(i,j),i=1,2),
     .                (Scord(i,j),i=1,3),(Lscrd(i,j),i=1,3),Kscrd(j)
   20                format(2a4,f18.11,2f18.13,8x,4i2,'##')
                  else
                     write(Ipunch,21) (Site(i,j),i=1,2),
     .                (Scord(i,j),i=1,3),(Lscrd(i,j),i=1,3),Kscrd(j)
   21                format(2a4,f18.11,f18.13,f18.11,8x,4i2,'##')
                  endif
               else
                  jd=T0site(j)
                  if(Kscrd(j).eq.0) then
                     write(Ipunch,25) (Site(i,j),i=1,2),
     .                (Scord(i,j),i=1,3),(Lscrd(i,j),i=1,3),Kscrd(j),
     .                (Scord(i,j),i=4,6),jd,(Lscrd(i,j),i=4,6)
   25                format(2a4,f18.11,2f18.13,' 6',6x,4i2,'##'/
     .                '...',5x,3f18.11,i8,3i2,2x,'##')
                  else
                     write(Ipunch,26) (Site(i,j),i=1,2),
     .                (Scord(i,j),i=1,3),(Lscrd(i,j),i=1,3),Kscrd(j),
     .                (Scord(i,j),i=4,6),jd,(Lscrd(i,j),i=4,6)
   26                format(2a4,f18.11,f18.13,f18.11,' 6',6x,4i2,'##'/
     .                '...',5x,3f18.11,i8,3i2,2x,'##')
                  endif
                  ncard = ncard + 1
               endif
               ncard = ncard + 1
               goto 100
            endif
         end do
  100 end do
 
      return
      end
