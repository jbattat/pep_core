      subroutine LVTBDY(lbdx,lbd,mbd,length)
 
      implicit none
c
c m.e.ash    oct 1967    subroutine lvtbdy
c
c arguments
      integer*2 lbdx(7),lbd(7),mbd(7)
      integer*4 length

c array dimensions
      include 'globdefs.inc'
c common
      include 'comdat.inc'

c local
      integer   i,len

      len = 6
      if(length.le.0) len = iabs(length)
      do i = 1,len
         lbdx(i) = 0
         if(Iabs1.gt.0) then
            if(mbd(i).gt.0) goto 50
         endif
         if(lbd(i).le.0) goto 100
         Lnotm   = 1
   50    lbdx(i) = 1
  100 end do
      len = length - len
      if(len.gt.0) call LVTPRM(lbdx(7),lbd(7),mbd(7),len)
      return
      end
