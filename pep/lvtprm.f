      subroutine LVTPRM(lprx,lpr,mpr,length)

      implicit none
c
c m.e.ash    oct 1967    subroutine lvtprm
c
c arguments
      integer*4 length
      integer*2 lprx(length),lpr(length),mpr(length)

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdat.inc'

c local
      integer   i,l,m

      l = 1
      m = 1
      do i = 1,length
         if(Iabs1.gt.0) then
            if(mpr(m).gt.0) then
               if(lpr(l).gt.0) then
                  if(mpr(m).ge.lpr(l)) then
                     if(mpr(m).eq.lpr(l)) then
                        m = m + 1
                     else
                        Lnotm = 1
                     endif
                     goto 50
                  endif
               endif
               lprx(i) = mpr(m)
               m = m + 1
               goto 100
            endif
         endif
         if(lpr(l).le.0) then
            lprx(i) = 0
            goto 100
         else
            Lnotm = 1
         endif
   50    lprx(i) = lpr(l)
         l = l + 1
  100 end do
      return
      end
