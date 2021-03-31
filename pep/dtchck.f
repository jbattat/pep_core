      subroutine DTCHCK(jddtm0, mumdt, mumdt1, dt, mdt)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ioff3, ioff5
 
c*** end of declarations inserted by spag
 
 
c subroutine dtchck  r.w.king and j.f.chandler  july 1980
c check for inconsitency between dt parmaeters on osblib or
c sne data sets, and current array sizes; and fix
c
      real*4    dt(600)
      integer*4 jddtm0
      integer*2 mdt(600), mumdt, mumdt1, imdt(6)
c
c initialize imdt corrections
      ioff3 = 0
      ioff5 = 0
 
      if( jddtm0 .le. 0 ) then
         if( mumdt1 .gt. 200 .and. mumdt1 .le. 300 ) then
            ioff3 = 100
            ioff5 = 200
            do i = 1, mumdt
               dt(i + 400)  = dt(i + 200)
               dt(i + 200)  = dt(i + 100)
               dt(i + 100)  = 0.
               mdt(i + 400) = mdt(i + 200)
               mdt(i + 200) = mdt(i + 100)
               mdt(i + 100) = 0
            end do
            mumdt1 = mumdt1 + 200
         endif
      endif
      return
 
      entry DTCKI(imdt)
c correct imdt's if dt arrays have been expanded
c
      if( ioff3 .ne. 0 ) then
         if( imdt(3) .gt. 0 ) imdt(3) = imdt(3) + ioff3
         if( imdt(4) .gt. 0 ) imdt(4) = imdt(4) + ioff3
         if( imdt(5) .gt. 0 ) imdt(5) = imdt(5) + ioff5
         if( imdt(6) .gt. 0 ) imdt(6) = imdt(6) + ioff5
      endif
 
      return
      end
