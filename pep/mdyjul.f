      subroutine MDYJUL(imonth, iday, iyear, itime, jds)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ic, jd, jds, nyr
 
c*** end of declarations inserted by spag
 
 
c
c t.forni  may   1968   subroutine mdyjul
c obtain month, day, year  for given julian day jds
c
 
c base year is  1600      (2305447)
      integer*2 mdn(13)/0, 31, 59, 90, 120, 151, 181, 212, 243, 273,
     .          304, 334, 365/
      integer*2 imonth, iday, iyear, itime
 
c jd= days since 0 january 1600
      jd  = jds - 2305447
      nyr = jd/365
      do while( .true. )
 
c ic= number of centuries since 0 january 1600
         ic = nyr/100
 
c days due to leap years
         iday = jd - nyr*365 - (nyr - 1)/4 + (nyr + 99)
     .          /100 - (nyr + 399)/400 - 1
         if( ic .eq. 0 ) then
            if( nyr .eq. 0 ) iday = iday + 1
         endif
         if( iday .gt. 0 ) then
 
c iyear=  (0  thru  99)   year of the century
            iyear = nyr - ic*100
            itime = ic - 3
            nyr   = iyear
            if( nyr .ne. 0 ) then
               if( mod(nyr,4) .ne. 0 ) go to 200
            else if( mod(ic,4) .ne. 0 ) then
               go to 200
            endif
            if( iday .lt. 60 ) then
            else if( iday .eq. 60 ) then
               go to 300
            else
               iday = iday - 1
            endif
         else
            nyr = nyr - 1
            go to 100
         endif
         go to 200
  100 end do
  200 do i = 2, 13
         if( iday .le. mdn(i) ) then
            imonth = i - 1
            iday   = iday - mdn(imonth)
            return
         endif
      end do
      call SUICID(' STOP IN MDYJUL   36', 5)
  300 imonth = 2
      iday   = 29
      return
      end
