      integer function JULDAY(imonth,iday,iyear)
 
      implicit none
 
c
c    determination of julian day number
c from given month,day and year since 1900 (gregorian calendar)
c valid from 1 AD onward
c
      integer*2 imonth,iday,iyear,iyr,iyr4,iyr100,iyr400
      integer*2 montot(12)/0,31,59,90,120,151,181,212,243,273,304,334/
c imonth=month           (between 1 and 12)
c iday  =day of month    (between 1 and 31)
c iyear =year since 1900 (negative before 1900)
c
      if(imonth.le.0 .or. imonth.gt.12)
     .    call SUICID(' MONTH INCORRECT, STOP IN JULDAY', 8)
      iyr = iyear+1900
      if(imonth.le.2) iyr = iyr - 1
      iyr4  = iyr/4
      iyr100=iyr/100
      iyr400=iyr/400
  100 JULDAY = (2415020 - 1899/4 + 1899/100 - 1899/400 + 365*iyear)
     .         + (montot(imonth) + iday + iyr4 - iyr100 + iyr400)
      return
      end
