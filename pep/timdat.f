      subroutine TIMDAT(date)
 
      implicit    none
c
c input: date - if the first 4 bytes are binary zero, then fill the
c               string with the current date as "mm/dd/yy"
c output: revision of values in /TIMSTF/ common as follows
c
c Dat1  = packed date for this reading and setting.
c Ireal1= real time (binary integer centiseconds) for this call.
c         Constrained to increase by at least 1 on each call.
c Itask = increment of task interval timer counted from when TIMDAT was
c         previously called (units of centiseconds).
c Idayr = day of year in a format (only calculated if date(1)
c         has value 0 when sent to subroutine timdat).
c
      character*8 date, tdate
      real*4      date4(2)
      equivalence (tdate,date4)
      integer*4   month, day, year, dayx(3), timx(3)
      equivalence (day,dayx(1)),(month,dayx(2)),(year,dayx(3))
      real*4      tarray(2), dtime
      integer*4   oreal/0/,
     1            daytot(12)/0,31,59,90,120,151,181,212,242,273,303,334/
      include 'timstf.inc'
 
c      call idate(month, day, year)
      call idate(dayx)
      call itime(timx)
      tdate=date
      if(date4(1).eq.0.0) write(date,100) month, day, mod(year,100)
  100 format(i2, '/', i2, '/', i2)
 
      Itask = dtime(tarray)*100.
 
      oreal = mod(oreal+1+Itask,8640000)
      Ireal1=(((timx(1)*60)+timx(2))*60+timx(3))*100
      if(oreal.gt.Ireal1 .and. oreal.lt.Ireal1+4320000) Ireal1=oreal
      oreal = Ireal1
 
      Dat1=daytot(month)+day
      if((year/4)*4.eq.year .and. month.gt.2) Dat1=Dat1+1
      if(date4(1).eq.0.0) write(Idayr,200) Dat1
  200 format(i3)
      Dat1=Dat1+year*1000
 
      return
      end
