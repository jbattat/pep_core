      subroutine TODAY(date_work)
 
      implicit    none
      character*20 date_work
      integer*4 month, day, year, dayx(3), timx(3), dat1, i
      integer*2 yrpak, datpak
      equivalence (day,dayx(1)),(month,dayx(2)),(year,dayx(3))
      integer*4 daytot(12)/0,31,59,90,120,151,181,212,242,273,303,334/
 
c      call idate(month, day, year)
      call idate(dayx)
      call itime(timx)
      if(year.gt.199) year=mod(year,100)
      dat1=daytot(month)+day
      if((year/4)*4.eq.year .and. month.gt.2) dat1=dat1+1
      datpak=(dat1/100)*4096 +mod(dat1/10,10)*256 +mod(dat1,10)*16 +15
      yrpak=mod(year/10,10)*16 +mod(year,10)
      write(date_work,100) month, day, year, yrpak, datpak, timx
  100 format(i2, '/', i2, '/', i2, 2a2, 2(i2,':'), i2)
      do i=1,20
         if((i.le.8 .or. i.ge.13) .and. 
     .    date_work(i:i).eq.' ') date_work(i:i)='0'
      end do
      return
      end
