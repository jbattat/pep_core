      subroutine TIMRIT(title, len)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      f
      integer   i, lem, len
 
c*** end of declarations inserted by spag
 
 
c        j.f.chandler  1977 mar
c        print title and elapsed times (real,task) in standard format
c        'title' is character string ending with at least 1 blank
c        'len' is length of 'title' in words
c
c        common
      include 'inodta.inc'
      include 'timstf.inc'
c
c local
      character*4 req(2)/'REQU','IRED'/, rlts(2)/'REAL','TASK'/,
     .          title(*)
      integer*4 ih(2), im(2), itot(2), itotc(2)
      real*4    sec(2)
      character*4 gotdat(2)/'GOT ','DATE'/
 
      call TIMSET(ih, im, sec, itot, gotdat)
      go to 100
 
      entry TIMRTC(title, len, itotc)
c print time elapsed since some instant other than the last
c call to timset -- saved time values passed in itotc
 
      call TIMSET(ih, im, sec, itot, gotdat)
 
c get time since saved instant
      itot(1) = Ireal0 - itotc(1)
      if(itot(1).lt.0) itot(1) = itot(1) + 8640000
      itot(2) = Itotsk - itotc(2)
      call CONVRT(ih(1), im(1), sec(1), itot(1))
      call CONVRT(ih(2), im(2), sec(2), itot(2))
 
  100 lem = min0(iabs(len), 31)
      call PAGCHK(60, 5, 0)
      write(Iout, 200) (title(i), i = 1, lem), req
  200 format('0', 33A4)
      write(Iout, 300) (ih(i), im(i), sec(i), rlts(i), i = 1, 2)
  300 format(i12, 'H', i3, 'M', f6.2, 'S ', a4, ' TIME')
      f = 999999.
      if(itot(1).gt.0) f = float(itot(2))/float(itot(1))
      write(Iout, 400) f
  400 format('  (TASK TIME)/(REAL TIME)=', f8.5)
      return
      end
