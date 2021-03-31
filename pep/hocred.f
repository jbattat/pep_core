      subroutine HOCRED(jd,fract)
      implicit none

c subroutine hocred - j.f.chandler - 2012/11/29
c read external file of corrections to lunar libration

c on entry: jd,fract have the desired date/time, or, if jd=0, then
c           just initialize
c on exit:  Ihoc=0 indicates there is a correction in HOCCOM common
c           Ihoc=1 indicates no correction

c format of external file
c  line 1 - free-form title

c  line 2
c Jd1hoc - ephemeris start epoch (pep convention)           i7
c Jd2hoc - ephemeris end epoch (pep convention)             i8
c  the ephemeris must include a block covering the end of day Jd2hoc,
c  even if the correction is zero

c  subsequent lines
c Thoca - time block start epoch (noon convention)          f10.2  
c Thocb - time block end epoch (noon convention)            f11.2  
c Hocx  - ad hoc correction to longitude libration (nrad)   f8.2, 8x
c Hocy  - ad hoc correction to latitude libration (nrad)    f8.2  

c  last line (if necessary)
c Jd2hoc + .5 and Jd2hoc + .5

c array dimensions
      include 'globdefs.inc'

c arguments
      integer*4 jd
      real*10 fract

c commons
      include 'comdateq.inc'
      include 'fcntrl.inc'
      include 'hoccom.inc'
      include 'inodta.inc'

c local variables
      character*72 hoched
      real*10 t

      Ihoc=1
      if(jd.eq.0) then
         if(Libhoc.le.0) return
         if(Itrwnd(Libhoc).gt.0) return
         read(Libhoc,100) hoched,Jd1hoc,Jd2hoc
  100    format(a72/i7,i8)
         Itrwnd(Libhoc)=1
         Thoca=Jd1hoc-1
         Thocb=Thoca
         if(Line.gt.55) call NEWPG
         write(Iout,120) hoched,Jd1hoc,Jd2hoc
  120    format('0APPLYING LIBRATIONS: ',a72,' ('i7,'-',i7')')
         Line=Line+2
         return
      endif

      if(jd.lt.Jd1hoc .or. jd.gt.Jd2hoc) return
      t=jd+fract-0.5
c start at the beginning if necessary
      if(t.lt.Thoca) then
         rewind Libhoc
         read(Libhoc,100) hoched
         Thoca=Jd1hoc-1
         Thocb=Thoca
      endif
c look ahead in the file if necessary
      do while(t.gt.Thocb)
         read(Libhoc,200) Thoca,Thocb,Hocx,Hocy
  200    format(f10.2,f11.2,9pf8.2,8x,9pf8.2)
      end do
      if(t.ge.Thoca .and. t.le.Thocb) then
         Ihoc=0
      endif
      return
      end
