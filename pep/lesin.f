      subroutine LEsin(nstop, lname, in0, les0)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ihr, imin, in0, ithrs0, ithrst, jd, JULDAY, les,
     .          nstop
 
c*** end of declarations inserted by spag
 
 
c
c m. e. ash       august 1975       subroutine lesin
c
c read les-8/9 thrust and attitude history cards
c
      integer*2 les0
 
      include 'funcon.inc'
      include 'inodta.inc'
 
      real*10 time0, time1, sec, etutc, fract
      integer*2 imonth, iday, iyear
      character*4 iblank/'    '/, name, lname
      real*4    thrst, er1pt, er2rl, er3yw
      real*10 of1pt, of2rl
c
c spool data from in to in0
      call PEPTIC(In, Iout, in0, 10,
     .            'THRUST AND ATTITUDE CARDS FOR LES-8/9   ', nstop, 1)
c
c initialization
      les    = les0
      time0  = 0._10
      ithrs0 = 0
c
c read thrust and attitude data
  100 read(in0, 200) name, imonth, iday, iyear, ihr, imin, sec, ithrst,
     .               thrst, er1pt, er2rl, er3yw, of1pt, of2rl
  200 format(a4, 5(1x,i2), 1x, f6.3, i3, f9.3, 3F8.4, 2F9.4)
c
c name   = satellite name (les8 or les9)
c imonth = month (1 to 12)
c iday   = day (1 to 31)
c iyear  = year (last two digits in 1900's, ge 100 in 2000's)
c ihr    = utc hour
c imin   = utc minute
c sec    = utc second
c
c ithrst = indicate if tangential thrusters are firing
c          0   no firing
c         +1   start or middle of firing (thrust is interpolated between
c              tabular points with +1, or last point for interpolation
c              can be -1)
c         -1   end of firing (can only occur after a +1 has occured,
c              thrust is nonzero before the -1 and zero after)
c
c thrst  = thrust (millipounds)
c          if thrust positive, thrusters on -z2 roll face are firing
c          with dv in +z2 roll direction
c              (west face on les9, dv with orbital velocity)
c              (east face on les8, dv against orbital velocity)
c          if thrust negative, thrusters on +z2 roll face are firing
c          with dv in -z2 roll direction
c              (east face on les9, dv against orbital velocity)
c              (west face on les8, dv with orbital velocity)
c note: les8 is turned upside down from les9
c
c er1pt  = attitude error about z1 pitch axis (deg)
c er2rl  = attitude error about z2 roll axis (deg)
c er3yw  = attitude error about z3 yaw axis (deg)
c
c of1pt  = offset angle about z1 pitch axis
c          (deg) (between -180 and +180)
c of2rl  = offset angle about z2 roll axis
c          (deg) (between -11.25 and +11.25)
c
c jd    = julian day number corresponding to imonth, iday, iyear
c fract = ephemeris (or coordinate) time fraction of day corresponding
c         to utc ihr, imin, sec
c
c           is this the right satellite
      if( name .eq. iblank ) then
 
         rewind in0
 
c end file les
         rewind les
         return
      else
         if( name .ne. lname ) then
            write(Iout, 220) name, lname
  220       format(1x, a4, ' DOES NOT MATCH ', a4,
     .             ', ERROR IN SUBROUTINE LESIN')
            nstop = nstop + 1
         endif
c
c convert from utc to ephemeris time
         jd    = JULDAY(imonth, iday, iyear)
         etutc = 46.184_10 + (iyear - 75)
         fract = (ihr*3600._10 + imin*60._10 + sec + etutc)/8.64E4_10
         i     = fract
         if( fract .lt. 0._10 ) i = i - 1
         jd    = jd + i
         fract = fract - i
         time1 = jd + fract
         if( time1 .le. time0 ) then
            write(Iout, 240) imonth, iday, iyear, ihr, imin, sec
  240       format(i3, '/', i2, '/', i2, i3, ':', i2, ':', f6.3,
     .             ' NOT INCREASING, ERRORIN SUBROUTINE LESIN')
            nstop = nstop + 1
         endif
c
c check for consistency of starting and stopping thrust
         if( ithrst .lt. 0 ) then
            if( ithrs0 .gt. 0 ) go to 300
         else if( ithrst .eq. 0 ) then
            if( ithrs0 .le. 0 ) go to 300
         else
            go to 300
         endif
         write(Iout, 250) imonth, iday, iyear, ihr, imin, sec, ithrst,
     .                    ithrs0
  250    format(i3, '/', i2, '/', i2, i3, ':', i2, ':', f6.3,
     .          ' ITHRST=', i3,
     .          ' NOT CONSISTENT WITH PREVIOUS ITHRST=', i3,
     .          ', ERROR IN SUBROUTINE LESIN')
         nstop = nstop + 1
      endif
c
c convert from millipounds to micropounds and from degrees to
c radians
  300 thrst = thrst*1.E3_10
      er1pt = er1pt*Convd
      er2rl = er2rl*Convd
      er3yw = er3yw*Convd
      of1pt = of1pt*Convd
      of2rl = of2rl*Convd
c
c write disk with thruster, attitude history
      time0  = time1
      ithrs0 = ithrst
      write(les) jd, fract, ithrst, thrst, er1pt, er2rl, er3yw, of1pt,
     .           of2rl
 
c recm=vbs, lrecl=52, blksize=524
      go to 100
      end
