      subroutine LESTIN(jd, fract, ther, offset, itherf)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      deljd, delt, delthr, deltt, frac
      integer   im1, ind, index, itherf, jd, jj, kth, nback, nread
 
c*** end of declarations inserted by spag
 
 
c
c ash/costello      september 1975      subroutine lestin
c
 
      real*4    ther(4)
      real*10 offset(2)
      real*10 fract, deloff
c
c  jd, fract = input julian day number and fraction of day
c
c  ther      = output thrust level and attitude errors and offsets
c  ther(1)   = signed thrust level (+ in positive z2 roll direction,
c              - in negative z2 roll direction) (micropounds)
c  ther(2)   = attitude error about z1 pitch axis (radians)
c  ther(3)   = attitude error about z2 roll axis (radians)
c  ther(4)   = attitude error about z3 yaw axis (radians)
c  offset(1) = offset angle about z1 pitch axis (radians)
c  offset(2) = offset angle about z2 roll axis (radians)
c
c  itherf    = 0   no such information exists
c  itherf    = 1   such information was read and interpolated from
c                  data set les
c
      include 'thratt.inc'
c
c           at tabular point i (i=1,2,3)
c  jdles(i)  = julian day number
c  frles(i)  = ephemeris time fraction of day
c  ithrst(i) = indicate if tangential thrusters are firing
c              0   no firing
c             +1   start or middle of firing (thrust is interpolated
c                  between tabular points with +1, or last point for
c                  interpolation can be -1)
c             -1   end of firing (can only occur after a +1 has occured,
c                  thrust is nonzero before the -1 and zero after)
c
c  thr(1,i)  = thrust (micropounds)
c              if thrust positive, thrusters on -z2 roll face are firing
c              with dv in +z2 roll direction
c                  (west face on les9, dv with orbital velocity)
c                  (east face on les8, dv against orbital velocity)
c              if thrust negative, thrusters on +z2 roll face are firing
c              with dv in -z2 roll direction
c                  (east face on les9, dv against orbital velocity)
c                  (west face on les8, dv with orbital velocity)
c  note: les8 is turned upside down from les9
c  thr(2,i)  = attitude error about z1 pitch axis (radians)
c  thr(3,i)  = attitude error about z2 roll axis (radians)
c  thr(4,i)  = attitude error about z3 yaw axis (radians)
c  off(1,i)  = offset angle about z1 pitch axis (radians)
c  off(2,i)  = offset angle about z2 roll axis (radians)
c  les       = data set reference number for thrust and attitude
c              information that was written in subroutine lesin of
c              input link from input cards after units were changed
c  nrecls      = record number that has been read on data set les
c               0   nothing read yet, get three records into storage
c              -1   end of file encountered, no more to be gotten
c                   from this tape
c              >0   we are in middle of data set, be careful not to
c                   back into load point
 
      integer   i, j, k, kndex
 
c
c is there data?
      ind    = 1
      itherf = 0
      if( Jendls .eq. 1 ) go to 200
      if( Les .le. 0 ) return
      if( Nrecls .lt. 0 ) return
      if( Nrecls .gt. 0 ) go to 200
      nread = 3
      ind   = 1
  100 kndex = 3 - nread
c
c read nread records from tape.
      do k = 1, nread
         kndex = kndex + 1
         index = kndex
         read(Les, end=600) Jdles(kndex), Frles(kndex), Ithrst(kndex)
     .                        , (Thr(j,kndex), j = 1, 4),
     .                        (Off(j,kndex), j = 1, 2)
         Nrecls = Nrecls + 1
      end do
c
c check to see if jd lies within interval read.
  200 do i     = ind, 3
         kndex = i
         if( Jdles(kndex) .ge. jd ) then
            if( Jdles(kndex) .eq. jd ) then
               if( Frles(kndex) .lt. fract ) go to 350
               if( Frles(kndex) .eq. fract ) then
c
c time on tape record equals input time
                  if( kndex .le. 1 ) kndex = kndex + 1
                  im1 = kndex - 1
                  go to 800
               endif
            endif
c
c time on tape exceeds input time.  if tape record is
c first record in storage, backspace.
            if( kndex .le. 1 ) go to 400
            im1 = kndex - 1
            if( Jdles(im1) .lt. jd ) go to 800
            if( Jdles(im1) .ne. jd ) go to 400
            if( Frles(im1) .le. fract ) go to 800
            go to 400
         endif
 
  350    if( Jendls .eq. 1 .and. index .eq. kndex ) go to 700
c
c input time exceeds time on tape.  check next tape value.
      end do
250   kndex = i
 
c
c shift indices and read two more records.
      nread     = 2
      ind       = 2
      Jdles(1)  = Jdles(3)
      Frles(1)  = Frles(3)
      Ithrst(1) = Ithrst(3)
      do k = 1, 4
         Thr(k, kndex) = Thr(k, 3)
      end do
      do k = 1, 2
         Off(k, 1) = Off(k, 3)
      end do
      go to 100
 
c
c backspacing
c
c if nrec greater than or equal to 5, backspace 5 and read 3.
c if nrec equal to 4, backspace 4 and read 3.
c if nrecls less than 4, don't backspace.  return 0 values.
  400 nback = 5
      if( Nrecls .lt. 4 .and. index .ge. Nrecls ) go to 700
      if( Nrecls .lt. 4 ) nback = Nrecls
      if( Nrecls .eq. 4 ) nback = 4
  500 do jj = 1, nback
         backspace Les
      end do
      Nrecls = Nrecls - nback
      Jendls = 0
      nread  = 3
      ind    = 1
      go to 100
 
c
c end of tape file reached.
  600 index  = index - 1
      Jendls = 1
      kndex  = index
      if( index .ne. 0 ) then
         if( Jdles(index) .lt. jd ) go to 700
         if( Jdles(index) .le. jd ) then
c
c day on tape record equals input julian day number.
            if( Frles(index) .lt. fract ) go to 700
         endif
c
c time onls tape record is greater than or equal to input time.
         if( index .gt. 1 ) then
            if( Jdles(kndex) .ge. jd ) then
               if( Jdles(kndex) .eq. jd ) then
                  if( Frles(kndex) .lt. fract ) go to 300
                  if( Frles(kndex) .eq. fract ) then
c
c time on tape record equals input time
                     if( kndex .le. 1 ) kndex = kndex + 1
                     im1 = kndex - 1
                     go to 800
                  endif
               endif
c
c time on tape exceeds input time.  if tape record is
c first record in storage, backspace.
               if( kndex .le. 1 ) go to 400
               im1 = kndex - 1
               if( Jdles(im1) .lt. jd ) go to 800
               if( Jdles(im1) .ne. jd ) go to 400
               if( Frles(im1) .le. fract ) go to 800
               go to 400
            endif
 
  300       if( Jendls .eq. 1 .and. index .eq. kndex ) go to 700
c
c input time exceeds time on tape.  check next tape value.
            go to 250
         end if
      endif
 
c
c backspace 4 records if nrecls g.e. 3.
      if( Nrecls .ge. 3 ) then
         nback = 4
         go to 500
      endif
c
c jd not in range of data points.  return 0 values.
  700 do j = 1, 4
         ther(j) = 0.0
      end do
      do j = 1, 2
         offset(j) = 0.0
      end do
      return
 
 
c check to see if tangential thrusters are firing.
  800 ther(1) = 0._10
      kth     = 2
      if( Ithrst(im1) .ne. 0 .and. Ithrst(kndex) .ne. 0 ) then
         if( Ithrst(im1) .ne. -1 .or. Ithrst(kndex) .ne. 1 ) kth = 1
      endif
 
c
c interpolate
c
c find linear interpolated values between data pts i-1 and i.
c determine time interval between data points.
      deljd = 0.0
      im1   = kndex - 1
      if( Jdles(im1) .ne. Jdles(kndex) ) deljd = Jdles(kndex)
     .    - Jdles(im1)
      deltt = deljd + (Frles(kndex) - Frles(im1))
c
c determine time interval between input data and data pt i-1.
      deljd = 0.0
      if( jd .gt. Jdles(im1) ) deljd = jd - Jdles(im1)
      delt = deljd + (fract - Frles(im1))
c
c determine fraction of interval corresponding to input data.
      frac = delt/deltt
      do k = kth, 4
         delthr  = Thr(k, kndex) - Thr(k, im1)
         ther(k) = Thr(k, im1) + frac*delthr
      end do
      do k = 1, 2
         deloff    = Off(k, kndex) - Off(k, im1)
         offset(k) = Off(k, im1) + frac*deloff
      end do
      itherf = 1
 
      return
      end
