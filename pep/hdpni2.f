      subroutine HDPNI2(int2, ndim, mdim, title, n)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, i1, i2, iatype, linrat, lts, mdim, n, ndim
 
c*** end of declarations inserted by spag
 
 
c
c subroutine hdpni2 paul macneil may, 1978
c print rne header information read by frmhed
c
c eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
c
c
c        ndim = number of array entries to be printed (max = 300)
c        mdim = first dimension of input array
c        caution: this subroutine must not set (modify) mdim or ndim
c
c        commons
      include 'fcntrl.inc'
      include 'inodta.inc'
 
      integer*2 int2(mdim, 1)
      integer*4 int4(mdim, 1)
      character*4 title(n)
      real*4    rl4(mdim, 1)
      real*10 rl8(mdim, 1)
 
      iatype = 1
      linrat = 20
      go to 100
c
c
c eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
      entry HDPNI4(int4, ndim, mdim, title, n)
      iatype = 2
      linrat = 10
      go to 100
c
c
c eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
      entry HDPNR4(rl4, ndim, mdim, title, n)
      iatype = 3
      linrat = 10
      go to 100
c
c
c eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
      entry HDPNR8(rl8, ndim, mdim, title, n)
      iatype = 4
      linrat = 5
c
c
c
c*start=700
  100 if(ndim.le.0 .or. mdim.le.0) then
         write(Iout, 150) ndim, mdim, title
  150    format(' NDIM=', i12, '  MDIM=', i12, '  TITLE='/(1x,30A4))
         call SUICID('BAD DIMENSIONS, STOP IN HDPNI2  ', 8)
      endif
 
c set up subtitle
      call PAGSET(title, -n)
      lts = 3 + (ndim - 1)/linrat
      if(Line.gt.50 .and. Line+lts.gt.59) call NEWPG
      call PAGHED(0)
 
      i1 = 1
      i2 = ndim
      if(ndim/linrat .gt. 58-Line) i2 = linrat*(59 - Line)
      do while( .true. )
 
         if(iatype.eq.1) then
c
c*start=1000
            write(Iout, 160) (int2(1,i), i = i1, i2)
  160       format(4(2x,5I6))
         else if(iatype.eq.2) then
c
            write(Iout, 180) (int4(1,i), i = i1, i2)
  180       format(2(3x,5I12))
         else if(iatype.eq.3) then
c
            write(Iout, 200) (rl4(1,i), i = i1, i2)
  200       format(2(1x,1p,5E13.5))
         else if(iatype.eq.4) then
c
            write(Iout, 220) (rl8(1,i), i = i1, i2)
  220       format(1x, 1p, 5D25.15)
         else
            return
         endif
c
c
         if(i2.ge.ndim) then
 
            Line = Line + (i2 - i1)/linrat + 1
            goto 300
         else
            call PAGE(0, 1)
            i1 = i2 + 1
            i2 = ndim
         endif
         end do
c
  300 return
      end
