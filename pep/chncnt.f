      subroutine CHNCNT(jd0, fract, nplnt, ncentr, mcentr, cond, icnd,
     .                  mcnd)
 
      implicit none
 
c
c m.e.ash   may 1974   subroutine chncnt
c change central body for initial conditions
c
c parameters
      integer*4 jd0,mcentr,mcnd
      real*10 fract, cond(6)
      integer*2 nplnt, ncentr, icnd

c array dimensions
      include 'globdefs.inc'

c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'param.inc'
      include 'prtcod.inc'

c local
      integer*2 kp(10), lcentr
      real*10 fractx,xcnt(6),xbod(6)
      integer*2 jtp0/0/, npln0/0/
      integer*4 i,iepoch,incnd,jd0x
c
c consistency checks
      if(icnd.ne.-1) call SUICID(' ICND.NE.-1, STOP IN CHNCNT ', 7)
      if(ncentr.le.0 .and. mcentr.le.0) goto 500
      if((ncentr.eq.10 .and. mcentr.ne.3) .or.
     .   (mcentr.eq.10 .and. ncentr.ne.3)) call SUICID(
     .   ' ONE OF NCENTR,MCENTR IS 10, OTHER NOT 3, STOP IN CHNCNT',14)
      if(ncentr.gt.0 .and. mcentr.gt.0 .and. ncentr.ne.3 .and.
     .   mcentr.ne.3) call SUICID(
     .    ' NCENTR.GT.0 AND MCENTR.GT.0, STOP IN CHNCNT',11)
      if(ncentr.gt.10 .or. mcentr.gt.10)
     .    call SUICID(' NCENTR AND/OR MCENTR .GT.10, STOP IN CHNCNT',11)
c
c determine central body coordinates
      do i = 1, 10
         kp(i) = -1
         if(ncentr.eq.i .or. mcentr.eq.i) kp(i) = 1
      end do
      if(kp(10).lt.1) then
         if(kp(3).eq.1) kp(10) = 1
         lcentr = ncentr
         if(mcentr.gt.0) lcentr = mcentr
      else
         lcentr = 10
      endif
      if(Ipert.le.0) call SUICID(
     .'IPERT=0, CANNOT CHANGE CENTRAL BODY. STOP IN CHNCNT ',13)
      call PRTRD1(1, npln0, lcentr, kp, 1, -1, -1)
      call PRTCRD(jd0, fract)
c velocity of body lcentr is numerically differentiated in prtcrd
c fix this?
      rewind Ipert
      Itrwnd(Ipert) = 0
c
c store central body coordinates
      if(lcentr.eq.3) then
         do i = 1, 6
            xcnt(i) = Xpert(i, 3) - Mass(10)*Xpert(i, 10)
         end do
      else
         do i = 1, 6
            xcnt(i) = Xpert(i, lcentr)
         end do
      endif
c
c change central body for initial conditions
      if(lcentr.lt.10) then
         if(mcentr.le.0) goto 200
      else if(ncentr.eq.10) then
         goto 200
      endif
      do i = 1, 6
         xbod(i) = cond(i) - xcnt(i)
      end do
      goto 300
  200 do i = 1, 6
         xbod(i) = cond(i) + xcnt(i)
      end do
 
c printout comments in midst of input stream
  300 write(Iout, 400) jd0, fract, ncentr, mcentr, lcentr, Ipert, cond,
     .                 xcnt, xbod
  400 format(' $', i8, f18.15, ' CENTRAL BODY', i3, ' CHANGED TO', i3,
     .       4x, 'ORIGINAL I.C., PLANET', i3,
     .       ' COORDINATES ON DATA SET', i3,
     .       ', NEW I.C. ARE'/(' $',1p,6D21.14))
      do i = 1, 6
         cond(i) = xbod(i)
      end do
c
c change initial conditions from cartesian coordinates
  500 ncentr = mcentr
      mcentr = -2
      incnd  = icnd
      icnd   = mcnd
      iepoch = Jct(13)+1
      if(iepoch.eq.2) then
         jd0x = 2451545
         fractx = 0.5_10
      else
         jd0x = 2433282
         fractx = 0.923_10
      endif
      call CHNCNE(nplnt,ncentr,cond,jd0x,fractx,0,jtp0,incnd,icnd,
     . iepoch)
 
      return
      end
