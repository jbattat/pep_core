      subroutine TRPLST(jd,fract,nvel,lcntl,label,x)
 
      implicit none
 
c        subr. trplst - j.f.chandler - 1999 nov
c        derived from subr. evtrp
c        This is a routine for listing the interpolated coordinates for
c        debugging purposes.
c        input parameters:
c          jd      integer value of the coord. time of the desired pt
c          fract   fractional part of the coordinate time of the point
c          nvel =0 positions
c                1 velocities
c               -2 accelerations
c          lcntl   control integer used by caller
c          label   indicates which is the calling routine
c          x       coordinate array (pos+vel or acc)
c
c arguments
      integer*4 jd,nvel,lcntl
      character*(*) label
      real*10 fract,x(6)
c
c        commons
      include 'fcntrl.inc'
      include 'inodta.inc'
c
c
      integer*4 i,jc,jo,nvela
      integer*2 jtest(4)/4,8,16,32/
 
      nvela=iabs(nvel)
c
c debug printout
      jc = Jct(6)/jtest(nvela + 1)
      if(mod(jc,2).ne.0) then
         if(Line.gt.56) call OBSPAG
         jo=3*nvela
         if(nvela.eq.2) jo = 0
         write(Iout,750) label,jd,fract,nvel,lcntl,
     .                    (x(i+jo),i = 1,3)
  750    format(1x,a,': JD.F=', i7, f13.12, ' NV=', i2,
     .          '  L=', i2, '  X=', 1p3d23.15)
         Line = Line + 1
      endif
c
      return
      end
