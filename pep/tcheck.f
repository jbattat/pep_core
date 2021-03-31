      subroutine TCHECK(name, t1, t0, t2, int, nstop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 hc, r1, range, t0, t1, t2
      integer   nstop
 
c*** end of declarations inserted by spag
 
 
c m.e.ash   oct 1968    subroutine tcheck
c check consistency of beginning,initial,ending times
c modified to use r*8 epochs - 1985 feb
c t0 is required to be consistent only if positive
      character*8 name
      integer*2 int
c
c common
      include 'inodta.inc'
 
      if(t0 .ne. 0._10) then
         if(t1 .gt. 0._10 .and. t2 .gt. 0._10) then
            if(t0 .lt. 0) return
            hc = 1._10
            if(int .lt. 8) hc = int/8._10
            if(int .le. 0) hc = 2._10**(int - 3)
            range = t2 - t1
            if(ABS(range) .ge. hc) then
               r1 = t0 - t1
               if(ABS(r1) .ge. hc) then
                  if(r1*(t2-t0) .le. 0._10) then
c
c error message for jd1,jd0,jd2 inconsistency
                     write(Iout, 10) name
   10                format(' ERROR FOR ', a8,
     .' BECAUSE VALUES OF JD1,JD0,JD2 NOT CONSISTENT FOR INTEGRATION WHI
     .CH WILL OCCUR')
                     if(Mout .gt. 0) write(Mout, 10) name
                     nstop = nstop + 1
                  endif
               endif
               return
            endif
         endif
c
c error message for jd1,jd2 inconsistency
         write(Iout, 50) name
   50    format(' ERROR FOR ', a8,
     .' BECAUSE VALUES OF JD1,JD2 NOT CONSISTENT FOR INTEGRATION WHICH W
     .ILL OCCUR')
         if(Mout .gt. 0) write(Mout, 50) name
         nstop = nstop + 1
      endif
 
      return
      end
