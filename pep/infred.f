      subroutine INFRED(s,sbcor,rfli,kreflt)
 
      implicit none
c
c s.margolis/r.mckinnis   may 1974   subroutine infred
c read a data set of infrared pressure at tabular time and eccentric
c anomaly, for interpolation by reflcc.
c
      real*10 s,szero,t,sbcor(6),rfli(3)
      real*10 time(2)
      real*4    press(100,3,2),up(100,2)
      integer*4 iumax(2),ja,jb,kreflt,m,m2

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'inodta.inc'
      include 'petuna.inc'
 
      data time(1)/0./,time(2)/0./
      data iumax/1,1/
      data up/200*0./,press/600*0./
      logical*4 uflag
      data szero/2440000.5/
      integer*4 jinf
c kkp(79) = data set for emitted radiation pressure (infrared)
c kkp(79) is read in the nmlst2 namelist as kk(79)
      jinf = Kkp(79)
 
c set uflag to calculate u
      uflag = .true.
 
c if albedo was used, then u  has been calculated.
      if(kreflt.eq.3) uflag = .false.
 
c t is the short form of time of present integration
      t = s - szero
      do while( t.ge.time(1) )
         if(t.le.time(2)) then
 
c the rest is common to both calculations.
            call REFLCC(t,sbcor,time,iumax,press,up,uflag,rfli)
            return
         else
c t is outside the time bracket (time(1), time(2)).
c must read in a new time.
c move right end of time interval to left
            time(1)  = time(2)
            iumax(1) = iumax(2)
            m = iumax(1)
            do ja = 1,m
               up(ja,1) = up(ja,2)
               do jb = 1,3
                  press(ja,jb,1) = press(ja,jb,2)
               end do
            end do
 
c get a new right end point
            read(jinf,end=300) time(2),m2,
     .       ((press(ja,jb,2),jb=1,3),up(ja,2),ja=1,m2)
 
c test time t for this new interval
            iumax(2) = m2
         endif
      end do
 
      write(Iout,100) s,t,time(1)
  100 format(1x,'INFRED TIME OUT OF BOUNDS',3D16.7)
      call SUICID('END IN INFRED   ',4)
 
  300 write(Iout,400) s,t,time(2)
  400 format(1x,'INFRED END OF FILE',3D16.7)
      call SUICID('END IN INFRED   ',4)
 
      end
