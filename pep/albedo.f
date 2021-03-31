      subroutine ALBEDO(s,sbcor,rfla,kreflt)
 
      implicit none
c
c s.margolis/r.mckinnis   may 1974   subroutine albedo
c read a data set of albedo pressure at tabular times and eccentric
c anomaly, for interpolation by reflcc.
c
      real*10 s,szero,t,sbcor(6),rfla(3)
      real*10 time(2)
      real*4    press(100,3,2),up(100,2)
      integer*4 iumax(2),ja,jb,kreflt,m1,m2

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'inodta.inc'
      include 'petuna.inc'
 
      data time(1)/0./,time(2)/0./
      data iumax/1,1/
      data up/200*0./,press/600*0./
      data szero/2440000.5/
      logical*4 uflag
      integer*4 jalb
c
c kkp(80) = data set for reflected radiation pressure (albedo)
c kkp(80) is read in the nmlst2 namelist as kk(80)
      jalb = Kkp(80)
 
c set uflag to calculate u
      uflag = .true.
 
c t is the short form of time of present integration
      t = s - szero
      do while( t.ge.time(1) )
         if(t.le.time(2)) then
 
c the rest is common to both calculations.
            call REFLCC(t,sbcor,time,iumax,press,up,uflag,rfla)
            return
         else
c t is outside the time bracket(time(1), time(2)).
c must read in a new time.
c move right end of time interval to left
            time(1)  = time(2)
            iumax(1) = iumax(2)
            m1 = iumax(1)
            do ja = 1,m1
               up(ja,1) = up(ja,2)
               do jb = 1,3
                  press(ja,jb,1) = press(ja,jb,2)
               end do
            end do
 
c get a new right end point
            read(jalb,end=300) time(2),m2,
     .       ((press(ja,jb,2),jb=1,3),up(ja,2),ja=1,m2)
 
c test time t for this new interval
            iumax(2) = m2
         endif
      end do
 
      write(Iout,100) s,t,time(1)
  100 format(1x,'ERROR IN ALBEDO TIME',3D16.7)
      call SUICID('END IN ALBEDO   ',4)
 
  300 write(Iout,400) s,t,time(2)
  400 format(1x,'ALBEDO END OF FILE',3D16.7)
      call SUICID('END IN ALBEDO   ',4)
 
      end
