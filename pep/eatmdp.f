      real function EATMDP(wetz, dryz, zen, dwetz, ddryz, zendot)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      ddryz, degthe, dryz, dwetz, v, wetz, zen, zendot
 
c*** end of declarations inserted by spag
 
 
c
c r.king  july 1974   function eatmdp
c
c     this function computes the change in phase delay rate caused by
c     the wet and dry components of the earth's neutral atmosphere.
c     this routine uses pfeiffer's anlytic formula to get the rate of
c     atmospheric delay as a function of zenith angle rate;  it is
c     it is therefore inconsistent with eatmdl, which uses chao's table.
c           wetz - wet zenith correction (sec.)
c           dryz - dry zenith correction (sec.)
c           zen  - zenith angle (rad)
c           dwetz- derivative of wet zenith correction (sec/sec)
c           ddryz- derivative of dry zenith correction (sec/sec)
c           zendot- derivative of zenith angle  (rad/sec)
      degthe = zen*57.2957795
      if((zen .gt. 0.) .and. (degthe .lt. 90.) ) then
         if( degthe .le. 50. ) then
            v = (2.14*sin(zen)/cos(zen)**2)/2.154
         else if( degthe .le. 75. ) then
            v = (3.964 + 2.04*(sin(zen)/cos(zen)**2-1.85))/2.154
         else if( degthe .le. 80. ) then
            v = (29.64 + 1.865*(sin(zen)/cos(zen)**2-14.41))/2.154
         else
            v = (63.71 + 1.442*(sin(zen)/cos(zen)**2-63.71))/2.154
         endif
         EATMDP = -(dryz + wetz)*zendot*v
      endif
      return
      end
