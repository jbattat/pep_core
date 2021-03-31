      subroutine PARTVL(der, n, kick)
 
      implicit none

      real*10 der(6,2)
      integer j, kick, n
 
c           subroutine partvl - j.f.chandler - 1980 jan 14
c           apply scale factors to obtain velocities as derivatives
c           w.r.t. observer's coordinate time.  see memo by jfc
c           "the light-time correction to partial derivatives in pep"
c           dated aug 28, 1978.  these partial derivatives of the
c           velocities are not exact (as are the velocities themselves)
c
c           der  - array of partials, either send/receive or reflect
c           n    - indicator of type of array
c                  n=2 -> apply dt3/dt1 to send-reflect leg
c                  n=1 -> apply dt2/dt1 to (only) leg
c           kick - type of observation series
c           kick =1 subroutine radar is calling program
c           kick =2 subroutine optic is calling program
c           kick =3 subroutine trnsit is calling program
c           kick =4 subroutine fermtr is calling program
 
      include 'fcntrl.inc'
      include 'partcm.inc'
      include 'rtrdvl.inc'
 
      if( Index .gt. 3 .and. Mouse .ge. n ) then
         if( kick .eq. 1 ) then
            if( Jct(67) .gt. 0 ) then
               do j = 4, 6
                  der(j, n) = der(j, n)*corfct(n)
               end do
            endif
         endif
      endif
 
      return
      end
