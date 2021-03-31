      subroutine ALALIZ
 
      implicit none
 
c
c r.reasenberg/d.white   jan 1974   subroutine alaliz
c form normal equations for maximum liklihood estimator or
c kalman-bucy filter from observation library tape
c
c
c common
      include 'fcntrl.inc'
      include 'obsdta.inc'
c
c setup least squares control constants
      call ANSET
      if( Iobcon .gt. 0 ) rewind Iobcon
c
c if only restoring solution, leave
      if((Iterat .le. 1) .and. (Ict(5) .gt. 1) ) return
c
c zero normal eqn
      call NRMSET(0)
c
c if not using obs to form norm eqn, leave
      if((Iterat .le. 1) .and. (Ict(5) .eq. 1) ) return
c
c form norm eqn from obslib tape and do series by series save
      call NRMICT
      return
      end
