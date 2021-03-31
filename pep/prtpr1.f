      subroutine PRTPR1(lgo, kick)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 DOT
      integer   i, ii, j, kick
 
c*** end of declarations inserted by spag
 
 
c
c           m. ash / w. decampli    aug. 1973   subroutine prtpr1
c              calculate euler angle planet rotation partials of
c           time delay and doppler observations of a spot (feature)
c           on a planet (right ascension, declination, angular velocity,
c           and rotational phase angle at epoch).
c              this subroutine is called by partl1.
c
c
      integer*2 lgo
      real*10 drot(3, 6), dd
      real*10 sitfct(6)/3*1._10, 3*8.64E4_10/
c
c common
      include 'ltrapobs.inc'
      integer*2 kind
      equivalence (kind, Numpar)
      include 'mnsprt.inc'
      include 'partcm.inc'
      include 'prtpr9.inc'
c-----------------------------------------------------------------------
c
      if(lgo.lt.6 .or. lgo.gt.9) then
c
c partials w.r.t. con(1-5) are all zero
c
         Deriv(kind, 1) = 0._10
         Deriv(kind, 2) = 0._10
         return
      else if(lgo.eq.6) then
c
c partial derivatives of the rotation matrix w.r.t. rotational
c phase angle con(6) at epoch con1(1).
c
         dd = 1._10
      else if(lgo.eq.7) then
c
c partial derivatives w.r.t. rotation angular frequency:
c twopi/con(7).
         dd = Stime
      else if(lgo.eq.8) then
c
c partial derivative w.r.t. minus declination of pole of
c rotation (-con(8)).
         do i = 1, 3
            drot(i, 1) = Salp*Rot(i, 3)
            drot(i, 2) = -Calp*Rot(i, 3)
            end do
         drot(1, 3) = Spsi*Cdel
         drot(2, 3) = Cpsi*Cdel
         drot(3, 3) = -Sdel
         go to 50
      else
c
c partial derivatives w.r.t. right ascension of pole of
c rotation (con(9)).
         do i = 1, 3
            drot(i, 1) = -Rot(i, 2)
            drot(i, 2) = Rot(i, 1)
            drot(i, 3) = 0._10
            end do
         go to 50
      endif
 
c for CON(6) and CON(7)
      do i = 1, 3
         drot(1, i) = Rot(2, i)*dd
         drot(2, i) = -Rot(1, i)*dd
         drot(3, i) = 0._10
         end do
c
c calculate time derivatives of the above partials
c
   50 if(Index.gt.3) then
         do i = 1, 3
            ii = i + 3
            drot(1, ii) = drot(2, i)*Omegsc
            drot(2, ii) = -drot(1, i)*Omegsc
            drot(3, ii) = 0._10
            end do
         if(lgo.eq.7) then
 
c extra term for frequency partial
            do i = 1, 3
               drot(1, i + 3) = drot(1, i + 3) + Rot(2, i)/sitfct(4)
               drot(2, i + 3) = drot(2, i + 3) - Rot(1, i)/sitfct(4)
               end do
         endif
      endif
c
c calculate partial derivative matrix of the equinox-equa-
c torial coordinates w.r.t. planet rotational parameters,
c con(6-9).
      do j = 1, Index
         derem(j, 1) = -DOT(drot(1,j), Yspcd(1,1))*sitfct(j)
         end do
      call PARTVL(derem, 1, kick)
      if(Mouse.eq.2) then
 
c copy to send array
         do j = 1, Index
            derem(j, 2) = derem(j, 1)
            end do
      endif
c
c partials of doppler and delay are calculated in cpartc
      call CPARTC(kick)
      return
      end
