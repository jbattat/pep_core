      subroutine REFLCC(t, sbcor, time, iumax, press, up, uflag, pvec)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 DOT
      real      a, amass, cosu, d, e, one, r2, ratio, rv, sinu,
     .          sqta, two, u, uone, v2, zero
      integer   iub, iut, ja, jb
 
c*** end of declarations inserted by spag
 
 
c
c s.margolis/r.mckinnis   may 1974   subroutine reflcc
c calculate eccectric anomaly from rectangular coordinates,
c use anomaly and time in tabular interpolation.
c
      real*10 t, sbcor(6), time(2), pvec(3)
      real*4    press(100, 3, 2), up(100, 2)
      real*4    pr(3, 2)
      integer*4 iumax(2)
      logical*4 uflag
      include 'funcon.inc'
      include 'inodta.inc'
c
c     this program is given tables of pressure vs eccentric anomaly
c     for two times.
c     it outputs an interpolated pressure vector based on sbcor and t.
c     t -- julian days after the time of the initial conditions
c          used to generate the tables.
c     sbcor  --  position and velocity of the spacecraft(au,days)
c     press(100,3,2) two tables of pressures, with iumax(2) points
c     and for time time(2).
c     up(100,2) -- the abscissas, eccentric anomaly, for press
c     uflag -- when .true. causes a new u to be calculated from sbcor
c     pvec -- the outputted pressure vector
c
      data iub/1/
      data zero, one, two/0.0E00, 1.0E00, 2.0E00/
      data amass/3.227E-07/
 
c do not recalculate u for infred if already for albedo
      if( uflag ) then
         r2 = DOT(sbcor, sbcor)
         rv = DOT(sbcor, sbcor(4))/(Gauss*sqrt(amass))
         v2 = DOT(sbcor(4), sbcor(4))/(Gauss**2*amass)
         d  = sqrt(r2)
 
c compute semi-major axis
         a    = one/(two/d - v2)
         sqta = sqrt(a)
         sinu = rv/sqta
         cosu = d*v2 - one
 
c compute eccentricity
         e    = sqrt(cosu*cosu + sinu*sinu)
         sinu = sinu/e
         cosu = cosu/e
 
c compute eccentric anomaly
         if( cosu .gt. one ) write(Iout, 50) time, cosu, sinu, e, a
   50    format(1x, 'ERROR IN REFLCC', 2D24.16, 3D18.10, '*****')
         u = ATAN2(sinu, cosu)
         if( u .lt. zero ) u = Twopi + u
      endif
 
c now look up u in the tables
      do ja = 1, 2
         iub = max0(1, min0(iub,iumax(ja)))
 
c find a mesh point which u is not below
         do while( u .lt. up(iub,ja) )
            if( iub .gt. 1 ) then
               iub = iub - 1
            else
               u   = u + Twopi
               iub = iumax(ja)
            endif
            end do
c find an upper mesh point for u
c if iub is iumax(ja), then upper mesh point is 1, by folding,
c and has been found.
         do while( iub .ne. iumax(ja) )
            iut = iub + 1
            if( u .lt. up(iut,ja) ) then
 
c between two points of array
               ratio = (u - up(iub,ja))/(up(iut,ja) - up(iub,ja))
               go to 100
            else
 
c u is still above the interval, so move it up another notch
               iub = iub + 1
            endif
            end do
c the linear interpolation formula has two cases
c at top of up array, fold over
         iut   = 1
         uone  = up(1, ja) + Twopi
         ratio = (u - up(iub,ja))/(uone - up(iub,ja))
 
c compute the interpolated pressure.  (iub,iut) has been set above.
  100    do jb = 1, 3
            pr(jb, ja) = press(iub, jb, ja)
     .                   + ratio*(press(iut,jb,ja) - press(iub,jb,ja))
            end do
         end do
 
c interpolate between the two times
      do jb = 1, 3
         pvec(jb) = pr(jb, 1) + ((t-time(1))/(time(2)-time(1)))
     .              *(pr(jb,2) - pr(jb,1))
         end do
      return
      end
