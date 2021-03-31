      subroutine PRCNT9(jd,fract,jdate,nutpr,jd9,fract9)

      implicit none


c*** start of declarations inserted by spag
      integer   i, j, jd, jd9

c*** end of declarations inserted by spag


c
c M.E.Ash   Sept 1972    subroutine PRCNT9
c The nutation-precession matrix is calculated.
c
c    JD    =Julian day number ("date")
c    FRACT = fraction of day from midnight
c    JDATE =0 NUTPR is output transformation matrix between
c             mean equinox & equator of epoch and mean equinox
c             & equator of date
c          =1 NUTPR is output transformation matrix between
c             mean equinox & equator of epoch and true equinox
c             & equator of date
c          =2 NUTPR is output transformation matrix between
c             mean equinox & equator of epoch and mean equinox
c             & ecliptic of date
c          =4 NUTPR is output transformation matrix between
c             mean equinox & equator of epoch and mean lunar
c             plane of date (x axis along intersection of mean
c             ecliptic and mean lunar plane of date)
c          =5 NUTPR is output transformation matrix between
c             mean equinox & equator of epoch and mean equinox
c             & true equator of date
c    JD9, FRACT9 = desired epoch
c
      integer*2 jdate
      real*10 fract, fract9, t9
      real*10 nutpr(3,3), prec(3,3), moblq, pc, ps, dpsi, deps
      real*10 t, sobliq, cobliq, oblq
      real*10 nut2(3,3)

c sine and cosine of mean obliquity in 1950.0
      real*10 smoblq/0.3978811865927521_10/
      real*10 cmoblq/0.91743695225096726346_10/

c sine and cosine of mean inclination of lunar plane on the ecliptic
      real*10 sincm/8.965339585272849E-2_10/
      real*10 cincm/0.99597302604642559968_10/

c given that SIN(inc/2.0_10)=4.4886967E-2_10
      real*10 t1, FUNCOF, ascm, cascm, sascm
c useful constants
      include 'funcon.inc'

      t = jd + (fract - 0.5_10)
      t9 = jd9 + (fract9 - 0.5_10)
c
c Determine transformation between mean ecliptic and equator
c of 1950.0 (JDATE=3).  Note that JDATE is not permitted to be 3 any
c longer when this subroutine is called -- it would be changed to 2
c and denoted by the values of JD and FRACT.
      if(jdate.eq.2 .and. ABS(t9-2433282.423_10).lt.1.E-5_10 .and.
     .    ABS(t-2433282.423_10).lt.1.E-5_10) then
         nutpr(1,1) = 1._10
         nutpr(1,2) = 0._10
         nutpr(1,3) = 0._10
         nutpr(2,1) = 0._10
         nutpr(2,2) = cmoblq
         nutpr(2,3) = smoblq
         nutpr(3,1) = 0._10
         nutpr(3,2) = -smoblq
         nutpr(3,3) = cmoblq
         return
      endif
c
c Determine precession matrix and mean obliquity
      call PRCES9(t,prec,moblq,t9)
      if(jdate.eq.5 .or. jdate.eq.1) goto 50
c
      if(jdate.eq.0) then
c Determine transformation matrix to mean equator of date
         do j = 1, 3
            do i = 1, 3
               nutpr(i,j) = prec(i,j)
               end do
            end do
         return
      endif
c
c Determine transformation matrix to ecliptic of date
c (if JDATE=2 or 4)
      cobliq = COS(moblq)
      sobliq = SIN(moblq)
      do j = 1, 3
         t1        = cobliq*prec(2,j) + sobliq*prec(3,j)
         prec(3,j) = -sobliq*prec(2,j) + cobliq*prec(3,j)
         prec(2,j) = t1
         do i = 1, 3
            nutpr(i,j) = prec(i,j)
            end do
         end do
      if(jdate.eq.2) return
c
c Determine transformation from mean equator of 1950.0 to
c mean lunar plane of date
      t1    = t - 2415020._10
      ascm  = FUNCOF(t1, 0.71995354167_10, -1,
     . -4.7094228332E-5_10, +0.432630E-14_10, +1.266E-22_10)*Twopi
      sascm = SIN(ascm)
      cascm = COS(ascm)
      nut2(1,1) = cascm
      nut2(1,2) = -sascm*cincm
      nut2(1,3) = sascm*sincm
      nut2(2,1) = sascm
      nut2(2,2) = cascm*cincm
      nut2(2,3) = -cascm*sincm
      nut2(3,1) = 0._10
      nut2(3,2) = sincm
      nut2(3,3) = cincm
      call PRODCT(nut2,prec,nutpr, -3,3,3)
      return
c
c Determine nutation
   50 call FUNARG(t)
      call NUTAT(t, dpsi, deps)
      dpsi = dpsi*Convds
      deps = deps*Convds

      oblq   = moblq + deps
      cobliq = COS(oblq)
      sobliq = SIN(oblq)
      pc     = cobliq*dpsi
      ps     = sobliq*dpsi
c
c is transformation to true equator and mean equinox of date
      if(jdate.eq.5) pc = 0._10
c
c calculate product of nutation and precession matrices
c
c nutation to second order
      nut2(1,1) = 1._10 - 0.5_10*dpsi**2
      nut2(1,2) = -pc
      nut2(1,3) = -ps
      nut2(2,1) = pc - deps*ps
      nut2(2,2) = 1._10 - 0.5_10*(deps**2 + pc**2)
      nut2(2,3) = -deps - 0.5_10*ps*pc
      nut2(3,1) = ps + deps*pc
      nut2(3,2) = deps - 0.5_10*ps*pc
      nut2(3,3) = 1._10 - 0.5_10*(deps**2 + ps**2)
      call PRODCT(nut2,prec,nutpr, -3,3,3)

      return
      end
