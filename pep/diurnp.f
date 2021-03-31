      subroutine DIURNP(jd, fract, thetg, dtheta)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 aa, c, d, d1, d3, dsq, dtheta, fract,
     .          term2, thetg, xm, xm2, ym, ym2
      integer   j, jd, m
 
c*** end of declarations inserted by spag
 
 
c
c
c     r.w.king        mar 1980        subroutine diurnp
c     calculate diurnal polar motion
c     copied with minor modifications from vlbi3 routine written origi-
c     nally by d.s.robertson and modified by c.ma, p.gasiz, and t.herrin
c
      logical   bothcs
      include 'funcon.inc'
      include 'nutprc.inc'
      real*4    dpsi(2), deps(2)
      equivalence (Nutat, dpsi), (Nutat(3), deps)
      real*10 farg(6), sum1(2,2), sum2(2), a(4,5), t1t2(2,7), t3(7)
      integer*2 nf(6,7)
      equivalence (xm, sum1(1,1)), (ym, sum1(2,1)), (xm2, sum1(1,2)),
     .            (ym2, sum1(2,2))
      data t1t2/-0.0012840_10, -0.0009149_10,
     *          -0.0012292_10, -0.0008659_10,
     *          -0.0065199_10, -0.0045924_10,
     *          -0.0028973_10, -0.0019980_10,
     *           0.0027652_10,  0.0019037_10,
     *           0.0059457_10,  0.0040934_10,
     *           0.0011784_10,  0.0008113_10/
      data t3/-0.00038502_10, -0.00041647_10, -0.00201046_10,
     *        -0.00093797_10,  0.00086146_10,  0.00185233_10,
     *         0.00042084_10/
      data nf/-1, 0, -2, 0, -2, 1,
     *         0, 0, -2, 0, -1, 1,
     *         0, 0, -2, 0, -2, 1,
     *         0, 0, -2, 2, -2, 1,
     *         0, 0,  0, 0,  0, 1,
     *         0, 0,  0, 0,  0, 1,
     *         0, 0,  0, 0, -1, 1/
      real*10 k, ks
      data c/8.068E40_10/, aa/8.042E40_10/, k/0.29_10/, ks/0.937_10/
      data a/296.104608_10,13.0649924465_10, .0006890_10, .000000295_10,
     .       358.475833_10, 0.9856002669_10,-.0000112_10,-.000000068_10,
     .        11.250889_10,13.2293504490_10,-.0002407_10,-.000000007_10,
     .       350.737486_10,12.1907491914_10,-.0001076_10, .000000039_10,
     .       259.183275_10,-0.0529539222_10, .0001557_10, .000000046_10/
 
      d1  = jd + fract - 0.5_10 - 2415020._10
      d   = d1*1E-4_10
      dsq = d*d
      d3  = dsq*d
      do j = 1, 5
         farg(j) = (a(1,j) + a(2,j)*d1 + a(3,j)*dsq + a(4,j)*d3)*Convd
         end do
      m = farg(5)/Twopi
      farg(5) = farg(5) + (1 - m)*Twopi
      farg(5) = MOD(farg(5), Twopi)
      farg(6) = thetg - Pc(1)
      bothcs  = .true.
      call DNPOLE(7, 6, t1t2, 2, nf, bothcs, farg, sum1)
      term2 = ((c-aa)/c)*(1._10 - k/ks)
c x,y are coordinates of angular momentum pole in bih lefthand syste
c
      Ywob   = Ywob - Ywob*term2 - ym2
      Xwob   = Xwob - Xwob*term2 + xm2
      dtheta = term2*Convds*(Xwob*COS(thetg) + Ywob*SIN(thetg))
      bothcs = .false.
      call DNPOLE(7, 5, t3, 1, nf, bothcs, farg, sum2)
      dtheta = (sum2(1)/TAN(Moblq+deps(1)))
      return
      end
