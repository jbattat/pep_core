      subroutine DLIBRA(jd,fract,librat,seli,nvel,beta,gamma)
 
      implicit none

c
c eckhardt/slade/king/ash  aug 1972   subroutine dlibra
c calculation of moon physical libration using eckhardt's second
c order model
c
c arguments
      integer*4 jd,nvel
      real*10 fract,seli,beta,gamma
      real*4 librat(2,3)

c common
      include 'funcon.inc'
 
c local
      real*10 td,julian,jd1900/2415020._10/
      real*4    dlibrt(6,2),ibet,igam
      real*10 sig(9),pea(9),tau(18)
      real*4    taua(15),taub(15),tauc(15),siga(9),sigb(9),sigc(9),
     .          peaa(9),peab(9),peac(9)
      real*10 ss(9),cp(9),st(18),ff,ssr(9),cpr(9),str(18)
      real*4    d,d2,f2,l,lp,l2,lp2,sn
      real*10 ldoub,argdub,sl,cl,den,denrt
      real*10 slong,soll,g,w,wp,per,gp
      real*10 t1,t2,t3,sl1,per1,soll1,gp1,sn1
      real*10 t,s,p
      real*4    tquand,tquanb,tquanr,lr,lpr,ffr,dr,l2r,lp2r,
     .          d2r,f2r
      real*4    arg,cd2,cf2,cf2d2,cl2,cld2,clf2,clf2d2,clp,gpr,
     .          gr,perr,sd2,sf2,sf2d2,sl2,slaf2,sld2,slf2,
     .          slf2d2,slongr,slp,snr,sollr,wpr,wr
      integer   i,j
      integer*4 newval/0/, nwval1/0/
c
c define constants
      data siga/-.26_10,-3.02_10,-10.61_10,-26.03_10,2.50_10,
     .     -102.63_10,.50_10,-.84_10,-.91_10/
      data peaa/-3.09_10,-10.81_10,.24_10,26.11_10,-1.97_10,-100.36_10,
     .     .53_10,-.74_10,-.39_10/
      data taua/-.47_10,1.58_10,.36_10,87.07_10,.22_10,-1.37_10,
     .     -.44_10,3.98_10,-3.31_10,-16.23_10,.22_10,.46_10,.91_10,
     .     9.54_10,-.43_10/
      data tauc/-.51_10,1.73_10,.40_10,96.12_10,.24_10,-1.37_10,
     .     -.40_10,4.33_10,-3.67_10,-17.55_10,.24_10,.43_10,1.02_10,
     .     10.45_10,-.47_10/
      data peac/-3.05_10,-10.85_10,.20_10,22.42_10,-1.89_10,-97.82_10,
     .     .48_10,-.73_10,-.36_10/
      data sigc/-.26_10,-2.99_10,-10.69_10,-22.35_10,2.43_10,
     .     -100.14_10,.45_10,-.83_10,-.87_10/
      data taub/-.47_10,1.57_10,.36_10,87.06_10,.22_10,-1.36_10,
     .     -.53_10,3.98_10,-3.32_10,-16.42_10,.22_10,.45_10,.91_10,
     .     9.54_10,-.43_10/
      data peab/-3.24_10,-11.46_10,.26_10,30.30_10,-2.13_10,-107.28_10,
     .     .59_10,-.79_10,-.42_10/
      data sigb/-.27_10,-3.17_10,-11.22_10,-30.22_10,2.69_10,
     .     -109.68_10,.56_10,-.89_10,-.98_10/

c portable single-precision arithmetic statement function
      real*4 SNG10
      real*10 x
      SNG10(x)=x
 
      julian = jd+fract-.5_10
      if(newval.ne.1) then
         newval = 1
         den    = gamma*1.E5_10-21.26_10
         denrt  = SQRT(den)
         tquand = (24.46_10/den-10.84_10/denrt+10.84_10/SQRT(denrt))
         tquanb = (1._10-(1.04_10/32.14_10)*(beta-.00063_10)/.00003_10)
         tau(1) = tquand*tquanb
         tau(2) = 7.5_10
         do j = 1,15
            tau(j+2) = taua(j)+(beta-.00063_10)
     .                   /.00003_10*(taub(j)-taua(j))
     .                  +(gamma-.00022_10)
     .                   /.00002_10*(tauc(j)-taua(j))
         end do
         tau(18) = 14.27_10
         do j = 1,9
            sig(j) = siga(j)+(beta-.00063_10)
     .       /.00003_10*(sigb(j)-siga(j))+(gamma-.00022_10)
     .       /.00002_10*(sigc(j)-siga(j))
            pea(j) = peaa(j)+(beta-.00063_10)
     .       /.00003_10*(peab(j)-peaa(j))+(gamma-.00022_10)
     .       /.00002_10*(peac(j)-peaa(j))
         end do
         seli = (5558.5_10+347.6_10*(beta-.00063_10)/.00003_10 -
     .    (gamma-.00022_10)/.00002_10)*Convds
      endif
 
c ***
      td = julian-jd1900
      t1 = td/10000._10
      t2 = t1*t1
      t3 = t2*t1
c
c explanatory supplement
c
      sl1   = (270.434164_10+1.31763965268E5_10*t1-8.5E-5_10*t2+
     .        3.9E-8_10*t3)*Convd
      slong = MOD(sl1,Twopi)
 
      per1 = (334.329556_10+1.114040803E3_10*t1-7.739E-4_10*t2+
     .       2.6E-7_10*t3)*Convd
      per  = MOD(per1,Twopi)
      sn1  = (259.183275_10-5.29539222E2_10*t1+1.557E-4_10*t2+
     .       5.00E-8_10*t3)*Convd
      sn   = MOD(sn1,Twopi)
 
      soll1 = (279.696678_10+9.856473354E3_10*t1+2.267E-5_10*t2)*Convd
      soll  = MOD(soll1,Twopi)
      gp1   = (358.475845_10+9.85600267E3_10*t1-1.12E-5_10*t2 -
     .        7.0E-8_10*t3)*Convd
      gp    = MOD(gp1,Twopi)
c
c arguments  g,w,wp
      g     = slong-per
      w     = per-sn
      wp    = soll-gp-sn
      ldoub = g
      argdub = 2._10*ldoub-2._10*(g+w)
      l   = g
      lp  = gp
      ff  = g+w
      d   = g-gp+w-wp
      l2  = 2.*l
      lp2 = 2.*lp
      d2  = 2.*d
      f2  = 2.*ff
c
c define
      sd2    = sin(d2)
      sf2d2  = sin(f2-d2)
      sf2    = sin(f2)
      slf2   = sin(l-f2)
      sld2   = sin(l-d2)
      sl     = SIN(ldoub)
      slf2d2 = sin(l+f2-d2)
      slaf2  = sin(l+f2)
      sl2    = sin(l2)
      cf2d2  = cos(f2-d2)
      cf2    = cos(f2)
      clp    = cos(lp)
      clf2   = cos(l-f2)
      cld2   = cos(l-d2)
      cl     = COS(ldoub)
      clf2d2 = cos(l+f2-d2)
      cl2    = cos(l2)
      slp    = sin(lp)
      cd2    = cos(d2)
      cf2d2  = cos(f2-d2)
 
c sig vector coefficients
      ss(1) = sd2
      ss(2) = sf2d2
      ss(3) = sf2
      ss(4) = slf2
      ss(5) = sld2
      ss(6) = sl
      ss(7) = slf2d2
      ss(8) = slaf2
      ss(9) = sl2
c
c pea vector coefficients
      cp(1) = cf2d2
      cp(2) = cf2
      cp(3) = clp
      cp(4) = clf2
      cp(5) = cld2
      cp(6) = cl
      cp(7) = clf2d2
      cp(8) = cos(l+f2)
      cp(9) = cl2
c
c tau vector coefficients
      st(1)  = SIN(argdub)
      st(2)  = sin(sn)
      st(3)  = sd2
      st(4)  = sf2d2
      st(5)  = sin(lp-f2+d2)
      st(6)  = slp
      st(7)  = sin(lp2)
      st(8)  = sin(l-lp-d)
      st(9)  = slf2
      st(10) = sld2
      st(11) = sin(l-d)
      st(12) = sl
      st(13) = sin(l+lp-d2)
      st(14) = sin(l2-lp2-d2)
      st(15) = sin(l2-lp-d2)
      st(16) = sin(l2-d2)
      st(17) = sin(l2)
      arg    = .53733431_10-1.0104982E-5_10*td+1.91E-14_10*td**2
      st(18) = sin(SNG10(Twopi*arg))
c
c calculate libration angles
      s = 0.
      p = 0.
      t = 0.
      do i = 1,9
         s = s+sig(i)*ss(i)
         p = p+pea(i)*cp(i)
      end do
      do i = 1,18
         t = t+tau(i)*st(i)
      end do
      librat(1,1) = Convds*t
      librat(1,2) = Convds*p
      librat(1,3) = (seli*t-s)*Convds
c
c calculate libration rates
      if(nvel.gt.0) then
         slongr= (13.1763965268_10 - 17.0E-13_10*t1 + 11.7E-20_10*t2)
     .    *Convd
         perr  = (.1114040803_10 - 15.478E-12_10*t1 - 7.8E-19_10*t2)
     .    *Convd
         snr   = (-.0529539222_10 + 3.114E-12_10*t1 + 15.0E-20_10*t2)
     .    *Convd
         gpr   = (.985600267_10 - 2.24E-13_10*t1 - 21.0E-20_10*t2)*Convd
         sollr = (.9856473354_10 - 4.534E-13_10*t1)*Convd
         gr    = slongr-perr
         wr    = perr-snr
         wpr   = sollr-gpr-snr
         lr    = gr
         lpr   = gpr
         ffr   = gr+wr
         dr    = gr-gpr+wr-wpr
         l2r   = 2.*lr
         lp2r  = 2.*lpr
         d2r   = 2.*dr
         f2r   = 2.*ffr
         ssr(1)= d2r*cd2
         ssr(2)= (f2r-d2r)*cf2d2
         ssr(3)= f2r*cf2
         ssr(4)= (lr-f2r)*clf2
         ssr(5)= (lr-d2r)*cld2
         ssr(6)= lr*cl
         ssr(7)= (lr+f2r-d2r)*clf2d2
         ssr(8)= (lr+f2r)*clf2
         ssr(9)= l2r*cl2
         cpr(1)= -(f2r-d2r)*sf2d2
         cpr(2)= -f2r*sf2
         cpr(3)= -lpr*slp
         cpr(4)= -(lr-f2r)*slf2
         cpr(5)= -(lr-d2r)*sld2
         cpr(6)= -lr*sl
         cpr(7)= -(lr+f2r-d2r)*slf2d2
         cpr(8)= -(lr+f2r)*slaf2
         cpr(9)= -l2r*sl2
         str(1)  = (l2r-f2r)*cos(l2-f2)
         str(2)  = snr*cos(sn)
         str(3)  = d2r*cd2
         str(4)  = (f2r-d2r)*cf2d2
         str(5)  = (lpr-f2r+d2r)*cos(lp-f2+d2)
         str(6)  = lpr*clp
         str(7)  = lp2r*cos(lp2)
         str(8)  = (lr-lpr-dr)*cos(l-lp-d)
         str(9)  = (lr-f2r)*clf2
         str(10) = (lr-d2r)*cld2
         str(11) = (lr-dr)*cos(l-d)
         str(12) = lr*cl
         str(13) = (lr+lpr-d2r)*cos(l+lp-d2)
         str(14) = (l2r-lp2r-d2r)*cos(l2-lp2-d2)
         str(15) = (l2r-lpr-d2r)*cos(l2-lp-d2)
         str(16) = (l2r-d2r)*cos(l2-d2)
         str(17) = l2r*cl2
         str(18) = (Twopi*(-1.0104982E-5_10+3.82E-14_10*td))
     .             *SIN(SNG10(Twopi*arg))
         s = 0.
         p = 0.
         t = 0.
         do i = 1,9
            s = s+sig(i)*ssr(i)
            p = p+pea(i)*cpr(i)
         end do
         do i = 1,18
            t = t+tau(i)*str(i)
         end do
         librat(2,1) = Convds*t
         librat(2,2) = Convds*p
         librat(2,3) = (seli*t-s)*Convds
      endif
      return
 
c calculate partial derivatives
      entry DLIBP1(nvel,dlibrt,ibet,igam)
      if(nwval1.ne.1) then
         nwval1 = 1
         tquanr = 1.E5*(-24.46/den**2+5.42/den/denrt -
     .            2.71/den/SQRT(denrt))
      endif
      s    = 0.
      p    = 0.
      ibet = 56.07678
      igam = -0.2424068
c
c partials w.r.t. beta
      do i = 1,9
         s = s+(sigb(i)-siga(i))/.00003_10*ss(i)
         p = p+(peab(i)-peaa(i))/.00003_10*cp(i)
      end do
      t = tquand*(-1.04_10/(32.14_10*.00003_10))*st(1)
      do i = 1,15
         t = t+(taub(i)-taua(i))/.00003_10*st(i+2)
      end do
      dlibrt(1,1) = Convds*t
      dlibrt(2,1) = Convds*p
      dlibrt(3,1) = Convds*s
      if(nvel.gt.0) then
         s = 0.
         p = 0.
         do i = 1,9
            s = s+(sigb(i)-siga(i))/.00003_10*ssr(i)
            p = p+(peab(i)-peaa(i))/.00003_10*cpr(i)
         end do
         t = tquand*(-1.04_10/(32.14_10*.00003_10))*str(1)
         do i = 1,15
            t = t+(taub(i)-taua(i))/.00003_10*str(i+2)
         end do
         dlibrt(4,1) = Convds*t
         dlibrt(5,1) = Convds*p
         dlibrt(6,1) = Convds*s
      endif
c
c partials w.r.t. gamma
      s = 0.
      p = 0.
      do i = 1,9
         s = s+(sigc(i)-siga(i))/.00002_10*ss(i)
         p = p+(peac(i)-peaa(i))/.00002_10*cp(i)
      end do
      t = tquanb*tquanr*st(1)
      do i = 1,15
         t = t+(tauc(i)-taua(i))/.00002_10*st(i+2)
      end do
      dlibrt(1,2) = Convds*t
      dlibrt(2,2) = Convds*p
      dlibrt(3,2) = Convds*s
      if(nvel.gt.0) then
         s = 0.
         p = 0.
         do i = 1,9
            s = s+(sigc(i)-siga(i))/.00002_10*ssr(i)
            p = p+(peac(i)-peaa(i))/.00002_10*cpr(i)
         end do
         t = tquanb*tquanr*str(1)
         do i = 1,15
            t = t+(tauc(i)-taua(i))/.00002_10*str(i+2)
         end do
         dlibrt(4,2) = Convds*t
         dlibrt(5,2) = Convds*p
         dlibrt(6,2) = Convds*s
      endif
c note that partials returned are tau, rho, isig, etc., not tau, rho
c i*(sig-tau)
c
      return
      end
