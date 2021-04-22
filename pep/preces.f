      subroutine PRECES(tt)

      implicit none

c
c ash/amuchastegui    sept 1968     subroutine preces
c calculate precession matrix and mean obliquity and also
c derivatives of precession w.r.t. time and precession,obliquity
c constants

c parameter
      real*10 tt

c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'leon.inc'
      real*10 dpprec(3, 3), dmprec(3, 3), ddprec(3, 3, 6)
      equivalence (dpprec(1,1), Dhprec(1,1,1)),
     .            (dmprec(1,1), Dhprec(1,1,2)),
     .            (ddprec(1,1,1), Dhprec(1,1,3))
      integer*4 lddprc(6)
      equivalence (lddprc, Lhprc), (lprc, Lhprc(7)),
     .            (loblq, Lhprc(8))
      include 'mnsprt.inc'
      include 'number.inc'
      include 'nutprc.inc'
c
c external function
      real*10 DOTN
c local
      real*10 ang,ca,cx,cy,sa,sx,sy,cmoblq,dhq1,dhq2,dxq1,dxq2,q1,
     .          q2,r,smoblq,t,tcntry,tddp,tf,tref,ttt
      integer*4 i,idegree,j,loblq,lprc
      real*10 z(6),cz(3),sz(3),dhcz(3),dxcz(3),dhsz(3),dxsz(3),
     .          dpranh(3),bias(3,3),prc0(3,3)
      real*10 ttrm(10),ocof(5),zcof(3,5),zcon,bang(3)
c
c        kindnp=0 precession and mean obliquity calculated at jed tt
c        kindnp=1 precession, mean obliquity and their derivatives
c               calculated at jed tt
c
      logical*4 j2000

      if(tt.le.0._10) goto 100
c
c compute precession
c r = julian days from reference epoch
      r = tt - tref
   10 t = r/tcntry
      ttrm(1) = Convds*t
      do i=2,idegree
         ttrm(i) = ttrm(i-1)*t
      end do
      call PRODCT(zcof, ttrm, z, 3, idegree, 1)
      Moblq = Oblq1 + DOTN(ocof, ttrm, idegree)
      if(Jct(21).eq.2) then
         z(1)=z(1)+zcon*Convds
         z(2)=z(2)-zcon*Convds
      endif
      do i = 1, 3
         cz(i) = COS(z(i))
         sz(i) = SIN(z(i))
      end do
c
c prec(j,k)= precession matrix  j,k=1,2,3
      q1 = cz(2)*cz(3)
      q2 = sz(2)*cz(3)
      Prec(1,1) = cz(1)*q1 - sz(1)*sz(2)
      Prec(1,2) = -sz(1)*q1 - cz(1)*sz(2)
      Prec(1,3) = -cz(2)*sz(3)
      Prec(2,1) = cz(1)*q2 + sz(1)*cz(2)
      Prec(2,2) = cz(1)*cz(2) - sz(1)*q2
      Prec(2,3) = -sz(2)*sz(3)
      Prec(3,1) = cz(1)*sz(3)
      Prec(3,2) = -sz(1)*sz(3)
      Prec(3,3) = cz(3)
c save basic precession matrix before applying bias or small rotations
      do i=1,3
         do j=1,3
            prc0(i,j)=Prec(i,j)
         end do
      end do
c include frame bias if applicable
      if(Jct(21).eq.2) call FULROT(Prec,bias)
c
c if just initializing, then all we needed was the precession matrix
      if(tt.eq.0._10) then
         do i=1,3
            do j=1,3
               Prec2000(i,j)=Prec(i,j)
            end do
         end do
         return
      endif
c
c calculate matrix from reference eq/eqn to mean ecliptic of date
c for moon rotation
      if(Klan.eq.17) then
         cmoblq = COS(Moblq)
         smoblq = SIN(Moblq)
         do j = 1, 3
            A(1,j) = Prec(1,j)
            A(2,j) = Prec(2,j)*cmoblq + Prec(3,j)*smoblq
            A(3,j) = Prec(3,j)*cmoblq - Prec(2,j)*smoblq
         end do
      endif
c  Note that the precession matrix used to calculate the Moon
c  rotation does not include the ad hoc small rotations used to
c  rotate Earth sites--this allows separate determination of Earth
c  and Moon orientation.
c        partial of precesion matrix with respect
c        to precession constant (seconds of arc per century)
      if(lprc.gt.0) then
         dhcz(1) = -sz(1)*t*Fact1
         dhsz(1) = cz(1)*t*Fact1
         dhcz(2) = -sz(2)*t*Fact1
         dhsz(2) = cz(2)*t*Fact1
         dhcz(3) = -sz(3)*t*Fact2
         dhsz(3) = cz(3)*t*Fact2
         dhq2    = dhsz(2)*cz(3) + sz(2)*dhcz(3)
         dhq1    = dhcz(2)*cz(3) + cz(2)*dhcz(3)
         dpprec(1,1)= dhcz(1)*q1 + cz(1)*dhq1 - dhsz(1)*sz(2)
     .                - sz(1)*dhsz(2)
         dpprec(1,2)= -dhsz(1)*q1 - sz(1)*dhq1 - dhcz(1)*sz(2)
     .                - cz(1)*dhsz(2)
         dpprec(1,3)= -dhcz(2)*sz(3) - cz(2)*dhsz(3)
         dpprec(2,1)= dhcz(1)*q2 + cz(1)*dhq2 + dhsz(1)*cz(2)
     .                + sz(1)*dhcz(2)
         dpprec(2,2)= dhcz(1)*cz(2) + cz(1)*dhcz(2) - dhsz(1)*q2
     .                - sz(1)*dhq2
         dpprec(2,3)= -dhsz(2)*sz(3) - sz(2)*dhsz(3)
         dpprec(3,1)= dhcz(1)*sz(3) + cz(1)*dhsz(3)
         dpprec(3,2)= -dhsz(1)*sz(3) - sz(1)*dhsz(3)
         dpprec(3,3)= dhcz(3)

         if(Jct(21).eq.2) call FULROT(dpprec,bias)
      endif
c
c partial of precession matrix with respect
c to mean obliquity of ecliptic (seconds of arc)
      if(loblq.gt.0) then
         dxcz(1) = sz(1)*t*Fact3
         dxsz(1) = -cz(1)*t*Fact3
         dxcz(2) = sz(2)*t*Fact3
         dxsz(2) = -cz(2)*t*Fact3
         dxcz(3) = -sz(3)*t*Fact4
         dxsz(3) = cz(3)*t*Fact4
         dxq1    = dxcz(2)*cz(3) + cz(2)*dxcz(3)
         dxq2    = dxsz(2)*cz(3) + sz(2)*dxcz(3)
         dmprec(1,1) = dxcz(1)*q1 + cz(1)*dxq1 - dxsz(1)*sz(2)
     .                  - sz(1)*dxsz(2)
         dmprec(1,2) = -dxsz(1)*q1 - sz(1)*dxq1 - dxcz(1)*sz(2)
     .                  - cz(1)*dxsz(2)
         dmprec(1,3) = -dxcz(2)*sz(3) - cz(2)*dxsz(3)
         dmprec(2,1) = dxcz(1)*q2 + cz(1)*dxq2 + dxsz(1)*cz(2)
     .                     + sz(1)*dxcz(2)
         dmprec(2,2) = dxcz(1)*cz(2) + cz(1)*dxcz(2) - dxsz(1)
     .                     *q2 - sz(1)*dxq2
         dmprec(2,3) = -dxsz(2)*sz(3) - sz(2)*dxsz(3)
         dmprec(3,1) = dxcz(1)*sz(3) + cz(1)*dxsz(3)
         dmprec(3,2) = -dxsz(1)*sz(3) - sz(1)*dxsz(3)
         dmprec(3,3) = dxcz(3)

         if(Jct(21).eq.2) call FULROT(dmprec,bias)
      endif
c
c modify precession matrix by small rotations
      if(Kddprc.gt.0) then
         tddp = t*Convds
         do i = 1,3
            dpranh(i) = Dprang(i)*tddp + Dprang(i+3)*Convds
         end do
c
c first calculate partials w.r.t. small rotations
         if(lddprc(1).gt.0) then
            do i = 1, 3
               ddprec(i,1,1) = 0._10
               ddprec(i,2,1) = -prc0(i,3)*tddp
               ddprec(i,3,1) = prc0(i,2)*tddp
            end do
            if(Jct(21).eq.2) call FULROT(ddprec(1,1,1),bias)
         endif
         if(lddprc(2).gt.0) then
            do i = 1, 3
               ddprec(i,1,2) = prc0(i,3)*tddp
               ddprec(i,2,2) = 0._10
               ddprec(i,3,2) = -prc0(i,1)*tddp
            end do
            if(Jct(21).eq.2) call FULROT(ddprec(1,1,2),bias)
         endif
         if(lddprc(3).gt.0) then
            do i = 1, 3
               ddprec(i,1,3) = -prc0(i,2)*tddp
               ddprec(i,2,3) = prc0(i,1)*tddp
               ddprec(i,3,3) = 0._10
            end do
            if(Jct(21).eq.2) call FULROT(ddprec(1,1,3),bias)
         endif
         if(lddprc(4).gt.0) then
            do i = 1, 3
               ddprec(i,1,4) = 0._10
               ddprec(i,2,4) = -prc0(i,3)*Convds
               ddprec(i,3,4) = prc0(i,2)*Convds
            end do
            if(Jct(21).eq.2) call FULROT(ddprec(1,1,4),bias)
         endif
         if(lddprc(5).gt.0) then
            do i = 1, 3
               ddprec(i,1,5) = prc0(i,3)*Convds
               ddprec(i,2,5) = 0._10
               ddprec(i,3,5) = -prc0(i,1)*Convds
            end do
            if(Jct(21).eq.2) call FULROT(ddprec(1,1,5),bias)
         endif
         if(lddprc(6).gt.0) then
            do i = 1, 3
               ddprec(i,1,6) = -prc0(i,2)*Convds
               ddprec(i,2,6) = prc0(i,1)*Convds
               ddprec(i,3,6) = 0._10
            end do
            if(Jct(21).eq.2) call FULROT(ddprec(1,1,6),bias)
         endif
c
c now multiply precession by small rotations
         if(Kindnp.le.0) call SMLROT(Prec,dpranh)
         if(lprc.gt.0) call SMLROT(dpprec,dpranh)
         if(loblq.gt.0) call SMLROT(dmprec,dpranh)
      endif
c
c dprec(j,k)= time derivative of precession matrix j,k=1,2,3
c   units of radian/sec
      if(Kindnp.gt.0) then
         ttt = Convds/tcntry/86400._10
         ttrm(6) = ttt
         do i=2,idegree
            ttt=ttt*t
            ttrm(i+5) = i*ttt
         end do
         Dmoblq = DOTN(ocof,ttrm(6),idegree)
         call PRODCT(zcof,ttrm(6),z(4),3,idegree,1)
         Dprec(1,3) = -prc0(2,3)*z(5) - q1*z(6)
         Dprec(2,3) = prc0(1,3)*z(5) - q2*z(6)
         Dprec(3,3) = -sz(3)*z(6)
         q1 = prc0(1,3)*z(6)
         q2 = prc0(2,3)*z(6)
         Dprec(1,1) = prc0(1,2)*z(4) + q1*cz(1) - prc0(2,1)*z(5)
         Dprec(1,2) = -prc0(1,1)*z(4) - q1*sz(1) - prc0(2,2)*z(5)
         Dprec(2,1) = prc0(1,1)*z(5) + q2*cz(1) + prc0(2,2)*z(4)
         Dprec(2,2) = prc0(1,2)*z(5) - q2*sz(1) - prc0(2,1)*z(4)
         q1 = cz(3)*z(6)
         Dprec(3,1) = prc0(3,2)*z(4) + cz(1)*q1
         Dprec(3,2) = -prc0(3,1)*z(4) - sz(1)*q1
         if(Jct(21).eq.2) call FULROT(Dprec,bias)
c
c multiply precession and its derivative by small rotations
         if(Kddprc.gt.0) then
            call SMLROT(Prec,dpranh)
            call SMLROT(Dprec,dpranh)
         endif
      endif

      if(mod(Jct(6)/128,2).ne.0) then
         if(Line.gt.56) call OBSPAG
         ttt = tt + 0.5_10
         write(Iout,20) ttt, Moblq, (z(i),i = 1,3)
   20    format(' PRECES:    JD.F=', f17.9, ' MOBLQ (RADIANS)=',
     .          f19.16, ' ANGLES (RADIANS)=', 1p, 3D15.7)
         Line = Line + 1
      endif
      if(mod(Jct(6)/512,2).ne.0) then
         if(Line.gt.54) call OBSPAG
         write(Iout,40) ((Prec(i,j),j=1,3),i = 1,3)
   40    format(' PRECES:    PREC=', (t18,3F20.16))
         Line = Line + 3
      endif
      return
c
c initialization
  100 j2000 = (Jct(13).gt.0)
      tref = 2433282.423_10
      if(Ercond(29).gt.0._10) Oblq1 = Ercond(29)*Convds
c
c see if old or new IAU precession quantities used

      if(Jct(21).lt.0) then
c use old expressions (tropical centuries, etc.)
         tcntry = 36524.21988_10
         if(Ercond(29).le.0._10) Oblq1 = 84404.8363494530058_10*Convds
         idegree=3
         zcof(1,1) = 2304.948_10 + Hxz1
         zcof(1,2) = 0.302_10
         zcof(1,3) = 0.0179_10
         zcof(2,1) = zcof(1,1)
         zcof(2,2) = 1.093_10
         zcof(2,3) = 0.0192_10
         zcof(3,1) = 2004.255_10 + Hxz2
         zcof(3,2) = -0.426_10
         zcof(3,3) = -0.0416_10
         ocof(1) = -46.8485418486769944_10
         ocof(2) = -3.18487539448290358E-3_10
         ocof(3) = 1.80988402545302038E-3_10
      else if(Jct(21).eq.0) then
c
c use IAU(1976) expressions
c (from Lieske et al., Astron. Astrophys. 58, 1-16, 1977)
         if(j2000) tref = 2451545._10
         tcntry = 36525._10
         tf     = (tref - 2451545._10)/tcntry
c at 1950.0 this gives oblq1= 84404.85522
         if(Ercond(29).le.0._10)
     .    Oblq1 = (84381.448_10 - 46.8150_10*tf - 0.00059_10*tf**2 +
     .    0.001813_10*tf**3)*Convds
         idegree=3
         ocof(1) = -46.8150_10 - 0.00117_10*tf + 0.005439*tf**2
         ocof(2) = -0.00059_10 + 0.005439_10*tf
         ocof(3) = 0.001813_10
c
c zeta=z(1),  z=z(2),  theta=z(3)
         zcof(1,1)= 2306.2181_10 + Hxz1 + 1.39656_10*tf -
     .    0.000139_10*tf**2
         zcof(1,2)= 0.30188_10 - 0.000345_10*tf
         zcof(1,3)= 0.017998_10
         zcof(2,1)= zcof(1,1)
         zcof(2,2)= 1.09468_10 + 0.000066_10*tf
         zcof(2,3)= 0.018203_10
         zcof(3,1)= 2004.3109_10 + Hxz2 - 0.85330_10*tf -
     .    0.000217_10*tf**2
         zcof(3,2)= -0.42665_10 - 0.000217_10*tf
         zcof(3,3)= -0.041833_10

      else
c use IAU 2006 precession model, including frame bias
c (from Capitaine et al., Astron. Astrophys. 412, 567-586, 2003)
         if(j2000) tref = 2451545._10
         tcntry = 36525._10
         tf     = (tref - 2451545._10)/tcntry
         if(Ercond(29).le.0._10) then
            if(j2000) then
               Oblq1 = 84381.406_10*Convds
            else
               Oblq1 = 84404.824187_10*Convds
            endif
         endif
         idegree=5
         ocof(1) = -46.836769_10 - 0.0003662_10*tf + 0.0060102*tf**2
     .    - 2.304E-6_10*tf**3 - 2.170E-7*tf**4
         ocof(2) = -0.0001831_10 + 0.0060102*tf - 3.456E-6_10*tf**2
     .    - 4.34E-7_10*tf**3
         ocof(3) = 0.00200340_10 - 2.304E-6_10*tf - 4.34E-7_10*tf**2
         ocof(4) = -5.76E-7_10 - 2.170E-7*tf
         ocof(5) = -4.34E-8_10
c
c zeta=z(1),  z=z(2),  theta=z(3)
         if(j2000) then
            zcof(1,1)= 2306.083227_10 + Hxz1
            zcof(1,2)= 0.2988499_10
            zcof(1,3)= 0.01801828_10
            zcof(1,4)= -5.971E-6_10
            zcof(1,5)= -3.173E-7_10
            zcof(2,1)= 2306.077181_10 + Hxz1
            zcof(2,2)= 1.0927348_10
            zcof(2,3)= 0.01826837_10
            zcof(2,4)= -2.8596E-5_10
            zcof(2,5)= -2.904E-7_10
            zcof(3,1)= 2004.191903_10 + Hxz2
            zcof(3,2)= -0.4294934_10
            zcof(3,3)= -0.04182264_10
            zcof(3,4)= -7.089E-6_10
            zcof(3,5)= -1.274E-7_10
            zcon = 2.650545_10
            bang(1) = -0.0146_10
            bang(2) = -0.016617_10
            bang(3) = -0.006819_10
         else
            zcof(1,1)= 2305.3878623659_10 + Hxz1
            zcof(1,2)=    0.2989656788_10
            zcof(1,3)=    0.0180250405_10
            zcof(1,4)=   -0.0000058714_10
            zcof(1,5)=   -0.0000003181_10
            zcof(2,1)= 2305.3809697081_10 + Hxz1
            zcof(2,2)=    1.0925747830_10
            zcof(2,3)=    0.0183260358_10
            zcof(2,4)=   -0.0000285924_10
            zcof(2,5)=   -0.0000002906_10
            zcof(3,1)= 2004.6213414132_10 + Hxz2
            zcof(3,2)=   -0.4293752342_10
            zcof(3,3)=   -0.0418246878_10
            zcof(3,4)=   -0.0000070989_10
            zcof(3,5)=   -0.0000001268_10
            zcon     =    2.6473102882_10
            bang(1)  =   -0.0146325090_10
            bang(2)  =   -0.0166917891_10
            bang(3)  =   -0.0065618068_10
         endif

c calculate bias matrix
         ang=bang(1)*Convds
         ca=COS(ang)
         sa=SIN(ang)
         ang=bang(2)*Convds
         cx=COS(ang)
         sx=SIN(ang)
         ang=bang(3)*Convds
         cy=COS(ang)
         sy=SIN(ang)
         bias(1,1)= cx*ca
         bias(1,2)= cx*sa
         bias(1,3)= -sx
         bias(2,1)= -cy*sa-sy*sx*ca
         bias(2,2)= cy*ca-sy*sx*sa
         bias(2,3)= -sy*cx
         bias(3,1)= -sy*sa+cy*sx*ca
         bias(3,2)= sy*ca+cy*sx*sa
         bias(3,3)= cy*cx

      endif
c
c fill common
      Coblq1 = COS(Oblq1)
      Soblq1 = SIN(Oblq1)
      if(j2000) then
         do i=1,3
            do j=1,i-1
               Prec2000(i,j)=0._10
               Prec2000(j,i)=0._10
            end do
            Prec2000(i,i)=1._10
         end do
         return
      endif
c
c get and save precession matrix to J2000
      r = 2451545._10 - tref
      goto 10
      end
