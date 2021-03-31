      subroutine LUNORB(jd,fract,nn)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i,int,int1,j,jd,k,mess,nacc,nder,nn,nn1
 
c*** end of declarations inserted by spag
 
 
 
c
c m.e.ash   sept 1965    subroutine lunorb
c determination of position,velocity,acceleration in browns mean
c lunar orbit and the partial derivatives of these quantities with
c respect to initial mean orbital elements
c
c        jd = julian day number (jed of noon on given day)
c      fract= fraction of day from midnight
c        nn = 0 position,velocity,acceleration determined
c        nn = 1 position,velocity,acceleration and partial derivatives
c               of position,velocity,accelertion determined
c        nn =-1 position,velocity,acceleration and partial derivatives
c               of position,velocity determined
      real*10 fract
c
c           ylun(i),i=1,9 position,velocity,acceleration in mean lunar
c                         orbit
c           dylun(i,j), i=1,9, j=1,6  partial derivatives of position,
c                         velocity,acceleration in mean luar orbit
c                         with respect to initial mean orbital elements
c
c        common
      include 'aacoff.inc'
      include 'funcon.inc'
      include 'orblun.inc'
      include 'orbstf.inc'
      include 'precmn.inc'
c
c miscellaneous specification statements
      real*10 bb(3,2),b(3,2),db(3,2),ddb(3,2),daa(3,3),
     .       ddaa(3,3),dc(3,2),ddc(3,2)
      real*10 pp1,pp2,pp3,pp4,pp5,pp6,pp7
      equivalence (bb(1,1),pp1), (bb(2,1),pp2), (bb(3,1),pp3),
     1 (bb(1,2),pp4), (bb(2,2),pp5), (bb(3,2),pp6), (ww(7),pp7)
      real*10 qq,qq1,qq2,qq3,qq4,qq5,qq6,qq7,qq8,qq9,
     .       secc,cecc,danom,ddanom,danom2,
     .       dper,ddper,pa1,pa2,pa3,spci,cpci,
     .       stuff(3),v(6),dv(2),ddv(2)
      equivalence (v(3),dv), (v(5),ddv)
c
c determination of lunar orbit quantities and AA matrix
      call ECLPRC(jd,fract,0)
c
c determination of eccentric anomaly (solution of kepler eq)
      secc = Anomx + E*(SIN(Anomx) + E22*SIN(Anomx*2._10))
  110 cecc = (Anomx - secc + E*SIN(secc))/(1._10 - E*COS(secc))
      secc = secc + cecc
      if(ABS(cecc).gt.1.E-14_10) goto 110
      cecc = COS(secc)
      secc = SIN(secc)
c
c determination of derivatives of mean anomaly
      danom = 2.2802713493961401E-1_10 +
     .         Tt*(2.404714643398E-13_10 + Tt*1.5655603390389E-20_10)
      ddanom = 2.404714643398E-13_10 + Tt*3.1311206780778E-20_10
      danom2 = danom**2
c
c determination of derivatives of perigee
      dper  = 2.868588295626071E-3_10 -
     .        Tt*(3.2449161453178E-13_10 + Tt*1.62315620435472E-20_10)
      ddper = -3.2449161453178E-13_10 - Tt*3.24631240870944E-20_10
      pa3   = dper*dper
      pa1   = Ascd*Ascd + pa3
      pa2   = Ascd*dper*2._10
c
c determination of position,velocity,acceleration in
c mean lunar orbital plane
      qq     = 1._10 - E*cecc
      qq1    = cecc - E
      qq2    = A/qq
      qq3    = -secc*qq2
      qq4    = Erta/qq
      qq5    = cecc*qq4
      qq6    = qq*qq
      qq7    = qq1/qq6
      qq8    = secc/qq6
      v(1)   = A*qq1
      v(2)   = Erta*secc
      dv(1)  = qq3*danom
      dv(2)  = qq5*danom
      ddv(1) = qq2*(-secc*ddanom - qq7*danom2)
      ddv(2) = qq4*(cecc*ddanom - qq8*danom2)
c
c determination of AA matrix derivatives
c (AA matrix calculated in ECLPRC)
      do j = 1, 3
         do i = 1, 3
            daa(i,j) = Daa1(1,i,j) + Tpr*(Daa1(2,i,j) + Tpr*(Daa1(3,i,j)
     .                 +Tpr*(Daa1(4,i,j)+Tpr*(Daa1(5,i,j)))))
            ddaa(i,j) = Ddaa1(1,i,j)
     .                 + Tpr*(Ddaa1(2,i,j) + Tpr*(Ddaa1(3,i,j)
     .                 +Tpr*Ddaa1(4,i,j)))
         end do
      end do
c
c determination of B matrix
      Casc   = COS(Asc)
      Sasc   = SIN(Asc)
      Cper   = COS(Per)
      Sper   = SIN(Per)
      spci   = Sper*Cinc
      cpci   = Cper*Cinc
      b(1,1) = Casc*Cper - Sasc*spci
      b(1,2) = -Casc*Sper - Sasc*cpci
      b(2,1) = Sasc*Cper + Casc*spci
      b(2,2) = -Sasc*Sper + Casc*cpci
      b(3,1) = Sper*Sinc
      b(3,2) = Cper*Sinc
      nn1    = 0
      nacc   = 3
      int1   = 1
c
c determination of time derivatives of B matrix
c
c in this region, the B matrix may be replaced by its partial derivative
c with respect to an orientation angle, with successive passes to
c compute all such partials
  300 db(1,1) = -b(2,1)*Ascd + b(1,2)*dper
      db(1,2) = -b(2,2)*Ascd - b(1,1)*dper
      db(2,1) = b(1,1)*Ascd + b(2,2)*dper
      db(2,2) = b(1,2)*Ascd - b(2,1)*dper
      db(3,1) = b(3,2)*dper
      db(3,2) = -b(3,1)*dper
      if(nn1.ge.0) then
         ddb(1,1)=-b(2,1)*Ascdd + b(1,2)*ddper - b(1,1)*pa1 - b(2,2)*pa2
         ddb(1,2)=-b(2,2)*Ascdd - b(1,1)*ddper - b(1,2)*pa1 + b(2,1)*pa2
         ddb(2,1)= b(1,1)*Ascdd + b(2,2)*ddper - b(2,1)*pa1 + b(1,2)*pa2
         ddb(2,2)= b(1,2)*Ascdd - b(2,1)*ddper - b(2,2)*pa1 - b(1,1)*pa2
         ddb(3,1)= b(3,2)*ddper - b(3,1)*pa3
         ddb(3,2)=-b(3,1)*ddper - b(3,2)*pa3
      endif
c
c determination of C matrix and its time derivatives
      do i = 1, 3
         do j = 1, 2
            C(i,j) = 0._10
            dc(i,j) = 0._10
            ddc(i,j) = 0._10
            do k = 1, 3
               C(i,j) = C(i,j) + Aa(k,i)*b(k,j)
               dc(i,j) = dc(i,j) + (Aa(k,i)*db(k,j) + daa(k,i)*b(k,j))
               if(nn1.ge.0) ddc(i,j) = ddc(i,j) + ((Aa(k,i)*ddb(k,j)+
     .                      ddaa(k,i)*b(k,j)) + 2._10*daa(k,i)*db(k,j))
            end do
         end do
      end do
c
c determination of position,velocity,acceleration for mean
c lunar orbit in coordinate system referred to mean equinox
c and equator of reference epoch
 
  400 do int=1,3
 
c matrix transformation
         do k = 1, nacc
            stuff(k) = 0._10
         end do
         do k = 1, 2
            stuff(1) = stuff(1) + C(int,k)*v(k)
            stuff(2) = stuff(2) + (C(int,k)*dv(k) + dc(int,k)*v(k))
            if(nn1.ge.0) stuff(3) = stuff(3) +
     .         C(int,k)*ddv(k) + 2._10*dc(int,k)*dv(k) + ddc(int,k)*v(k)
         end do

         if(int1.eq.1) then
c store coordinates
            do i = 1, 3
               j = 3*(i - 1) + int
               Ylun(j) = stuff(i)
            end do
         else if(int1.eq.2) then
c store partials w.r.t. a, e, or anom
            do i = 1, nacc
               j = 3*(i - 1) + int
               Dylun(j,nder) = stuff(i)
            end do
         else
c store partials w.r.t. orientation angles
            do i = 1, nacc
               j = 3*(i - 1) + int
               do k = 3, 5
                  Dylun(j,k) = Dylun(j,k) + Dintl(nder,k-2)*stuff(i)
               end do
            end do
         endif
      end do
 
      if(int1.eq.1) then

c just doing motion, no partials
         if(nn.eq.0) return
         if(nn.lt.0) nacc = 2
         nn1 = nn
         if(Luncon.ne.1) goto 500
         return
      else if(int1.eq.2) then

c doing partials w.r.t. a, e, or anom. go on to the next
         if(nder.lt.2) then

c partial derivatives with respect to eccentricity
            nder  = 2
            pp7   = Ert/qq
            pp2   = cecc*cecc
            pp3   = pp6/qq
            v(1)  = -A + qq3*secc
            v(2)  = secc*(qq5 - Ertae)
            dv(1) = qq3*(cecc + pp5)/qq
            dv(2) = qq2*(-Erte*cecc + pp7*(pp2-pp3))
            if(nn1.gt.0) then
               pp1    = A*qq9*danom
               ddv(2) = ddanom*dv(2)
     .                  + pp1*secc*(Erte - pp7*(4._10*cecc-pp4))
               pp1    = pp1/qq
               ddv(1) = ddanom*dv(1)
     .                  + pp1*(-qq1*qq1 - Ert2 + Ert23*pp3)
            endif
            dv(1) = danom*dv(1)
            dv(2) = danom*dv(2)
            goto 400
         else if(nder.eq.2) then
c
c partial derivatives w.r.t. initial mean anomaly
            nder  = 6
            v(1)  = qq3
            v(2)  = qq5
            dv(1) = -qq9*qq1*A
            dv(2) = -qq9*Erta*secc
            if(nn1.gt.0) then
               ddv(1) = -qq2*qq7*ddanom + pp1*secc*(1._10 + qq6)
               ddv(2) = -qq4*qq8*ddanom + pp1*Ert*(-cecc + pp4)
            endif
            goto 400
         else

c finished partials that use the same B matrix with partial derivatives
c of v,dv,ddv

c restore values of v,dv,ddv
            do i = 1, 6
               v(i) = Ww(i)
            end do
 
c save values of B
            do i = 1, 3
               do j = 1, 2
                  bb(i,j) = b(i,j)
               end do
            end do
 
c initilize dylun(i,j), j=3,4,5, to zero
            mess = 3*nacc
            do j = 3, 5
               do i = 1, mess
                  Dylun(i,j) = 0._10
               end do
            end do
            int1 = 3
c
c partial derivatives w.r.t. initial inclination, ascending
c node, perigee relative to mean equinox and equator of ref. epoch
c calculated by substituting the derivatives of the B matrix for
c the original B matrix

c partial derivative w.r.t. inclination
            nder   = 1
            b(1,1) = bb(3,1)*Sasc
            b(1,2) = bb(3,2)*Sasc
            b(2,1) = -bb(3,1)*Casc
            b(2,2) = -bb(3,2)*Casc
            b(3,1) = spci
            b(3,2) = cpci
            goto 300
         endif
      else
         if(nder.lt.2) then

c partial derivative w.r.t. ascending node
            nder = 2
            do i = 1, 2
               b(1,i) = -bb(2,i)
               b(2,i) = bb(1,i)
               b(3,i) = 0._10
            end do
            goto 300
         else if(nder.eq.2) then

c partial derivative w.r.t. perigee
            nder = 3
            do i = 1, 3
               b(i,1) = bb(i,2)
               b(i,2) = -bb(i,1)
            end do
            goto 300
         else
c
c adjustment of partial derivatives w.r.t. ascending node and
c perigee
            do i = 1, mess
               Dylun(i,4) = Dylun(i,4) - Dylun(i,5)
               Dylun(i,5) = Dylun(i,5) - Dylun(i,6)
            end do
            return
         endif
      endif

c
c * * * entry lunbro * * *
c
      entry LUNBRO
c
c determination of partial derivatives of position,velocity,
c (if nn1=1,-1), and acceleration (if nn1=1) in mean lunar
c orbit
  500 int1 = 2
 
c save values of v,dv,ddv
      do i = 1, 6
         Ww(i) = v(i)
      end do
c
c partial derivatives with respect to semi-major axis
      nder  = 1
      qq9   = Mot32/qq
      pp1   = qq9*secc
      pp2   = qq9*cecc*Ert
      pp7   = Tt - T0
      pp3   = pp7*pp1
      pp4   = pp7*pp2
      pp5   = qq1/qq
      pp6   = secc*secc
      pp7   = pp7*qq9
      v(1)  = qq1 + pp3
      v(2)  = Ert*secc - pp4
      dv(1) = (secc - pp7*pp5)/qq
      dv(2) = (cecc + pp3/qq)/qq*Ert
      qq9   = danom/qq6/qq
      if(nn1.gt.0) then
         pp4    = E3/qq
         qq6    = qq1*pp4
         pp4    = pp6*pp4
         ddv(1) = -ddanom*dv(1)
     .            + qq9*(Mot3*qq1 - danom*(qq1+pp3*(1._10+qq6)))
         ddv(2) = ddanom*dv(2)
     .            + qq9*Ert*(Mot3*secc - danom*(secc+pp7*(pp4-cecc)))
      endif
      dv(1) = pp1 - danom*dv(1)
      dv(2) = -pp2 + danom*dv(2)
      goto 400
 
      end
