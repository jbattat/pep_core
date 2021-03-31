      subroutine ROTAB(a,b,t)
 
      implicit none
c
c        r. goldstein   july 1977
c        calculates the expansions needed to obtain the inclination
c           and precession angles
c        ref:  reasenberg and king, rotation of mars, jgr
c              memo by rdr dated 6/22/77
c
      include 'rotcom.inc'
      real*10 omeg3,cosi,sini,m
      equivalence (omeg3,Trig(25)),(cosi,Trig(15)),(sini,Trig(14)),
     .            (Trig(29),m)
c
c        units of input parameters
c           t:  days (jd.fr) from con1(1)
c           i:  radians
c           omeg3:  radians/day
c
      real*10 skm(6),ckm(6),skmq(6),ckmq(6),skmk(6),ink(6)
      real*10 a,b,c1,c2,fk,q,s1,s2,s3,s4,t,tani
      integer   k,k1,ks
c
c initial calculations
      tani = sini/cosi
      q    = Qqdot(1) + Qqdot(2)*t
      c1   = Mmot*cosi/omeg3
      c2   = Mmot/cosi/omeg3
c
c build up arrays containing sin(km+q),cos(km+q),sin(km)/k,
c sin(km),cos(km), and 1/k
c
      skm(1)  = SIN(m)
      ckm(1)  = COS(m)
      ckmq(1) = COS(m + q)
      skmq(1) = SIN(m + q)
      skmk(1) = SIN(m)
      ink(1)  = 1._10
      do k = 2, 6
         k1     = k - 1
         fk     = k
         ink(k) = 1._10/fk
         ckm(k) = ckm(k1)*ckm(1) - skm(k1)*skm(1)
         skm(k) = skm(k1)*ckm(1) + ckm(k1)*skm(1)
         ckmq(k) = ckmq(k1)*ckm(1) - skmq(k1)*skm(1)
         skmq(k) = skmq(k1)*ckm(1) + ckmq(k1)*skm(1)
         skmk(k) = skm(k)*ink(k)
      end do
c
c build up sums of all series needed k=0 not included
c
      s1 = 0._10
      s2 = 0._10
      s3 = 0._10
      s4 = 0._10
      do k = 1, 6
         ks = k + 1
         s1 = Ecof(ks,2)*ckmq(k)*(ink(k) + c1) + s1
         s2 = Ecof(ks,1)*ckm(k) + s2
         s3 = Ecof(ks,1)*skmk(k) + s3
         s4 = Ecof(ks,2)*skmq(k)*(ink(k) + c2) + s4
      end do
 
      b = tani*(s1 - c1*(s2+Ecof(1,1)))
 
      a = -(Ecof(1,1)*Mmot*t + s3 - s4)
 
      return
      end
