      subroutine PHSCOR(ibsi,xni,xa,xb,rhai,rhbi)
 
      implicit none

c        subroutine phscor - j.f.chandler, 1977 aug
 
c  Compute approximate effect of solar phase angle on apparent position
c  of occulted or eclipsed body in mutual event.  Correction is applied
c  to that body's position.  The correction is fr*rhb*(1-cos(alpha)),
c  where fr is about .5 (fr depends just on the projected minimum
c  separation and the radius of the other body).

c arguments
      real*10 xni(5),xa(5),xb(5),rhai,rhbi
      integer*4 ibsi
c    ibsi=    index of obscured body coordinates in xscs (saved as ibs)
c    xni=     vector from earth observer to that body (saved as xn)
c              and xn(4),xn(5) are the distance squared and distance
c    xa,xb=   vectors from proper center to obscurer, obscured
c    rhai,rhbi=respective radii (saved as rha,rhb)
 
c      partials computations added 1978 aug

c array dimensions
      include 'globdefs.inc'

c common
      include 'comdateq.inc'
      include 'coord.inc'
 
c phase correction quantities
      common/EQNPHS/ Xpp(3),Rprp,Dfrf,Xrho(3),Xbsv(3),P(3),V(3),
     .        Dxcr(3),Drp(3),Rbrn,Xnx,Rmin,Gamma,Sasq,
     .        Salph,Calph,Dgam(2),Xn(5),Rha,Rhb,
     .        Ias,Ibs,Irh,Jfr
      real*10 Xpp,Rprp,Dfrf,Xrho,Xbsv,P,V,Dxcr,Drp,Rbrn,Xnx,
     .       Rmin,Gamma,Sasq,Salph,Calph,Dgam,Xn,Rha,Rhb
      integer*4 Ias,Ibs,Irh,Jfr
      real*10 rnsq,rn
      equivalence (Xn(4),rnsq),(Xn(5),rn)
      include 'partcm.inc'
      include 'trnocc.inc'

c external functions
      real*10 DOT

c local
      real*10 bb,beta,boa,bsq,dffdr,fr,geom,omdf,pn,
     .      psq,pv,qqa,ra,rb,rbsq,vai,vbi,vn,vsq
      integer   i,ivzb,ivzx,ix,ixa,jbod

c pointers into der.. and ivze for bodies 'a' and 'b'
      integer*2 ixx(2)/37,57/,ivix(2)/3,5/
 
c coefficients for quadratic function for fr
      real*10 c(6)/ -8.883631E-2_10, .1785059_10, -8.996493E-2_10,
     .          -.1340744_10, .1385038_10, .3917556_10/
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      do i=1,5
         Xn(i) = xni(i)
      end do
      Rha = rhai
      Rhb = rhbi
      if(Rhb.gt.0._10) then
         Ibs  = ibsi
         Ias  = 1 - Ibs
 
c get position of obscured body relative to sun
         do i = 1,3
            Xbsv(i) = Xscsun(i,Ibs)
         end do
 
c get projection factor
         if(Ibs.eq.0) then
 
c eclipse, position given is from sun
            ra   = SQRT(DOT(xa,xa))
            rbsq = DOT(xb,xb)
            rb   = SQRT(rbsq)
            bsq  = rbsq
            boa  = rb/ra
         else
 
c occultation, xa and xb include radii
            bsq  = xb(4)
            boa  = xb(5)/xa(5)
            rbsq = DOT(Xbsv,Xbsv)
            rb   = SQRT(rbsq)
         endif
         bb = boa*Rha/Rhb
 
c get relative position and velocity
         do i = 1, 3
            P(i) = boa*xa(i) - xb(i)
            vai  = Xscsun(i + 3,Ias)
            vbi  = Xscsun(i + 3,Ibs)
            if(Ibs.ne.0) then
 
c get velocity rel. to earth
               vai = vai - Xem(i + 3,1)
               vbi = vbi - Xem(i + 3,1)
            endif
            V(i) = boa*vai - vbi
         end do
         pn   = DOT(P,Xn)/rn
         vn   = DOT(V,Xn)/rn
         pv   = DOT(P,V)
         psq  = DOT(P,P)
         vsq  = DOT(V,V)
         Rmin = SQRT(psq - pn**2 - (pv-pn*vn)**2/(vsq-vn**2))/Rhb
 
c make sure obscuration can occur
         if(Rmin - 1._10.le.bb) then
            if(bb - 1.2_10*Rmin.gt.1.38_10) then
 
c use approximate constant
               fr  = 0.41_10
               Jfr = 2
            else
 
c use quadratic function
               Jfr = 1
               fr  = (c(1)*Rmin + c(2)*bb + c(4))*Rmin + (c(3)*bb +
     .               c(5))*bb + c(6)
            endif
            Xnx   = DOT(Xbsv,Xn)
            beta  = Xnx/rnsq
            Rbrn  = rb*rn
            geom  = SQRT((Rbrn-Xnx)/(Rbrn+Xnx))
            Gamma = Rhb*geom*fr/rb
 
c calculate apparent position offset
            do i = 1,3
               Xrho(i) = Gamma*(beta*Xn(i) - Xbsv(i))
            end do
 
c correct position of obscured body
            do i = 1,3
               Xscsun(i,Ibs) = Xscsun(i,Ibs) + Xrho(i)
            end do
            go to 100
         endif
      endif
c
c the correction is zero
      Jfr = 3
      fr  = 0._10
      do i = 1,3
         Xrho(i) = 0._10
      end do
      go to 100
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry PHSCRS
 
c set up for partials computations
      if(Jfr.ne.3) then
         ixa   = ixx(Ias+1)
         ix    = ixx(Ibs+1)
         Calph = Xnx/Rbrn
         Sasq  = 1._10 - Calph**2
         do i = 1,3
            Dxcr(i) = ((1._10-2._10*beta)*Xn(i) + Xbsv(i))/rnsq
            Xpp(i)  = ((Xbsv(i)/rbsq+Xn(i)/rnsq)*Calph - (Xbsv(i)+Xn(i))
     .                /Rbrn)/Sasq - Xbsv(i)/rbsq
         end do
         if(Jfr.eq.1) then
            dffdr = (2._10*c(1)*Rmin + c(2)*bb + c(4))/Rhb/fr
            omdf  = (Rmin*Rhb*dffdr)/bsq
            Dfrf  = dffdr*bsq/Rmin/Rhb
            do i = 1,3
               Xpp(i) = Xpp(i) + omdf*Xscsun(i,Ibs) - Dfrf*Xpr(i,Ibs+1)
            end do
 
c set up for radii partials
            Dgam(Ias+1) = (c(2)*Rmin + 2._10*c(3)*bb + c(5))/Rhb/fr
            Dgam(Ibs+1) = -dffdr*Rmin - Dgam(Ias+1)*bb + 1._10/Rhb
         endif
      endif
      go to 100
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry PHSCRH(jbod)
c set up partial with respect to radius (called from partl1)
c jbod=1 if body is the observed body
c jbod=2 if body is the second body
      Irh = jbod
      return
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      entry PHSCRP
c make corrections to partials of obscured body
c called from cpartc
c derpl and dersc contain partials of observed and other body,
c respectively (or else garbage to be ignored).
      if(Jfr.ne.3) then
         ivzx = ivix(Ibs+1)
         ivzb = Ivze(ivzx)
         if(ivzb.gt.0) then
 
c partials array contains info
            qqa  = DOT(Dxcr,Derpr(ix,1,1))
            Rprp = DOT(Xpp,Derpr(ix,1,1))
         else if(Irh.lt.0) then
 
c short cut for zero dot products
            qqa  = 0._10
            Rprp = 0._10
         else
 
c partial is w.r.t. a body radius
            if(Jfr.ne.1) go to 100
 
c no non-zero partials  yet (i.e. ivze(.)=0)
            Rprp = Dgam(Irh)/aukm
         endif
         if(Jfr.eq.1 .and. Ivze(ivix(Ias+1)).gt.0) Rprp = Rprp -
     .       Dfrf*DOT(Xpr(1,Ias+1),Derpr(ixa,1,1))
         do i = 1,3
            Drp(i) = Rprp*Xrho(i)
            if(ivzb.gt.0) Drp(i) = Drp(i)
     .                                 + Gamma*(qqa*Xn(i) + (beta-1._10)
     .                                 *Derpr(i+ix-1,1,1))
            Derpr(i+ix-1,1,1) = Derpr(i + ix - 1,1,1) + Drp(i)
         end do
         Ivze(ivzx) = 1
      endif
  100 Irh = -999
      return
      end
