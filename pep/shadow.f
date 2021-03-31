      subroutine SHADOW
 
      implicit none

c
c beebe,king,reasonberg,preston      6/71  subroutine shadow
c
 
c computation of fraction (lambda) of solar disk seen by
c spacecraft
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'funcon.inc'
      include 'lothrf.inc'
      include 'param.inc'
      include 'prtcod.inc'
      include 'sbthng.inc'

c external functions
      real*10 DOT

c
c lpsx set from lps/prtpin/ which is set from ncentr
c
c quantities internal to the routine
      integer*4 i
      real*10 day/8.6400E4_10/
      real*10 aultda,rsbx,sbcorx(3),area3,tlag,area1,area2,ari,rs,rp,
     .       sep,phi,r1,r2,thet,hgt,sepp(3),ubcor(3),pcvel(3)
c
c          shadow computation - geometric model
c          lambda = 1 - no shadow
c          lambda = 0 - no sunlight
c          0.lt.lambda.lt.1 - partial shadow
c
c        no consideration is given to the change of direction associated
c        with partial shadow.
      Lambda = 1.0_10
      if(Rc.le.0.0_10) call SUICID(
     .' SUN MUST BE PERTURBING FORCE TO USE RADIATION PRESSURE, STOP SHA
     .DOW',17)
      if(Rb.gt.Rc) then
 
c adjust sbcor for movement of planet since light passed it
         aultda = Aultsc/day
         do i = 1, 3
            ubcor(i) = Bcor(i)/Rb
            pcvel(i) = Xpert(i + 3,Lpsx)
         end do
c
c        tlag is the light propagation time (days) from the central
c        body to the s/c.
c        sbcorx is the position of the s/c wrt the central body---
c        s/c at current time / central body earlier by tlag
c
         tlag = aultda*DOT(Sbcor,ubcor)
         do i = 1, 3
            sbcorx(i) = Sbcor(i) + tlag*pcvel(i)
         end do
         call CROSS(sbcorx,ubcor,sepp)
 
c change rsbx to its projection along bcor
         rsbx = DOT(sbcorx,ubcor)
c rs, rp are apparent (from probe) radii of sun and planet
c sep is apparent separation of their centers
         rs  = Sunrad/(Rb*Aultsc*Ltvel)
         rp  = Pcrad/(rsbx*Aultsc*Ltvel)
         sep = SQRT(sepp(1)**2 + sepp(2)**2 + sepp(3)**2)/rsbx
         if(rs+rp.gt.sep) then
            if(rp-rs.ge.sep) then
 
c sun is completely eclipsed by planet
               Lambda = 0.0_10
            else if(sep.le.rs-rp) then
 
c planet lies within sun's disc
               Lambda = (rs**2 - rp**2)/rs**2
            else
 
c set r1 = smaller disc, r2 = larger
               if(rs.gt.rp) then
                  r1 = rp
                  r2 = rs
               else
                  r1 = rs
                  r2 = rp
               endif
 
c phi = 1/2 angle subtended in disc 1 by arc of intersection
               phi = ACOS((r1*r1+sep*sep-r2*r2)/(2.0_10*r1*sep))
               if(r2/r1.gt.5.0_10) then
 
c one disc much bigger - treat boundary as a straight line
                  hgt   = SQRT(r1**2 - (sep-r2)**2)
                  area2 = hgt*(sep - r2)
                  area3 = 0.0_10
               else
c thet = 1/2 angle subtended in disc 2 by arc of intersection
c hgt = 1/2 linear distance between ends of arc of intersection
                  hgt   = r1*SIN(phi)
                  thet  = ASIN(hgt/r2)
                  area2 = sep*hgt
                  area3 = thet*r2**2
               endif
               area1 = (Pi - phi)*r1**2
 
c ari = area of non-overlapped portion of small disc
               ari   = area1 + area2 - area3
               area1 = Pi*rs**2
               if(rs.gt.rp) then
 
c planet is small disc
                  area2  = Pi*rp**2
                  Lambda = (area1 + ari - area2)/area1
               else
 
c sun is small disc
                  Lambda = ari/area1
               endif
            endif
         endif
      endif
c
c write(6,101)lambda
c101  format(' shadow... lambda= ',f7.3)
      return
      end
