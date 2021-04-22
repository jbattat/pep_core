      subroutine rotelt(nplnt,cond,belpt,mat,dir,opt)
      implicit none
c subroutine rotelt - J.F.Chandler - 2012 Dec 6
c Transform elliptic orbital elements or rotation state vector
c from one frame to another
c   Input:
c NPLNT - planet number of object: orbit if positive, rotation if neg.
c COND  - array of elements a,e,i,asc,per,anom (angles in degrees)
c         or rotation state psi,theta,phi, and rates (in radians, rad/d)
c BELPT - rotation matrix from orbital plane's frame to reference frame
c         if opt=1, otherwise to be computed here from COND
c         or matrix to go from body-fixed frame to reference frame
c MAT   - rotation matrix between reference frames. Sense depends on DIR
c DIR   - if 3, MAT transforms coordinates from old frame to new
c         if -3, from new frame to old
c OPT   - if 0, compute BELPT here
c         if 1, assume BELPT already filled in
c   Output:
c COND  - contains transformed i,asc,per (inc assumed positive)
c         or transformed rotation state (theta assumed positive)
c BELPT - contains transformed rotation matrix

      real*10 cond(6),belpt(3,3),mat(3,3)
      integer*4 dir,opt
      integer*2 nplnt

c common
      include '../pep/funcon.inc'

c local
      real*10 brot(3,3),cinc,sinc,casc,sasc,cper,sper,cpci,spci,sinc2,
     . dbelpt(3,3),dbrot(3,3),inc,asc,per,dinc,dasc,dper,eulelt(6),
     . peroff
      equivalence (asc,eulelt(1)),(inc,eulelt(2)),(per,eulelt(3)),
     . (dasc,eulelt(4)),(dinc,eulelt(5)),(dper,eulelt(6))
      integer*4 i,j

      if(opt.eq.0) then
         if(nplnt.gt.0) then
            inc= cond(3)*Convd
            asc= cond(4)*Convd
            per= cond(5)*Convd
         else
            inc= cond(2)
            asc= cond(1)
            per= cond(3)
         endif
         cinc = COS(inc)
         sinc = SIN(inc)
         casc = COS(asc)
         sasc = SIN(asc)
         sper = SIN(per)
         cper = COS(per)
         spci = sper*cinc
         cpci = cper*cinc
         belpt(1,1)= casc*cper - sasc*spci
         belpt(1,2)=-casc*sper - sasc*cpci
         belpt(1,3)= sasc*sinc
         belpt(2,1)= sasc*cper + casc*spci
         belpt(2,2)= casc*cpci - sasc*sper
         belpt(2,3)=-casc*sinc
         belpt(3,1)= sper*sinc
         belpt(3,2)= cper*sinc
         belpt(3,3)= cinc
      endif
      call PRODCT(mat,belpt,brot,dir,3,3)
      inc=ACOS(brot(3,3))
      asc=ATAN2(brot(1,3),-brot(2,3))
      per=ATAN2(brot(3,1),brot(3,2))
      if(nplnt.gt.0) then
         cond(3)=inc/Convd
         cond(4)=asc/Convd
         cond(5)=per/Convd
         if(cond(4).lt.0._10) cond(4)=cond(4)+360._10
         if(cond(5).lt.0._10) cond(5)=cond(5)+360._10
      else
         dbelpt(1,1)=-cond(4)*belpt(2,1)+cond(5)*sper*belpt(1,3)
     .    +cond(6)*belpt(1,2)
         dbelpt(1,2)=-cond(4)*belpt(2,2)+cond(5)*cper*belpt(1,3)
     .    -cond(6)*belpt(1,1)
         dbelpt(1,3)=-cond(4)*belpt(2,3)+cond(5)*sasc*cinc
         dbelpt(2,1)= cond(4)*belpt(1,1)+cond(5)*sper*belpt(2,3)
     .    +cond(6)*belpt(2,2)
         dbelpt(2,2)= cond(4)*belpt(1,2)+cond(5)*cper*belpt(2,3)
     .    -cond(6)*belpt(2,1)
         dbelpt(2,3)= cond(4)*belpt(1,3)-cond(5)*casc*cinc
         dbelpt(3,1)= cond(5)*spci+cond(6)*belpt(3,2)
         dbelpt(3,2)= cond(5)*cpci-cond(6)*belpt(3,1)
         dbelpt(3,3)=-cond(5)*sinc
         call PRODCT(mat,dbelpt,dbrot,dir,3,3)
         cinc=brot(3,3)
         sinc2=1._10-cinc**2
         sinc=SQRT(sinc2)
         dinc=-dbrot(3,3)/sinc
         dasc=(-dbrot(1,3)*brot(2,3)+dbrot(2,3)*brot(1,3))/sinc2
         dper=(dbrot(3,1)*brot(3,2)-dbrot(3,2)*brot(3,1))/sinc2
         peroff=(cond(3)-per)/Twopi
         per=per+AINT(peroff+SIGN(0.5_10,peroff))*Twopi
         do i=1,6
            cond(i)=eulelt(i)
         end do
      endif
      do i=1,3
         do j=1,3
            belpt(i,j)=brot(i,j)
         end do
      end do
      return
      end
