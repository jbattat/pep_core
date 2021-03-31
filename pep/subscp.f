      subroutine SUBSCP(idopob)
 
      implicit none
c
c        subscp computes the latitude and longitude of the
c        sub-spacecraft point and the height of the spacecraft.
c        the longitude is put in dstf(2), the lat in dstf(3) and
c        the height in savb(48,1)
c        -180.lt.longitude.le.180 (measured east),
c        -90.le.latitude.le.90
c        height from center of planet in kilometers
c        this subroutine is turned on by ict(55)=1

c arguments
      integer*4 idopob

c array dimensions
      include 'globdefs.inc'

c        common
      include 'comdat.inc'
      real*10 aukm
      equivalence (Comcon(129),aukm)
      include 'coord.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'mnsprt.inc'
      include 'number.inc'
      include 'obscrd.inc'
      include 'sbcom.inc'
      include 'yvect.inc'
c
c external functions
      real*10 DOT

c variables internal to this subroutine
      real*10 dt,dum,rysp,xsbsv(6),frctsv
      integer i,j,jdsve,mnspt1,nso
c-----------------------------------------------------------------------
c
      if(idopob.eq.0) then
c-----------------------------------------------------------------------
c
c instant. doppler
         do j = 1,3
            Xspcd(j,1) = Xsb(j)
         end do
         jdsve  = Jdx
         frctsv = Fract
      else
c phase delay doppler
c
         if(idopob.gt.0) goto 100
c
c second pass of phase delay
         do j = 1,3
            xsbsv(j) = (xsbsv(j) + Xsb(j))/2.0_10
         end do
         dt = (frctsv - Fract + (jdsve-Jdx))/2._10
         call TIMINC(Jdx,Fract,jdsve,frctsv,dt)
c
c xsbsv(j) now contains the average spacecraft position in au
c with respect ot the planet and jdsve.frctsv is the avg time
c
         do j = 1,3
            Xspcd(j,1) = xsbsv(j)
         end do
      endif
c-----------------------------------------------------------------------
c
c
c expand 'save' vector
      if(Numsav.lt.48) then
         nso    = Numsav + 1
         Numsav = 48
         do i=nso,Numsav
            Save(i)=0._10
            if(i.ge.41 .and. i.le.46) Save(i)=-1._10
         end do
      endif
c
c obtain rotation matrix
      mnspt1 = 0
      call SPOTCD(jdsve,frctsv,mnspt1,-1,Ncp0,dum,dum,0)
 
c rotate spacecraft vector
      call PRODCT(Rot,Xspcd,Yspcd,3,3,1)
 
c rotate l.o.s. vector
      if(Jct(58).gt.0) then
         call PRODCT(Rot,Xsitep,Save(49),3,3,1)
         do i = 49,51
            Save(i) = -Save(i)
         end do
         if(Numsav.le.51) Numsav = 51
      endif
 
c compute lat, long, and height
      rysp = SQRT(DOT(Yspcd,Yspcd))
      if(Yspcd(1,1).ne.0.0_10 .or. Yspcd(2,1).ne.0.0_10) then
         Dstf(2) = ATAN2(Yspcd(2,1),Yspcd(1,1))/Convd
      else
         Dstf(2) = 0.0_10
      endif
      Dstf(3) = ASIN(Yspcd(3,1)/rysp)/Convd
 
c compute height in kilometers
      Save(48) = rysp*aukm
c
c first pass of phase delay
  100 do j = 1,3
         xsbsv(j) = Xsb(j)
      end do
      jdsve  = Jdx
      frctsv = Fract
      return
      end
