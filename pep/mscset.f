      subroutine MSCSET(ncentr)

c m.ash   jan 1993    setup central body masscon quantities

      implicit none
      integer*2 ncentr

c array dimensions
      include 'globdefs.inc'

c common
      include 'funcon.inc'
      include 'inodta.inc'
      include 'mascon.inc'
      include 'mascn1.inc'
      include 'param.inc'

      integer*4 i,j
      real*10 clong,slong,clat,slat,aukm

c everything zero initially
      Ms0msc = 0._10
      do i=1,3
         X0msc(i) = 0._10
         do j=1,100
            Xmsc(i,j) = 0._10
         end do
      end do
      Fcnmsc = 1._10
      Nmmsc1 = 0

c decide if masscons are included in force model of this numerical integration
      if (Nummsc.le.0) return
      if (ncentr.le.0) return
      if (ncentr.ne.Nplmsc) return

c setup central body fixed masscon cartesian coordinates in kilometers
      Nmmsc1 = Nummsc
      do j=1,Nmmsc1
         clong  = Lngmsc(j)*Convd
         slong  = SIN(clong)
         clong  = COS(clong)
         clat   = Latmsc(j)*Convd
         slat   = SIN(clat)
         clat   = COS(clat)
         Xmsc(1,j)  = Rkmmsc(j)*clong*clat
         Xmsc(2,j)  = Rkmmsc(j)*slong*clat
         Xmsc(3,j)  = Rkmmsc(j)*slat
         Ms0msc    = Ms0msc + Masmsc(j)
         do i=1,3
            X0msc(i)  = X0msc(i) + Masmsc(j)*Xmsc(i,j)
         end do
      end do
      R0msc   = SQRT(x0msc(1)**2 + x0msc(2)**2 + x0msc(3)**2)
      Fcnmsc = 1._10 - Ms0msc

c print out masscon information in fractions of cental body mass, km
      write(Iout,55) Nplmsc,Nmmsc1,Ms0msc,Fcnmsc,R0msc,x0msc
   55 format(/' CENTRAL BODY',i3,' HAS',i4,' MASCONS INCLUDED IN FORCE',
     .' MODEL OF SATELLITE MOTION NUMERICAL INTEGRATION '/
     .' TOTAL FRACTION OF MASCON MASS RELATIVE TO CENTRAL BODY MASS ='
     . 1pd20.13/' TO COMPENSATE CENTRAL BODY FORCE REDUCED BY',1pd20.13,
     .' OFFSET OF MASCON CENTER OF MASS FROM BODY CENTER OF MASS =',
     . 1pd14.7,' KM  (X,Y,Z = ',1p,3d12.5,')')

c convert masscon body fixed cartesian coordinates to astronomical units
      aukm  = Aultsc*Ltvel
      do j=1,Nmmsc1 
         do i=1,3
            Xmsc(i,j)  = Xmsc(i,j)/aukm
         end do
      end do

      return
      end
