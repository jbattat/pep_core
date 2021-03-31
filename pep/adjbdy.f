      subroutine ADJBDY(lbdy,conbdy,fctbdy,wrds,wrdl1,wrdl2,name,
     .                  nplnt,jd0,nbdy,ncentr,icnd)
 
      implicit none
 
c
c ash/forni  october 1967  subroutine adjbdy
c adjust earth-moon barycenter, earth rotation, moon, moon rotation,
c planet initial conditions and parameters
c
      character*8 wrds(nbdy),name
      character*4 wrdl1,wrdl2
      real*10 conbdy(nbdy),fctbdy(nbdy)
      integer*2 lbdy(nbdy),nplnt,ncentr,icnd
c           quantities set up in adjust
c           lbdy  = control constants for least squares adjustment
c           conbdy= initial conditions, constants for body
c           wrdl1 = one of char*4 'LEM(' , 'LMN(' , 'LPL('
c           wrdl2 = one of char*4  ', 1)' , ', 2)' , ETC ',16)' , ')   '
c           fctbdy= conversion factors
c           char*8 wrds/'a   (km)e       inc (dg)  etc'/
c           name  = eight character planet name
c           nplnt = planet number
c           jd0   = initial times
c           nbdy  = number of parameters to be printed
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'adjstf.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
c
c internal to subroutine adjbdy
      real*10 dom
      integer i,ict66,int,j,jd0,m,nbdy,maxang
      character*4 astrik(3)/'*   ','&   ','    '/
      character*8 earth/' EARTH  '/,neme,wrdssd
      character*4 wrdsss(2)
      equivalence (wrdssd,wrdsss)
      character*4 wrdss2/'(KM)'/
 
      neme = name
      call LINCHK
      do j = 1,nbdy
c
c calculate adjustment
         if(j.gt.6) then
            if(lbdy(j).le.0) go to 200
            m = lbdy(j) + 6
            if(nplnt.eq.3) then
               neme = earth
               if(m.eq.30) neme = name
            else if(nplnt.gt.0 .and. nplnt.le.30) then
               ict66 = Ict(66)
               if(mod(ict66,2).ne.1 .and. m.eq.13) then
                  dom = Twopi/conbdy(13)
                  call ADJAST(dom,fctbdy(m))
                  Nwv = Twopi/Nwv
                  Adj = Nwv - conbdy(13)
                  Sig = Sig*conbdy(13)/dom
                  go to 50
               endif
            endif
         else
            if(lbdy(j).le.0) go to 100
            m = j
         endif
         call ADJAST(conbdy(m),fctbdy(m))
         if(nplnt.le.30 .and. nplnt.gt.0 .and. nplnt.ne.3 .and.
     .      nplnt.ne.10) then
            if(m.eq.12 .or. (m.ge.14 .and. m.le.17))
     .         Nwv = MOD(Nwv,360._10)
         endif
c
c printout adjustment to initial conditions
   50    if(j.gt.6) then
c
c printout adjustment to parameters
            write(Iout,60) astrik(Ntype),N,wrdl1,j,wrdl2,lbdy(j),
     .                      wrds(m),neme,nplnt,conbdy(m),Adj,Nwv,
     .                      Sig,Fract
   60       format(1x,a1,i4,'. ',a4,i2,a4,' =',i2,1x,a8,
     .             ' OF ',a8,i3,8x,1pd22.15,1pd16.8,1pd22.15,
     .             1pd10.3,0pf8.3)
            if(Jout.gt.0) write(Jout,60) astrik(Ntype),N,wrdl1,
     .                        j,wrdl2,lbdy(j),wrds(m),neme,
     .                        nplnt,conbdy(m),Adj,Nwv,Sig,Fract
         else
            write(Iout,80) astrik(Ntype),N,wrdl1,j,wrdl2,lbdy(j),
     .                      wrds(j),name,nplnt,jd0,conbdy(j),Adj,
     .                      Nwv,Sig,Fract
   80       format(1x,a1,i4,'. ',a4,i2,a4,' =',i2,1x,a8,
     .             ' OF ',a8,i3,i8,1pd22.15,1pd16.8,1pd22.15,
     .             1pd10.3,0pf8.3)
            if(Jout.gt.0) write(Jout,80) astrik(Ntype),N,wrdl1,
     .                        j,wrdl2,lbdy(j),wrds(j),name,nplnt,
     .                        jd0,conbdy(j),Adj,Nwv,Sig,Fract
c
c save adjustment to satellite semi-major axis in kilometers
            if(nplnt.gt.0 .and. ncentr.gt.0 .and. j.le.1) then
               Iboda = Iboda+1
               cbod4(1,Iboda) = astrik(Ntype)
               Jboda(1,Iboda) = N
               cbod4(2,Iboda) = wrdl1
               Jboda(2,Iboda) = j
               cbod4(3,Iboda) = wrdl2
               Jboda(3,Iboda) = lbdy(j)
               wrdssd    = wrds(j)
               wrdsss(2) = wrdss2
               cbod8(1,Iboda) = wrdssd
               cbod8(2,Iboda) = name
               Jboda(4,Iboda) = nplnt
               Jboda(5,Iboda) = jd0
               Boda(1,Iboda)  = conbdy(j)*Aukm
               Boda(2,Iboda)  = Adj*Aukm
               Boda(3,Iboda)  = Nwv*Aukm
               Boda(4,Iboda)  = Sig*Aukm
               Boda(5,Iboda)  = Fract
            endif
         endif
c
c store new value into old value
         if(Keepit) conbdy(m) = Nwv
  100 end do
c
c fix adjusted eccentricity and angles
  200 if(nplnt.gt.0) then
         if(nplnt.gt.0) then
            if(icnd.lt.0 .or. icnd.gt.1) return
         endif
         if(lbdy(2).gt.0 .and. lbdy(5).gt.0 .and. lbdy(6).gt.0 .and.
     .          conbdy(2).lt.0._10) then
            conbdy(2) = -conbdy(2)
            conbdy(5) = conbdy(5) + 180._10
            conbdy(6) = conbdy(6) + 180._10
         endif
         maxang=6
         if(conbdy(1).lt.0._10) maxang=5
         do i=4,maxang
            if(lbdy(i).gt.0) then
               int = conbdy(i)/360._10
               if(conbdy(i).lt.0._10) int = int - 1
               conbdy(i) = conbdy(i) - int*360
            endif
         end do
 
c rotation initial conditions could be fixed here
      endif
 
      return
      end
