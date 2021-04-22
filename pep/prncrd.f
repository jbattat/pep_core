      subroutine PRNCRD(nstop)
 
      implicit none

c ash/forni  october 1967  subroutine prncrd
c printout for observing sites and spots on observed bodies

c parameters
      integer*4 nstop

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'namtim.inc'
      include 'sptcrd.inc'
      include 'stcord.inc'
c
c internal to subroutine prncrd only
      integer   i, j, k, kount, mersps, merspt, mresps, nplchk, nplnum
      character*8 ksmark,planam,sptch1,
     .         astrik/'********'/, pound/'########'/,
     .          blank/'        '/, earth/' EARTH  '/, star/'  STAR  '/
      equivalence (ksmark, sptch1)
      integer*4 nsv
c
c printout for observing sites
      if(Numsit.gt.0) then
         kount = 0
         if(Line.gt.52) call NEWPG
         nsv=0
         do j=1,Numsit
            do i=4,6
               if(Lscrd(i,j).gt.0 .or. Scord(i,j).ne.0._10) then
                  nsv=1
                  goto 5
               endif
            end do
         end do
    5    continue
         do j=1,Numsit
            ksmark = blank
            if(Kscrd(j).eq.-1) then
               kount  = kount + 1
               ksmark = astrik
            endif
            call PAGCHK(59,1,0)
            if(j.le.1 .or. Line.le.2) then
               if(nsv.eq.0) then
                  write(Iout,20)
               else
                  write(Iout,25)
               endif
               Line = Line + 4
            endif
            if(nsv.eq.0) then
               write(Iout,40) j,sitd(j),(Lscrd(i,j),i=1,3),
     .          (Scord(i,j),i=1,3),Kscrd(j),ksmark
            else
               write(Iout,45) j,sitd(j),(Lscrd(i,j),i=1,6),T0site(j),
     .          (Scord(i,j),i=1,6),Kscrd(j),ksmark
            endif
   20       format('0OBSERVING SITES WITH ADJUSTABLE COORDINATES'//
     .       8x,'SITE',3x,'LSCRD',3x,
     .       'RADIUS (KM)  LONGITUDE (DEG)  LATITUDE (DEG) KSCRD')
   40       format(i4,'.',1x,a8,3i2,f15.9,f16.10,f16.10,i3,a1)
c  1. SSSSSSSS 1 1 1 111.111111111  22.2222222222  33.3333333333  3*
   25       format('0OBSERVING SITES WITH ADJUSTABLE COORDINATES'//
     .       8x,'SITE',5x,'LSCRD',6x,'T0',9x,
     .'RADIUS (KM)  LONGITUDE (DEG)  LATITUDE (DEG) VERT(MM/Y) WEST(MM/Y
     .) NRTH(MM/Y) KSCRD')
   45       format(i4,'.',1x,a8,6i2,f10.1,f15.9,f16.10,f16.10,3f11.5,
     .       i3,a1)
         end do
         if(kount.gt.0) then
            write(Iout,60) kount
   60       format(i5,' CYLINDRICAL COORDINATES(EQUAT. RADIUS(KM),',
     .             'LONGITUDE(DEG),Z(KM)) *')
            Line = Line + 1
         endif
      else
         call PAGCHK(60,2,0)
         write(Iout,100)
  100    format(
     .'0THERE IS NO INPUT DATA FOR OBSERVING SITES WITH ADJUSTABLE COORD
     .INATES')
      endif
c
c printout for spots on observed bodies
c initialize
      nsv=0
      do j=1,Numspt
         do i=4,6
            if(Lspcrd(i,j).gt.0 .or. Spcord(i,j).ne.0._10) then
               nsv=1
               goto 105
            endif
         end do
      end do
  105 continue

      merspt = 0
      mersps = 0
      mresps = 0
      if(Numspt.gt.0) then
 
c start checking each individual spot
         nplchk = -2
         if(Line.gt.52) call NEWPG
         do j = 1,Numspt
 
c is this a star
            if(Nsplnt(j).lt.0) then
               nplnum = 9999
               planam = star
 
c is this moon,earth or sun
            else if(Nsplnt(j).eq.10) then
               k = 17
               nplnum = -1
               planam = Aplnt(k)
            else if(Nsplnt(j).eq.3) then
               planam = earth
               nplnum = -2
            else if(Nsplnt(j).eq.0) then
               k = 18
               nplnum = 0
               planam = Aplnt(k)
            else
 
c find k such that nsplnt(j)= nplnt(k)
               do k = 1,Numpln
                  if(Nsplnt(j).eq.Nplnt(k)) then
                     nplnum = k
                     planam = Aplnt(k)
                     go to 120
                  endif
                  end do
               nplnum = 9998
               do i = 1,3
                  if(Lspcrd(i,j).gt.0) then
                     mersps = mersps + 1
                     planam = astrik
                     go to 120
                  endif
                  end do
               merspt = merspt + 1
               planam = pound
            endif
 
c check order of spots
  120       sptch1 = blank
            if(nplnum.lt.nplchk) then
               sptch1 = astrik
               mresps = mresps + 1
            endif
 
c printout
            call PAGCHK(59,1,0)
            if(j.le.1 .or. Line.le.2) then
               if(nsv.eq.0) then
                  write(Iout,130)
  130             format('0OBSERVED SPOTS ON OTHER BODIES'//
     .                6x,'SPOT  PLANET  NO. LSPCRD',2x,
     .                'RADIUS (KM)  LONGITUDE (DEG) LATITUDE (DEG)')
               else
                  write(Iout,135)
  135             format('0OBSERVED SPOTS ON OTHER BODIES'//
     .             6x,'SPOT  PLANET  NO.    LSPCRD',6x,'T0',7x,
     .             'RADIUS (KM)  LONGITUDE (DEG) LATITUDE (DEG)',
     .            ' VERT(MM/Y) WEST(MM/Y) NRTH(MM/Y)')
               endif
               Line = Line + 4
            endif
            if(nsv.eq.0) then
               write(Iout,140) j,Spot(j),planam,Nsplnt(j),
     .          (Lspcrd(i,j),i = 1,3),(Spcord(i,j),i = 1,3),sptch1
  140          format(i4,'.',1x,1A4,1x,a8,i4,1x,3I2,f14.8,
     .             f16.10,f15.10,1x,a4)
            else
               write(Iout,145) j,Spot(j),planam,Nsplnt(j),
     .          (Lspcrd(i,j),i=1,6),T0spot(j),(Spcord(i,j),i=1,6),sptch1
  145          format(i4,'.',1x,1A4,1x,a8,i4,1x,6I2,f10.1,f14.8,
     .             f16.10,f15.10,3f11.5,1x,a4)
            endif
            if(Mout.gt.0 .and. planam.eq.astrik)
     .          write (Mout,160) j,Spot(j),planam,Nsplnt(j),
     .                            (Lspcrd(i,j),i = 1,3),sptch1
  160       format(i4,'.',1x,1A4,1x,a8,i4,1x,3I2,1x,a4)
 
c set for planet order check
            nplchk = nplnum
            end do
         if(merspt.ne.0) then
            call PAGCHK(60,1,0)
            write(Iout,180) merspt
  180       format(11x,'########',i4,' PLANETS WHICH ARE NOT ',
     .             'INPUT, BUT SPOTS NOT ADJUSTED, SO WARNING ONLY')
         endif
         if(mersps.ne.0) then
            call PAGCHK(60,1,0)
            write(Iout,200) mersps
  200       format(11x,'********',i3,
     .     ' PLANETS WHICH ARE NOT INPUT AND SPOTS ARE ADJUSTED, ERROR')
            if(Mout.gt.0) write(Mout,200) mersps
         endif
         if(mresps.ne.0) then
            call PAGCHK(60,1,0)
            write(Iout,220) mresps
  220       format(40x,i4,' SPOTS OUT OF ORDER, ERROR ****')
            if(Mout.gt.0) write(Mout,220) mresps
         endif
         nstop = nstop + mersps + mresps
      else
         call PAGCHK(60,2,0)
         write(Iout,250)
  250    format('0THERE IS NO INPUT DATA FOR SPOTS ON OBSERVED BODIES')
      endif
      return
      end
