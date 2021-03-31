      subroutine FRMSPT(rst,mspot,mplnt,ntop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i
 
c*** end of declarations inserted by spag
 
 
c
c m.e. ash - february 1970 - subroutine frmspt
c spot logic routine for forming normal equations
c from saved normal equations
c
c
c arguments
      real*10 rst(1000)
      integer*4 mspot,ntop
      integer*2 mplnt

c array dimensions
      include 'globdefs.inc'

c commons
      include 'anctrl.inc'
      include 'restor.inc'
      include 'sptcdm.inc'
      include 'sptcrd.inc'
 
      do while( mspot.lt.Mumspt )
 
         if(mplnt.ge.0) then
 
            if(mplnt.ne.Msplnt(mspot+1)) return
            mspot = mspot + 1
         else
            mspot = mspot + 1
            if(Msplnt(mspot).ge.0) goto 100
         endif
 
         if(Numspt.gt.0) then
            do i = 1, Numspt
               if(mplnt.ge.0 .or. Nsplnt(i).ge.0) then
                  if(Nsplnt(i).ne.mplnt) goto 20
               endif
               if(Spot(i).eq.Spot1(mspot)) then
                  Nrst = Lspt1(i) - 1
                  call FRMBDY(rst,Lspcrd(1,i),Mspcrd(1,mspot),-3,
     .                        ntop)
                  goto 100
               endif
   20       end do
         endif
 
         do i = 1, 3
            if(Mspcrd(i,mspot).gt.0) Nsav = Nsav + 1
         end do
  100 end do
 
      return
      end
