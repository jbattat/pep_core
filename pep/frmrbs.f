      subroutine FRMRBS(rst,mrbias,mplnt,ntop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i
 
c*** end of declarations inserted by spag
 
 
c
c m.e. ash - february 1970 - subroutine frmrbs
c radar biases logic routine for forming normal
c equations from saved normal equations
c
c

c arguments
      real*10 rst(1000)
      integer*4 mrbias,ntop
      integer*2 mplnt

c array dimensions
      include 'globdefs.inc'

c common
      include 'anctrl.inc'
      include 'rbiasm.inc'
      include 'rdbias.inc'
      include 'restor.inc'
 
      do while( mrbias.lt.Mumrbs )
 
         if(mplnt.ge.0) then
 
            if(mplnt.ne.Mplrbs(mrbias+1)) return
            mrbias = mrbias + 1
         else
            mrbias = mrbias + 1
            if(Mplrbs(mrbias).ge.0) goto 100
         endif
 
         if(Numrbs.gt.0) then
            do i = 1, Numrbs
               if(mplnt.ge.0 .or. Nplrbs(i).ge.0) then
                  if(Nplrbs(i).ne.mplnt) goto 20
               endif
               if(Rdbsit(1,i).eq.Rdbst1(1,mrbias)) then
                  if(Rdbsit(2,i).eq.Rdbst1(2,mrbias)) then
                     if(Rdbser(i).eq.Rdbsr1(mrbias)) then
                        Nrst = Lrbs1(i) - 1
                        call FRMBDY(rst,Lrbs(1,i),Mrbs(1,mrbias),-2,
     .                              ntop)
                        goto 100
                     endif
                  endif
               endif
   20       end do
         endif
 
         do i = 1, 2
            if(Mrbs(i,mrbias).gt.0) Nsav = Nsav + 1
         end do
  100 end do
 
      return
      end
