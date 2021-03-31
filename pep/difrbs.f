      subroutine DIFRBS(mrbias, mplnt, ndif, cdf)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i
 
c*** end of declarations inserted by spag
 
 
c subr. difrbs - j.f.chandler - 1980 may
c compare input radar biases with saved (see difnom)
c
c arguments
      integer*2 mplnt
      integer*4 mrbias,ndif(2)
      logical*4 cdf(5)

c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'anctrl.inc'
      include 'rbiasm.inc'
      include 'rdbias.inc'
      include 'restor.inc'
      include 'wrkcomrs.inc'
 
      character*6 rbsn
c
c set up name text
      i = mplnt
      call EBCDIX(i, rbsn, 5, 2)
 
      do while( mrbias.lt.Mumrbs )
         if(mplnt.ne.Mplrbs(mrbias+1)) return
         mrbias = mrbias + 1
 
         do i = 1, Numrbs
            if(Nplrbs(i).eq.mplnt .and. Rdbsit(1,i).eq.Rdbst1(1,mrbias)
     .       .and. Rdbsit(2,i).eq.Rdbst1(2,mrbias)
     .       .and. Rdbser(i).eq.Rdbsr1(mrbias)) then
               rbsn(1:4)= Rdbser(i)
               Nsav     = Lrbs1(i) - 1
               call DIFPM4(Rbias(1,i), Rbias1(1,mrbias),
     .                     Lrbs(1,i), 2, ndif, cdf, Rdbsit(1,i), rbsn)
            endif
            end do
         end do
 
      return
      end
