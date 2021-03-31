      subroutine ADJRBS(nrbias,nplnt,name)
 
      implicit none
c
c ash/forni  october 1967  subroutine adjrbs
c adjust planetary radar observation biases
c
c arguments
      integer   nrbias
      character*8 name
      integer*2 nplnt
c           nrbias= bias counter
c           nplnt = planet number
c           name  = eight character planet name
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'adjstf.inc'
      include 'inodta.inc'
      include 'rdbias.inc'
c
c internal to subroutine adjrbs
      real*4    snwv
      integer   i
      real*10 old
      real*10 fctrbs/1._10/
      character*10 wrds(2)/'TMDLY BIAS','DOPLR BIAS'/
      character*2 astrik(3)/'* ','& ','  '/
 
      do while( nrbias.lt.Numrbs )
 
         if(nplnt.ge.0) then
 
            if(nplnt.ne.Nplrbs(nrbias+1)) return
            nrbias = nrbias + 1
         else
            nrbias = nrbias + 1
            if(Nplrbs(nrbias).ge.0) go to 100
         endif
 
         call LINCHK
         do i = 1,2
            if(Lrbs(i,nrbias).gt.0) then
               old = Rbias(i,nrbias)
               call ADJAST(old,fctrbs)
               snwv = Nwv
               write(Iout,10) astrik(Ntype),N,i,nrbias,
     .                         Lrbs(i,nrbias),wrds(i),name,
     .                         Rdbsit(1,nrbias),Rdbser(nrbias),
     .                         Rdbsit(2,nrbias),Rbias(i,nrbias),
     .                         Adj,snwv,Sig,Fract
   10          format(1x,a1,i4,'. LRBS(',i1,',',i2,') = ',i1,
     .                1x,a10,2x,a8,1x,a4,1x,a4,1x,a4,1x,
     .                1pe12.5,5x,1pd16.8,5x,1pe12.5,5x,1pd10.3,
     .                0pf8.3)
               if(Jout.gt.0) write(Jout,10) astrik(Ntype),N,i,
     .            nrbias,Lrbs(i,nrbias),wrds(i),name,
     .            Rdbsit(1,nrbias),Rdbser(nrbias),Rdbsit(2,nrbias),
     .            Rbias(i,nrbias),Adj,snwv,Sig,Fract
               if(Keepit) Rbias(i,nrbias) = snwv
            endif
         end do
 
         call FRSTAT(1,2,'OBS BIAS')
  100 end do
 
      return
      end
