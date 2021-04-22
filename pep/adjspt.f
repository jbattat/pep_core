      subroutine ADJSPT(nspot,nplnt,name)
 
      implicit none
c
c ash/forni  october 1967  subroutine adjspt
c adjust spot coordinates on given body
c
c parameters
      integer*4 nspot
      integer*2 nplnt
      character*8 name
c           nspot= spot counter
c           nplnt= planet number
c           name = eight character planet name
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'adjstf.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'param.inc'
      include 'sptcrd.inc'
c
c internal to subroutine adjspt
      real*10 fctspt(6)
      character*12 wrdc(6)/' RADIUS   KM',' LONGITUDEDG',' LATITUDE DG',
     .                     ' VERT.V MM/Y',' WEST.V MM/Y',' NRTH.V MM/Y'/
      character*2 astrik(3)/'* ','& ','  '/
      integer   i,int,j
 
      fctspt(1) = Aultsc*Ltvel
      fctspt(2) = Aultsc/Convd
      fctspt(3) = fctspt(2)
      do i=4,6
         fctspt(i)=1._10
      end do

      do while(nspot.lt.Numspt)
         if(nplnt.ge.0) then
            if(nplnt.ne.Nsplnt(nspot+1)) return
            nspot = nspot + 1
         else
            nspot = nspot + 1
            if(Nsplnt(nspot).ge.0) goto 100
         endif
         call LINCHK
         do i = 1,6
            if(Lspcrd(i,nspot).gt.0) then
               call ADJAST(Spcord(i,nspot),fctspt(i))
               write(Iout,10) astrik(Ntype),N,i,nspot,
     .                         Lspcrd(i,nspot),wrdc(i),
     .                         Spot(nspot),name,Spcord(i,nspot),
     .                         Adj, Nwv, Sig, Fract
   10          format(1x,a1,i4,'. LSPCD(', i1, ',', i2, ')=', i2,
     .                a12, ' OF ', a4, ' ON ', a8, 1pd22.15, 1pd16.8,
     .                1pd22.15, 1pd10.3, 0pf8.3)
               if(Jout.gt.0) write(Jout,10) astrik(Ntype),N,i,
     .            nspot, Lspcrd(i,nspot),wrdc(i),
     .            Spot(nspot),name,Spcord(i,nspot),Adj,Nwv,Sig,
     .            Fract
               if(Keepit) Spcord(i,nspot) = Nwv
            endif
         end do
c fix spot coordinates if latitude goes beyond a pole,
c or radius goes negative, or longitude goes beyond +/- 360
         if(Lspcrd(2,nspot).gt.0 .and. Lspcrd(3,nspot).gt.0) then
            if(ABS(Spcord(3,nspot)).gt.90._10) then
               int=Spcord(3,nspot)/360._10
               if(Spcord(3,nspot).lt.0._10) int=int-1
               Spcord(3,nspot)=Spcord(3,nspot)-360*int
               if(Spcord(3,nspot).gt.270._10) then
                  Spcord(3,nspot)=Spcord(3,nspot)-360._10
               else if(Spcord(3,nspot).gt.90._10) then
                  Spcord(3,nspot)=180._10-Spcord(3,nspot)
                  Spcord(2,nspot)=Spcord(2,nspot)+180._10
               endif
            endif
            if(Lspcrd(1,nspot).gt.0 .and. Spcord(1,nspot).lt.0._10) then
               Spcord(2,nspot)=Spcord(2,nspot)+180._10
               Spcord(3,nspot)=-Spcord(3,nspot)
               Spcord(1,nspot)=-Spcord(1,nspot)
            endif
            if(ABS(Spcord(2,nspot)).ge.360._10)
     .       Spcord(2,nspot)=MOD(Spcord(2,nspot),360._10)
         endif
 
         call FRSTAT(1,2,'SPOT CRD')
  100 end do
 
      return
      end
