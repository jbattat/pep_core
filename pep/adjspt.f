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
      real*10 fctspt(3)
      character*4    wrds(3,3)/' RAD', 'IUS ', '  KM', ' LON', 'GITU',
     .          'DEDG', ' LAT', 'ITUD', 'E DG'/
      character*2 astrik(3)/'* ','& ','  '/
      integer   i, j
 
      fctspt(1) = Aultsc*Ltvel
      fctspt(2) = Aultsc/Convd
      fctspt(3) = fctspt(2)
 
      do while( nspot.lt.Numspt )
         if(nplnt.ge.0) then
            if(nplnt.ne.Nsplnt(nspot+1)) return
            nspot = nspot + 1
         else
            nspot = nspot + 1
            if(Nsplnt(nspot).ge.0) goto 100
         endif
         call LINCHK
         do i = 1, 3
            if(Lspcrd(i,nspot).gt.0) then
               call ADJAST(Spcord(i,nspot),fctspt(i))
               write(Iout,10) astrik(Ntype),N,i,nspot,
     .                         Lspcrd(i,nspot),(wrds(j,i),j = 1,3),
     .                         Spot(nspot),name,Spcord(i,nspot),
     .                         Adj, Nwv, Sig, Fract
   10          format(1x,a1,i4,'. LSPCD(', i1, ',', i2, ')=', i2,
     .                3A4, ' OF ', a4, ' ON ', a8, 1pd22.15, 1pd16.8,
     .                1pd22.15, 1pd10.3, 0pf8.3)
               if(Jout.gt.0) write(Jout,10) astrik(Ntype),N,i,
     .            nspot, Lspcrd(i,nspot),(wrds(j,i),j = 1,3),
     .            Spot(nspot),name,Spcord(i,nspot),Adj,Nwv,Sig,
     .            Fract
               if(Keepit) Spcord(i,nspot) = Nwv
            endif
         end do
 
         call FRSTAT(1,2,'SPOT CRD')
  100 end do
 
      return
      end
