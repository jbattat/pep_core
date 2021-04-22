      subroutine SPTPRM(nplnt,k,names,iskale,xnom)
 
      implicit none
c
c d. white  april 1973  subroutine sptprm
c
c generates names for spot coordinates
c
c parameters
      character*8    names(2,1)
      real*10 iskale(1),xnom(1)
      integer*4 k, nplnt

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'sptcrd.inc'
c
c external
      real*10 SKALE
c
c locals
      integer   i,nspot
      character*4 blanks(4)/4*'    '/
      character*8 coord(6)/'RAD','LONG','LAT','UP','WEST','NORTH'/
c
c loop thru spots
c now do loop over all nspot.  no longer need 'data nspot/0/'
      do nspot = 1, Numspt
         if(nplnt.eq.Nsplnt(nspot)) then
c
c loop thru coordinates
            do i = 1,6
               if(Lspcrd(i,nspot).ne.0) then
                  k = k + 1
                  call MVC(Spot(nspot),1,4,names(1,k),1)
                  call MVC(blanks,1,4,names(1,k),5)
                  names(2,k) = coord(i)
                  if(i.le.3) then
                     iskale(k)=SKALE(71 + i/2)
                  else
                     iskale(k)=1._10
                  endif 
                  xnom(k)= Spcord(i,nspot)/iskale(k)
               endif
            end do
         endif
      end do
      return
      end
