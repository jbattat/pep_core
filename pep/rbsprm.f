      subroutine RBSPRM(nplnt,k,names,iskale,xnom)
 
      implicit none
c
c d. white  april 1973  subroutine rbsprm
c
c generates names for radar biases
c
c parameters
      character*8  names(2,1)
      real*10 iskale(1),xnom(1)
      integer*4 k, nplnt

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'rdbias.inc'
c
c locals
      integer   i, n, nrbias
      character*1 rb(3) /'T', 'D', 'S'/
      character*1 minus / '-' /
      character*8 temp
      character*1 tem(8)
      equivalence (tem(1),temp)
c
c loop thru observations series
c now do loop over all nrbias.  no longer need 'data nrbias/0/'
      if(Numrbs.gt.0) then
         n = 100000 + iabs(nplnt)
         call EBCDI(n,temp,6)
         if(nplnt.lt.0) tem(5) = minus
         do nrbias = 1, Numrbs
            if(nplnt.eq.Nplrbs(nrbias)) then
c
c loop thru biases
               do i = 1, 2
                  if(Lrbs(i,nrbias).ne.0) then
                     k = k + 1
                     xnom(k) = Rbias(i,nrbias)
                     call MVC(Rdbsit(1,nrbias),1,8,names(1,k),1)
                     call MVC(Rdbser(nrbias),1,4,temp,1)
                     tem(7) = rb(i)
                     tem(8) = rb(i + 1)
                     names(2,k) = temp
                  endif
               end do
            endif
         end do
      endif
      return
      end
