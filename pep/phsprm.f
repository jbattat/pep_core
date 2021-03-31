      subroutine PHSPRM(nplnt,k,names,iskale,xnom)
 
      implicit none
c
c d. white  april 1973  subroutine phsprm
c
c generates names for optical phase corrections
c name format implemented by m. ratner may 1979.
c
c parameters
      character*8 names(2,1)
      real*10 iskale(1),xnom(1)
      integer*4 k,nplnt

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'phase.inc'
c
c locals
      integer*4 i,iabsnp,j1,j10,n,nphase
      character*1 minus / '-' /
      character*8 word2 / 'PL00PHC0' /
      character*1 nums(10)/'1','2','3','4','5','6','7','8','9','0'/
      character*8 temp
      character*1 tem(8)
      equivalence (tem,temp)
c
c loop thru optical observations series
      if(Numphs.gt.0) then
         do nphase = 1,Numphs
            if(nplnt.eq.Nplphs(nphase)) then
               n = Ncphs(nphase)
               if(n.gt.0) then
c
c loop thru phases
                  do i = 1,n
                     if(Lphs(i,nphase).ne.0) then
                        k = k + 1
                        xnom(k) = Aphase(i,nphase)
c construct name from 1st site,series,nplnt, & phase coef index.
c for example, 'hayssera  pl04pc1' or 'ovro9900  pl-4pc2'.
                        call MVC(Phsit(nphase),1,4,names(1,k),1)
                        call MVC(Phser(nphase),1,4,names(1,k),5)
                        temp   = word2
                        iabsnp = iabs(nplnt)
                        j1     = mod(iabsnp,10)
                        if(j1.eq.0) j1 = 10
                        j10 = iabsnp/10
                        if(j10.eq.0) j10 = 10
                        tem(3) = nums(j10)
                        tem(4) = nums(j1)
                        if(nplnt.lt.0) tem(3) = minus
                        tem(8) = nums(i)
                        names(2,k) = temp
                     endif
                  end do
               endif
            endif
         end do
      endif
      return
      end
