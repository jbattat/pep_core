      subroutine GRDPUN(lgrid,transf,nsize,klam,nshape,scontl,
     .                  ncard)
 
      implicit none

c subr. grdpun - j.f.chandler - 1980 november
c punch grid model l-vector and parameters

c array dimensions
      include 'globdefs.inc'
c
c parameters
      real*4    transf(u_stdsz,nsize,1000/u_stdsz),scontl(nsize,9)
      integer*2 lgrid(nsize,1000)
      integer*4 nsize,klam,nshape,ncard
c
c common
      include 'inodta.inc'
c
c local
      integer   i,j,jj,jjs,jt,k,ngdpts,npt2
      real*4    gdpts,gtmp(4)
      equivalence (gdpts,ngdpts)
 
      if(nshape.ne.2) call SUICID('ILLEGAL NSHAPE, STOP IN GRDPUN  ',8)
 
      gdpts = scontl(klam,9)
      write(Ipunch,100) (lgrid(klam,i),i = 1,ngdpts)
  100 format(' LGRID=',(t9,20(i2,',')))
      ncard = ncard + (ngdpts - 1)/20 + 1
 
c punch parameter values
      jj   = 0
      jt   = 0
      jjs  = 1
      npt2 = (ngdpts+u_stdsz-1)/u_stdsz
      do j = 1,npt2
         do k = 1,u_stdsz
            jj = jj + 1
            if(jj.gt.ngdpts) goto 200
            jt = jt + 1
            if(jt.gt.4) then
               write(Ipunch,110) jjs,(gtmp(i),i = 1,4)
  110          format(' GRID(',i3,')=',4(1pe14.7,','))
               ncard = ncard + 1
               jjs   = jj
               jt    = 1
            endif
            gtmp(jt) = transf(k,klam,j)
         end do
      end do
 
c write out last buffer
  200 write(Ipunch,110) jjs,(gtmp(i),i = 1,jt)
      ncard = ncard + 1
      return
      end
