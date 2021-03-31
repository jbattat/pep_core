      subroutine DIFSTR(ndif,cdf)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   j, k, n
 
c*** end of declarations inserted by spag
 
 
c subr. difstr - j.f.chandler - 1983 feb
c compare input sky corrections with saved (see difnom)
c
c arguments
      integer*4 ndif(2)
      logical*4 cdf(5)

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'anctrl.inc'
      include 'restor.inc'
      include 'skystf.inc'
      include 'skystm.inc'
      include 'wrkcomrs.inc'
 
      if(Mumstr.gt.0 .and. Numstr.gt.0) then
         do k = 1, Mumstr
 
            do j = 1, Numstr
               if(Ctlgnm(j).eq.Ctlgn1(k)) then
                  Nsav = Lsky1(j) - 1
                  n    = Nskycf(j)
                  if(Mskyc(k).gt.n) n = Mskyc(k)
                  call DIFBDY(Skycf(1,j),Skycf1(1,k),Lskycf(1,j),-n,
     .                        ndif, cdf, Ctlgn1(k),'COEFF.')
                  goto 50
               endif
            end do
 
   50    end do
      endif
 
      return
      end
