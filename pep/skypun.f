      subroutine SKYPUN(ncard)
 
      implicit none

c           subr. skypun - j.f.chandler - 1983 feb
c           punch adjusted sky corrections (star catalog errors)
c           update 'ncard': count of cards punched
      integer*4 ncard

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'inodta.inc'
      include 'skystf.inc'

c local
      integer*4 i,j,k,n,ns

      if(Numstr.le.0) return
      do j = 1,Numstr
         ns = Nskycf(j)
         if(ns.gt.0) then
            do k = 1,ns
               if(Lskycf(k,j).gt.0) then
                  n = (ns + 3)/4
                  write(Ipunch,10) Ctlgnm(j),n,
     .                              (Lskycf(i,j),i = 1,ns)
   10             format('*SKYCORR'/' CTLG=''',a8,''',N=',i4,
     .                   '  LSKY='/8(1x,4(i1,',')))
                  ncard = ncard + 3 + (ns - 1)/32
                  write(Ipunch,20) (Skycf(i,j),i = 1,ns)
   20             format(' SKY=',(t6,3(1pd21.14,',')))
                  ncard = ncard + 1 + (ns - 1)/3
                  goto 100
               endif
            end do
         endif
  100 end do
 
      return
      end
