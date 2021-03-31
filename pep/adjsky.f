      subroutine ADJSKY
 
      implicit none

c subr. adjsky - j.f.chandler - 1983 feb
c adjust sky correction coefficients
c
c array dimensions
      include 'globdefs.inc'
c
c common
      include 'adjstf.inc'
      include 'inodta.inc'
      include 'skystf.inc'

c local 
      character*1 astrik(3)/'*','&',' '/
      real*10 fctsky/1._10/
      character*4 cssn(2)/' COS',' SIN'/,RADEC(2)/'*RA','*DEC'/
      integer   i,ics,j,k,n4,nsky
 
      if(Numstr.le.0) return
      do k = 1,Numstr
         nsky = Nskycf(k)
         if(nsky.gt.0) then
            call LINCHK
            j   = 1
            ics = 1
            do i = 1,nsky
               if(Lskycf(i,k).gt.0) then
                  call ADJAST(Skycf(i,k),fctsky)
                  n4 = (i + 3)/4
                  write(Iout,10) astrik(Ntype),N,i,k,Lskycf(i,k),
     .                            Ctlgnm(k),cssn(ics),n4,radec(j),
     .                            Skycf(i,k),Adj,Nwv,Sig,Fract
   10             format(1x,a1,i4,'. LSKY(',i2,',',i2,')=',i2,
     .                   1x,a8,a4,i3,a4,' (ARCSEC)   ',1pd22.15,
     .                   d16.8,d22.15,d10.3,0pf8.3)
               endif
               ics = 3 - ics
               if(ics.eq.1) j = 3 - j
            end do
            call FRSTAT(1,2,Ctlgnm(k))
         endif
      end do
      call FRSTAT(2,2,'STR CTLG')
 
      return
      end
