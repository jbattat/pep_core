      subroutine PRNSTR
 
      implicit none

c subr. prnstr - j.f.chandler - 1983 jan
c print input star catalog error models, if any

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'skystf.inc'
 
c local
      integer*4 i,i0,j,k,lex,n
      character*4 rpar/')   '/
 
      if(Numstr.gt.0) then
         if(Line.gt.51) call NEWPG
         call PAGSET('SKY CORRECTION STAR CATALOG ERROR MODELS',-10)
         call PAGHED(0)
         do j = 1,Numstr
            n = Nskycf(j)/4
            call PAGCHK(60,2,1)
            write(Iout,20) Ctlgnm(j),n
   20       format('0CATALOG ',a8,' HAS COEFFICIENTS UP TO ORDER',
     .             i3,':')
            i0 = 0
            do i = 1,n
               call PAGCHK(60,1,1)
               write(Iout,30) i,(Skycf(i0+k,j),k = 1,4)
   30          format(' SKY(',i2,')=',1p,4D22.15)
               i0 = i0 + 4
            end do
            lex = 1 + (n - 1)/10
            call PAGCHK(60,lex,1)
            write(Iout,40) (Lskycf(i,j),i = 1,i0)
   40       format(' LSKY=',(t10,10(1x,4I2)))
         end do
         return
      else
         call PAGCHK(60,2,0)
         write(Iout,50)
   50    format('0THERE ARE NO SKY CORRECTION STAR CATALOG ERROR MODELS'
     .        )
         return
      endif
      end
