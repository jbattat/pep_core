      subroutine DIURN(imonth,fract,em,emdot)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      actlng, diff, em, emday, emdot, emnite, fourpi, fract,
     .          tempor, twopi
 
c*** end of declarations inserted by spag
 
 
c
c d.robertson/l.friedman   dec 1969   subroutine diurn
c to give a diurnal model for ionospheric elctron density
c the model is a rectified sine wave - peak at noon;const at nite
c *****  data statements for emday,emnite serve as the input (12 values)
c
      dimension emnite(12),emday(12),em(2),emdot(2),actlng(2)
      integer*2 imonth
      dimension fract(2)
 
      include 'sitcrd.inc'
      data twopi/6.283185308/
 
      integer   i, k
 
 
c electron density in  elec/cm**3  are put in at monthly values
      data emday/12*1.E6/
      data emnite/12*1.E5/
      diff = emday(imonth) - emnite(imonth)
      do i = 1, 2
         actlng(i) = twopi*fract(i) - Longr(i)
         fourpi    = 2._10*twopi
         if(abs(actlng(i)).gt.fourpi) then
            write(6,20) (k,actlng(k),fract(k),Longr(k),k=1,2)
   20       format(2x,'ANGULAR DIST.FROM UT MIDNIGHT INCORRECT IN DIURN'
     .       /(' SITE',I2,' OFF,FR,LONG=',3F10.3))
            return
         else
            do while( actlng(i).le.0.0 )
               actlng(i) = actlng(i) + twopi
            end do
            do while( actlng(i).gt.twopi )
               actlng(i) = actlng(i) - twopi
            end do
            if((actlng(i).gt.(0.25*twopi)) .and.
     .         (actlng(i).lt.(0.75*twopi))) then
               tempor   = actlng(i) - 0.25*twopi
               em(i)    = emnite(imonth) + diff*sin(tempor)
               emdot(i) = diff*cos(tempor)*twopi
               emdot(i) = emdot(i)/86400.
            else
               em(i)    = emnite(imonth)
               emdot(i) = 0.0
            endif
         endif
      end do
      return
      end
