      subroutine NEWHED(jmats,jtape,key,iptr,iprnt,sumsq)
 
      implicit none

c         d. white  april 1973  subroutine newhed
c         reads header records of apriori coeff matrix in format
c         of saved normal eqn
c
c         parameters
      integer*4 jmats,jtape,key,iprnt
      integer*2 iptr(1)
      real*10 sumsq
c jmats - unit number for saved a priori information
c jtape - not used
c key   - not used
c iprnt - code to govern printout: if 0 suppress print
c sumsq - weighted sum-squared residuals encoded in a priori (return)
c
c array dimensions
      include 'globdefs.inc'
c
c         common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'mcnfrm.inc'
c
c local variables
      character*80 title
      integer*4 i
c
c read title record
      read(jmats) title
      Itrwnd(jmats) = 1

c read error statistics
      read(jmats) i,sumsq

c read iptr
      read(jmats) Mparam,(iptr(i),i = 1,Mparam)
 
      if(iprnt.le.0) return
 
c print heading for apriori
      if(Line.ge.50) then
         write(Iout,50) Iterat,Heding,Date,Npage
   50    format('1ADD APRIORI COEFF.  ITERAT=',i3,t42,18A4,1x,2A4,
     .          ' PAGE',i5)
         Npage = Npage + 1
         Line  = 1
      endif
      write(Iout,100) jmats,title
  100 format('-APRIORI DATA SET',i3,' HAS TITLE=',A80/)
      Line = Line + 4
      return
      end
