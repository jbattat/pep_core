      subroutine MNRD1(lice)
 
      implicit none
c
c m.e.ash   march 1967    subroutine mnrd1
c first five records of moon tape are read
c
      integer*2 lice
c lice =0 printout of data on first two records of moon tape
c lice =1 no such printout

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'comdateq.inc'
      include 'comnut.inc'
      include 'pemctl.inc'
      include 'tapdta.inc'
      include 'trpcom.inc'
      include 'yvectrd1.inc'

c local
      integer*4 intmx,itp
      integer*2 np10/10/, klmn/-2/
      real*10 frm1(2)
c
c read moon constants from disk
c test if coordinates available from n-body tape and,
c if so, set up controls as if from individual tape
      call XXRDBD(np10,itp,klmn,Jdm,intmx,Idirmn,Km,Mnintx,Nkimn,
     . Kimn,0,Mcom)
c
c read first two records of moon peripheral data set
      if(Jtest.eq.0) call XXRD1(lice,np10,itp,klmn,Jdm1,Jdm2,
     . Iparm,i_mxplprt+1,intmx,Idirmn,Km,Mnintx,frm1,Nkimn,Kimn)
 
      Mnint = Mnintx
      if(Idirmn.lt.0) Mnint = -Mnintx
      T0svmn = Jdxx0
c
c count partials and reset jd0 for next iteration
      call XXRDCK(Lparm,Kimn,klmn,1)
c
c read first three data records of moon peripheral data set
      if(Jtest.le.0) call MNRED1(0)
      Inter = Idirmn
      Dnut  = Mnint
      return
      end
