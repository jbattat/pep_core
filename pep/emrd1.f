      subroutine EMRD1(lice)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, itp
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   march 1967    subroutine emrd1
c first five records of earth-moon barycenter tape are read
c
      integer*2 lice
c lice =0 printout of data on first two records of embary tape
c lice =1 no such printout
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'comdateq.inc'
      include 'empcnd.inc'
      include 'funcon.inc'
      include 'namtim.inc'
      include 'param.inc'
      include 'pemctl.inc'
      include 'tapdta.inc'
      include 'trpcom.inc'
      include 'yvectrd1.inc'
 
      integer*2 np3/3/, klem/-3/
c
c read embary constants from disk
c test if coordinates available from n-body tape and,
c if so, set up controls as if from individual tape
      call XXRDBD(np3,itp,klem,Jdem,Intemx,Idirem,Kem,Emintx,
     . Nkiem,Kiem,0,Ecom)
c
c read first two records of earth-moon barycenter peripheral
c data set
      if(Jtest.eq.0) call XXRD1(lice,np3,itp,klem,Jdem1,Jdem2,
     . Iparem,i_mxplprt+1,Intemx,Idirem,Kem,Emintx,Frem1,Nkiem,Kiem)
 
      Emint = Emintx
      if(Idirem.lt.0) Emint = -Emintx
      Emint5 = Emint*5._10
      T0svem = Jdxx0
c
c count partials and reset jd0 for next iteration
      call XXRDCK(Lparem,Kiem,klem,1)
 
c see if any i.c. partials needed
      if(Jdem0.ne.0) then
c
c initial condition partials to be gotten from elliptic
c orbit partials
         if(Kiem(1).eq.0) then
            Kiem(1) = -2
            i = 1
            do while( .true. )
               call IMITL(Gauss,Mass(3),Econd,1,i)
               if(i.gt.1) goto 100
               i = 2
            end do
         endif
      endif
c
c read first three data records of earth-moon barycenter
c peripheral data set
  100 if(Jtest.le.0) call EMRED1(0)
      return
      end
