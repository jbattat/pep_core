      subroutine SSBRD1(npl,klax,lice)

      implicit none

c
c j.f.chandler   jan 2019    subroutine ssbrd1
c first five records of planet tape are read, to be used in overriding
c n-body for computing solar-system barycenter offset
c
c parameters 
      integer*2 npl,klax,lice
c npl  = planet to be initialized
c klax = index of planet in list of input planets
c lice =0 printout of data on first two records of planet tape
c lice =1 no such printout
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'pemctl.inc'
      include 'tapdta.inc'
      include 'trpcom.inc'

      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'lcntrl.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'param.inc'
      include 'stats.inc'
      include 'yvectrd1.inc'

c local variables
      real*10 frp1(2),pcomc(12)
      integer*4 i,itp,j,jdum
      integer*2 neg999/-999/,kssb(u_nmprm)

c
c read body constants from disk
c preclude taking coordinates from n-body tape
      call XXRDBD(neg999,itp,klax,jdum,jdum,jdum,kssb,
     . Bintxz(npl+3),Nkissb(npl),Kissb(1,npl),0,pcomc)
      Iplss(npl)=itp
c
c read first two records of planet peripheral data set
      call XXRD1(lice,npl,itp,klax,Jdss1(npl),Jdss2(npl),
     . Iparss(npl),i_mxplprt+1,Intss(npl),Idirz(npl+3),kssb,
     . Bintxz(npl+3),frp1,Nkissb(npl),Kissb(1,npl))

      Intss5(npl) = 5*Intss(npl)
      Bintz(npl+3)= Bintxz(npl+3)
      if(Idirz(npl+3).lt.0) Bintz(npl+3) = -Bintxz(npl+3)
      T0svss(npl) = Jdxx0

c
c count partials and reset jd0 for next iteration
      call XXRDCK(Lparss(npl),Kissb(1,npl),klax,1)
c
c read first three data records of planet data set
      call SSRED1(0,npl,itp)
      return
      end
