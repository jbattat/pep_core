      subroutine TRPNIT
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ip
 
c*** end of declarations inserted by spag
 
 
c subroutine trpnit - j.f.chandler - 1979 nov
c initialize interpolation routines
c
c commons
      include 'bddta.inc'
      include 'trpcom.inc'
 
      do i = 1, 10
         Nvels(i) = 9999
         Ntb1s(i) = 9999
         Nbtrp(i) = 0
      end do
c
c interval quantities
c
c embary, moon, planet set up in EMRD1, MNRD1, PLRD1
 
c satellite/probe set up in SBRD1, SCRD1
 
c rotation set up in RTRD1
 
c solar system barycenter, earth, mercury for ssbc
      do i = 8, 9
         Bint(i)  = Intb(2)
         Bintx(i) = Bint(i)
         Idirb(i) = Ibdsgn
         if(Ibdsgn.lt.0) Bint(i)=-Bint(i)
      end do
      Bint(10)=Intb(1)
      Bintx(10)=Bint(10)
      Idirb(10)=Ibdsgn
      if(Ibdsgn.lt.0) Bint(10)=-Bint(10)
 
      return
      end
