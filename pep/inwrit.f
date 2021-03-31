      subroutine INWRIT
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, ntki
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   march 1970   write input quantities onto disk
c
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'france.inc'
      include 'namtim.inc'
      include 'plndta.inc'
 
c
c write body quantities onto disk
      do j=1,u_mxpl+4
         ntki = Ndumki(j)
         write(Iplcon) (Dumcon(i,j),i = 1,12),
     .                 (Dumeps(i,j),i = 1,6),(Kkk(i,j),i = 1,100),
     .                 Jd1(j),Jd2(j),Int(j),Intp1(j),Intp2(j),
     .                 Ihrp(j),Iminp(j),Secp(j),
     .                 (Kkp(i,j),i = 1,100),Nplnt(j-4),
     .                 (Dtcon(i,j),i = 1,30),Ndumki(j),
     .                 (Kdumi(i,j),i = 1,ntki)
      end do
 
c end file iplcon
      rewind Iplcon
 
      return
      end
