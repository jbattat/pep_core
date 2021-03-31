      subroutine ELIPT(nn, time)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 Elpte, time
      integer   nn, nv
 
c*** end of declarations inserted by spag
 
 
c     m.e.ash   april 1965   subroutine elipt
c     optimized 1977 j.f.chandler
c        code moved to jnitl/jlipt  1977 dec  j.f.chandler
c     subroutine elipt determines position and velocity in the encke
c         elliptic orbit and the partial derivatives of position and
c         velocity in this orbit with respect to the orbital elements
c             nn=0 position and velocity determined
c          nn=-1,1 position,velocity,and partial derivatives determined
      include 'ellips.inc'
 
c elliptic orbit quantities (encke labeled common)
      common/ENCKE/ Elpte(31)
 
      nv = 1
      if(nn.ne.0) nv = 7
      call JLIPT(time,Elpte,nv,Ylpt,Rylpt,Rylpt2,Rylpt3,Dylpt)
      return
      end
