      subroutine EMIPT(nn,time,j)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   j,nn,nv
      real*10 time
 
c*** end of declarations inserted by spag
 
 
c     ash/becker  may 1968  subroutine emipt
c        code moved to jnitl/jlipt  1977 dec  j.f.chandler
c     subroutine emipt determines position and velocity in the encke
c         elliptic orbit and the partial derivatives of position and
c         velocity in this orbit with respect to the orbital elements
c     for target bodies
c             nn=0 position and velocity determined
c          nn=-1,1 position,velocity,and partial derivatives determined
c
c common
      include 'emcke.inc'
      include 'emmips.inc'
 
      nv = 1
      if(nn.ne.0) nv = 7
      call JLIPT(time,Elptg(1,j),nv,Ympt(1,j),Rympt(j),Rympt2(j),
     . Rympt3(j),Dympt(1,1,j))
      return
      end
