      subroutine ROTSET(cons, a, d, n, j, i0, psi0, phi0, p, mu, pdot)
 
      implicit none
 
      real*10 cons(30), a, d, n, j, i0, psi0, phi0, p, mu, pdot
c
c        r. goldstein   july 1977
c        fills the cons array as follows
c           cons(1)=alpha0 in radians
c           cons(2)=sin(alpha0)
c           cons(3)=cos(alpha0)
c           cons(4)=delta0 in radians
c           cons(5)=sin(delta0)
c           cons(6)=cos(delta0)
c           cons(7)=n in radians
c           cons(8)=sin(n)
c           cons(9)=cos(n)
c           cons(10)=j in radians
c           cons(11)=sin(j)
c           cons(12)=cos(j)
c           cons(13)=i0 in radians
c           cons(14)=sin(i0)
c           cons(15)=cos(i0)
c           cons(16)=psi0 in radians
c           cons(17)=sin(psi0)
c           cons(18)=cos(psi0)
c           cons(19)=phi0 in radians
c           cons(20)=sin(phi0)
c           cons(21)=cos(phi0)
c           cons(22)=alpha0-n
c           cons(23)=sin(alpha0-n)
c           cons(24)=cos(alpha0-n)
c           cons(25)=omeg3
c           cons(26)=p
c           cons(27)=mu
c           cons(28)=-pdot/(2*p0)
c
      include 'funcon.inc'
c
c convd=pi/180
c
      cons(1)  = a*Convd
      cons(2)  = SIN(cons(1))
      cons(3)  = COS(cons(1))
      cons(4)  = d*Convd
      cons(5)  = SIN(cons(4))
      cons(6)  = COS(cons(4))
      cons(7)  = n*Convd
      cons(8)  = SIN(cons(7))
      cons(9)  = COS(cons(7))
      cons(10) = j*Convd
      cons(11) = SIN(cons(10))
      cons(12) = COS(cons(10))
      cons(13) = i0*Convd
      cons(14) = SIN(cons(13))
      cons(15) = COS(cons(13))
      cons(16) = psi0*Convd
      cons(17) = SIN(cons(16))
      cons(18) = COS(cons(16))
      cons(19) = phi0*Convd
      cons(20) = SIN(cons(19))
      cons(21) = COS(cons(19))
      cons(22) = a - n
      cons(23) = SIN(cons(22)*Convd)
      cons(24) = COS(cons(22)*Convd)
      cons(25) = Twopi/p
      cons(26) = p
      cons(27) = mu
      cons(28) = -pdot/(2._10*p)
 
      return
      end
