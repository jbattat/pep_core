      subroutine REFRCT(rfc,ist)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 f1, fsq, rfc, rstp, tlon
      integer   i, ist, jh, n
 
c*** end of declarations inserted by spag
 
 
c       subroutine refrct  - j.f.chandler  -  1979 february
c       compute the refraction constant (n-1) at send or receive site
c       setup must be performed via a call to refrc1 after eshape
c
c       ist=1  receive site
c       ist=2  send site
c       rfc    value of refraction constant returned
c
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obscrd.inc'
      include 'sitcrd.inc'
c
c parameters of earth climate
c t=t1+t2*cos(l)+t3*cos(l)*diurn+t4*l*annu+t5*h
      real*4    tcn(5)/249.7257,51.27402,7.9365,11.40036,-9.759547/
 
c dependence of refraction on f, p, t
      real*4    rcn(6)/77.4935,4.85766,0.463989,4.0041E-6,-1.1778E-8,
     .          1.0735E4/
 
c local variables
      real*4    tf(2),ta(2),td(2),th(2),ts,t,p
      integer*2 ihrs(2)
      real*4    antrm(12)/-1.,-.8660254,-.5,0.,.5,.8660254,1.,
     .          .8660254, .5, 0., -.5, -.8660254/
 
c default frequency
      real*10 f0/6E14_10/
 
      jh = Ihr - ihrs(ist)
      if(jh.lt.1) jh = jh + 24
      jh  = (jh + 1)/2
      ts  = tf(ist) + ta(ist)*antrm(Imonth) + td(ist)*antrm(jh)
      t   = ts + th(ist)
      p   = 1013.*(t/ts)**3.5
      rfc = (rstp*p*(1.+p*(rcn(4)+rcn(5)*t)) + rcn(6))/t*1E-6_10
      if(Ict(24).eq.0) return
      if(Line.gt.57) call OBSPAG
      write(Iout,100) ist,t,p,rfc
  100 format(5x,'SITE', i2, ': T=', f10.5, ', P=', f10.4, ', N=1+',
     .       1pd13.6)
      Line = Line + 1
      return
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c setup
      entry REFRC1(n)
c n=1  do just receive site
c n=2  also do send site
c
      f1 = Freq
      if(Freq.le.0._10) f1 = f0
      fsq  = (f1*1E-15_10)**2
      rstp = rcn(1) + fsq*(rcn(2) + fsq*rcn(3))
      do i = 1, n
         tf(i) = tcn(1) + tcn(2)*Cnrm(i)
         th(i) = tcn(5)*Shgt(i)
         ta(i) = tcn(4)*ATAN2(Snrm(i),Cnrm(i))
         td(i) = tcn(3)*Cnrm(i)
         tlon  = Coords(2,i)/15._10
         if(tlon.lt.0._10) tlon = tlon + 24._10
         ihrs(i) = tlon + .5_10
      end do
      return
      end
