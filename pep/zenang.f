      subroutine ZENANG(izctl,x,u,r,s,su,sr,zen,zendot)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 DOT, dum, dz, r, sr, zs
      integer   ic
 
c*** end of declarations inserted by spag
 
 
      integer*2 izctl
c
c        r.b. goldstein  r.w. king  april 1978
c        routine to calculate the zenith angle and/or zenith
c        angle rate of a source
c        the routine needs to be called separately for each
c        site.
c
c        izctl:
c             1 bit=1: calculate za
c             2 bit=1: calculate zar
c             4 bit=1: calculate sitnrm (otherwise it is assumed to exis
c
c        note: the zenith angle is calculated using the site
c              normal, whereas the rate is calculated using xsite/r.
c              these two are different in that the site normal is
c              calculated on the basis of an elipsoidal shape.
c
c
      real*10 x(6),u(3),s(3),su(3)
c        x(6):  coordinates of the observed source w.r.t. the observing
c        u(3):  unit vector of x
c        r   :  magnitude of x
c        s(3):  planetocentric site coordinates
c        sr  :  magnitude of s
c        su(3):  coordinates of site normal
c
      real*4    zen, zendot
c zen:  zenith angle (radians)
c zendot:  zenith angle rates (rad/sec)
c
c quantities internal to this routine
      real*10 du(3),dsu(3)
c
c
      ic = izctl
c
c get unit normal if necessary
      if(mod(ic/4,2).ne.0) call UVECTR(3,s,sr,su,dum)
c
c zenith angle
      if(mod(ic,4).ne.0) then
         zs = DOT(u,su)
         if(mod(ic,2).ne.0) zen = ACOS(-zs)
c
c
c zenith angle rate
c
         if(mod(ic/2,2).ne.0) then
c
c calculate derivative of unit vectors and derivatives of
c site vectors
            call UVECTR(4,x,r,u,du)
            call UVECTR(4,s,sr,su,dsu)
 
            dz     = DOT(su,du) + DOT(u,dsu)
            zendot = dz/SQRT(1._10 - zs**2)
         endif
      endif
 
      return
      end
