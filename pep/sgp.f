      subroutine SGP
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 a, allng, atem, axnl, aynl, cos2u, cosi, cosn,
     .          coso, d, dm, dnode, domega, elng, EXANM,
     .          omegal
      real*10 p1, pl, sin2u, sini, sinn, sino, sis, snodes, t, tlng, 
     .          ts, ttscsi
      integer   i
 
c*** end of declarations inserted by spag
 
 
c
c w.taylor   may 1974   subroutine sgp
c determine osculating cartesian coordinates from mean elements
c
c common
      include 'ebloc.inc'
      include 'funcon.inc'

c local
      real*10 er2km/6378.160_10/,
     . ed2ks/7.382130E-2_10/,
     . mu/11467.874226_10/,
     . aj2/.108248E-2_10/,
     . aj3/-2.562E-6_10/
c
c compute semi-major axis from mean motion
      a  = (mu/(Tno*Tno))**.33333333333333333333_10
      d  = -1.5_10*aj2*(1.0_10/a)**2/SQRT(1.0_10 - Eo*Eo)
     .     **3*(1.0_10 - 1.5_10*SIN(Aio)**2)
      Ao = a*(1.0_10 + .33333333333333333333_10*d*(1.0_10-d))
c
c update elements to new epoch
      p1     = 1.5*aj2*(1./(Ao*(1.-Eo*Eo)))**2*Tno
      Anodot = -p1*COS(Aio)
      Omgadt = p1*(2. - 2.5*SIN(Aio)**2)
      t      = Time - Tepoch
      domega = Omgadt*t
      dnode  = Anodot*t
      dm     = Tno*t
      Alm    = MOD(Al + dm + domega + dnode + Twopi,Twopi)
      Omegam = MOD(Omegao + domega + Twopi,Twopi)
      Anodem = MOD(Anodeo + dnode + Twopi,Twopi)
 
      cosi = COS(Aio)
      sini = SIN(Aio)
      cosn = COS(Anodem)
      sinn = SIN(Anodem)
      coso = COS(Omegam)
      sino = SIN(Omegam)
 
c compute and apply long periodic terms.
      tlng   = (aj3*sini)/(aj2*Ao*(1.-Eo*Eo))
      axnl   = Eo*coso
      aynl   = Eo*sino - .5*tlng
      elng   = SQRT(axnl*axnl + aynl*aynl)
      omegal = ATAN2(aynl,axnl)
      atem   = Alm - .25*tlng*axnl*(3. + 5.*cosi)/(1. + cosi)
      allng  = MOD(atem + Twopi,Twopi)
 
      Ave    = MOD(Twopi + EXANM(allng-omegal-Anodem,elng),Twopi)
      Rmag   = Ao*(1. - elng*COS(Ave))
      pl     = Ao*(1. - elng*elng)
      Rvdt   = SQRT(mu*pl)/Rmag
      Rmgdt  = SQRT(mu*Ao)*elng*SIN(Ave)/Rmag
      Atrueu = 2.*ATAN(SQRT((1.+elng)/(1.-elng))*TAN(.5*Ave))
     .         + omegal
c
c compute and apply short periodic terms.
c
      ts     = .25*aj2/(pl*pl)
      sin2u  = SIN(2.*Atrueu)
      cos2u  = COS(2.*Atrueu)
      Rmag   = Rmag + ts*cos2u*pl*sini*sini
      atem   = Atrueu - .5*ts*(6. - 7.*sini*sini)*sin2u
      Atrueu = MOD(atem + Twopi,Twopi)
      ttscsi = 3.*ts*cosi
 
      sis    = Aio + ttscsi*sini*cos2u
      snodes = Anodem + ttscsi*sin2u
      sinn   = SIN(snodes)
      cosn   = COS(snodes)
      sini   = SIN(sis)
      cosi   = COS(sis)
 
      Atruev = MOD(Twopi + Atrueu - Omegam,Twopi)
      Amm    = MOD(Twopi + Alm - Omegam - Anodem,Twopi)
      if(Amm.lt.0.0) Amm = Amm + Twopi
 
      U(1) = COS(Atrueu)*cosn - SIN(Atrueu)*sinn*cosi
      U(2) = COS(Atrueu)*sinn + SIN(Atrueu)*cosn*cosi
      U(3) = SIN(Atrueu)*sini
 
      V(1) = -SIN(Atrueu)*cosn - COS(Atrueu)*sinn*cosi
      V(2) = -SIN(Atrueu)*sinn + COS(Atrueu)*cosn*cosi
      V(3) = COS(Atrueu)*sini
 
      W(1) = sini*sinn
      W(2) = -sini*cosn
      W(3) = cosi
 
      do i = 1, 3
         R(i) = Rmag*U(i)
 
         Rdot(i) = Rmgdt*U(i) + Rvdt*V(i)
      end do
c
c convert from earth radii and days to kilometers and seconds
      do i = 1, 3
         R(i)    = R(i)*er2km
         Rdot(i) = Rdot(i)*ed2ks
      end do
 
      return
      end
