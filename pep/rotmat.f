      subroutine ROTMAT(ktype,rot,s)
 
      implicit none

c f amuchastegui/m ash - april 1969 - subroutine rotmat
c
 
c arguments
      real*10 rot(3,3),s
      integer*4 ktype
c
c        s =   integration time variable (julian date+.5)
c        ktype =-1 determine rotation matrix for integrated body nplnt
c        ktype = 0 determine rotation matrix for central body ncentr
c        ktype = 1,2,3,4 determine rotation matrix for target body
c
c        rot = output rotation matrix from body fixed coordinates to
c              those refered to mean equinox & equator of the
c              reference epoch
c
c        if partials wrt target body rotation parameters should become
c        necessary, /rtpars/ can be expanded as follows:
c     common /rtpars/ pcnrot(4,3,3),ptrot(4,3,3,i_mxtrg),icrotp(4,2),
c    1 itrotp(4,2,i_mxtrg),icrot1,itrot1(i_mxtrg)
c
c array dimensions
      include 'globdefs.inc'

c commons
      include 'cnthar.inc'
      include 'funcon.inc'
      include 'intstf.inc'
      include 'mnrtlb.inc'
      include 'nutces.inc'
      include 'petuna.inc'
      include 'rtpars.inc'
      include 'sbstuf.inc'
      include 'trghar.inc'
c
c local partials array
      real*10 drot(3,3,4)
      real*10 alph,calph,cdelt,csidtm,ctut1,delt,dera,fracts,omar,
     . omeg,psi,psi0,salph,sdelt,sidtm,sidtm0,ssidtm,trot
      integer   i,icp,ido,ix,j,jcrot,ktess,ll

c external functions
      real*10 ETUTF
c
c determine if it is the integrated body, central body, or
c target body
      if(ktype.lt.0) then
c
c rotation matrix for integrated body
         ktess = 0
         jcrot = 0
         if(Nplnt.ne.3) then
            alph  = (90._10 + Con(9))*Convd
            salph = SIN(alph)
            calph = COS(alph)
            delt  = (90._10 - Con(8))*Convd
            sdelt = SIN(delt)
            cdelt = COS(delt)
            goto 300
         endif
      else if(ktype.eq.0) then
c
c rotation matrix for central body
         jcrot = Icrot
         ktess = Nctess
         if(Ncentr.ne.3) then
c
c determine rotation matrix for moon
            if(Ncentr.eq.10) goto 100
c
c determine rotation matrix for planet
            if(Alphc(2).eq.0.0_10) then
               calph = Calphc
               salph = Salphc
            else
               alph  = Alphc(1) + Alphc(2)*(s - Epochc)
               calph = COS(alph)
               salph = SIN(alph)
            endif
            if(Deltc(2).eq.0.0_10) then
               cdelt = Cdeltc
               sdelt = Sdeltc
            else
               delt  = Deltc(1) + Deltc(2)*(s - Epochc)
               cdelt = COS(delt)
               sdelt = SIN(delt)
            endif
            if(ktess.le.0) goto 300
            omeg = Omegc
            psi0 = Psic0
            trot = s - Epochc
            goto 200
         endif
      else
c
c
c this is a target body
         ktess = Nttess(ktype)
         trot  = s - Epocht(ktype)
         jcrot = 0
         if(Ntrg(ktype).ne.3) then
            if(Ntrg(ktype).eq.10) goto 100
            if(Alpht(ktype,2).eq.0.0_10) then
               calph = Calpht(ktype)
               salph = Salpht(ktype)
            else
               alph  = Alpht(ktype,1) + Alpht(ktype,2)*trot
               calph = COS(alph)
               salph = SIN(alph)
            endif
            if(Deltt(ktype,2).eq.0.0_10) then
               cdelt = Cdeltt(ktype)
               sdelt = Sdeltt(ktype)
            else
               delt  = Deltt(ktype,1) + Deltt(ktype,2)*trot
               cdelt = COS(delt)
               sdelt = SIN(delt)
            endif
            if(ktess.le.0) goto 300
            omeg = Omegt(ktype)
            psi0 = Psit0(ktype)
            goto 200
         endif
      endif
c
c determine rotation matrix for earth
      call PRCES(s - 0.5_10)
      call PRENUT
      do i = 1,3
         rot(3,i) = Nutprc(3,i)
      end do
      if(ktess.gt.0) then
         ctut1 = ETUTF(Jd,Fract)
         call SIDTIM(Jd,ctut1,sidtm0,fracts,dera)
c use approximation of projected nutation in longitude as the offset of
c true from mean sidereal time, even with new models
         sidtm = sidtm0+Pc+fracts*(Fract*8.64E4_10-ctut1)
         csidtm = COS(sidtm)
         ssidtm = SIN(sidtm)
         do i = 1,3
            rot(1,i) = csidtm*Nutprc(1,i) + ssidtm*Nutprc(2,i)
            rot(2,i) = -ssidtm*Nutprc(1,i) + csidtm*Nutprc(2,i)
 
c earth wobble is ignored for now
         end do
      endif
      return
  100 call ECLPRC(Jd,Fract,0)
      call MONLIB(0,0,0)
      do j = 1,3
         do i = 1,3
            rot(i,j) = Mrotlb(i,j)
         end do
      end do
      return
  200 psi  = MOD(psi0 + omeg*trot,Twopi)
      Cpsi = COS(psi)
      Spsi = SIN(psi)
  300 rot(3,1) = salph*sdelt
      rot(3,2) = -calph*sdelt
      rot(3,3) = cdelt
      if(jcrot.gt.0) then
         do ix = 1,jcrot
            icp = Icrotp(ix,2)
            if(icp.eq.3) then
 
c partial w/r/t con(8)= dec. of pole
               drot(3,1,ix) = -cdelt*salph*Convd
               drot(3,2,ix) = cdelt*calph*Convd
               drot(3,3,ix) = sdelt*Convd
            else if(icp.eq.4) then
 
c partial w/r/t con(9)= r.a. of pole
               drot(3,1,ix) = calph*sdelt*Convd
               drot(3,2,ix) = salph*sdelt*Convd
               drot(3,3,ix) = 0._10
            else
 
c partial w/r/t con(6)=psic0 & con(7)=rotation rate
               drot(3,1,ix) = 0._10
               drot(3,2,ix) = 0._10
               drot(3,3,ix) = 0._10
            endif
         end do
      endif
      if(ktess.gt.0) then
         rot(1,1) = calph*Cpsi - salph*Spsi*cdelt
         rot(1,2) = salph*Cpsi + calph*Spsi*cdelt
         rot(1,3) = Spsi*sdelt
         rot(2,1) = -calph*Spsi - salph*Cpsi*cdelt
         rot(2,2) = -salph*Spsi + calph*Cpsi*cdelt
         rot(2,3) = Cpsi*sdelt
         if(jcrot.le.0) return
         ido = 0
         do ix = 1,jcrot
            icp = Icrotp(ix,2)
            if(icp.eq.2) then
 
c partial w/r/t con(7)
               omar = -trot*omeg**2/Twopi
               if(ido.gt.0) then
                  omar = omar/Convd
                  do i = 1,2
                     do j = 1,3
                        drot(i,j,ix) = omar*drot(i,j,ido)
                     end do
                  end do
                  goto 350
               endif
            else if(icp.eq.3) then
 
c partial w/r/t con(8)
               drot(1,1,ix) = Convd*(-sdelt*salph*Spsi)
               drot(1,2,ix) = Convd*(sdelt*calph*Spsi)
               drot(1,3,ix) = Convd*(-cdelt*Spsi)
               drot(2,1,ix) = Convd*(-sdelt*Cpsi*salph)
               drot(2,2,ix) = Convd*(sdelt*calph*Cpsi)
               drot(2,3,ix) = Convd*(-cdelt*Cpsi)
               goto 350
            else if(icp.eq.4) then
 
c partial w/r/t con(9)
               drot(1,1,ix) = Convd*(-salph*Cpsi - calph*cdelt*Spsi)
               drot(1,2,ix) = Convd*(Cpsi*calph - cdelt*salph*Spsi)
               drot(1,3,ix) = 0._10
               drot(2,1,ix) = Convd*(salph*Spsi - calph*Cpsi*cdelt)
               drot(2,2,ix) = Convd*(-calph*Spsi - salph*cdelt*Cpsi)
               drot(2,3,ix) = 0._10
               goto 350
            else
 
c partial w/r/t con(6)
               omar = Convd
               ido  = ix
            endif
            drot(1,1,ix) = omar*(-Spsi*calph - Cpsi*cdelt*salph)
            drot(1,2,ix) = omar*(-Spsi*salph + Cpsi*cdelt*calph)
            drot(1,3,ix) = omar*(Cpsi*sdelt)
            drot(2,1,ix) = omar*(-Cpsi*calph + Spsi*cdelt*salph)
            drot(2,2,ix) = omar*(-salph*Cpsi - calph*cdelt*Spsi)
            drot(2,3,ix) = omar*(-Spsi*sdelt)
  350    end do
      endif
c
c alter partial derivatives
c alter partial derivative matrix for gravity gradient code
c in sbfn1
      if(jcrot.gt.0) then
         ll = 3*jcrot
         call PRODCT(rot,drot,Pcnrot,-3,3,ll)
      endif
 
      return
      end
