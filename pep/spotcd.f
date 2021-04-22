      subroutine SPOTCD(jda,frect,mnspt1,nvel,npl,ctutsx,atutsx,n)
 
      implicit none
c
c m.ash   oct 1971   subroutine spotcd
c calculate observed spot coordinates
c
c arguments
      real*10 frect,ctutsx,atutsx
      integer*4 jda,mnspt1,nvel,n
      integer*2 npl

c array dimensions
      include 'globdefs.inc'

c commons
      include 'comdateq.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'mnsprt.inc'
      real*10 dmrtlb(3,3)
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 sbepoch,sbstate(6)
      equivalence (Save(56),sbepoch),(Save(57),sbstate(1))
      include 'param.inc'
      include 'prtpr9.inc'

c external functions
      real*10 A1UT1

c local
      real*10 alp,csidtm,ctvcor,del,fract,psi,sidtm,sidtm0,ssidtm,
     . temp(6),utsec,ututsx
      integer*4 i,ict66,ii,ioff,iprot,itstb,j,jdm,nvp
      real*10 jd2000/2451545._10/
c
c           modifications for planetary rotations by w. decampli  8/73
c
c           set up printout flags
      itstb = 4
      ioff  = 0
      nvp   = 0
      ict66 = Ict(66)
      iprot = mod(ict66/256,2)
 
c any ctvary correction should not be applied to rotation
      if(Spcdx(1,n).ne.-5._10) then
         ctvcor = Ctvary*(jda - Prm97 - 0.5_10 + frect)**2
         call TIMINC(jda,frect,jdm,fract,-ctvcor)
      else
c but it does apply to a moving spacecraft
         jdm=jda
         fract=frect
      endif
 
c calculate rotation matrix for earth
      if(npl.eq.3) then
         if(mnspt1.le.0) then
            mnspt1 = 1
            call MNREED(jdm)
            if(jdm.le.0) return
            Kindnp = 0
            ututsx = atutsx - A1UT1(jdm,fract + ctutsx/Secday)
            call SIDTIM(jdm,ctutsx-ututsx,sidtm0,Sidvel,Dera)
            call PRCNUT(jdm,fract)
            do i = 1,3
               Rot(3,i) = Nutpr(3,i)
            end do
            sidtm0 = sidtm0 + Dgst
         end if
         utsec  = fract*Secday - ctutsx + ututsx
         sidtm  = sidtm0 + Sidvel*utsec
         csidtm = COS(sidtm)
         ssidtm = SIN(sidtm)
         do i = 1,3
            Rot(1,i) = csidtm*Nutpr(1,i) + ssidtm*Nutpr(2,i)
            Rot(2,i) = -ssidtm*Nutpr(1,i) + csidtm*Nutpr(2,i)
         end do
 
c wobble not included in earth rotation matrix
         if(nvel.ge.0) call PRODCT(Rot,Yspcd(1,n),Xspcd(1,n),-3,3,1)
         goto 200
c
      else if(npl.eq.10) then
c calculate rotation matrix for moon
         goto 600
c
      else if(npl.lt.0) then
c See if any proper motion
         Stime=(jdm-Jdps0)+fract
         if(Nplsr.gt.0 .and.(Alphc(2).ne.0._10 .or. Deltc(2).ne.0._10))
     .    then
            alp  = Alphc(1) + Alphc(2)*Stime
            Calp = COS(alp)
            Salp = SIN(alp)
            del  = Deltc(1) + Deltc(2)*Stime
            Cdel = COS(del)
            Sdel = SIN(del)
            Yspcd(1,n)=Cdel*Calp
            Yspcd(2,n)=Cdel*Salp
            Yspcd(3,n)=Sdel
            Dydphi(1,n)=-Sdel*Calp
            Dydphi(2,n)=-Sdel*Salp
            Dydphi(3,n)=Cdel
         endif

c star rotation matrix is the identity matrix
         do i = 1,3
            Xspcd(i,n) = Yspcd(i,n)
         end do
         goto 400
 
      else if(npl.eq.0) then
c
c calculate rotation matrix for sun
         goto 700
      else
 
         if(Spcdx(1,n).eq.-5._10) then
c coordinates of spot are supplied externally in the obslib
c interpolate from the supplied reference epoch
            Stime = (jdm-jd2000) + (fract-0.5_10)
            Stime = Stime*secday - 3.90e-4_10 - sbepoch
            do i=1,3
               Xspcd(i,n) = (sbstate(i) + sbstate(i+3)*Stime)/Ltvel
               Xspcd(i+3,n) = sbstate(i+3)/Ltvel
            end do
            if(Jct(13).ne.1) then
c convert J2000 coordinates into our reference frame
               do i=1,6
                  temp(i)=Xspcd(i,n)
               end do
               call PRODCT(Prec2000,temp,Xspcd(1,n),-3,3,2)
            endif
            goto 400
c
c calculate rotation matrix for planet
         else if(mod(ict66,2).eq.0 .or. npl.ne.4) then
 
            Stime = (jdm - Pcom(1)) + (fract - 0.5_10)
            if(mnspt1.le.0) then
               mnspt1 = 1
               if(Alphc(2).eq.0._10) then
                  Calp = Calphc
                  Salp = Salphc
               else
                  alp  = Alphc(1) + Alphc(2)*Stime
                  Calp = COS(alp)
                  Salp = SIN(alp)
               end if
               if(Deltc(2).eq.0._10) then
                  Cdel = Cdeltc
                  Sdel = Sdeltc
               else
                  del  = Deltc(1) + Deltc(2)*Stime
                  Cdel = COS(del)
                  Sdel = SIN(del)
               end if
               Rot(3,1) = Salp*Sdel
               Rot(3,2) = -Calp*Sdel
               Rot(3,3) = Cdel
            end if
            psi  = mod(Psic0+Omegm*Stime,Twopi)
            Cpsi = COS(psi)
            Spsi = SIN(psi)
            Rot(1,1) = Calp*Cpsi - Salp*Spsi*Cdel
            Rot(1,2) = Salp*Cpsi + Calp*Spsi*Cdel
            Rot(1,3) = Spsi*Sdel
            Rot(2,1) = -Calp*Spsi - Salp*Cpsi*Cdel
            Rot(2,2) = -Salp*Spsi + Calp*Cpsi*Cdel
            Rot(2,3) = Cpsi*Sdel
c
c*  start=700
c calculate spot position
            if(nvel.ge.0) call PRODCT(Rot,Yspcd(1,n),Xspcd(1,n),
     .          -3,3,1)
         else
 
c nvel must be passed to spotpl for generality
            call SPOTPL(jdm,fract,n,ict66,0,nvel)
         end if
         goto 200
      end if
c
  200 if(iprot.eq.1) then
         if(Line.gt.55) call OBSPAG
         write(Iout,250) ((Rot(i,j),j=1,3),i = 1,3)
  250    format(' PLANET ROTATION MATRIX'/(1p,3D22.12))
         Line = Line + 4
      end if
      goto 400
 
      entry SPOTCV(nvel,npl,n)
c
c*  start=1000
c set up printout flags
      itstb = 8
      ioff  = 3
      nvp   = 1
c
c determine partial derivatives of spot position
c with respect to local radius
      if(Lspot(1,n).gt.0) then
         do i = 1,3
            Dxspcd(i,1,n) = Xspcd(i,n)/Rspot(n)
         end do
      end if

      if(npl.ge.0) then
c
c with respect to spot longitude
         if(Lspot(2,n).gt.0) then
            do i = 1,3
               Dxspcd(i,2,n) = -Rot(1,i)*Yspcd(2,n) + Rot(2,i)
     .                           *Yspcd(1,n)
            end do
         endif
c
c with respect to spot latitude
         if(Lspot(3,n).gt.0) then
            call PRODCT(Rot,Dydphi(1,n),Dxspcd(1,3,n),-3,3,1)
         end if
c
c with respect to spot upward velocity
         if(Lspot(4,n).gt.0) then
            call PRODCT(Rot,Dydv(1,1,n),Dxspcd(1,4,n),-3,3,1)
            do i=1,3
               Dxspcd(i,4,n)=Dxspcd(i,4,n)*Dtspt(n)/Aultsc
            end do
         endif
c
c with respect to spot westward velocity
         if(Lspot(5,n).gt.0) then
            call PRODCT(Rot,Dydv(1,2,n),Dxspcd(1,5,n),-3,3,1)
            do i=1,3
               Dxspcd(i,5,n)=Dxspcd(i,5,n)*Dtspt(n)/Aultsc
            end do
         endif
c
c with respect to spot northward velocity
         if(Lspot(6,n).gt.0) then
            call PRODCT(Rot,Dydv(1,3,n),Dxspcd(1,6,n),-3,3,1)
            do i=1,3
               Dxspcd(i,6,n)=Dxspcd(i,6,n)*Dtspt(n)/Aultsc
            end do
         endif

      else
c
c star coordinates
c
c with respect to spot longitude
         if(Lspot(2,n).gt.0) then
            Dxspcd(1,2,n) = -Yspcd(2,n)
            Dxspcd(2,2,n) = Yspcd(1,n)
            Dxspcd(3,2,n) = 0._10
            do i = 4,6
               Dxspcd(i,2,n) = 0._10
            end do
         endif
c     
c with respect to spot latitude
         if(Lspot(3,n).gt.0) then
            do i = 1,3
               ii = i + 3
               Dxspcd(ii,3,n) = 0._10
               Dxspcd(i,3,n)  = Dydphi(i,n)
            end do
         endif
      endif
c
c*  start=1200
c spot velocity calculations
      if(nvel.le.0) return
      if(npl.lt.0) call SUICID(
     .               ' CANNOT CALCULATE STAR VELOCITY, STOP IN SPOTCD '
     .               ,12)
c
c calculate time derivative of earth rotation matrix
      if(npl.eq.3) then
         Omegsc = Sidvel
c
c calculate time derivative of moon rotation matrix
      else if(npl.eq.10) then
         goto 600
 
c change if real calculations are to be done
c
c calculate time derivative of planet rotation matrix
      else if(npl.eq.0) then
         goto 700
      else if(Spcdx(1,n).eq.-5._10) then
         goto 400
      else
         Omegsc = Omegm/Secday
      end if
      do i = 1,3
         dmrtlb(1,i) = Rot(2,i)*Omegsc
         dmrtlb(2,i) = -Rot(1,i)*Omegsc
         dmrtlb(3,i) = 0._10
      end do
c
c calculate spot velocity
      call PRODCT(dmrtlb,Yspcd(1,n),Xspcd(4,n),-3,3,1)
c
c*  start=1400
c determine partial derivatives of spot velocity
c with respect to local radius
      if(Lspot(1,n).gt.0) then
         do i = 4,6
            Dxspcd(i,1,n) = Xspcd(i,n)/Rspot(n)
         end do
      end if
c
c with respect to spot longitude
      if(Lspot(2,n).gt.0) then
         do i = 4,6
            j = i - 3
            Dxspcd(i,2,n) = -dmrtlb(1,j)*Yspcd(2,n) + dmrtlb(2,j)
     .                        *Yspcd(1,n)
         end do
      end if
c
c with respect to spot latitude
      if(Lspot(3,n).gt.0)
     .    call PRODCT(dmrtlb,Dydphi(1,n),Dxspcd(4,3,n),-3,3,1)
c
c*  start=2000
c print appropriate coordinates
  400 if(mod(Jct(6)/itstb,2).ne.0) then
         if(Line.gt.58) call OBSPAG
         write(Iout,450) jdm,fract,nvp,n,
     .                    (Xspcd(i+ioff,n),i = 1,3)
  450    format(' SPOTCD:    JD.F=',i7,f13.12,' NV=',i2,'  N=',
     .          i2,'  X=',1p,3D23.15)
         Line = Line + 1
      end if
 
c*  start=9900
      return
 
  600 call SUICID(' SATELLITE BASED OBSERVATIONS OF SPOT ON MOON '//
     .            'NOT PROGRAMMED, STOP IN SPOTCD',19)
 
  700 call SUICID(' OBSERVATIONS OF A SUN SPOT NOT PROGRAMMED, STOP '//
     .            'IN SPOTCD  ',15)
 
      end
