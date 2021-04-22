      subroutine PARTL(kick)
 
      implicit none

c     m.e.ash    july 1969     subroutine partl
c     calculate radar,optical,transit,interferometer partial derivatives
c     made to utilize saved partial derivatives, oct 1967

c arguments
      integer*4 kick
c           kick =1 subroutine radar is calling program
c           kick =2 subroutine optic is calling program
c           kick =3 subroutine trnsit is calling program
c           kick =4 subroutine fermtr is calling program
c           partl calculates partial derivatives of observations with
c           respect to parameters to be adjusted

c array dimensions
      include 'globdefs.inc'

c commons
      include 'bddta.inc'
      include 'cmcke.inc'
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 dutrec
      equivalence (Angdum(10),dutrec)
      include 'emmips.inc'
      include 'eqnphs.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'lfix.inc'
      include 'ltrap.inc'
      integer*2 kind
      equivalence (kind,Numpar)
      include 'mnsprt.inc'
      include 'mtrap.inc'
c
c           the following vectors control calculation of the following
c                                 partial derivatives
c           lprm(100),mprm(100)   solar system parameters in /param/
c           lem(30),mem(30)       earth-moon barycenter initial
c                                 conditions and earth parameters from
c                                 /empcnd/
c           lpl(30),mpl(30)       planet initial conditions and
c                                 parameters from /empcnd/ for given
c                                 value of klan
c           the m-vectors control which partial derivatives exist in
c           real*10 vired(296,2) vector
c           the l-vectors control which partial derivatives will be put
c           into real*8 deriv(296,2) vector.
c           every partial in vired exists in deriv.
c           every partial in deriv which exists in vired (except for
c           planet shape parameters) is transferred from vired to deriv.
c           every partial in deriv which is not transferred from vired
c           is calculated.
c           kind is counter for deriv, mind is counter for vired.
c           vired(mind,i) = partial derivitive number mind of
c                           measurement i
c           deriv(kind,i) = partial derivative number kind of
c                           measurement i
c           i=1    time delay or right ascension
c           i=2    doppler shift or declination
      include 'namtim.inc'
      include 'number.inc'
      include 'obscrd.inc'
      real*10 ctrecf,freqtr(2),xfactr
      equivalence (Dstf,ctrecf),(Dstf(7),freqtr)
      equivalence (Save(33),xfactr)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'partcm.inc'
      integer*4 kmrt
      equivalence (Kprt,kmrt)
      include 'pemctl.inc'
      include 'radcrd.inc'
      include 'rotdtamr.inc'
      include 'rtrdvl.inc'
      include 'sbcom.inc'
      include 'sbdta.inc'
      include 'scdta.inc'
      include 'sitcrd.inc'
c unifying site partials
c order: rad lon lat up wes nor
      real*10 sitpar(6,2,6)
      equivalence (sitpar,Pstrad)
      include 'stats.inc'
      include 'tapdta.inc'
      include 'tidal.inc'
      include 'trnocc.inc'
      include 'yvect.inc'

c external functions
      real*10 DOT

c local variables
      real*10 sitfct(6)/3*1._10,3*8.64E4_10/
      integer*2 neg3/-3/, neg2/-2/
      integer*2 lprml
      integer*2 lemic
      real*10 yss(3),dnmra,dnmdec,dnmd2,aut,ctcor,t1,t2,
     .          tau10,tau20,tau30,term,timlpt,unrec,unsend,uvrec,
     .          uvsend,x1,xffact
      integer   i,i1,iflag,j,jetest,k,kpl,l,lembry,
     .          lmon,lmrt,lplnt,lplss(9),lrad,lsprb,lsprc,m,m1
 
      Mnstf = Mnau
 
c nk1 =-1 for first observation of series set in cmpar2
      if(Nk1.lt.0) call ZFILL(Lold,2*296)
 
c set der** status indicators   (could be generalized)
      do i = 1, 8
         Ivus(i) = 0
      end do
      if(kick.ne.3 .or. Ibtrn.ne.9) Ivus(1) = 1
      if(Nplnt0.eq.10 .or. Ncp0.eq.10) Ivus(2) = 1
      if(Nplnt0.eq.10 .and. Nspot.gt.0) Ivus(6) = 1
      if(Klan.gt.0 .and. Klan.le.u_mxpl) Ivus(3) = 1
 
c if(klanb.gt.0) ivus(4)=1
      if(Nps1.gt.0) Ivus(5) = 1
 
c initially, all is garbage
      do i = 1, 8
         Ivze(i) = -1
      end do
      if(Klan.ne.u_mxpl+1 .or. kick.eq.3) Mnstf = Mnau*Mass(10)
      kind  = 2
      Mind  = 2
      Mouse = 1
      Index = 3
      if(kick.eq.2) then
         if(Jct(39).lt.0) Index = 6
      else if(kick.eq.3) then
         if(Ibtrn.gt.6) Index = 6
      else
         if(Nice.ge.0) Index = 6
         if(Ncodf.ne.19 .and. nddiff.ne.-1) Mouse = 2
         Ivze(1)=0
         do k = 1, 2
            do j = 1, 2
               do i = 1, Index
                  Derpr(i,j,k) = 0._10
               end do
            end do
         end do
         if(Nice.ge.0) then
            Dopfct(1) = Dt3dt2/Omv2
            Dopfct(2) = 1._10/Opv3
 
c maybe unit vector derivatives should already be done
            do j = 1, Mouse
               call UVECTR(4,Xsitep(1,j),Rsitp(j),Xsitp0(1,j),
     .                     Dxsit0(1,j))
            end do
         endif
      endif
c
c copy partials quickly
      if(Ict(4).lt.-2 .or. (Ict(4).eq.-2.and.Lnotm.ne.1)) then
         Numpar = Munpar
         do i = 3, Munpar
            Deriv(i,1) = Vired(i,1)
            Deriv(i,2) = Vired(i,Mun2)
            Lold(i)     = 1
         end do
 
c copy all 6 dt pointers
         call MVC(Imdt1,1,12,Ildt1,1)
         if(Nk1.lt.0) Nk1 = 0
         return
      endif
c
c
c-----------------------------------------------------------------------
c*  start=1000
c
c           partial derivatives w.r.t. observing site coordinates
c
      m = 0
 
  100 m = m + 1
c m=1 receiving site
c m=2 sending site
      l = 0
      call PCOPS(m,'SCRD', Iabs1)
 
  110 iflag = 0
      call PCOPY(l,6,iflag,1,Lscrd(1,m),Mscrd(1,m))
      if(iflag.le.0) then
c
c*  start=1100
         do m1 = 1, Mouse
            if(m1.ne.m .and. Nsite1.ne.Nsite2) then
               do i = 1, Index
                  derem(i,m1) = 0._10
               end do
            else
               do i = 1, Index
                  derem(i,m1) = sitpar(i,m1,l)*sitfct(i)
               end do
 
c partial of diurnal term of ct-at
               ctcor = 0._10
               if(Ncodf.eq.18)
     .          ctcor = -DOT(sitpar(1,m1,l),Xemlsc(4,1))*Aultsc*Dpphdt
            endif
         end do
         call PARTVL(derem,2,kick)
c
c*  start=1400
         if(Ncodf.gt.20) then
            do j = 1, 2
               do i = 1, Index
                  Derpr(i,j,2) = derem(i,j)
               end do
            end do
         endif
         Lprspt  = 0
         Ivze(1) = 1
         call CPARTC(kick)
         if(Nice.le.0) Deriv(kind,1) = Deriv(kind,1) + ctcor
      endif
      if(l.lt.6) goto 110
      if(m.lt.Mouse .and. Nsite1.ne.Nsite2) goto 100
c
c-----------------------------------------------------------------------
c
c*  start=2000
c partial derivatives w.r.t. equinox-equator-latitude
c corrections for optical observation series
      if(kick.eq.3) then
c
c error message for transit observations
         if(Neqnox.le.0) goto 10
         call SUICID(
     .' EQ-EQ-DECL BIASES NOT ALLOWED FOR TRANSIT OBS, STOP PARTL  ',15)
      else if(kick.eq.2) then

c setup for pure rotation biases
         if(Sitf(1)(1:4).eq.'@REF') then
            dnmra=Xsitp0(1,2)**2+Xsitp0(2,2)**2
            if(dnmra.ne.0._10) then
               dnmdec=1._10/SQRT(dnmra)
               dnmra=1._10/dnmra/15._10
            else
               dnmdec=0._10
            endif
            dnmd2=Xsitp0(3,2)*dnmdec/DOT(Xsitp0(1,2),Xsitp0(1,2))
            call CROSS(Xsitp0(1,1),Xsitp0(1,2),yss)
            do i=1,3
               yss(i)=yss(i)*dnmd2
            end do
            yss(1)=yss(1)-Xsitp0(2,1)*dnmdec
            yss(2)=yss(2)+Xsitp0(1,1)*dnmdec
         endif
 
c we must increment counters even if neqnox=0
         do i = 1, 3
            if(Leqn(i).gt.0) then
               kind = kind + 1
               if(Iabs1.gt.0) then
                  if(Meqn(i).gt.0) Mind = Mind + 1
               endif
               if(Sitf(1)(1:4).eq.'@REF') then

c partials for pure rotation biases
                  Deriv(kind,1)=-Xsitp0(3,1)*Xsitp0(i,2)*dnmra
                  Deriv(kind,2)=yss(i)
                  if(i.eq.3) Deriv(kind,1)=(Xsitp0(1,1)*Xsitp0(1,2)+
     .             Xsitp0(2,1)*Xsitp0(2,2))+dnmra
               else

c partials for equiniox/equator/latitude biases
                  if(i.eq.1) then
                     Deriv(kind,1) = 1._10
                     Deriv(kind,2) = 0._10
                     if(Ncodf.gt.20) then
c for differential obs. the origin is arbitrary
c use equinox bias for plate scale (1+denox)*nominal
                        Deriv(kind,1) = Angdum(8)
                        Deriv(kind,2) = Angdum(9)
                     endif
                  else if(i.eq.2) then
                     Deriv(kind,1) = -Calph*Tdelt/15._10
                     Deriv(kind,2) = Salph
                  else
                     Deriv(kind,1) = 0._10
                     Deriv(kind,2) = 1._10
                  endif
               endif
            endif
         end do
         goto 10
      endif
c
c*  start=2100
c  partial derivatives w.r.t.radio or interferometer clock
c  corrections
c         jct(59)=0  eqnx are terms added to observable
c         jct(59)=1  eqnx are true clock terms
      do l = 1, 3
         if(Leqn(l).gt.0) then
            kind = kind + 1
            if(Iabs1.gt.0) then
               if(Meqn(l).gt.0) Mind = Mind + 1
            endif
            if(l.eq.1) then
               Deriv(kind,1) = 1._10
               Deriv(kind,2) = 0._10
            else if(l.eq.2) then
               Deriv(kind,1) = dutrec
               Deriv(kind,2) = -1._10
               if(kick.eq.4) Deriv(kind,2) = 1._10
            else
               Deriv(kind,1) = dutrec**2/2._10
               Deriv(kind,2) = -dutrec
               if(kick.eq.4) Deriv(kind,2) = dutrec
            endif
            if(Ncodf.eq.18) then
               Deriv(kind,1) = Deriv(kind,1)*Dpphdt
               Deriv(kind,2) = Deriv(kind,2)*Dpphdt
            endif
            if(kick.eq.4 .and. Jct(59).ne.0) then
               if(Nice.ge.0) call SUICID(
     .' CANNOT CALCULATE CLOCK PARAMETER PARTIAL FOR RATE OBSERVABLE, ST
     .OP IN PARTL', 19)
               do i = 1, 3
                  if(Nplnt0.gt.0) derem(i,1)
     .                = Xsitep(i + 3,1)*Deriv(kind,1)/Aultsc
                  if(Nplnt0.lt.0) derem(i,1)
     .                = Xsite(i + 3,1)*Deriv(kind,1)/Aultsc
                  derem(i,2) = 0._10
                  if(nddiff.ge.1) then
                     if(Nplnt2.gt.0) Derpr(i,1,2)
     .                   = Ysitep(i + 3,1)*Deriv(kind,1)/Aultsc
                     if(Nplnt2.lt.0) Derpr(i,1,2) = derem(i,1)
                  endif
               end do
               Ivze(1) = 1
               call CPARTC(kick)
               xffact = 1._10
               if(nintrf.ne.-1) xffact = (freqtr(2)-freqtr(1))*xfactr
               if(l.eq.1) then
                  x1 = 1._10
                  if(nintrf.eq.0) x1 = 0._10
                  Deriv(kind,1) = Deriv(kind,1) + xffact*x1
               else if(l.eq.2) then
                  Deriv(kind,1) = Deriv(kind,1) + xffact*dutrec
               else
                  Deriv(kind,1) = Deriv(kind,1) + xffact*dutrec**2/2._10
               endif
            endif
         endif
      end do
c
c*   start=2300
c partials with respect to star catalog error models
   10 call SKYPAR(kick)
c
c-----------------------------------------------------------------------
c
c*  start=2500
c  partial derivatives with respect to solar system parameters
c
      lmon   = 8
      Kmon   = Lparm
      lembry = 8
      Kembry = Lparem
      lplnt  = 8
      Kplnt  = Lparp
      lsprb  = 8
      Ksprb  = Lparsb
      lsprc  = 8
      Ksprc  = Lparsc
      lmrt   = 8
      Kmrt   = Lparmr
      do i=1,9
         lplss(i)=8
         Kplss(i)=Lparss(i)
      end do
      m      = 1
      l      = 0
      call PCOPS(m,'PRM ', Iabs1)
  200 iflag = 1
      call PCOPY(l,100,iflag,1,Lprm,Mprm)
      if(iflag.gt.0) goto 1200
      lprml = Lprm(l)
c
c see if observed body is moon or probe in cislunar space
      if(lprml.le.50) then
         if(kick.ne.3 .or. Ibtrn.eq.1 .or. Ibtrn.eq.9) then
            if(Nplnt0.eq.10 .or. Ncp0.eq.3 .or. Ncp0.eq.10) goto 300
         endif
c
c partial derivatives of earth-moon barycenter position,vel.
         iflag = -1
         call PBDPRM(Nkiem,Kiem,lembry,Kembry,lprml,iflag)
         if(iflag.gt.0) goto 300
         call CPARTL(1,1,kick)
         goto 400
c
c*  start=3300
c
c           partials with respect to parameters 51-100 which do not
c           affect motion of bodies
c
      else if(lprml.eq.51) then
c partial derivatives w.r.t. astronomical unit in
c light-seconds multiplied by the current value of the
c astronomical unit in light-seconds
         if(kick.eq.2 .or. kick.eq.3) goto 1080
         do i = 1, Mouse
            do j = 1, Index
               derem(j,i) = 0._10
               if(Nplnt0.eq.10 .or. Ncp0.eq.10) then
                  if(Mdstsc.le.0._10) derem(j,i) = -Xm(j,i)
               else if(Ncp0.ne.3) then
                  derem(j,i) = Xem(j,i)
                  if(Mdstsc.le.0._10) derem(j,i) = derem(j,i)
     .                - Mass(10)*Xm(j,i)
                  if(Nswcns.eq.1 .and. kick.eq.1) derem(j,i)
     .                = derem(j,i) - Xslcns(j,i)
                  if(Klan.gt.0) derem(j,i) = derem(j,i) - Xp(j)
               endif
               if(Klanb.gt.0) derem(j,i) = derem(j,i) - Xsb(j)*Cmfct
            end do
         end do
         Ivze(1) = 1
         goto 1000
c
c suppress RELDEL, RELDOP, LTVARY, and interplanetary plasma partials
c for observations of the Moon
      else if(lprml.lt.62 .and. Nplnt0.eq.10) then
         goto 1080
c
      else if(lprml.eq.53) then
c
c*  start=3400
         if(kick.eq.2 .or. kick.eq.3) goto 1080
         Deriv(kind,1) = Raddum(1)*(1._10 + Gamapm)*0.5_10
         Deriv(kind,2) = Raddum(6)*Freq*(1._10 + Gamapm)*0.5_10
      else if(lprml.ge.60 .and. lprml.le.63) then
         if(kick.eq.2 .or. kick.eq.3) goto 1080
         lrad = lprml-58
         Deriv(kind,1) = Raddum(lrad)
         Deriv(kind,2) = Raddum(lrad + 5)*Freq
c
c*  start=3500
c partial w.r.t. ctvary
      else if(lprml.eq.72) then
         if(kick.ne.1) call SUICID(
     .'CTVARY PARTIAL ONLY IMPLEMENTED FOR RADAR LINK, STOP PARTL  ',
     . 15)
         Deriv(kind,1) = 0._10
         if(Nice.le.0) then
            t1     = DOT(Xsitp0(1,2),Xsitep(4,2))
            t2     = DOT(Xsitp0(1,1),Xsitep(4,1))
            tau30  = Jd - Prm97 - 0.5_10 + ctrecf
            tau20  = tau30 - Dstf(4)/Secday
            tau10  = tau30 - Tmdly1/Secday
            unrec  = -DOT(Xsitp0,Xsite(4,1))
            unsend = DOT(Xsitp0(1,2),Xsite(4,2))
            uvrec  = -DOT(Xsitp0,Xspcd(4,1))
            uvsend = DOT(Xsitp0(1,2),Xspcd(4,1))
 
c doppler ctvary partial not implemented
            Deriv(kind,1) = (tau30**2*(t1/Opv3+(t2+unrec)*Dt3dt2/Omv2)
     .                     + tau20**2*(uvsend/Opv3-uvrec*Dt3dt2/Omv2)
     .                     - tau10**2*unsend/Opv3)*Secday
     .                     - Tmdly1*(tau10 + tau30)
         endif
         Deriv(kind,2) = 0._10
c
c partial w.r.t. ct-at scale factor
      else if(lprml.eq.81) then
         if(kick.ne.1 .or. Ncodf.eq.3) goto 1080
         if(Ncodf.eq.18) then
c pulsar timing
c (ought to exclude site vel. from dot(xessbc...) here, since
c earth rotation is in ut, but error is only 1E-6_10)
            if(Nice.le.0) Deriv(kind,1) = Dpphdt*Ctatv(1)
     .          *(DOT(Xessbc(4,1),Xsitp0) - 1._10)
         else
 
c radar/radio ranging
            if(Nice.le.0) Deriv(kind,1) = Ctatv(2)
     .          - Ctatv(1)*(Omv1 + DOT(Xsitp0(1,1),Xsite(4,1)))
     .         /Omv2*Opv2/(Opv3 - DOT(Xsitp0(1,2),Xsite(4,2)))
         endif
         Deriv(kind,2) = 0._10
c
c*  start=3600
c partials w.r.t sun radius for transit observation
      else if(lprml.eq.95) then
         if(kick.ne.3 .or. Ibtrn.ne.4) goto 1080
         if(Nice.le.0) Deriv(kind,1)= (Re(1)-Rmxc(1))/Dfdt3(1)/Aukm
         if(Nice.ge.0) Deriv(kind,2)= (Re(2)-Rmxc(2))/Dfdt3(2)/Aukm
      else
         call SUICID(
     .' PARTIAL DERIVATIVES W.R.T. PARAMETERS 52,54-59,64-71,73-100 CANN
     .OT BE CALCULATED, STOP IN PARTL', 24)
      endif
      goto 1100
c
c zero derem if no earth partials because earth is logical center
  300 do i = 1, Mouse
         do j = 1, Index
            derem(j,i) = 0._10
         end do
      end do
      Ivze(1) = 0
      if(Ncp0.eq.3 .and. Klanb.gt.0) goto 500
c
c*  start=2600
c increment by partial derivatives of moon position,velocity
  400 iflag = -1
      call PBDPRM(Nkimn,Kimn,lmon,Kmon,lprml,iflag)
      if(iflag.le.0) then
 
         call CPARTL(2,3,kick)
 
c partial not found on moon tape, decide if result is zero
      else if(Klan.eq.u_mxpl+1) then
         if(Ivze(1).gt.0) Ivze(1) = -1
         goto 1080
      endif
c
c increment by partial derivatives of moon libration
      if(Ivus(6).eq.1) then
         iflag = -1
         call PBDPRM(Nkimr,Kimr,lmrt,kmrt,lprml,iflag)
         if(iflag.le.0) call CPARTL(6,3,kick)
      endif
c
c increment by partial derivatives of planet position,velocity
      if(Klan.gt.0 .and. Klan.le.u_mxpl) then
         iflag = -1
         call PBDPRM(Nkipl,Kipl,lplnt,Kplnt,lprml,iflag)
         if(iflag.le.0) then
            call CPARTL(3,1,kick)
            goto 600
         endif
      endif
 
c planet partial not found, assume zero
  500 do j = 1, Index
         Derpl(j) = 0._10
      end do
      Ivze(3) = 0
  600 if(Klan.gt.0 .and. Klan.le.u_mxpl) then
         if(lprml.gt.10 .and. lprml.le.30) then
 
c check for satellites of this planet
            do i = 1, Numpcm
               kpl = Kplcm(i)
               if(lprml.eq.Nplnt(kpl)) then
 
c correct for offset of planet from c of m
                  do j = 1, Index
                     Derpl(j) = Derpl(j) - Xpcm(j,i)
                  end do
                  if(lprml.eq.Nplnt0 .and. Cmfct.ne.1._10) then
                     do j=1,Index
                        Derpl(j)=Derpl(j)-Xsb(j)
                     end do
                  endif
                  Ivze(3) = 1
                  goto 700
               endif
            end do
         endif
      endif
c
c*  start=2800
c get partials for secondary body
  700 if(Nps1.le.0) goto 900
      if(Klans1.le.0) then
         if(Klan.le.0 .or. Klan.gt.u_mxpl) goto 900
         if(Nps1.eq.Nplnt(Klan)) goto 800
      endif
      iflag = -1
      call PBDPRM(Nkisc,Kisc,lsprc,Ksprc,lprml,iflag)
      if(iflag.le.0) call CPARTL(5,1,kick)
      if(Ncs1.le.0) goto 900
      if(Klan.le.0 .or. Klan.gt.u_mxpl) goto 900
      if(Ncs1.ne.Nplnt(Klan)) goto 900
 
c increment partials by central body
      if(Ivze(3).le.0) goto 900
      if(Ivze(5).eq.1) then
         do i = 1, Mouse
            do j = 1, Index
               Dersc(j,i) = Dersc(j,i) + Derpl(j)
            end do
         end do
         goto 900
      endif
 
c dersc contained garbage
  800 do j = 1, Index
         Dersc(j,1) = Derpl(j)
      end do
      Ivze(5) = Ivze(3)
c
c*  start=2900
c increment by partial derivitives of probe pos.and vel.
  900 if(Klanb.gt.0) then
         iflag = -1
         call PBDPRM(Nkisb,Kisb,lsprb,Ksprb,lprml,iflag)
         if(iflag.le.0) then
            call CPARTL(4,1,kick)
            do j = 1, Index
               Derpl(j) = Derpl(j) + Dersb(j)*Cmfct
            end do
            Ivze(3) = 1
         endif
      endif
 
      if(kick.ne.3) then
 
c for anything but transits, get partial relative to earth
         if(Ivze(3).gt.0) then
            do i = 1, Mouse
               do j = 1, Index
                  derem(j,i) = derem(j,i) - Derpl(j)
               end do
            end do
         endif
      endif
c     ivze(1)=1
c
c*  start=3000
c for all solar-system parameters that affect motion
c include indirect effect on barycenter if one-way delay from a star
      if((kick.eq.1 .or. kick.eq.4) .and.
     . Nplnt0.eq.-4 .and. Nswcns.eq.1 .and. Ncodf.le.20) then
c indirect effect through the dependence of planet coordinates
c we can ignore the effects on the terrestrial planets
         do i=5,8
            if(Ssbkl(i).gt.0) then
               iflag = -1
               call PBDPRM(Nkissb(i),Kissb(1,i),lplss(i),
     .          Kplss(i),lprml,iflag)
               if(iflag.le.0) then
                  call CPARTL(10+i,3,kick)
                  Ivze(1)=1
               endif
            endif
         end do
      endif

      if(lprml.lt.10) then
c
c partials w.r.t masses 1-9
c include direct effect on barycenter if one-way delay from a star
         if((kick.eq.1 .or. kick.eq.4) .and.
     .    Nplnt0.eq.-4 .and. Nswcns.eq.1 .and. Ncodf.le.20) then
c direct dependence of barycenter offset on masses
            do j = 1, 3
               derem(j,1) = derem(j,1) + 
     .          (Xslcns(j,1)-Xpcm(j,lprml))/Mascnt
            end do
            Ivze(1) = 1
         endif
c further increment if partial derivatives are with respect
c to mass(10) = (mass of moon)/(mass of earth+moon)
      else if(lprml.eq.10) then
         if(Nplnt0.ne.10) then
            if(Klanb.gt.0) then
               if(Ncp0.eq.3 .or. Ncp0.eq.10) goto 1000
            endif
            do i = 1, Mouse
               do j = 1, Index
                  derem(j,i) = derem(j,i) - Xm(j,i)*Mnau
               end do
            end do
            Ivze(1) = 1
         endif
c
      else if(lprml.le.30) then
c partials with respect to masses 11-30
c include asteroid effect on barycenter if one-way delay from a star
         if((kick.eq.1 .or. kick.eq.4) .and.
     .    Nplnt0.eq.-4 .and. Nswcns.eq.1 .and. Ncodf.le.20) then
            do i=1,Numpcm
               if(Kplcm(i).gt.0 .and. Nplnt(Kplcm(i)).eq.lprml) then
                  do j = 1, 3
                     derem(j,1) = derem(j,1) + 
     .                (Xslcns(j,1)-Xpcm(j,i))/Mascnt
                  end do
                  Ivze(1) = 1
               endif
            end do
         endif
 
      else if(lprml.gt.30) then
c*  start=3100
c partials with respect to parameters 31-50 (which affect
c motions of bodies)
      endif
c
 1000 if(Ncodf.gt.20 .and. (Nplnt2.eq.Nplnt0 .or. Nplnt2.eq.Ncp0)) then
         do j = 1, Mouse
            do i = 1, Index
               Derpr(i,j,2) = Derpr(i,j,1)
            end do
         end do
      endif
c
c*  start=3200
c calculate partial of observation
      call CPARTC(kick)
      if(lprml.gt.30 .and. lprml.le.50) then
 
c do scale adjustment for parameters 31-50
         if(Nice.le.0) Deriv(kind,1) = Deriv(kind,1)/Hippo(lprml-30)
         if(Nice.ge.0) Deriv(kind,2) = Deriv(kind,2)/Hippo(lprml-30)
      endif
      if(kick.le.1) then
         if(lprml.eq.42) then
            if(Nice.le.0) Deriv(kind,1) = Deriv(kind,1)
     .          + 0.5_10*Raddum(1)*Reldel
            if(Nice.ge.0) Deriv(kind,2) = Deriv(kind,2)
     .          + 0.5_10*Raddum(6)*Freq*Reldel
         else if(lprml.eq.44) then
            term = 0.5_10*Reldel*Hippo(12)
            if(Nice.le.0) Deriv(kind,1) = Deriv(kind,1)
     .          + Raddum(1)*term
            if(Nice.ge.0) Deriv(kind,2) = Deriv(kind,2)
     .          + Raddum(6)*Freq*term
         endif
      endif
      goto 1100
 
c zero partials and loop
 1080 Deriv(kind,1) = 0._10
      Deriv(kind,2) = 0._10
 
c*  start=3900
 1100 if(l.lt.100) goto 200
c
c-----------------------------------------------------------------------
c*  start=4000
c
c           partial derivatives w.r.t. earth-moon barycenter initial
c           conditions and earth parameters
c
 1200 do i = 1, Index
         Derpl(i) = 0._10
      end do
      Ivze(3) = 0
      jetest  = 0
      Kembry  = 1
      lembry  = 2
      Kmon    = Lparm
      lmon    = 8
      kmrt    = Lparmr
      lmrt    = 8
      Kplnt   = Lparp
      lplnt   = 8
      m = 7
      l = 0
      call PCOPS(m,'EM  ', Iabs1)
 1300 do while( .true. )
c
c*  start=4100
c search for partial w.r.t. embary init.cond. from embary tape
         iflag = 0
         call PCOPY(l,6,iflag,1,Lem,Mem)
         if(iflag.gt.0) then
c
c search for partial w.r.t. earth parameter from embary tape
            if(l.le.6) then
               Kembry = Lparem
               lembry = 8
            endif
            iflag = 1
            call PCOPY(l,u_nmbod,iflag,1,Lem,Mem)
            if(iflag.gt.0) goto 1800
            iflag = 1
            call PBDPRM(Nkiem,Kiem,lembry,Kembry,Lem(l),iflag)
            if(iflag.le.0) goto 1500
c
c*  start=4300
c           calculate partial w.r.t. earth parameter not on embary tape
c
c           love-number scale factors and time lag for earth tides
            if(Jct(10).gt.0) then
               if(kick.ne.1 .and. kick.ne.4) goto 1600
               if(Lem(l).eq.3) then
                  do j = 1, Mouse
                     do i = 1, Index
                        Derpr(i,j,1) = Dxdhe(i,j)
                     end do
                  end do
               else if(Lem(l).ne.4) then
                  if(Lem(l).ne.5) goto 1320
                  do j = 1, Mouse
                     do i = 1, Index
                        Derpr(i,j,1) = Dxdte(i,j)
                     end do
                  end do
               else
                  do j = 1, Mouse
                     do i = 1, Index
                        Derpr(i,j,1) = Dxdle(i,j)
                     end do
                  end do
               endif
c tide calculation was performed in lt.sec and lt.sec/sec
c must convert to au and au/day to match integrated partials
               do j = 1, Mouse
                  do i = 1, Index
                     if(i.le.3) then
                        Derpr(i,j,1) = Derpr(i,j,1)/Aultsc
                     else
                        Derpr(i,j,1) = Derpr(i,j,1)/Aultvl
                     endif
                  end do
               end do
               call PARTVL(derem,2,kick)
               if(Ncodf.gt.20) then
                  do j = 1, Mouse
                     do i = 1, Index
                        Derpr(i,j,2) = Derpr(i,j,1)
                     end do
                  end do
               endif
               Ivze(1) = 1
               call CPARTC(kick)
               goto 1700
            endif
c
c*  start=4700
c relativity parameter
 1320       if(Ncp0.eq.3 .or. Ncp0.eq.10) goto 1600
            call PRTREL(Lem(l),Nkiem,Kiem,Lparem,1,kick,0)
            goto 1700
         else
            if(kick.ne.3 .or. Nps1.le.0) then
               if(Klan.eq.u_mxpl+1) then
c
c initial condition partial interpolated from moon tape
                  if(Lemmn(l).le.0) goto 1600
                  lemic = 300 + l
                  iflag = -1
                  call PBDPRM(Nkimn,Kimn,lmon,Kmon,lemic,iflag)
                  if(iflag.gt.0) goto 1600
                  if(Ivus(6).eq.0) then
                     call CPARTL(2,2,kick)
                  else
                     call CPARTL(2,3,kick)
c increment by partial derivatives of moon libration
                     iflag = -1
                     call PBDPRM(Nkimr,Kimr,lmrt,kmrt,lemic,iflag)
                     if(iflag.le.0) call CPARTL(6,3,kick)
                     call CPARTC(kick)
                  endif
                  goto 1400
               else if(Ncp0.eq.3) then
                  goto 1600
               endif
            endif
            if(Kiem(1).ge.-1) then
c
c initial condition partial interpolated from embary tape
               call PBDIC(Kiem,lembry,Kembry,l,'EM', neg3, 0)
               goto 1500
            else
 
c initial condition partial to be gotten from elliptic orbit formula
               if(jetest.le.0) then
                  jetest = 1
                  timlpt = Jd - Jdem0
                  timlpt = timlpt + ctrecf
                  call EMIPT(1,timlpt,1)
                  if(Mouse.gt.1) then
                     timlpt = timlpt - Dstf(4)/Secdy2
                     call EMIPT(1,timlpt,2)
                  endif
               endif
               do i = 1, Mouse
                  do j = 1, Index
                     derem(j,i) = Dympt(j,l,i)
                  end do
               end do
               call PARTVL(derem,2,kick)
               do i = 1, Mouse
                  do j = 1, Index
                     Derpr(j,i,2) = derem(j,i)
                  end do
               end do
               Ivze(1) = 1
               call CPARTC(kick)
            endif
         endif
 1400 end do
c
c*  start=4200
c calculate partial of coordinates
 1500 call CPARTL(1,1,kick)
 
      lemic = 300 + l
c
c search for cross partial w.r.t. embary i.c. or parameter from planet tape
      if(Jct(14).gt.0) then
         if(Klan.gt.0 .and. Klan.le.u_mxpl) then
            iflag = -1
            call PBDPRM(Nkipl,Kipl,lplnt,Kplnt,lemic,iflag)
            if(iflag.le.0) call CPARTL(3,3,kick)
         endif
      endif
c
c add contribution from integrated moon partial, if any
      iflag = -1
      call PBDPRM(Nkimn,Kimn,lmon,Kmon,lemic,iflag)
      if(iflag.le.0) then
         call CPARTL(2,3,kick)
         if(Ivus(6).gt.0) then
c add contribution from moon libration partials, if any
            iflag = -1
            call PBDPRM(Nkimr,Kimr,lmrt,kmrt,lemic,iflag)
            if(iflag.le.0) call CPARTL(6,3,kick)
         endif
      endif
c
c calculate partial of observation
      call CPARTC(kick)
      goto 1700
 
 1600 Deriv(kind,1) = 0._10
      Deriv(kind,2) = 0._10
c
c*  start=4900
 1700 if(l.lt.u_nmbod) goto 1300
c
c-----------------------------------------------------------------------
c
c           partial derivatives w.r.t. earth rotation initial conditions
c           and parameters
c
 1800 call ERTPAR(kick)
c
c-----------------------------------------------------------------------
c
c           partial derivatives w.r.t. earth gravitational potential
c           harmonic coefficients for observations of probe with
c           earth as central body
c
      if(Klanb.gt.0 .and. Ncp0.eq.3) then
         Ksprb  = Lparsb
 
c zonal harmonics
         call HPARTL(kick,Nczone,Lczhar,Mczhar,Ncp0,31,8,Iabs1)
 
c tesseral cosine harmonics
         call HPARTM(Nctess,Lcchar,Mcchar,41)
 
c tesseral sine harmonics
         call HPARTM(Nctess,Lcshar,Mcshar,51)
      endif
c
c-----------------------------------------------------------------------
c*  start=5000
c
c           partial derivatives w.r.t. moon initial conditions and
c           and parameters
      do i = 1, Index
         Derpl(i) = 0._10
      end do
      Ivze(3) = 0
      Kembry  = Lparem
      lembry  = 8
      Kmon    = 1
      lmon    = 2
      kmrt    = Lparmr
      lmrt    = 8
      Kplnt   = Lparp
      lplnt   = 8
      m = 7
      l = 0
      call PCOPS(m,'MN  ', Iabs1)
c
 1900 iflag = 0
      call PCOPY(l,6,iflag,1,Lmn,Mmn)
      if(iflag.le.0) then
         lprml=1000+l
         if(Klan.ne.u_mxpl+1 .and. Jct(14).gt.0) then
c
c search for partial w.r.t. moon initial condition from embary tape
            iflag = -1
            call PBDPRM(Nkiem,Kiem,lembry,Kembry,lprml,iflag)
            if(iflag.le.0) call CPARTL(1,1,kick)
c
c search for partial w.r.t. moon initial condition from planet tape
            if(Klan.gt.0 .and. Klan.le.u_mxpl) then
               iflag = -1
               call PBDPRM(Nkipl,Kipl,lplnt,Kplnt,lprml,iflag)
               if(iflag.le.0) call CPARTL(3,3,kick)
            endif
         endif
c
c search for partial w.r.t. moon initial condition from moon tape
         call PBDIC(Kimn,lmon,Kmon,l,'MN', neg2, 0)
         goto 1950
      endif
c
c partials w.r.t. moon parameters
      if(l.le.6) then
         Kmon = Lparm
         lmon = 8
      endif
      iflag = 1
      call PCOPY(l,u_nmbod,iflag,1,Lmn,Mmn)
      if(iflag.gt.0) goto 2100
c search for partial w.r.t. moon parameter from moon tape
      iflag = 1
      call PBDPRM(Nkimn,Kimn,lmon,Kmon,Lmn(l),iflag)
      if(iflag.le.0) then
c*  start=5200
c calculate partial of observation
         lprml=1006+Lmn(l)
         goto 1950
      else
c
c*  start=5300
c calculate partial w.r.t. moon parameter not on moon tape
         if(kick.eq.2 .or. kick.eq.3 .or. kick.eq.4) then
            if(Lmn(l).le.1) goto 1980
         else if(Lmn(l).eq.1) then
 
c partials w.r.t. moon radius (km) times a.u. (km)
            if(Nplnt0.ne.10) goto 1980
            Deriv(kind,1) = -2._10*Aultsc
            Deriv(kind,2) = 0._10
            goto 2000
         endif
c
c ad hoc precession about eliptic
         if(Lmn(l).eq.6) then
            if(kick.eq.3) goto 1980
            i1 = 1
            do j = 1, Mouse
               if(j.eq.i1) then
                  call PRODCT(Dgdprc,Xm(1,i1),Dermn(1,i1),3,3,
     .                        Index/2)
                  aut = Mnstf*Tgdp(i1)
                  do i = 1, 3
                     if(Index.gt.3) Dermn(i + 3,i1)
     .                   = aut*Dermn(i + 3,i1) + Mnstf*Dermn(i,i1)
                     Dermn(i,i1) = aut*Dermn(i,i1)
                  end do
               endif
               do i = 1, Index
                  Derpr(i,j,1) = -Dermn(i,i1)
               end do
               if(Klan.ne.u_mxpl+1) i1 = i1 + 1
            end do
            Ivze(1) = 1
            Ivze(2) = 1
 
         else
c love-number scale factors and time lag for lunar tides
            if(Jct(11).le.0) goto 1970
            if(kick.ne.1 .and. kick.ne.4) goto 1980
            if(Nplnt0.ne.10 .or. Nspot.le.0) goto 1980
            if(Lmn(l).eq.3) then
               do j = 1, Mouse
                  do i = 1, Index
                     Derpr(i,j,1) = -Dxdhm(i,j)
                  end do
               end do
             else if(Lmn(l).eq.4) then
               do j = 1, Mouse
                  do i = 1, Index
                     Derpr(i,j,1) = -Dxdlm(i,j)
                  end do
               end do
            else if(Lmn(l).eq.5) then
               do j = 1, Mouse
                  do i = 1, Index
                     Derpr(i,j,1) = -Dxdtm(i,j)
                  end do
               end do
            else
               goto 1970
            endif
c tide calculation was performed in lt.sec and lt.sec/sec
c must convert to au and au/day to match integrated partials
            do j = 1, Mouse
               do i = 1, Index
                  if(i.le.3) then
                     Derpr(i,j,1) = Derpr(i,j,1)/Aultsc
                  else
                     Derpr(i,j,1) = Derpr(i,j,1)/Aultvl
                  endif
               end do
            end do
         endif
         do j = 1, Mouse
            call PARTVL(Derpr(1,j,1),1,kick)
         end do
         if(Ncodf.gt.20) then
            do j = 1, Mouse
               do i = 1, Index
                  Derpr(i,j,2) = Derpr(i,j,1)
               end do
            end do
         endif
         Ivze(1) = 1
         call CPARTC(kick)
      endif
      goto 2000

c increment by moon orbit partials
 1950 call CPARTL(2,3,kick)
      if(Ivus(6).gt.0) then
c increment by partial derivatives of moon libration for spot observed
         iflag = -1
         call PBDPRM(Nkimr,Kimr,lmrt,kmrt,lprml,iflag)
         if(iflag.le.0) call CPARTL(6,3,kick)
      endif
      call CPARTC(kick)
      goto 2000
      
c do relativity factor or print error message
 1970 call PRTREL(Lmn(l),Nkimn,Kimn,Lparm,2,kick,0)
      goto 2000
 
c zero partials and loop
 1980 Deriv(kind,1) = 0._10
      Deriv(kind,2) = 0._10
c
c*  start=5900
 2000 if(l.lt.u_nmbod) goto 1900
c
c
c-----------------------------------------------------------------------
c
c           partial derivatives w.r.t. moon rotation initial conditions
c           and parameters
c
 2100 call MRTPAR(kick)
c
c-----------------------------------------------------------------------
c
c           partial derivatives w.r.t. moon gravitational potential
c           harmonic coefficients for observations of probe with
c           moon as central body
c
      if(Klanb.gt.0) then
         if(Ncp0.eq.10) then
            Ksprb  = Lparsb
 
c zonal harmonics
            call HPARTL(kick,Nczone,Lczhar,Mczhar,Ncp0,31,8,Iabs1)
 
c tesseral cosine harmonics
            call HPARTM(Nctess,Lcchar,Mcchar,41)
 
c tesseral sine harmonics
            call HPARTM(Nctess,Lcshar,Mcshar,51)
         endif
      endif
c
c-----------------------------------------------------------------------
c*  start=9000
c
c           continue calculating partial derivatives
c*  start=9990
      call PARTL1(kick)
      if(Nk1.lt.0) Nk1 = 0
      return
      end
