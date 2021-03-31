      subroutine HPARTL(kicki,length,lhar,mhar,nqlnti,ntype,nstrti,mti)
 
      implicit none
 
 
c
c m.e.ash   june 1969    subroutine hpartl
c calculate partial derivative w.r.t. gravitational potential
c harmonic coefficients for artificial space probe observations
c or observations of the moon
c

c parameters
      integer kicki,ntype,nstrti,mti
      integer*2 length,lhar(100),mhar(100),nqlnti
c kicki - calling routine: 1=radar, 2=optic, 3=trnsit, 4=interf
c length- size of L vector
c lhar  - input L vector
c mhar  - L vector from old observation library
c nqlnti- object whose integration is to be searched for partials
c ntype - code for type of harmonics: 31=zonal, 41=tesseral cosine,
c         51=tesseral sine, 46+56 not implemented
c nstrti- starting index into ki vector
c mti   - flag indicating presence of old obslib if >0

c array dimensions
      include 'globdefs.inc'

c commons
      include 'partcm.inc'
      include 'pemctl.inc'
      include 'plndta.inc'
      include 'rotdta.inc'
      include 'sbdta.inc'
      include 'tapdta.inc'
c
c local
      integer i,iflag,iswtch,khar,kharz,kick,kjar,kmona,kmonz,l,
     . l1,l1m,l2,lim,ll,lmon,ltg0,mt,n,n1,n1m,n2,n2m,nr,nstart
      integer*2 ltest,nqlnt,nplnt
      include 'maxkidat.inc'
      integer*2 nkh,kh(maxki)
c
c save variables for later calls to HPARTM
      kick  = kicki
      nqlnt = nqlnti
      nstart= nstrti
      mt    = mti
 
c
c kh control vector set up from probe tape controls or
c rotation tape controls if nqlnt pos. or neg., respectively
c
c this setup is done even if length=0, since it serves three successive calls
c
      if(nqlnt.lt.0) then
         nplnt = -nqlnt
         if(nplnt.eq.3) then
            khar = Kert
            nkh  = Nkier
            do i = 1,nkh
               kh(i) = Kier(i)
            end do
         else
            khar = Kprt
            nkh  = Nkipr
            do i = 1,nkh
               kh(i) = Kipr(i)
            end do
            if(nplnt.eq.10) then
               Kmon = Lparm
               lmon = 8
            endif
         endif
      else
         nplnt = nqlnt
         khar  = Ksprb
         nkh   = Nkisb
         do i = 1,nkh
            kh(i) = Kisb(i)
         end do
      endif
      ltg0 = nplnt*100
 
c after first call to hpartl, the other two calls need no setup
c shorter entry 'hpartm' bypasses copying, etc.
      entry HPARTM(length,lhar,mhar,ntype)
      if(length.le.0) return
      lim = length
      l1  = 0
      l2  = 0
      l   = 0
      n   = 0
      l1m = 0
      call PCOPS(l,'HAR ',mt)
  100 iswtch = -l2
  200 iflag  = 0
      call PCOPY(l,lim,iflag,iswtch,lhar,mhar)
      if(iflag.gt.0) goto 500
c
c have we reached proper point in kh control vector
c for partial on probe or rotation tape
      if(l1.gt.0) goto 400
      ltest = ltg0 + ntype
      iflag = -1
      call PBDPRM(nkh,kh,nstart,khar,ltest,iflag)
      if(iflag.le.0) then
c
c we have reached the desired harmonic partials
         l1     = 1
         nr     = 0
         if(ntype.le.31) then
c
c zonal harmonics setup
            n1     = kh(nstart) - 1
            n2     = kh(nstart+1) - 1
            nstart = nstart+2
c
c tesseral harmonic set up
         else if(ntype.eq.46 .or. ntype .eq. 56) then
c
c resonant tesseral harmonics
crs   nr=kh(nstart)
crs   n1= kh(nstart+1)
crs   n2= kh(nstart+2)
c get actual ranges for resonant tesserals
crs   n1=((n1-1)*n1)/2+nr-1
crs   n2=((n2-1)*n2)/2+nr-1
crs   nstart=nstart+2
            call SUICID(
     .        'RESONANT TESSERAL HARMONIC PARTIALS NOT IMPLEMENTED YET '
     .        , 14)
         else
            n1     = (kh(nstart)*(kh(nstart)-1))/2 + kh(nstart+1)-1
            n2     = (kh(nstart+2)*(kh(nstart+2)-1))/2 + kh(nstart+3)-1
            nstart = nstart+4
         endif
c
c           interpolate for partial from integration tape
         kjar  = khar - n1
         kharz = kjar + n2
         goto 400
      else
 
c desired harmonic partials not on tape. use saved partials only
         l2 = 1
      endif
  300 if(iswtch.le.0) call SUICID(
     .' HARMONIC PARTIAL NOT ON PROBE OR ROTATION TAPE, STOP IN HPARTL '
     ., 16)
 
c partial not on tape, copy saved partial after all
      iswtch= -1
      goto 200

  400 if(l.lt.n1 .or. l.gt.n2) goto 300
      ll = l
crs   if(nr.eq.0) goto 155
crs   lp=8*(l-nr)+9
crs   lpp=sqrt(lp+0.1)
crs   if(lp.ne.lpp**2) goto 113
crs   ll=lpp/2
      khar = kjar + ll
      if(nqlnt.lt.0) then
 
         if(nplnt.eq.3) then
            Kert = khar
         else
            Kprt = khar
         endif
         call CPARTL(6,3,kick)
 
c if partial is w.r.t. j2 of moon, increment by  orbit partials
c         if(ntype.eq.31 .and. Imn .gt. 0 .and. lhar(1) .gt. 0)
c     .       call HPART2(kick)
         if(nplnt.eq.10 .and. Imn.gt.0) then
            if(l1m.eq.0) then
               iflag = -1
               call PBDPRM(Nkimn,Kimn,lmon,Kmon,ltest,iflag)
               if(iflag.le.0) then
c
c we have reached the desired harmonic partials
                  l1m = 1
                  if(ntype.le.31) then
c zonal harmonics
                     n1m = Kimn(lmon)-1
                     n2m = Kimn(lmon+1)-1
                     lmon = lmon+2
c tesseral harmonic
                  else
                     n1m = (Kimn(lmon)*(Kimn(lmon)-1))/2
     .                + Kimn(lmon+1)-1
                     n2m = (Kimn(lmon+2)*(Kimn(lmon+2)-1))/2
     .                + Kimn(lmon+3)-1
                     lmon = lmon+4
                  endif
                  kmona = Kmon-n1m
                  kmonz = kmona+n2m
               else
c desired harmonic partials not on tape. skip this contribution
                  n1m=0
                  n2m=0
                  kmonz=Kmon
               endif
            endif
c if desired partial is on orbit integration, interpolate and add
            if(l.ge.n1m .and. l.le.n2m) then
               Kmon=kmona+l
               call CPARTL(2,3,kick)
            endif
         endif

         call CPARTC(kick)
      else
         Ksprb = khar
         call CPARTL(4,2,kick)
      endif
c
c loop on l
      if(l.lt.length) goto 100
c
c flush out ksprb,etc. if there was interpolation from probe tape
  500 if(l1.gt.0) then
         khar   = kharz
         if(nqlnt.ge.0) then
            Ksprb = kharz
         else if(nqlnt.eq.-3) then
            Kert = kharz
         else
            Kprt = kharz
         endif
      endif
      if(l1m.gt.0) Kmon=kmonz
      return
      end
