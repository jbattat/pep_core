      subroutine PLPCRD(jd,fract)
 
      implicit none
 
 
c  j.f.chandler - 1976 sep
c    revised 1977 feb
c  read in planet tape for partials during integration of satellite.
c  supply interpolated positions, velocities, and partials wrt the
c  initial conditions or other parameters of the central body's (or
c  a target body's) orbit
c     planet tape is read either forward or backward in time
c  also performs the same functions for embary in a moon integration
c  called by morfn or sbfn
c
c entry point PLPCRD refers to the central body, implying kt=0
c entry point PLPCRDT selects a body via kti, copied into kt
c
c  initial call from emprd1 or plprd1 to entry RDPTAP or RDPTAPT,
c  with the same distinctions of kt
c
c arguments
      integer jd,kt,kti
      real*10 fract
c jd,fract - time needed for interpolation
c kt/kti - code for needed body: 0=>central body, other=>target(kt)

c array dimensions
      include 'globdefs.inc'
c common
      include 'emcke.inc'
      include 'emmips.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'orblun.inc'
      include 'prtcod.inc'
      include 'prtpin.inc'
      include 'tapdtplp.inc'
      include 'yvectplp.inc'

c external functions
      real*10 DOT
c local 
      real*10 fl,t2m
      integer i,iq,ivl,j,jj,k,k1,k2,k3,k4,ll,la,lb,m,mm,mtab,
     . n,n1,n2,ntab
      character*8 erms(12)/' BAD REC','ORDS ON ','DATA SET','***, REC',
     .          'ORD NO.=', '******, ', 'STOP IN ', 'PLPCRD  ',
     .          'VELOCITY', ' NOT ON ', ' IVEL=**', ' <** ON '/

c local status variables for the collection of bodies, not kept in common
      integer jdcst(0:i_mxtrg),nplrec(0:i_mxtrg),ntp(0:i_mxtrg)

      kt=0
      goto 10

      entry PLPCRDT(jd,fract,kti)
      kt=kti
   10 continue
c
c use elliptic approximation for i.c. partials if no tape provided
      if(Jplntg(kt).le.0) goto 700
      if(jd.le.Jdt1(kt) .or. jd.ge.Jdt2(kt)) goto 800
  100 n = (jd - Jdt(2,kt))*Idirt(kt)
      if(n.lt.0) then
c
c backspace logic
         n = (-n - 1)/Intt5(kt) + 4
         do i = 1, n
            nplrec(kt) = nplrec(kt) - 1
            backspace Jplntg(kt)
         end do
         goto 200
      else if(n.eq.0) then
         goto 600
      else
         m = n/Intt5(kt)
 
c test if any twiddling is needed
         if(m.lt.1) goto 600
         n = m - 3
         if(n.lt.0) then
c
c shift storage by one or two records
            Jdt(1,kt)    = Jdt(m+1,kt)
            Ipvelt(1,kt) = Ipvelt(m+1,kt)
            Tfract(1,kt) = Tfract(m+1,kt)
            Ibadt(1,kt)  = Ibadt(m+1,kt)
            if(m.eq.1) then
               Jdt(2,kt)    = Jdt(3,kt)
               Ipvelt(2,kt) = Ipvelt(3,kt)
               Tfract(2,kt) = Tfract(3,kt)
               Ibadt(2,kt)  = Ibadt(3,kt)
            endif
            mm = m*Nptspr(kt)
            n2 = (3 - m - Ibadt(3,kt))*Nptspr(kt)
            n1 = 1 + Ibadt(1,kt)*Nptspr(kt)
            do k = n1, n2
               n = k + mm
               do j = 1, Ipart(kt)
                  do i = 1, Lmvlt(kt)
                     Trgbd(i,k,j,kt) = Trgbd(i,n,j,kt)
                  end do
               end do
            end do
            mm = 4 - m
            goto 300
         else if(n.ne.0) then
            do i = 1, n
               read(Jplntg(kt), end=800, err=110)
               goto 130
  110          read(Jplntg(kt))
 
c second read of error record might not be needed if system changes
               call PAGCHK(60, 1, 0)
               write(Iout, 120) nplrec(kt), Jplntg(kt)
  120          format(' **** ERROR ON RECORD', i6,
     .                ' SKIPPED ON PLANET DATA SET', i3, '  IN PLPCRD')
  130          nplrec(kt) = nplrec(kt) + 1
            end do
         endif
         goto 200
      endif
c
c
c read planet data into storage
      entry RDPTAP(jd,kti)
      kt=kti

      nplrec(kt) = 2
      Lmvlt(kt)  = Limvel
      if(kt.eq.0 .and. Ler.lt.0) Lmvlt(kt) = 6
      if(Lmvlt(kt).gt.6) call SUICID('LIMVEL>6, STOP IN PLPCRD',6)
  200 mm = 1
  300 do ll  = mm, 3
         lb = ll*Nptspr(kt)
         la = lb-Nptspr(kt)+1
         Ibadt(ll,kt)  = 0
         Ipvelt(ll,kt) = 0
         Jdt(ll,kt)    = 0
         Tfract(ll,kt) = 0
         nplrec(kt)= nplrec(kt) + 1
         read(Jplntg(kt), err=350, end=800) Jdt(ll,kt),Tfract(ll,kt),
     .        ivl,(((Trgbd(i,j,k,kt),i=1,ivl),k=1,Ipart(kt)),j=la,lb)
         Ipvelt(ll,kt) = ivl
         goto 400
  350    Ibadt(ll,kt) = 1
         read(Jplntg(kt))
  400 end do
c
c reconstructs dates of start of bad records
      if(Ibadt(1,kt).gt.0 .and. Jdt(1,kt).le.0) then
         if(Ibadt(2,kt).le.0 .or. Jdt(2,kt).gt.0) then
            Jdt(1,kt) = Jdt(2,kt) - ISIGN(Intt5(kt),Idirt(kt))
            Jdt(3,kt) = Jdt(2,kt) + ISIGN(Intt5(kt),Idirt(kt))
            goto 500
         else if(Ibadt(3,kt).gt.0 .and. Jdt(3,kt).le.0) then
            if(jd.le.0) goto 900
            Jdt(1,kt) = jdcst(kt)+(Intt5(kt)*(nplrec(kt)-5)*Idirt(kt))
         else
            Jdt(2,kt) = Jdt(3,kt)-ISIGN(Intt5(kt),Idirt(kt))
            Jdt(1,kt) = Jdt(2,kt)-ISIGN(Intt5(kt),Idirt(kt))
            goto 500
         endif
      endif
      Jdt(2,kt) = Jdt(1,kt)+ISIGN(Intt5(kt),Idirt(kt))
      Jdt(3,kt) = Jdt(2,kt)+ISIGN(Intt5(kt),Idirt(kt))
  500 ntp(kt)   = -99
      if(jd.gt.0) goto 100
      jdcst(kt) = Jdt(1,kt)
      erms(1)=' BAD REC'
      erms(2)='ORDS ON '
      return
c
c interpolate to requested (jd,fract)
  600 n    = jd - Jdt(2,kt)
      fl   = (n + (fract-Tfract(2,kt)))/Tint(kt)
      ntab = fl
      if(fl.lt.0._10) ntab=ntab-1
      P(1) = fl - ntab
      ntab = ntab + 1
      mtab = ntab + 1
      P(3) = 1._10 - P(1)
      P(2) = P(1)**2
      P(4) = P(3)**2
      ivl  = Lmvlt(kt)
      if(ntab.ne.ntp(kt)) then
         k1 = mtab
 
c see if requested time is ok
         iq = Ibadt(2,kt) + Ibadt(3,kt)
         if(k1.lt.6) iq = iq + Ibadt(1,kt)
         if(iq.gt.0) goto 900
         k3 = 1
         if(mtab.ne.ntp(kt)) then
            if(ntab.eq.(ntp(kt)+1)) then
               k1 = k1 + 1
               k3 = 2
            else
               k2 = 2
               goto 650
            endif
         endif
         k2 = 1
         k4 = 3 - k3
c
c shift y-vectors
         do i = 1, 5
            do j = 1, ivl
               do jj = 1, Nqt(kt)
                  Yplp(i,j,k4,jj,kt) = Yplp(i,j,k3,jj,kt)
               end do
            end do
         end do
  650    ntp(kt) = ntab
         if(ivl.gt.Ipvelt(2,kt)) then
            erms(1)=erms(11)
            erms(2)=erms(12)
            call EBCDI(Ipvelt(2,kt),erms(1)(7:8),2)
            call EBCDI(ivl,erms(2)(3:4),2)
            goto 900
         endif
c ntab,mtab,k1 are all relative to start of 2nd record
         k1=k1+Nptspr(kt)-5
         call YPLPCD(Trgbd(1,k1,1,kt),Yplp(1,1,k3,1,kt),Nqt(kt),ivl,k2,
     .    Krt(1,kt))
      endif
      call PLPTRP(kt)

c copy into perturbing planet block for central or target body
c 'further quantities' not needed if integrated body is a satellite
      n=Nplpt(kt)
      if(n.le.10) then
         do i = 1,ivl
            Xpert(i,n) = Ytp(i,1,kt)
         end do
      endif
      if(n.le.9) then
         if(Lps.le.0) then
            Rpert2(n) = DOT(Xpert(1,n), Xpert(1,n))
            Rpert(n)  = SQRT(Rpert2(n))
            Rpert3(n) = Rpert2(n)*Rpert(n)
            do i = 1, 3
               Xpert3(i,n) = Xpert(i,n)/Rpert3(n)
            end do
         endif
      endif

      return

c elliptic approximation if necessary
  700 continue
      if(Nqt(kt).le.1) return
      if(Nplpt(kt).ne.10) then
         fl = jd-T0mpt(kt)
         fl = fl+fract
         call EMIPT(1,fl,kt)
         do i = 1, 6
            do jj = 1, 6
               Dtcor(i,jj,kt) = Dympt(i,jj,kt)
            end do
         end do
c elliptic approx. for partial w.r.t. planet's mass
         if(Mppt(kt).gt.0) then
            jj=Mppt(kt)-1
            t2m=Elptg(7,kt)*Gauss*Gauss*0.5_10/Elptg(12,kt)
            do i=1,3
               Dtcor(i,jj,kt)=Ympt(i+3,kt)*t2m
            end do
         endif
      else
         call LUNORB(jd,fract,-1)
         do i = 1, 6
            do jj = 1, 6
               Dtcor(i,jj,kt) = Dylun(i,jj)
            end do
         end do
      endif

      return
c
c error stops
c normally 'bad records on data set'
c 'date not on data set...' encode date into message
  800 call EBCDI(jd,erms,8)
      erms(2) = erms(10)
 
c encode record no. into message
  900 call EBCDI(nplrec(kt),erms(6),6)
 
c encode data set number into message
      call EBCDI(Jplntg(kt),erms(4),3)
      call SUICID(erms,16)
      stop
      end
