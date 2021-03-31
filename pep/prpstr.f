      subroutine PRPSTR(ipct, i, kick, kobj, cor)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   id, igo, ij, ip, ip1
      real*10 xfactr
 
c*** end of declarations inserted by spag
 
 
      integer*4 ipct, i, kick, kobj
      real*4    cor(2)
c
c        r.b. goldstein  april 1978
c
c           kick= 1  prpstr called for radar-link observable
c           kick= 4  prpstr called for fermtr-link observable
c
c        functions: for values of ipct
c             -2:  store the correction in fcal0(i)
c             -1:  store the correction in fcal(i)
c              0:  if kick=1, form doppler corrections from fcal0
c                  and fcal and store in doppler slot
c                  if kick=4, form accumulated phase vlbi corrections
c                   from fcal and fcal0 and store in delay slot
c              1:  store corrections in cal(i) for delay
c              2:  store corrections in cal(i) for rate
c              3:  shift fcal to fcal0 for contiguous
c                  observations
c
      include 'difnct.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'number.inc'
      include 'obscrd.inc'
 
      real*10 freqtr(2)
      equivalence(Dstf(6), xfactr), (Dstf(7), freqtr)
 
      integer*2 nintrf
      equivalence(Kob(12), nintrf)
 
      include 'prpgat.inc'
 
      logical*4 init/.false./
      real*4    sgn, cor2
      integer*2 ni/3/, inch(10)/6, 8, 10, 7*0/
 
      real*10 frqtr0(2)
 
 
c setup controls
      if(  .not. (init) ) then
         ij   = Jct(3)
         ip1  = mod(ij/32, 2)
         init = .true.
      endif
      igo = ipct + 3
c
c
c-----------radar-link observables--------------------------------------
c
      if( kick .eq. 1 ) then
         cor2 = cor(2)
         if( igo .eq. 1 ) go to 200
         if( igo .eq. 2 ) then
c
c ipct= -1  store delay correction in fcal
            Fcal(i, kobj) = cor(1) + cor2
            return
         else if( igo .eq. 3 ) then
            go to 300
         else if( igo .eq. 4 ) then
         else if( igo .eq. 5 ) then
         else if( igo .eq. 6 ) then
            go to 500
         else
            go to 100
         endif
         go to 400
      endif
c
c
c-----------interferometry-link observables-----------------------------
c
  100 if( kick .ne. 4 ) call SUICID(
     .' CANNOT CALCULATE PROPAGATION CORRECTION FOR OPTICAL OBSERVABLE '
     ., 16)
      if( Nsite2 .eq. 0 ) cor(2) = 0.
      cor2 = -cor(2)
      if( igo .eq. 2 ) then
         Fcal(i, kobj) = cor(1) + cor2
         return
      else if( igo .eq. 3 ) then
         go to 300
      else if( igo .eq. 4 ) then
      else if( igo .eq. 5 ) then
      else if( igo .eq. 6 ) then
         go to 500
      else
         go to 200
      endif
      go to 400
 
c
c ipct= -2  store delay correction at initial time
  200 Fcal0(i, kobj) = cor(1) + cor2
      if( kick .eq. 4 .and. nintrf .gt. 0 ) Fcal0(i, kobj) = 0.
      return
c
c           ipct= 0  form accumulated phase correction and store it in
c                    delay slot for interferometry calculated correction
c                    or doppler slot for radio tracking (kick=1)
c                    for internal corrections only (i<10 or i=19)
c
c           change sign of charged particle corrections (since xmtr freq
c           are needed, this call after completion of both objects, and
  300 id = i
      if( id .le. 10 .or. id .eq. 19 ) then
         if( kick .eq. 1 ) id = id + 1
         sgn = 1.
         do ip = 1, ni
            if( id .eq. inch(ip) ) sgn = -1.
         end do
         if( kick .eq. 1 ) then
            Cal(id) = sgn*(Fcal0(i,kobj) - Fcal(i,kobj))/Tcs
         else
            if( Nk1 .le. 0 ) then
               frqtr0(1) = freqtr(1)
               frqtr0(2) = freqtr(2)
            endif
            if( kobj .eq. 1 ) Cal(id)
     .          = sgn*(Fcal(i,1)*freqtr(1) - Fcal0(i,1)*frqtr0(1))
     .          *xfactr
            if( kobj .eq. 2 ) Cal(id) = Cal(id)
     .          - sgn*(Fcal(i,2)*freqtr(2) - Fcal0(i,2)*frqtr0(2))
     .          *xfactr
         endif
         if( ip1 .eq. 1 ) then
            if( Line .gt. 57 ) call OBSPAG
            write(Iout, 320) kobj, i, Fcal0(i, kobj), Fcal(i, kobj),
     .                       id, Cal(id)
  320       format(' PRPSTR: KOBJ=', i2, '  FCAL(', i2, ')=', 1p,
     .             2E18.8, '  CAL(', i2, ')=', e18.8)
            Line = Line + 1
         endif
      endif
      return
c
c ipct= 1,2  store delay, rate in cal
  400 if( kobj .eq. 1 ) Cal(i) = cor(1) + cor2
      if( kobj .eq. 2 ) Cal(i) = Cal(i) - (cor(1) + cor2)
      return
c
c ipct=3  shift fcal to fcal0 for contiguous obs.
  500 do i = 1, 19, 2
         Fcal0(i, kobj) = Fcal(i, kobj)
      end do
 
      return
      end
