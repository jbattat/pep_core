      subroutine SBEXP(sb,t,mtabx,sbfmt)
 
      implicit none

c
c        r.reasenberg  4/72  subroutine sbexp
c     extra printout for probe integration
c     updated 1977  j.f.chandler
c        additional options added feb. 1979 by rdr and rbg
c        updated jan 1980 by rdr and kcl
c        close-approach print added 1989  jfc
c
c array dimensions
      include 'globdefs.inc'
c
c parameters
      real*10 sb(6,i_mxeqn,5),t
      character*(*) sbfmt
      integer   mtabx
c
c common
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'intstf.inc'
      include 'namtimq.inc'
      include 'param.inc'
      include 'petuna.inc'
      character*8 aname
      equivalence (Name,aname)
      include 'sbembr.inc'
      include 'sbrot.inc'
      include 'sbstuf.inc'
      include 'sbthng.inc'

c external functions
      real*10 DOT

c local
      integer   jedz,jedzi
      real*10 djedz,djedzi
      real*10 temp(6)
      character*8 cntrnm/' CENTER '/,refbnm(3)/' EARTH  ','  SUN   ',
     .                 ' EMBARY '/
      character*4 lola(2)/'LONG',' LAT'/
      character*4
     .    month(12)/'JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG',
     1 'SEP','OCT','NOV','DEC'/
      character*6 cepch(2)/'1950.0','2000.0'/
      integer*2 icmo,icda,icyr,icce
      real*10 akm,audkms,aukm,clalph,clat,cldelt,clgcal,
     .          clong,persec,q,tsv
      integer   i,i3,ijd2,ijd3,ikp5,indcl,itaps5,itpclg,
     .          itpcls,j,jx,kkp3,kkp4,kkp6,kt,l,
     .          linexp,linusd,linust
 
      integer*4 indx/3/
c
c        flags     if true then
c        iflg      iout=kout
c        xflg      print 'x' in km and date
c        aflg      print 'x' in au and date
c        dflg      print partials
c        kflg      print kepler elements
c        tflg      print apsidal time
c        bflg      print body vectors
c        1         earth     s/c
c        2         earth     center
c        3         sun       s/c
c        4         sun       center
c        5        target     s/c
c        6         t(1)
c        7         t(2)
c        8         t(3)
c        9         t(4)
c        lflg      print landmark data
c        1         lat. and long of s/c w.r.t. central body
c        2         motion w.r.t. 'stationary' central body
c        3         lat and long. of s/c wrt central body at apsidal time
c        clflg     print close approaches: body vector + angles
c        1         earth-site figure of merit
c        2         earth
c        3-14      targets
c
      logical*1 iflg/.false./,xflg/.false./,dflg/.false./,
     .          kflg/.false./
      logical*1 bflg(9)/9*.false./
      logical*1 tflg/.false./
      logical*1 aflg/.false./
      logical*1 lflg(4)/4*.false./
      logical*4 clflg(14),incrsv(14),decrs
      integer*4 kpage/5000/, nexp/0/
      real*10 juldas,juldat,julsav/0._10/
      real*10 tcond(6)
      real*10 phiaps,motpi,taps(6)/6*0._10/,pnom,tt,ttt
      real*10 sbvp(3),sbvpc(3),clatd,alph,clongd,clonge,omegac,
     .          incd,tl,clat1,clon1,ince
      integer*4 kline/1000/
      real*4    peps(2)/2*0.05/
      integer*4 imin,jmin,naps/0/,jdclo
      character*3 ihrs,jhrs
      character*8 clname
      real*4    clcor(5),cldif(5),clcors(5,14),cldifs(5,14),
     .          clcomn(5),cldist,clt1,clt2,clt3,clt4,clfct1,
     .          clfct2,clgain,clszen,clczen,clclat,clslat,clsz0,
     .          clsz1,clgfct,cldt,cldts,clsz0p,clgfcp,clcorx
      equivalence (clcor(4),cldist)
c antenna gain (db) = g0 if zenith angle < lim0
c                     taper to g1 at zenith angle = lim1
c                     zero at greater angles
c noise: constant for ZA < lim0', increasing by 10**(-dg) out to lim1
c could make this model a set of input parameters
      character*8 clsite/'ARECIBO'/
c                       lat    g0   lim0  g1  lim1  scale    lim0'  dg
c
      real*4 clgmod(8)/18.23, 73.6, 15., 72.3, 20., 4.0E-17, 18.5, -0.2/
c     real*4 clgmod(6)/18.23, 71.,  11., 35.,  20., 2.0E-17/  (old)
 
      linusd = linexp
 
c print page head if special print needs more lines than are left
      if(iflg) then
         if(Line+linexp.ge.60) then
            call NEWPGT(Iout, Npage, jedz)
            Line = 2
         endif
      else if(kline+linexp.ge.60) then
         call NEWPGT(Kout, kpage, jedz)
         kline = 2
      endif
 
      juldat = t - 0.5_10
      juldas = juldat - djedz
      if(Kkp(3).ne.0) then
         if(aflg) write(Kout, sbfmt) juldas,
     .                    (sb(i,1,mtabx), i = 1, 6)
         if(xflg) then
            do i = 1, 3
               i3 = i + 3
               temp(i)  = sb(i, 1, mtabx)*aukm
               temp(i3) = sb(i3, 1, mtabx)*audkms
               end do
            write(Kout, 1100) temp
 1100 format(14x, ' X=', 1pd21.14, '(KM)', 11x, 'Y=', d21.14, '(KM)',
     .       11x, 'Z=', d21.14, '(KM)'/ 10x, ' DX/DT=', d21.14,
     .       '(KM/SEC)   DY/DT=', d21.14, '(KM/SEC)   DZ/DT=', d21.14,
     .       '(KM/SEC)')
         endif
         if(dflg) then
 
c write out partials
            jx = 1
            do while( .true. )
               jx = jx + 1
               write(Kout, 300) jx, (sb(i,jx,mtabx), i = 1, 6)
  300          format(i11, 1x, 1p, 6D20.13)
               if(jx.ge.Iparp) go to 50
               end do
         endif
   50    if(kflg) then
           call CHNCNC(Goose, sb(1,1,mtabx), tcond)
           akm = tcond(1)*aukm
           q   = akm*(1._10 - tcond(2))
           write(Kout, 400) cepch(Jct(13)+1), tcond(1), akm, tcond(2),
     .      q, (tcond(i),i=3,6)
  400 format('0ORBITAL ELEMENTS ', a6, '  A=', 1pd21.14, ' =', d21.14,
     .       '(KM)', 3x, 'E=', d21.14, 3x, 'Q=', d21.14, '(KM)'/
     .       ' INC=', d21.14, 5x, 'ASC=', d21.14, 5x, 'PER=', d21.14,
     .       5x, 'ANOM=', d21.14)
         endif
c
c process lflg requests
c
         if(Kkp(5).ne.0) then
            if(lflg(1)) then
               clong = ATAN2(Cslng(1), Cclng(1))/Convd
               clat = ASIN(Cslat)/Convd
               write(Kout, 800) aname, lola, cntrnm, clong, Cclng(1),
     .                          Cslng(1), clat, Cclat, Cslat
  800 format('0', a8, 8x, 2(16x,a4,'ITUDE',12x,'COS',12x,'SIN')/
     .       ' FROM ', a8, 9x, 2(5x,f20.15,2F15.10))
            endif
 
            if(lflg(2)) then
 
c find velocity perpendicular to radius
               alph = DOT(sb(1,1,mtabx), sb(4,1,mtabx))/Rsb2
               do j = 1, 3
                  sbvp(j) = sb(j + 3, 1, mtabx) - alph*sb(j, 1, mtabx)
                  end do
 
c transform to central body coords
               call PRODCT(Cntrot, sbvp, sbvpc, 3, 3, 1)
 
c find angle rates
               clatd  = sbvpc(3)/(Rsb*Cclat*Convd)
               alph   = sbvpc(1)*Cclng(1) - sbvpc(2)*Cslng(1)
               clongd = -alph/(Rsb*Cclat*Convd)
               clonge = clongd - omegac/Convd
 
c find inclinations
               incd = ATAN(clatd/(clongd*Cclat))/Convd
               ince = ATAN(clatd/(clonge*Cclat))/Convd
 
c print it out
               write(Kout, 900) clatd, clongd, incd, clonge, ince
  900 format('0LATDOT(DEG/DAY), LONGDOT(DEG/DAY)           ',
     .       2(f14.6,5x), 10x, 'INC(DEG) ', 2x, f14.6/
     .       '            WRT ROTATING PLANET               ',
     .       18x, f14.6, 26x, f14.6)
 
            endif
         endif
 
         if(tflg) then
            linusd = linusd - linust
 
c write out time of periapsis
            if(tcond(2).lt.1._10) then
              phiaps = MOD(tcond(6) + 90._10, 180._10) - 90._10
              phiaps = phiaps/360._10
              if(ABS(phiaps).le.peps(1)) then
                motpi   = Goose/tcond(1)/SQRT(tcond(1))/Twopi
                pnom    = 1._10/motpi
                taps(1) = phiaps/motpi
 
c save old apsidal time
                if(ABS(julsav-juldas).ge.pnom/3._10) then
                  naps    = naps + 1
                  julsav  = juldas
                  ijd3    = ijd2
                  taps(3) = taps(2)
                endif
                ijd2    = juldas
                taps(2) = -taps(1) + MOD(juldas, 1._10)
                call D2DMS(ihrs, imin, tt, ABS(taps(1))*24._10, 4)
                call D2DMS(jhrs, jmin, ttt, taps(2)*24._10, 4)
                write(Kout, 500) taps(1), ihrs, imin, tt, pnom,
     .                           taps(2), jhrs, jmin, ttt
  500 format('0TIME FROM APSIDAL EPOCH =', 1pd26.15, '  ABS TIME =(',
     .       a3, 'H', i3, 'M', 0pf8.4, 'S )', '   NOMINAL PERIOD= ',
     .       1pd28.15/' APSIDAL EPOCH = JD +', 1pd26.15, ' (', a3, 'H',
     .       i3, 'M', 0pf8.4, 'S )')
                linusd = linusd + 3
                if(lflg(3)) then
                  tl    = -taps(1)
                  clat1 = clat + tl*clatd
                  clon1 = clong + tl*clongd
                  write(Kout, 1000) clat1, clon1
 1000 format('0LAT AND LONG AT APSES(DEG)', 10x, 2(5x,f14.6))
                  linusd = linusd + 2
                endif
                taps(4) = ABS((ijd2-ijd3) + (taps(2)-taps(3)))
                if(naps.ge.2) then
                  taps(5) = taps(4)/pnom
                  itaps5  = 2._10*taps(5) + 0.5_10
                  tt      = itaps5
                  taps(6) = taps(4)*2._10/tt
                  persec  = taps(6)*86400._10
                  write(Kout, 600) taps(4), taps(5), taps(6), persec
  600 format(' TIME BETWEEN PRINTED APSES=', 1pd23.15, 7x, '*/P=',
     .       f8.4, 7x, 'PERIOD=', d21.14, ' =', d21.14, '(SEC)')
                  linusd = linusd + 1
                endif
              endif
            endif
         endif
      endif
      if(Kkp(4).ne.0) then
 
c earth     s/c
         if(bflg(1)) write(Kout, 700) aname, refbnm(l), Rpb(3),
     .                       (Pbcor(j,3), j = 1, indx)
  700 format('0VECTOR FROM ', a8, ' TO ', a8, 10x, 'RANGE =', f21.12/
     .       5x, 1p, 3D23.15, 5x, 3D17.9)
 
c earth     center
         if(bflg(2)) write(Kout, 700) cntrnm, refbnm(l), Rpc(3),
     .                       (Pccor(j,3), j = 1, indx)
 
c sun       s/c
         if(bflg(3)) write(Kout, 700) refbnm(2), aname, Rb,
     .                       (Bcor(j), j = 1, indx)
 
c sun       center
         if(bflg(4)) write(Kout, 700) refbnm(2), cntrnm, Rc,
     .                       (Ccor(j), j = 1, indx)
         if(bflg(5)) then
            do i = 1, 4
               if(bflg(i+5)) then
                  kt = Ktrg(i)
 
c target          s/c
                  write(Kout, 700) aname, Aplnt(kt), Rtb(i),
     .                             (Tbcor(j,i), j = 1, indx)
               endif
               end do
         endif
      endif
      go to 1300
 
 
 
      entry SBEXPS(djedzi, jedzi)
c
c     entry sbexps
c        called from sbout
c
c        kkp(3) = state
c        state is a packed-bits number
c          1  print cartesian state
c          2  print all partials
c          4  print apsidal time if near
c          8  print elliptic elements  standard frame
c n.i.
c         16  print elliptic elements  central body of date
c         32  print elliptic elements  plane of the sky
c
c        kkp(4) = bodys
c        bodys is a packed-bits number
c          1  print earth - s/c vector
c          2  print earth - central body vector
c          4  print sun - s/c vector
c          8  print sun - central body vector
c         16  print target - s/c vector
c
c        kkp(5) = lmark
c        lmark is a packed-bits number
c          1  print latitude and longitude of s/c w.r.t. central body
c          2  print projected motion of s/c on  central body
c          4  print latitude and long. of s/c wrt central body at apsida
c              time if near
c
c        kkp(6) = close
c        close is a packed-bits number
c          1     print distance, snr, zenith angle at good approach
c          2     print vector of s/c w.r.t. earth at close approach
c          3-14  print vector of s/c w.r.t. target at close approach
c
      jedz = jedzi
      djedz = djedzi
c
      iflg = .false.
      xflg = .false.
      aflg = .false.
      dflg = .false.
      kflg = .false.
      tflg = .false.
      do i = 1, 9
         bflg(i) = .false.
      end do
      do i = 1, 4
         lflg(i) = .false.
      end do
      if(Iout.eq.Kout) iflg = .true.
      l = 3
      if(Kp(40).ge.0) l = 1
c
c look at state requests
      kkp3 = Kkp(3)
      if(mod(kkp3,2).eq.1) xflg = .true.
      if(mod(kkp3,2).eq.1) aflg = .true.
      kkp3 = kkp3/2
      if(mod(kkp3,2).eq.1) dflg = .true.
      kkp3 = kkp3/2
      if(mod(kkp3,2).eq.1) tflg = .true.
      kkp3 = kkp3/2
      if(mod(kkp3,2).eq.1) kflg = .true.
 
c change state requests
      if(iflg) then
         if(Kp(99).le.0) then
            if(Kp(98)**2.le.1) aflg = .false.
            if(Kp(98).eq.1) dflg    = .false.
         endif
      endif
      if(Iparp.le.1) dflg = .false.
      if( .not. kflg) tflg  = .false.
c
c look at body requests
      kkp4 = Kkp(4)
      do j = 1, 5
         if(mod(kkp4,2).eq.1) bflg(j) = .true.
         kkp4 = kkp4/2
      end do
 
c change body requests
      if(Nplnt.eq.3) bflg(1)  = .false.
      if(Ncentr.eq.3) bflg(2) = .false.
      if(Kp(33).lt.0) bflg(1) = .false.
      if(Kp(33).lt.0) bflg(2) = .false.
      if(Ncentr.le.0) bflg(3) = .false.
      if(Ncentr.le.0) bflg(4) = .false.
      if(Nplnt.gt.30) then
         if(Kp(30+Ncentr).lt.0) bflg(3) = .false.
         if(Kp(30+Ncentr).lt.0) bflg(4) = .false.
      endif
      if(Numtar.le.0) bflg(5) = .false.
      if(bflg(5)) then
         do j = 1, 4
            bflg(j + 5) = (j.le.Numtar)
         end do
      endif
c
c look at landmark requests
      if(Ncentr.gt.0) then
         ikp5 = Kkp(5)
         do i = 1, 4
            if(mod(ikp5,2).eq.1) lflg(i) = .true.
            ikp5 = ikp5/2
         end do
      endif
      if(.not.lflg(2) .or. .not.tflg) lflg(3) = .false.
c
c look at close-approach requests
      kkp6 = Kkp(6)
      do i = 1, 14
         clflg(i) = (mod(kkp6,2).eq.1)
         if(i-2.gt.Numtar) clflg(i) = .false.
         clcors(4, i) = -1E30
         incrsv(i)    = .true.
         kkp6 = kkp6/2
      end do
      cldts = 0._10
      tsv   = 0._10
      if(Con(1).le.0._10) clflg(1) = .false.
      if(Kp(33).lt.0 .or. Ncentr.eq.3) clflg(1) = .false.
      if(Kp(33).lt.0 .or. Ncentr.eq.3) clflg(2) = .false.
      if(clflg(1)) then
         clsz0  = SIN(clgmod(3)*Convd)
         clsz1  = SIN(clgmod(5)*Convd)
         clsz0p = SIN(clgmod(7)*Convd)
         clgfct = (clgmod(4) - clgmod(2))/(clsz1 - clsz0)
         clgfcp = clgmod(8)/(clsz1 - clsz0p)
         clclat = COS(clgmod(1)*Convd)
         clslat = SIN(clgmod(1)*Convd)
         clcors(5, 1) = -1E30
         clgcal = clgmod(6)*Con(1)**2
      endif
c
c find number of lines needed by special print
      linexp = 2
      if(xflg) linexp = linexp + 2
      if(dflg) linexp = linexp + Iparp - 1
      if(kflg) linexp = linexp + 3
      if(tflg) linexp = linexp + 4
      if(aflg) linexp = linexp + 1
      do j = 1, 4
         if(bflg(j)) linexp = linexp + 3
         if(bflg(5) .and. bflg(j+5)) linexp = linexp + 3
         if(lflg(j)) linexp = linexp + 3
      end do
      linust = 4
      if(lflg(3)) linust = linust + 3
 
      write(Iout, 100) linexp, iflg, xflg, aflg, dflg, tflg, kflg,
     .                 bflg, lflg, clflg
  100 format(' EXTENDED PRINT SETUP  LINEXP=', i5, 10x, 'FLAGS',
     .       ' I,X,A,D,T,K =', 6L2, ' BFLG =', 4L2, l3, 1x,
     .       4L2/' LFLG =', 4L2, ' CLFLG =', 14L2)
      aukm   = Aultsc*Ltvel
      audkms = aukm/86400._10
      Line   = Line + 3
      if(Ncentr.ne.0) then
         do i = 1,u_mxpl
            if(Ncentr.eq.Nqlnt(i)) then
               cntrnm = Aplnt(i)
               omegac = Twopi/Pcond(13, i)
               go to 200
            endif
         end do
      else
         cntrnm = refbnm(2)
         omegac = 1._10
      endif
 
c if not found, cntrnm stays 'center'
  200 return
 
c
c find time of day h, m, s  set line counter  write end of print
 1300 nexp = nexp + 1
      call D2DMS(ihrs, imin, tt, MOD(t,1._10)*24._10, 4)
      write(Kout, 1200) nexp, juldas, ihrs, imin, tt
 1200 format(' END    EXTENDED PRINT NUMBER =', i6, '    MJED =',
     .       f18.10, '(', a3, 'H', i3, 'M', f8.4, 'S )'/)
      if(iflg) Line = Line + linusd
      if( .not. iflg) kline = kline + linusd
      return
 
 
 
 
      entry SBEXPT(t)
c
c print close approaches, if any.  called from sbout.
c
      juldas = t - 0.5_10 - djedz
      cldt   = t - tsv
      linusd = 0
      do i = 1, 14
         if(clflg(i)) then
 
c usual quantity to minimize is clcor(4)=cldist
            indcl = 4
            if(i.gt.2) then
 
c target body distance and vector
               kt     = Ktrg(i - 2)
               cldist = Rtb(i - 2)
               clcor(1) = -Tbcor(1, i - 2)
               clcor(2) = -Tbcor(2, i - 2)
               clcor(3) = -Tbcor(3, i - 2)
               clname   = Aplnt(kt)
            else
 
c earth distance and vector
               cldist   = Rpb(3)
               clcor(1) = -Pbcor(1, 3)
               clcor(2) = -Pbcor(2, 3)
               clcor(3) = -Pbcor(3, 3)
               clname   = refbnm(l)
               if(i.le.1) then
 
c quantity to minimize is clcor(5) = -log(signal/noise)
                  indcl  = 5
                  clname = clsite
 
c compute zenith angle
                  clcorx = sqrt(clcor(1)**2 + clcor(2)**2)
                  clczen = (clcorx*clclat + clcor(3)*clslat)/cldist
                  clszen = sqrt(1. - clczen*clczen)
 
c compute antenna gain in db
                  clgain = clgmod(2)
                  if(clszen.gt.clsz0) clgain = clgmod(2)
     .                + (clszen - clsz0)*clgfct
                  if(clszen.gt.clsz0p) clgain = clgain
     .                + (clszen - clsz0p)*clgfcp
                  if(clszen.gt.clsz1 .or. clczen.lt.0.) clgain = 0.
                  itpclg = 0
                  if(clszen.gt.clsz0) itpclg = 1
                  if(clszen.gt.clsz0) itpclg = 2
                  if(clszen.gt.clsz1 .or. clczen.lt.0.) itpclg = 3
 
c save -log(gain**2/r**4) as clcor(5)
                  clcor(5) = -clgain*(2./10.) + 4.*ALOG10(cldist)
               endif
            endif
            decrs = clcor(indcl).lt.clcors(indcl, i)
 
c form differences from previous tabular point
            do j = 1, indcl
               cldif(j) = clcor(j) - clcors(j, i)
            end do
 
c see if just passed minimum
            if( .not. (decrs .or. incrsv(i))) then
 
c print page head if necessary
               if(iflg) then
                  if(Line+linusd.ge.56) then
                     call NEWPGT(Iout, Npage, jedz)
                     Line   = 2
                     linusd = 0
                  endif
               else if(kline+linusd.ge.56) then
                  call NEWPGT(Kout, kpage, jedz)
                  kline  = 2
                  linusd = 0
               endif
 
c use previous tabular epoch if crossed a cusp
               clt1   = 0.
               clt2   = 0.
               clt3   = 1.
               clfct1 = 0.
               clfct2 = 0.
               if(i.ne.1 .or. itpclg.eq.itpcls) then
 
c get interpolation coefficients for minimum
                 clt1   = cldifs(indcl, i)*cldt**2
                 clt2   = cldif(indcl)*cldts**2
                 clt3   = cldifs(indcl, i)*cldt - cldif(indcl)*cldts
                 clt4   = (clt1 + clt2)/(4.*(cldts+cldt)*clt3**2)
                 clfct1 = clt4*(clt1 - clt2 - 2.*cldif(indcl)
     .                    *cldts*cldt)/cldts
                 clfct2 = clt4*(clt1 - clt2 + 2.*cldifs(indcl,i)
     .                    *cldts*cldt)/cldt
              endif
              do j = 1, indcl
                 clcomn(j) = clcors(j, i) + cldifs(j, i)*clfct1
     .                       + cldif(j)*clfct2
              end do
              juldat = (tsv - 0.5_10 - djedz) + 0.5*(clt1 + clt2)/clt3
              jdclo=juldat + djedz + 0.5_10
              call MDYJUL(icmo,icda,icyr,icce,jdclo)
              icyr=icyr+100*(19+icce)
              clalph = ATAN2(clcomn(2), clcomn(1))/Convd
              if(clalph.lt.0.) clalph = clalph + 360.
              cldelt = asin(clcomn(3)/clcomn(4))/Convd
              if(i.gt.1) then
                write(Kout, 1310) juldat, icyr, month(icmo), icda,
     .                           clname, clcomn(4), clalph, cldelt
 1310           format('0CLOSE APPROACH AT', f10.3,
     .                  ' (',I4,1X,A3,I3,') TO ',A8, 1pd11.3,
     .                  ' AU,  RA,DEC=', 0p, 2F7.2, 2x, a8, f9.2,f5.1)
                linusd = linusd + 2
              else
                clgain = clgcal*10.**(-clcomn(5))
                if(clgain.ge.1.) then
                  clcorx = sqrt(clcomn(1)**2+clcomn(2)**2)
                  clczen = (clcorx*clclat+clcomn(3)*clslat)/clcomn(4)
                  if(abs(clczen).gt.1.) clczen = SIGN(1., clczen)
                  clszen = acos(clczen)/Convd
                  write(Kout, 1310) juldat, icyr, month(icmo), icda,
     .                          clname, clcomn(4), clalph, cldelt,
     .                          'SNR,ZEN=', clgain, clszen
                  linusd = linusd + 2
                  endif
                endif
            endif
            incrsv(i) =  .not. decrs
 
c save current differences and coordinates
            do j = 1, indcl
               cldifs(j, i) = cldif(j)
               clcors(j, i) = clcor(j)
            end do
         endif
      end do
      cldts  = cldt
      tsv    = t
      itpcls = itpclg
      if(linusd.eq.0) return
      linusd = linusd + 2
      goto 1300
c
c
      end
