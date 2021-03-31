      subroutine STATON(i8or9,jd,fract,hmx)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real*10 amin, aplong, ep, epw, etutc, fet, fract, frl, hmx, phi,
     .          sec, sthrst, sum, sum1, sumest, sumwst, tdiv, 
     .          thrust, thrusx
      real*10 timep, tot, totfr, z3, z3save, z4
      integer   i, i8or9, idd, ieast, if, ilev, iset, j, jd, jdl, 
     .          jj, l, nbeg, nbeg1, nend, nswt
 
c*** end of declarations inserted by spag
 
 
c
c     subroutine staton        june, 1976              l.weiner
c     this subroutine will emulate the stationkeeping logic
c     used in les8 and les9. this program simulates an 8th order
c     filter with input of 4 solar transit times from the calling
c     program lesthp and returns the thrust level and duration.
c     thrust is quantized to begin and end on tabular points
c     of integration output.  thrust is adjusted accordingly.
c     nb:  program may only be used with forward integration due
c          to the memory requirements of the filter.
c
c
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'sscon.inc'
      real*4    aicoe
      real*4    a1/3.350E-2/, a2/3.506E-4/, a3/5.09E-6/
      character*4    dirct(2)/'EAST', 'WEST'/
      integer*2 imn, idy, iyr, itm, ihr, imin
 
c sun perigee    1977    jan. 3.03889
      integer*4 jdsper/2443147/
      real*10 frsper/0.03889/
      real*10 ap(4),pfun(4),xs(4)

      real*10 PHI1
      real*4 arg
      PHI1(arg) = a1*sin(arg) + (a2*sin(2.*arg) + a3*sin(3.*arg))
c
c     --calling parameters--
c     i8or9 = 8 for les8
c             9 for les9
c     jd    = the current date
c     fract = the current time in pep
c     hmx   = the interval for constant step size integration
c
c     --input parameters--
c     jdtrns(i),frtrns(i) = julian day number and ephemeris time
c                           fraction of day of sun transit
c     itrns(i) = +1  90 degree sun transit
c              = -1 270 degree sun transit
c     thlev =  thrust level     (micropounds)
c
c     --output parameters--
c     jdfire(i,j),frfire(i,j) = julian day number and fraction of
c                               day for station-keeping thrusting
c                           j = 1 start, j = 2 end
c                           i =1,2,3,4 successive half orbit firing time
c     itrns(5) = number of transit points accumulated (0,1,2,3, or 4)
c                when =4 reset to 0
c     sthrst  = thrust level and direction
c               + for positive roll face thrusters   (micropounds)
c               - for negative roll face thrusters   (micropounds)
c     itrns(6) = 0   no thrusting
c              = 4   station keeping thrusting, prev thr finished
c              = 5  station keeping thrusting, prev thr unfinished
c
c
c     --other variables--
c      irout = 0 first time staton is called, 1 oth.
c      timep = time since sun perigee
c      thrust = thrust duration (given by the statonkp algorithm)
c      thlev = actual thrust level used in the satellite for staton
c      thrusx = adj thrust duration so that thrust begins and ends on
c               tabular points
c      vthrst = adjust thrust level so thrusx*vthrst=thlev*thrust
c
c     difference between ephemeris and universal time
c       staton works in ephemermis time
      i = jd - 2443145
      j = i/365
      if(i.lt.0) j = j - 1
      etutc = (48.184_10 + j)/86400._10
c etutc=48.184 in 1977
c
c inclination coefficient
      aicoe = Ainc*Ainc/4. + (Ainc**4)/24.
      write(6,100) (Itrns(i),Jdtrns(i),Frtrns(i),i = 1,4),
     .              Itrns(5),Itrns(6)
  100 format(/' SUN TRANSITS', 4(i3,i8,f13.10),2I2)
      if(Itrns(5).ne.4) return
 
      if(Irout.ne.1) then
         write(7,150) Sslong,Damp,Coast,Satdr,Ainc,Aomeg
  150    format(' SSLONG,DAMP,COAST,SATDR,AINC,AOMEG:'/4D11.5, 2E11.5)
         sumest = 0.
         sumwst = 0.
      endif
c
c input to filter
      sum = 0.
      do l = 1, 4
         timep  = (Jdtrns(l) - jdsper) + (Frtrns(l) - frsper)
         tdiv   = timep/365.27_10
         timep  = MOD(tdiv,1._10)
         arg    = Twopi*timep
         phi    = PHI1(arg) + (aicoe/a1)*PHI1(2.*arg + 2.*Aomeg)
         aplong = -(Frtrns(l) - etutc)*360._10 + 90._10
         if(Itrns(l).eq.-1) aplong = aplong + 180._10
         aplong = MOD(aplong,360._10)
         sum1   = (aplong + phi/Convd - Sslong)
         sum1   = MOD(sum1,360._10)
         if(sum1.gt.90._10) sum1  = sum1 - 180._10
         if(sum1.gt.90._10) sum1  = sum1 - 180._10
         if(sum1.lt.-90._10) sum1 = sum1 + 180._10
         if(sum1.lt.-90._10) sum1 = sum1 + 180._10
         sum     = sum + sum1
         xs(l)   = sum1
         ap(l)   = aplong
         pfun(l) = phi/Convd
      end do
      sum = MOD(sum,360._10)
c
c z3 is the new drift, z4 the drift rate
      z3 = sum/4._10
      z4 = .5*(z3 - z3save)
      Itrns(5) = 0
      iset     = Itrns(6)
      nswt     = 0
      if(iset.ne.0) then
         if(jd.lt.Jdfire(iset,2)) nswt = 1
         if(jd.eq.Jdfire(iset,2) .and. fract.le.Frfire(iset,2))
     .       nswt = 1
      endif
      ep = z4 + z3/Damp
      write(6,200) (ap(l),pfun(l),l = 1,4),z3,z4,ep
  200 format(' APL,PHI', 4(3x,f8.3,f8.3),3x,'Z3 =', f8.4, '  Z4 =',
     .       f8.4, '  EP =', f8.4)
      write(6,300) (xs(l),l = 1,4)
  300 format(' X3', 4F19.6)
      write(7,400) z3,z4,Date(1),Date(2)
  400 format(1x,'Z3 =', f19.6, 10x, 'Z4 =', f19.6, 10x, 2A4)
      z3save = z3
      if(Irout.ne.0) then
c
c sign of thrust negative to go east, positive to go west
c for les9.  reverse signs for les8.
c
         if(ep.lt.0) then
            if(z4.lt.Coast) then
               thrust = MIN(.5_10,ABS(z3)/(2.*Satdr))
               sthrst = -Thlev
               sumest = sumest + 4.*thrust
               goto 500
            endif
         else if(ep.ne.0) then
            if(z4.gt.-Coast) then
               thrust = MIN(.5_10,ABS(z3)/(2.*Satdr))
               sthrst = Thlev
               sumwst = sumwst + 4.*thrust
               goto 500
            endif
         endif
      endif
      Itrns(6) = 0
      Irout    = 1
      if(nswt.eq.1) Itrns(6) = iset
      return
  500 if(i8or9.eq.8) sthrst = -sthrst
c
c if prev thrusting has not finished, save last thrusting times
c in the first entries of the jd and frfire arrays.
      if(nswt.eq.0) then
         nbeg     = 1
         nend     = 4
         Itrns(6) = 4
      else
         Itrns(6)     = 5
         Jdfire(1,1) = Jdfire(iset,1)
         Jdfire(1,2) = Jdfire(iset,2)
         Frfire(1,1) = Frfire(iset,1)
         Frfire(1,2) = Frfire(iset,2)
         Vthrst(1)    = Vthrst(iset)
         nbeg         = 2
         nend         = 5
      endif
c
c     change thrust duration and level to match integral number
c     of tabular points.  thrust must span at least 2 tabular intervals
c     in order for sbout and lesthp to know whether we are in the
c     start-up integration procedure.
      ilev = thrust/hmx
      if(thrust.gt.0._10) ilev = ilev + 1
      if(ilev.eq.1) ilev = 2
      thrusx = MIN(.5_10,ilev*hmx)
      if(thrusx.gt.0._10) sthrst = sthrst*thrust/thrusx
c
c
c return the firing times. assume the closer of fire1 or fire2
c from last sun transit as the commencement of thrusting
      Jdfire(nbeg,1) = Jdtrns(4)
      fet = Frtrns(4)
      if(fet.lt.Fire1 .or. fet.ge.Fire2) Frfire(nbeg,1)  = Fire1
      if(fet.ge.Fire1 .and. fet.lt.Fire2) Frfire(nbeg,1) = Fire2
      if(fet.ge.Fire2) Jdfire(nbeg,1) = Jdfire(nbeg,1) + 1
      Frfire(nbeg,2) = Frfire(nbeg,1) + thrusx
      if = Frfire(nbeg,2)
      Frfire(nbeg,2) = Frfire(nbeg,2) - if
      Jdfire(nbeg,2) = Jdfire(nbeg,1) + if
      Vthrst(nbeg)    = sthrst
c
c generate new start and end of thrusting times
      nbeg1 = nbeg + 1
      do j = 1, 2
         do l = nbeg1, nend
            Frfire(l,j) = Frfire(l - 1,j) + .5
            if = Frfire(l,j)
            Frfire(l,j) = Frfire(l,j) - if
            Jdfire(l,j) = Jdfire(l - 1,j) + if
            Vthrst(l)    = sthrst
         end do
      end do
c difference between fire start and end is thrusx at vthrst level
c at end of output statement, have thrust as commanded by logic
c at the actual (thlev) thrust level
      write(6,600) (Jdfire(i,1),Frfire(i,1),Jdfire(i,2),Frfire(i,2)
     .              , Vthrst(i),thrust,Thlev,i = 1,nend)
  600 format(' STATION KEEPING THRUSTING', i8, 0pf13.10, i8, 0pf13.10,
     .       1pd17.9, 9x, '(', 0pf10.8, ' AT', 0pf8.2, ' ULB)',
     .       /(26x,i8,0pf13.10,i8,0pf13.10,1pd17.9,9x,'(',0pf10.8,
     .       ' AT',0pf8.2,' ULB)'))
c
c output total east and west thrust, as well as net thrust
      epw   = sumest - sumwst
      ieast = 1
      if(epw.lt.0._10) ieast = 2
      epw   = ABS(epw)
      idd   = epw
      frl   = epw - idd
      totfr = frl*86400.
      ihr   = totfr/3600.
      tot   = totfr - ihr*3600
      amin  = tot/60.
      write(6,700) sumest,sumwst,dirct(ieast),idd,ihr,amin,Thlev
  700 format(f12.6,' DAYS EAST THRUST', f14.6, ' DAYS WEST THRUST',
     .       6x, 'NET THRUST ', a4, ':', i6, 'DY', i4, 'HR', f6.2,
     .       'MIN   AT', f8.2, ' ULB')
c
c punch out firing times to make thruster log which can be
c altered if does not agree with actual satellite actions
c and then used as input to pep
      do l = nbeg, nend
         do j = 1, 2
            frl = Frfire(l,j) - etutc
            jdl = Jdfire(l,j)
            if(frl.le.0.) then
               frl = frl + 1._10
               jdl = jdl - 1
            endif
            call MDYJUL(imn,idy,iyr,itm,jdl)
c
c convert frl to ihr,imin,sec
            totfr = frl*86400.
            ihr   = totfr/3600.
            tot   = totfr - ihr*3600
            imin  = tot/60.
            sec   = tot - 60*imin
            jj    = j
            if(j.eq.2) jj = -1
            write(8,720) i8or9,imn,idy,iyr,ihr,imin,sec,jj,
     .                    Vthrst(l)
  720       format('LES', i1, 1x, 2(i2,'/'), i2, 1x, 2(i2,':'), f6.3,
     .             i3, f9.3)
         end do
      end do
      return
      end
