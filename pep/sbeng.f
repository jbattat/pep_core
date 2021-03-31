      subroutine SBENG(ssav, sdel, s)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      real      ss
      integer   i, j, LEG, nm9, npvo, nsc, nvb
 
c*** end of declarations inserted by spag
 
 
c
c           r.reasenberg january 1972   subroutine sbeng
c           t. forni march 1976 fixed and added mariner6&7
c           spacecraft name used in &nmlst2 and sbeng must be identical
c           spacecraft engineering data inputed in this routine
c           returns sbmass, sbarea, ssav, sdel
c
      real*10 ssav, sdel, s
 
c common
      include 'engstf.inc'
c sbarea(1) - radiation pressure area
c (2) - principal drag area
c (3) - secondary drag area
      include 'inodta.inc'
c        internal to this routine
c
c        units    area in cm*cm   mass in grams
c
c        mariner 9 specific data
      data nm9/5/
      real*4    area9/10.361E4/
      real*4    sref9(6)/1102., 1107.2, 1269.4, 1271., 1316.3, 2000./
      real*10 mass9(5)/2*1.0E6, .6E6, .560084E6, .55189E6/
c mariner 9 specific data       end
c
c viking specific data
      data nvb/4/
      real      svb(05)/2949.375, 2951.14, 2980.354, 2982.073, 4000./
      real      mvb(4)/2.2E6, 2.192E6, 1.135E6, 1.048E6/
      real      avb(4)/4*2.6607E5/
c viking specific data          end
c
c pvo specific data
c
      data npvo/5/
      real*4    spvo(6)/3700., 3846., 3900., 4195., 4321., 4500./
      real*4    mpvo(5)/5.78E5, 3.64E5, 3.53E5, 3.48E5, 3.44E5/
      real*4    apvo(2)/5.067E4, 3.097E4/
c apvo is area seen from 2 directions
c
      real*10 sd/2.44E6_10/
      character*8 scidsv
c
c for determining which spacecraft
c nidv= dimension of arrays nscv & scidv which should be equal
c nscv array is pointers for the  go to statement
      character*8 scid
      integer*4 nidv/16/
      integer*2 nscv(16)/1, 2, 3, 4, 4, 5, 6, 7, 8, 9, 10, 11, 11, 12,
     .          13, 13/
      character*8 scidv(16)/' PLANET ', 'MARINR71', 'MRNR10  ',
     .          'LES8    ', 'LES9    ', 'MARINER5', 'MARNR4  ',
     .          'VKNG1   ', 'VKNG2   ', 'LAGEOS  ', 'MRNR2   ',
     .          'MARINER6', 'MARINER7', 'PVO     ', 'NAVSTAR', 'NS'/
c gps satellites have names navstarn (n between 1 & 9) or nsn
c (n between 1 and 99). extra logic is needed to identify these char
c
      ss = s - sd
 
      if(nsc.eq.1) go to 500
      if(nsc.eq.2) then
      else if(nsc.eq.3) then
c* start=200
c mrnr10
         Sbarea(1) = 1.087E5_10
         Sbarea(2) = Sbarea(1)
         Sbmass    = 5.E5_10
c* start=300
c les8  les9
         go to 500
      else if(nsc.eq.4) then
         go to 500
      else if(nsc.eq.5) then
c* start=400
c mariner5
         Sbarea(1) = 6.60519E4_10
         Sbarea(2) = Sbarea(1)
         Sbmass    = 2.4571E5_10
         go to 500
      else if(nsc.eq.6) then
c* start=500
c marnr4
         Sbarea(1) = 11.0E4_10
         Sbarea(2) = Sbarea(1)
         Sbmass    = 2.5880E5_10
         go to 500
      else if(nsc.eq.7) then
c* start=600
c vkng1  (a,b)
         Sbarea(1) = 2.6607E5_10
         Sbarea(2) = Sbarea(1)
         if(ss.gt.2949.375) then
 
c orbit  (vkng1-b)
            do i = 1, nvb
               if(ss.le.svb(i+1)) then
                  Sbarea(1) = avb(i)
                  Sbarea(2) = Sbarea(1)
                  Sbmass    = mvb(i)
                  ssav = sd + (svb(i) + svb(i+1))/2.
                  sdel = (svb(i+1) - svb(i))/2.
                  go to 600
               endif
            end do
            write(Iout, 100) svb(1), svb(nvb + 1), ss
            go to 200
         else
 
c cruise (vkng1-a)
            Sbmass = 3.3E6_10
            ssav   = 2442749.375_10
            sdel   = 200._10
            go to 600
         endif
      else if(nsc.eq.8) then
c* start=700
c vkng2  (d,e)
         Sbmass    = 1.E6_10
         Sbarea(1) = 2.6607E5_10
         Sbarea(2) = Sbarea(1)
         go to 500
      else if(nsc.eq.9) then
c* start=900
c lageos
         Sbarea(1) = 1.52E3_10
         Sbarea(2) = Sbarea(1)
         Sbmass    = 6.80E5_10
         go to 500
      else if(nsc.eq.10) then
c* start=1000
c mrnr2
         Sbarea(1) = 3.83E4_10
         Sbarea(2) = Sbarea(1)
         Sbmass    = 1.9822E5_10
         go to 500
      else if(nsc.eq.11) then
c* start=1100
c mariner6 and mariner7 are identical spacecrafts (r.reasenberg)
         Sbarea(1) = 9.0E4_10
         Sbarea(2) = Sbarea(1)
         Sbmass    = 0.384E6_10
         go to 400
      else if(nsc.eq.12) then
c* start=800
c pvo
         do i = 1, npvo
            if(ss.le.spvo(i+1)) then
               Sbarea(1) = apvo(2)
               Sbarea(2) = apvo(1)
               Sbarea(3) = apvo(2)
               Sbmass    = mpvo(i)
               ssav = sd + (spvo(i) + spvo(i+1))/2.
               sdel = (spvo(i+1) - spvo(i))/2.
               go to 600
            endif
         end do
         write(Iout, 100) spvo(1), spvo(npvo + 1), ss
         go to 200
      else if(nsc.eq.13) then
         go to 400
      else
         call SUICID(' SBENG CALL WITH UNKNOWN SPACECRAFT ', 9)
      endif
c
c        table look up format
c        n.s/c-id-string is the number of elements in the
c        area and mass vectors.  it is one less than the number of
c        elements in the time vector.
c
c* start=100
c        marinr71
      do j = 1, nm9
         if(ss.le.sref9(j+1) .and. ss.gt.sref9(j)) go to 300
      end do
      write(Iout, 100) sref9(1), sref9(nm9 + 1), ss
  100 format(' SREF = ', 2F10.2, '  SS = ', f12.5)
  200 call SUICID(' SBENG OUT OF RANGE ', 5)
  300 Sbmass    = mass9(j)
      Sbarea(1) = area9
      Sbarea(2) = Sbarea(1)
      ssav      = sd + 0.5*(sref9(j) + sref9(j+1))
      sdel      = 0.5*(sref9(j+1) - sref9(j))
      go to 600
c* start=1200
c gps satellites  -  navstar1-navstar8
c set m=280 kg and a to make a/m=2.146E-8_10 km**2/kg (carr spherical m
  400 Sbarea(1) = 6.0088E4_10
      Sbarea(2) = Sbarea(1)
      Sbmass    = 2.8E5_10
c go to 750              when additions made
c* start=2000
  500 ssav = s
      sdel = 10000.
  600 call PAGCHK(60, 3, 0)
      write(Iout, 700) scidsv, nsc, Sbarea(1), Sbmass, s, ssav, sdel
  700 format(' AT END OF CALL TO SBENG FOR ', a8, i4, ' AREA,MASS=',
     .       2D20.10 / ' TIMES (S,SSAV,SDEL)=', 2(f16.3,4x), f10.3)
      if(nsc.eq.12) write(Iout, 800) Sbarea(2), Sbarea(3)
  800 format(' OTHER TWO AREAS ', 2d22.10)
      return
c
c set up call
c
c* start=3000
      entry SBZENG(scid)
      scidsv = scid
      nsc    = 0
      if( LEG(7,1,scid,1,scidv(15)).ne.0 .and.
     .    LEG(2,1,scid,1,scidv(16)).ne.0 ) then
         do j = 1, nidv
            if(scid.eq.scidv(j)) then
               nsc = nscv(j)
               Sbarea(3) = 0._10
               go to 900
            endif
         end do
 
c name not found
         call PAGCHK(60, 3 + (nidv-1)/5, 0)
         write(Iout, 1000) scid, nidv, (scidv(i), i = 1, nidv)
      else
         nsc = 13
         Sbarea(3) = 0._10
      endif
  900 return
 1000 format('0SPACECRAFT NAME ', a8, ' NOT KNOWN IN SBENG', '  THE',
     .       i5, ' KNOWN NAMES ARE'/ (5x,10A10))
      end
