      subroutine NAMPRM(numpar, names, iskale, xnom)
 
      implicit none

c
c d. white  april 1973  subroutine namprm
c
c         generates an array of names for all parameters adjusted in
c         this run.  most names have two 8 character parts.  formats are
c         1.  observing site coordinates
c             sitename rad, long, lat, up, west, or north
c         2.  equator - equinox corrections
c             sitesrie denox, dequat, or dlat
c         2a. sky corrections (star catalog error models)
c             ctlgname nn
c         3.  solar system parameters
c             mass     bodyname
c             prmterxx
c             relfct, gmvary, sunhar, etc.
c         4.  et-ut2 dates
c             et-ut2   jd
c         5.  spot coordinates
c             spot     rad, long, or lat
c         6.  radar biases
c             rsitssit srienp(td or ds)
c             e.g., 14ds14ds mm7104ds
c         7.  phase corrections
c             sitesrie  plxxphcx
c             e.g., 6usn6956 pl01phc3
c         8.  body ic's
c             bodyname a (or e, inc, asc, etc)
c             bodyname x (or y, x', etc)
c             bodyname r (or ra, flaz, etc)
c             bodyname cond(x)
c         9.  body parameters
c             bodyname con(x)
c         10. body harmonics
c             bodyname j(x), c(x,x), or s(x,x)
c         11. body rotation/shape parameters
c             first character of bodyname is minus sign, except for
c             planets, which are called 'shpmrcry', etc.
c         12. pulsar parameters
c             plsrxxxx con(x)
c
c         parameters
      character*8 names(2, 1)
      real*10 iskale(1), xnom(1)
      integer*4 numpar

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'aprtbf.inc'
      include 'dtparm.inc'
      include 'empcnd.inc'
      include 'eqenox.inc'
      include 'ethhar.inc'
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'monhar.inc'
      include 'namtim.inc'
      include 'param.inc'
      include 'plnhar.inc'
      integer*4 ngdpts(4)
      equivalence (Scontl(1,9), ngdpts)
      include 'prmnms.inc'
      include 'psrstf.inc'
      include 'scoef4.inc'
      include 'skystf.inc'
      include 'stcord.inc'
c
c external
      real*10 SKALE, SKALEINIT
c
c locals
      integer*4 k
      character*8 coord(6)/'RAD','LONG','LAT','UP','WEST','NORTH'/
      character*8 correc(3)/'DENOX', 'DEQUAT', 'DLAT'/
      character*8 massn/'MASS'/
      character*8 temp
      character*1 tem(8)
      equivalence (tem(1), temp)
      character*8  bodies(10)/'MERCURY', 'VENUS', 'EARTH', 'MARS',
     .          'JUPITER', 'SATURN', 'URANUS', 'NEPTUNE', 'PLUTO',
     .          'MOON'/,
     .    shpbod(10)/'SHPMRCRY','SHPVENUS','SHPEARTH','SHPMARS ',
     .          'SHPJUPTR', 'SHPSATRN', 'SHPURANS',
     .          'SHPNEPTN', 'SHPPLUTO', 'SHPMOON '/
      character*8 erotat /'EROTAT'/, mrotat /'MROTAT'/
      character*8 embary /'EMBARY'/
      character*1 minus /'-'/
      character*8 blanks /'        '/
      character*8 dtwrd(4)/'ET-UT2','A1-UT1','XWOB','YWOB'/
      character*8 zed /'Z'/
      character*8 body/'BODY    '/
      character*8 plsr/'PLSR    '/, con1/'CON( )'/, con2/'CON(  )'/
      real*10 skaledummy
      integer*4 dummy,i,ii,j,jj,l,labs,m,missng,nd,nplic,nsky
c
c initialize
      k = 0
      skaledummy = SKALEINIT(dummy)
c
c observing site coordinates
      do j = 1,Numsit
         do i = 1,6
            if(Lscrd(i,j).ne.0) then
               k = k + 1
               call MVC(Site(1,j), 1, 8, names(1,k), 1)
               names(2, k) = coord(i)
               if(i.le.3) then
                  iskale(k)=SKALE(71 + i/2)
               else
                  iskale(k)=1._10
               endif
               if(Kscrd(j).lt.0 .and. i.eq.3) then
                  names(2,k)= zed
                  iskale(k)= SKALE(71)
               endif
               xnom(k)= Scord(i,j)/iskale(k)
            endif
         end do
      end do
c
c equator - equinox corrections
      do j = 1, Numeqn
         do i = 1, 3
            if(Leqn(i,j).ne.0) then
               k = k + 1
               call MVC(Eqnsit(j), 1, 4, names(1,k), 1)
               call MVC(Eqnser(j), 1, 4, names(1,k), 5)
               names(2, k) = correc(i)
               xnom(k)     = deqnx(j, i)
            endif
         end do
      end do
c
c sky corrections
      do j = 1, Numstr
         nsky = Nskycf(j)
         do i = 1, nsky
            if(Lskycf(i,j).ne.0) then
               k = k + 1
               names(1, k) = Ctlgnm(j)
               names(2, k) = blanks
               nd = 1
               if(i.gt.9) nd  = 2
               if(i.gt.99) nd = 3
               call EBCDI(i, names(2,k), nd)
               xnom(k) = Skycf(i, j)
            endif
         end do
      end do
c
c solar system parameters
      missng = 0
      do i = 1, 100
         if(Lprm(i).eq.0) goto 200
         k = k + 1
         l = Lprm(i)
         xnom(k) = prmter(l)
         if(l.gt.30) then
c
c other params
            names(2, k) = blanks
            iskale(k)   = SKALE(l - 30)
            xnom(k)     = xnom(k)/iskale(k)
            do j = 1, Prmsmx
               if(l.eq.Iprms(j)) then
c
c special name
                  names(1, k) = Prms(j)
                  goto 100
               endif
            end do
c
c no special name
            temp = Qprmtr
            call EBCDI(l, tem(7), 2)
            names(1, k) = temp
         else
c
c masses
            names(1, k) = massn
            if(l.le.10) names(2, k) = bodies(l)
            if(l.eq.3) names(2, k)  = embary
            if(l.gt.10) then
               do j = 1, Numpln
                  if(l.eq.Nplnt(j)) then
                     if(Aplnt(j).eq.blanks) goto 10
                     names(2, k) = Aplnt(j)
                     goto 100
                  endif
               end do
               missng = 1
   10          temp   = body
               call EBCDI(l, tem(5), 2)
               names(2, k) = temp
            endif
         endif
  100 end do
  200 if(missng.gt.0) call SUICID(
     .   'ADJUSTING MASS OF NON-INPUT BODY, WARNING IN NAMPRM ',-13)
c
c embary and earth parameters
      nplic = 203
      call BDYPRM(Lem, Econd, bodies(3), nplic, k, names, iskale, xnom)
c
c earth rotation
      nplic = -203
      call BDYPRM(Ler, Ercond, erotat, nplic, k, names, iskale, xnom)
c
c dates for et-ut2, a1-ut1, and/or wobble
      jj = 0
      if(Jddt0.lt.0) jj = 1
      m = 1
      if(Jddt0.le.1) m = 2
      if(Jddt0.lt.0) m = 3
      do while( .true. )
         do ii = 1, Numdt
            i = ii + jj*200
            if(Ldt(i).ne.0) then
               k = k + 1
               names(1, k) = dtwrd(m)
               temp = blanks
               call EBCDI(Jddt(ii), temp, 7)
               names(2, k) = temp
               xnom(k)     = Dt(i)
            endif
         end do
         if(Jddt0.gt.0 .or. jj.ge.2) goto 300
         jj = jj + 1
         m  = m + 1
      end do
c
c earth harmonics
  300 if(Nezone.gt.1) call HRMPRM(Lezhar, 1, bodies(3), k, names,
     .                          iskale, 1, Nezone, xnom, Ezhar)
      if(Netess.gt.1) then
         call HRMPRM(Lechar, 1, bodies(3), k, names, iskale, 2, Netess,
     .               xnom, Echar)
         call HRMPRM(Leshar, 1, bodies(3), k, names, iskale, 3, Netess,
     .               xnom, Eshar)
      endif
c
c moon parameters
      nplic = 210
      call BDYPRM(Lmn, Mcond, bodies(10), nplic, k, names, iskale, xnom)
c
c moon rotation
      nplic = -210
      call BDYPRM(Lmr, Mrcond, mrotat, nplic, k, names, iskale, xnom)
c
c moon harmonics
      if(Nmzone.gt.1) call HRMPRM(Lmzhar, 1, bodies(10), k, names,
     .                          iskale, 1, Nmzone, xnom, Mzhar)
      if(Nmtess.gt.1) then
         call HRMPRM(Lmchar, 1, bodies(10), k, names, iskale, 2, Nmtess,
     .               xnom, Mchar)
         call HRMPRM(Lmshar, 1, bodies(10), k, names, iskale, 3, Nmtess,
     .               xnom, Mshar)
      endif
c
c moon spot coordinates
      call SPTPRM(10, k, names, iskale, xnom)
c
c moon radar biases
      call RBSPRM(10, k, names, iskale, xnom)
c
c moon optical phase corrections
      call PHSPRM(10, k, names, iskale, xnom)
c
c sun spot coordinates
      call SPTPRM(0, k, names, iskale, xnom)
c
c sun radar biases
      call RBSPRM(0, k, names, iskale, xnom)
c
c sun optical phase corrections
      call PHSPRM(0, k, names, iskale, xnom)
c
c main loop for planet parameters
      do i = 1, Numpln
         if(Nplnt(i).eq.0) goto 400
c
c set up planet name in temp
         l     = Nplnt(i)
         nplic = ISIGN((Icnd(i)+2)*100, l) + l
         temp  = Aplnt(i)
         labs  = iabs(l)
         if(labs.gt.10) then
            if(labs.gt.30 .or. temp.eq.blanks) then
               temp = body
               call EBCDI(labs, tem(5), 2)
            endif
            if(l.lt.0) tem(1) = minus
         else
            temp = bodies(labs)
            if(l.lt.0) temp = shpbod(labs)
         endif
c
c planet ic's and con's
         call BDYPRM(Lpl(1,i), Pcond(1,i), temp, nplic, k, names,
     .               iskale, xnom)
c
c planet gravitational potential harmonics (nplnt positive) or
c planet shape (nplnt negative)
         do j = 1, 4
            if(Nplhar(j).eq.l) then
               if(l.lt.0 .and. Nshape(j).gt.0) then
 
c planet shape models other than spherical harmonic
                  call SHPPRM(j,4,temp,k,names,Nshape(j),ngdpts(j),xnom)
               else
                  if(Npzone(j).gt.1)
     .                call HRMPRM(Lpzhar(j,1), 4, temp, k, names,
     .                iskale, 1, Npzone(j), xnom, Pzhar(j,1))
                  if(Nptess(j).gt.1) then
                     call HRMPRM(Lpchar(j,1), 4, temp, k, names, iskale,
     .                           2, Nptess(j), xnom, Pchar(j,1))
                     call HRMPRM(Lpshar(j,1), 4, temp, k, names, iskale,
     .                           3, Nptess(j), xnom, Pshar(j,1))
                  endif
               endif
            endif
         end do
c
c no spots, biases, or phases for planet rotation
         if(l.gt.0) then
c
c planet spot coordinates
            call SPTPRM(l, k, names, iskale, xnom)
c
c planet radar biases
            call RBSPRM(l, k, names, iskale, xnom)
c
c planet optical phase corrections
            call PHSPRM(l, k, names, iskale, xnom)
         endif
      end do
c
c star coordinates
  400 call SPTPRM(-4, k, names, iskale, xnom)
c
c star radar bias
      call RBSPRM(-4, k, names, iskale, xnom)
c
c star atmosphere correction
      call PHSPRM(-4, k, names, iskale, xnom)
c
c pulsar parameters
      do j = 1, Numpsr
         call MVC(Sptpsr(j), 1, 4, plsr, 5)
         temp = con1
         nd   = 1
         do i=1,u_nmpsr
            l = Lpsrcn(i, j)
            if(l.le.0) goto 450
            k = k + 1
            names(1, k) = plsr
            if(l.gt.9) then
               temp = con2
               nd   = 2
            endif
            call EBCDI(l, tem(5), nd)
            names(2, k) = temp
            xnom(k)     = Psrcn(l, j)
         end do
  450 end do
c
c write out parameter names in ibuf1
      if(Ibuf1.gt.0) then
         write(Ibuf1) (names(1,i), names(2,i), i = 1, k)
         write(Ibuf1) (xnom(i), iskale(i), i = 1, k)
         endfile Ibuf1
         rewind Ibuf1
      endif
c
c write names, nominal values, and scale factors on lout
      if(Lout.gt.0) then
         do j = 1, k
            write(Intern, 460) j, (names(l,j), l = 1, 2), xnom(j),
     .                         iskale(j)
  460       format(i5, 2x, a8, 1x, a8, 1p, 2D24.16)
         end do
         endfile Intern
         rewind Intern
      endif
c
c set numpar = k (should = nparam)
      numpar = k
      return
      end
