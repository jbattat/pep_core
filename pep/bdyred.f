      subroutine BDYRED(in0,nstop,jdpad,kall)
 
      implicit none

c     m.ash   jan 1972    subroutine bdyred
c     move quantities for earth,moon,planets,space probe motion and
c     rotation into appropriate locations (from &nmlst2 locations)

c arguments
      integer in0,nstop,jdpad,kall
c in0  - fortran unit for editing input stream
c nstop- cumulative error count
c jdpad- adjustment, if any, for all jd1-jd2 integration ranges
c        meaning of kall:
c          -1 - just initialize
c           0 - initialize and read in all nmlst2 namelists
c           1 - already initialized, just read in one namelist
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'empcnd.inc'
      real*10 erad,eflat,mrad,beta,gamma
      equivalence (Econd(7), erad), (Econd(8), eflat), (Mcond(7), mrad)
      equivalence (Mrcond(9), beta), (Mrcond(10), gamma)
      integer*4 zempcn/9600/   !r8=4800,r10=9600
      include 'ethhar.inc'
      integer*4 zemhar/46750/   !r8=25974,r10=46750
      include 'france.inc'
      integer*4 zfranc/26680/   !r8=19800,r10=26680
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'monhar.inc'
c monhar has the same length as ethhar and so uses 'zemhar' for zeroing
      include 'namtim.inc'
      integer*2 klan,npln3,npln10
      equivalence (klan,Numpln),(npln3,Nplnt(-3)),(npln10,Nplnt(-2))
      include 'param.inc'
      include 'plndta.inc'
      include 'plnhar.inc'
      integer*4 zplnhr/190/   !r8=190,r10=190
      include 'psrstf.inc'
      integer*4 zpsrst/4582/   !r8=2582,r10=4582
      include 'scoef4.inc'
      real*4 transf(u_stdsz,4,1000/u_stdsz)
      equivalence (Pzhar(1,1),transf(1,1,1))
      integer*2 klam
      equivalence (klam, Nmphar)
      integer*4 zscof4/40000/   !r8=24000,r10=40000
      include 'smlbdy.inc'
      integer*4 zsmlbd/24002/   !r4=7202,r10=24002
 
c earth,moon,planet control and data constants (&nmlst2)
      common/WRKCOM/ Bconst(u_nmbod),Bcon1(12),Btcon(30),Secb,Intb1,
     . Intb2,Beps(6),Bname,Jdb1,Jdb0,Jdb2,Itape,Limb(u_nmbod),I4dumm,
     . Nshp,Scntrl(9),Bz(u_mxzon-1),Bc((u_mxtes*(u_mxtes+1))/2-1),
     . Bs((u_mxtes*(u_mxtes+1))/2-1),Lbz(u_mxzon-1),
     . Lbc((u_mxtes*(u_mxtes+1))/2-1),Lbs((u_mxtes*(u_mxtes+1))/2-1),
     . Denptr,Kimb(u_nmprm),Imb,Nimb,Ncentb,Nzone,Ntess,Ihrb,Iminb,
     . Kkimb(100),Icndb,Numki,Ki(99),Small,Pulsar
      character*8 Bname
      real*10 Bconst,Bcon1,Btcon,Secb,Bz,Bc,Bs
      real*4 Beps,Scntrl
      integer*4 Intb1,Intb2,Jdb1,Jdb0,Jdb2,Itape,I4dumm,Nshp
      integer*2 Limb,Lbz,Lbc,Lbs,Denptr,Kimb,Imb,Nimb,Ncentb,
     .  Nzone,Ntess,Ihrb,Iminb,Kkimb,Icndb,Numki,Ki
      logical*4 Small,Pulsar
      real*4 shape(1000)
      equivalence (Scntrl(9),ngdpts)
      equivalence (Bz,shape)

c external functions
      integer*4 NSCAN

c local variables
      character*8 blank/'        '/
      integer i,ishp,j,jplnt,k,kkln,n1,ngdpts

c
c skip initialization if already done
      if(kall.le.0) then
         do j = -3, u_mxpl
            Npcent(j) = 0
            Icnd(j)   = 0
            Jdpl0(j)  = 0
            Inttyp(j) = 3
         end do
         do j = 1, u_mxpl
            Nplnt(j)  = 0
            Aplnt(j)  = blank
         end do
c special built-in cases for earth and moon
         Nplnt(-3)=3
         Nplnt(-2)=10
         Npcent(-2)=3
         Nplnt(-1)=-3
         Nplnt(0)=-10
c
c numpln = number of input planets,asteroids,satellites and
c            artificial space probes (exclusive of earth-moon
c            barycenter,moon,earth rotation,moon rotation)
c
c jdem0   = initial time for earth-moon barycenter integration
c jdmn0   = initial time for moon integration
c jder0   = initial time for earth rotation integration
c jdmr0   = initial time for moon rotation integration
c jdpl0(j)= initial time for planet nplnt(j) integration, j=1,numpln
c            initial time is julian day number, with numerical
c            integration commencing from midnight beginning of day of
c            calendar date corresponding to julian day number. if 0,
c            no integration. if initial conditions are adjusted, there
c            will be integration on second and subsequent least squares
c            analysis iterations because initial time will be set to
c            initial time on tape during comparison of theory and
c            observation. if negative, checkpoint restart mode.
c            see explanation of jdbdy0 above or explanation of jd0 in
c            subroutine bodred where these quantiies are read.
c
c Aplnt(j)=8 character planet name, j=1,u_mxpl=16
c Aplnt(-3)=' EMBARY '
c Aplnt(-2)='  MOON  '
c Aplnt(-1)=' EROTAT '
c Aplnt(0) =' MROTAT '
c Aplnt(17)='  MOON  '
c Aplnt(18)='  SUN   '
c earth, moon, sun names are setup in block data routine for main link.
c
c nplnt(j)= planet number, j=1,numpln
c           1 mercury
c           2 venus
c           4 mars
c           5 jupiter
c           6 saturn
c           7 uranus
c           8 neptune
c           9 pluto
c           11-30  other natural planets,asteroids,satellites
c           31,... artificial space probes
c           -1,-2,-4,...,-9,-11,... rotation for planet iabs(nplnt)
c           3,10,-3,-10 do not exist in nplnt(j),j=1,numpln, but instead
c             are built into nplnt(j),j=-3,0,  imagined
c             to be associated with earth-moon barycenter, moon,
c             earth rotation, moon rotation
c           0 does not exist in nplnt(j),j=1,numpln, imagined to be
c             associated with sun
c
c npcent(j)=central body for planet nplnt(j), j=1,numpln
c          -1 sun with integration in satellite-probe link
c           0  sun with integration in planet link if nplnt.le.30
c           3  earth (not earth-moon barycenter)
c           10 moon
c           other numbers as above
c
c
c           initialize earth,moon,planet init.cond.and parameters
         call ZFILL(Econd, zempcn)
         erad  = 6378.166_10
         eflat = 1._10/298.3_10
         mrad  = 1738._10
c
c erad    = equatorial radius of earth in kilometers
c eflat   = degree of flattening of earth
c mrad    = radius of moon in kilometers
c
c moment of inertia ratios for lunar libration
         beta  = .630E-3_10
         gamma = .228E-3_10
c
c           the remainder of earth-moon barycenter,moon,planet control
c           and data constants (not to be adjusted in least squares
c           analysis) are stored in peripheral data set iplcon
c
c
c           initialize earth gravitational potential harmonic coeff.
         call ZFILL(Ezhar, zemhar)
         call HARNTL(npln3, Ezhar, Echar, Eshar, Nezone, Netess)
c
c initialize moon gravitational potential harmonic coeff.
         call ZFILL(Mzhar, zemhar)
         call HARNTL(npln10, Mzhar, Mchar, Mshar, Nmzone, Nmtess)
c
c initialize coefs,l-vectors, and controls for planet grav.pot.
c or shape
         call ZFILL(Nshape, zplnhr)
         call ZFILL(Pzhar, zscof4)
         do j = 1, 4
            Szero(j) = .true.
         end do
c
c initialize quantities to go into iplcon peripheral data set
         call ZFILL(Dumcon, zfranc)
         do j = 1,20
            if(j.le.4) Kkk(88,j) = 3
         end do
c
c clear limited asteroid and pulsar quantities
         call ZFILL(Scond, zsmlbd)
         call ZFILL(Psrcn, zpsrst)
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c read earth,moon,planet control and data constants (&nmlst2)
         klan   = 0
         klam   = 0
         Numsml = 0
         Numpsr = 0
         if(kall.eq.-1) return
      endif
  100 call BODRED(in0,nstop,jdpad)
      if(Small) goto 280
      if(Pulsar) then
c
c segregate pulsars
         if(Numpsr.lt.u_mxpsr) then
            Numpsr = Numpsr + 1
            i = NSCAN(Bname, 5, ' ')
            if(i.lt.0) i = 4
            call MVC(Bname, i + 1, 4, Sptpsr(Numpsr), 1)
            Jdpsr0(Numpsr) = Jdb0
            do i = 1, u_nmpsr
               Psrcn(i, Numpsr)  = Bconst(i + 6)
               Lpsrcn(i, Numpsr) = Limb(i)
            end do
            Plspr(Numpsr)  = Bcon1(2)
            Ntypsr(Numpsr) = Nshp
         else
            nstop = nstop + 1
            write(Iout, 110) u_mxpsr
  110       format(' ***** MORE THAN', i4, ' PULSARS, ERROR IN BDYRED')
            if(Mout.gt.0) write(Mout, 110) u_mxpsr
         endif
         goto 400
      endif
c
c earth rotation constants
      if(Nimb.eq.-3) then
         jplnt = 3
         goto 250
 
      else if(Nimb.eq.-10) then
c moon rotation constants
         jplnt = 4
         goto 250
 
      else if(Nimb.eq.0) then
c end of planet list
         goto 500
c
c earth-moon barycenter constants
      else if(Nimb.eq.3) then
         jplnt = 1
         if(Nzone.gt.u_mxzon)
     .    call SUICID('NEZONE TOO LARGE, STOP IN BDYRED',8)
         if(Nzone.gt.0) Nezone = Nzone
         n1 = Nzone - 1
         do i = 1, n1
            if(Bz(i).ne.0._10 .or. Lbz(i).gt.0) goto 120
         end do
         goto 125
  120    do i = 1, n1
            Ezhar(i)  = Bz(i)
            Lezhar(i) = Lbz(i)
         end do
  125    if(Ntess.gt.u_mxtes)
     .    call SUICID('NETESS TOO LARGE, STOP IN BDYRED',8)
         if(Ntess.gt.0) Netess = Ntess
         if(Ntess.le.1) goto 250
         n1 = (Ntess*(Ntess+1))/2 - 1
         do i = 1, n1
            if(Bc(i).ne.0._10 .or. Bs(i).ne.0._10 .or.
     .         Lbc(i).gt.0 .or. Lbs(i).gt.0) goto 130
         end do
         goto 250
  130    do i = 1, n1
            Echar(i)  = Bc(i)
            Eshar(i)  = Bs(i)
            Lechar(i) = Lbc(i)
            Leshar(i) = Lbs(i)
         end do
         goto 250
c
c moon constants
      else if(Nimb.eq.10) then
         jplnt = 2
         do i = 1, 6
            if(Limb(i).gt.0 .and. Mdstsc.lt.0._10) Mdstsc = 0._10
         end do
         if(Jdb0.gt.0 .and. Mdstsc.lt.0._10) Mdstsc = 0._10
         if(Nzone.gt.u_mxzon)
     .    call SUICID('NMZONE TOO LARGE, STOP IN BDYRED',8)
         if(Nzone.gt.0) Nmzone = Nzone
         n1 = Nzone - 1
         do i = 1, n1
            if(Bz(i).ne.0._10 .or. Lbz(i).gt.0) goto 135
         end do
         goto 140
  135    do i = 1, n1
            Mzhar(i)  = Bz(i)
            Lmzhar(i) = Lbz(i)
         end do
  140    if(Ntess.gt.u_mxtes)
     .    call SUICID('NMTESS TOO LARGE, STOP IN BDYRED',8)
         if(Ntess.gt.0) Nmtess = Ntess
         if(Ntess.le.1) goto 250
         n1 = (Ntess*(Ntess+1))/2 - 1
         do i = 1, n1
            if(Bc(i).ne.0._10 .or. Bs(i).ne.0._10 .or.
     .         Lbc(i).gt.0 .or. Lbs(i).gt.0) goto 150
         end do
         goto 250
  150    do i = 1, n1
            Mchar(i)  = Bc(i)
            Mshar(i)  = Bs(i)
            Lmchar(i) = Lbc(i)
            Lmshar(i) = Lbs(i)
         end do
         goto 250
      endif
c
c planet rotation goes to same storage moves as planet motion
c planet constants
      if(klan.ge.u_mxpl) then
         nstop = nstop + 1
         write(Iout, 210) u_mxpl
  210    format('0 * * * MORE THAN', i3,
     .          ' INPUT PLANETS, ERROR IN BDYRED')
         goto 400
      endif
      klan = klan + 1
      jplnt = klan + 4
      Nplnt(klan)  = Nimb
      Npcent(klan) = Ncentb
      Aplnt(klan)  = Bname
  250 kkln = jplnt - 4
      Icnd(kkln) = Icndb
      do i = 1,30
         Pcond(i,kkln)  = Bconst(i)
         Dtcon(i,jplnt) = Btcon(i)
         Lpl(i,kkln)    = Limb(i)
      end do
      do i = 1,12
         Dumcon(i,jplnt) = Bcon1(i)
      end do
      do i = 1,6
         Dumeps(i,jplnt) = Beps(i)
      end do
      Jd1(jplnt)   = Jdb1
      Jdpl0(kkln)  = Jdb0
      Jd2(jplnt)   = Jdb2
      Inttyp(kkln) = Kimb(88)
      Int(jplnt)   = Imb
      Intp1(jplnt) = Intb1
      Intp2(jplnt) = Intb2
      Ihrp(jplnt)  = Ihrb
      Iminp(jplnt) = Iminb
      Secp(jplnt)  = Secb
      do i = 1,100
         Kkp(i,jplnt) = Kkimb(i)
         Kkk(i,jplnt) = Kimb(i)
      end do
      Ndumki(jplnt) = Numki
      do i = 1,Numki
         Kdumi(i,jplnt) = Ki(i)
      end do
      if(Itape.ge.0 .or. kkln.gt.0) Iplnt(kkln) = Itape
      if(kkln.gt.0 .and.
     . (Nzone.gt.1 .or. Ntess.gt.1 .or. Nshp.gt.0)) then
         if(klam.ge.4) then
            write(Iout,260)
  260       format(
     .'0* * * MORE THAN 4 PLANETS WITH GRAVITATIONAL OR SHAPE COEFS. ERR
     .OR IN BDYRED')
            nstop = nstop + 1
            goto 400
         endif
         klam = klam + 1
         Nplhar(klam) = Nimb
         Nshape(klam) = Nshp
         do i = 1,9
            Scontl(klam,i) = Scntrl(i)
         end do
         if(Nshp.le.0) then
            if(Nzone.gt.25) call SUICID(
     .          'NPZONE.GT.25, STOP IN BDYRED',7)
            if(Nzone.gt.0) Npzone(klam) = Nzone
            n1 = Nzone - 1
            do i = 1,n1
               Pzhar(klam,i)  = Bz(i)
               Lpzhar(klam,i) = Lbz(i)
            end do
            if(Ntess.gt.20) call SUICID(
     .          'NPTESS.GT.20, STOP IN BDYRED',7)
            if(Ntess.gt.0) then
               Nptess(klam) = Ntess
               n1 = (Ntess*(Ntess+1))/2 - 1
               do i = 1,n1
                  Pchar(klam,i)  = Bc(i)
                  Pshar(klam,i)  = Bs(i)
                  Lpchar(klam,i) = Lbc(i)
                  Lpshar(klam,i) = Lbs(i)
               end do
            endif
c planet shape for nshp gt 0
c special coding needed to transfer shape coefficients
c and their l-vectors from overlay arrays in
c /scoef/ into /scoef4/
         else if(Nshp.eq.1) then
 
c shape model is 2-dimensional fourier series
            Npzone(klam) = 122
            Nptess(klam) = 0
 
c transfer coefficients
            do i = 1,122
               Pzhar(klam,i)  = Bz(i)
               Lpzhar(klam,i) = Lbz(i)
            end do
         else if(Nshp.eq.2) then
c     nshp=2, local shape model-altitude grid
c     must transfer real*4 coefs from /scoef/ to
c     real*10 slots in /scoef4/ using array transf(_,4,___)
c     to maintain allignment of different planet models
c     to allow one planet to use r*4 grid model, while another
c     planet used r*8 fourier or sph. harmonic model
            ishp = 0
            do i = 1,1000/u_stdsz
               do k = 1,u_stdsz
                  ishp = ishp + 1
                  transf(k,klam,i) = shape(ishp)
                  Lpzhar(klam,ishp)= Lbz(ishp)
               end do
            end do
 
c check for non-zero elements of grid and l-vector
            do i = 1,ngdpts
               if(shape(i).ne.0. .or. Lbz(i).ne.0) then
                  Szero(klam) = .false.
                  goto 400
               endif
            end do
         endif
      endif
      goto 400
c
c segregate asteroids
  280 if(Numsml.lt.u_mxsml) then
         Numsml = Numsml + 1
         Nbsml(Numsml)  = Nimb
         Denpts(Numsml) = Denptr
         do i = 1,7
            Scond(i,Numsml) = Bconst(i)
         end do
         Jdsml0(Numsml) = Bcon1(1) + 0.50001_10
      else
         nstop = nstop + 1
         write(Iout,300) u_mxsml
  300    format(' ***** MORE THAN',i4,' ASTEROIDS, ERROR IN BDYRED')
         if(Mout.gt.0) write(Mout,300) u_mxsml
      endif
 
c coding for additional shape models goes here
  400 if(kall.eq.0) goto 100
c
c set value of mdstsc (moon distance unit in light seconds
c for compar link) if it has not been set yet
  500 if(Mdstsc.lt.0._10) Mdstsc = 2.127521368E-2_10
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
      return
      end
