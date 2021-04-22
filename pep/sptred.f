      subroutine SPTRED(in0,nstop,init)
 
      implicit none

c
c m.e.ash/k.m.becker     august 1967    subroutine sptred
c names and coordinates for spots on other bodies are initialized
c and read
c
c arguments
      integer*4 in0,nstop
      logical   init
c
c modified for *command july 1978  r.b. goldstein
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'empcnd.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'namtim.inc'
      include 'sptcrd.inc'
c
c for k=1,...,numspt we have
c spot(k)  = 4 character observed spot name
c nsplnt(k)= planet number of body on which spot lies (if.ge.0)
c nsplnt(k)= -1,-2,-3 for various types of earth spot coordinates
c            (changed to 3 before return to calling program)
c nsplnt(k)= -4 spot is a star or quasar
c nsplnt(k)= -5 spot is a pulsar
c spcord(1,k) =radius of spot relative to center of body (kilometers)
c spcord(2,k) =east longitude of spot (degrees).  on moon or planet
c              longitude is angle in right hand coordinate system,
c              whereas on earth longitude is angle in left hand
c              coordinate system for observing sites.
c              for earth spots it is right handed coordinate system
c spcord(3,k) =latitude of spot (degrees).
c lspcrd(i,k) =0 spcord(i,j) not adjusted (i=1,3)
c lspcrd(i,k) =1 spcord(i,j) adjusted in least squares analysis (i=1,3)
c kspt        =flag indicating type of coordinates
c            if nsplnt=3 and kspt<0, same logic as for nsplnt<0
c            otherwise if kspt=-1, convert cylindrical to spherical
c            otherwise if kspt=-4, convert cartesian to spherical
c            otherwise if kspt=-5 spot is moving and its position
c              relative to planet must be supplied in the obslib,
c              in which case, spcord(1) is set to -5 on return
c
c                  v. spot coordinates           (sptred)
c  for each spot
c card 1a
c  columns
c    1- 4  spot name                                         (1a4)
c    5- 7  planet number                                     (i3)
c    8- 8  not used                                          (blank)
c    9-24  first  spot coordinate                            (f16.9)
c            or  9-26   (f18.11) ##
c   25-40  second spot coordinate                            (f16.9)
c            or 27-44   (f18.11) ##
c   41-56  third  spot coordinate                            (f16.9)
c            or 45-62   (f18.11) ##
c   57-58  not used in old format                            (blank)
c            64 velocity flag: if '6', another card follows  (1x,a1) ##
c   59-64  not used                                          (blank)
c            or 65-70 ##
c   65-70  lsp=1,0  adjust or not the spot coordinates       (3i2)
c            or 71-76   (3i2) ##
c   71-72  kspt  type of coordinates                         (i2)
c            or 77-78   (i2) ##
c
c if the second spot card is omitted, then the velocity is assumed to
c be zero, and the components are not to be adjusted
c ## if columns 79-80 are '##", use the alternate (new) format
c
c card 1b - implemented only in new format
c note: velocities for pulsars are implemented separately, not here
c    1- 4  not used, but must be non-blank
c    5- 8  not used
c    9-26  upward spot velocity                            (f18.11) ##
c   27-44  westward spot velocity                          (f18.11) ##
c   45-62  northward spot velocity                         (f18.11) ##
c   63-70  reference epoch                                 (f8.0)   ##
c   71-76  ls= 1  adjust each of three spot velocities     (3i2)    ##
c          ls= 0  do not adjust 
c   79-80  should be '##', but will be ignored
c
c shared work space in input link
      common/WRKCOM/ Spcort,Spc,T0sp,T0,Spott,Lsp,Lspcrt,Nspltt
      real*10 Spcort(6,u_mxspt),T0sp(u_mxspt)
      real*10 Spc(6),T0
      character*4 Spott(u_mxspt)
      integer*2 Lsp(6),npl
      integer*2 Lspcrt(6,u_mxspt),Nspltt(u_mxspt)

c external function
      real*10 DOT

c local variables
      real*10 ss(4),cc(4),erad,eflat,latr,qq,q2,q3,q4,rc,rs
      integer*4 nplcop,ks
      character*4 blank/'    '/,sp
      integer   i, j, l, n, ns
      character*80 card
      real*10 jd2000/2451545._10/
      character*1 vflg
 
c standard spots (all moon)
      character*4 spt0(14)/'SUR1','SUR3','SUR5','SUR6','SUR7','APL2',
     .          'AP11', 'AP14', 'AP15', 'AL12', 'AL14', 'AL15', 'AL16',
     .          'AL17'/
      real*10 sptc0(3,14)/ 1735.4740_10,316.6760_10,  -2.50200_10,
     1                    1736.1060_10,336.682510_10,-3.0550_10,
     2                    1735.1140_10, 23.2170_10,   1.4060_10,
     3                    1736.4390_10,358.632290_10, 0.4590_10,
     4                    1739.0510_10,348.5630_10, -40.9750_10,
     5                    1735.639_10,  23.4602_10,   0.6707_10,
     6                    1735.501_10,  23.41198_10,  0.69231_10,
     7                    1736.359_10, -17.53876_10, -3.62446_10,
     8                    1735.5087_10,  3.56822_10, 26.15404_10,
     9                    1736.012_10, -23.48504_10, -2.99051_10,
     A                    1736.363_10, -17.53798_10, -3.62409_10,
     B                    1735.5087_10,  3.56965_10, 26.15478_10,
     C                    1737.447_10,  15.43625_10, -8.95541_10,
     D                    1734.814_10,  30.70843_10, 20.20966_10/
c
c initialize spot data
      do j = 1, u_mxspt
         Spot(j)   = blank
         Spott(j)  = blank
         Nsplnt(j) = 0
         Nspltt(j) = 0
         T0spot(j) = 0.0_10
         T0sp(j)   = 0.0_10
         do i = 1, 6
            Spcord(i,j) = 0.0_10
            Spcort(i,j) = 0.0_10
            Lspcrd(i,j) = 0
            Lspcrt(i,j) = 0
         end do
      end do
c
c set up standard spot names and planets
      l = 0
      Numspt = 14
      do j = 1, Numspt
         Spott(j)  = spt0(j)
         Nspltt(j) = 10
         do i = 1, 3
            Spcort(i,j) = sptc0(i,j)
         end do
      end do
 
      if(.not. init) then
c
c spool spot cards from in to in0 with a-format printout
         call PEPTIC(In,Iout,in0,5,'SPOT OR STAR CARDS  ', nstop, 1)
         do while( .true. )
c
c read spot data
            read(in0,10) card
   10       format(a80)
            if(card(79:80).eq."##") then
               read(card,15) sp,npl,(Spc(i),i=1,3),vflg,
     .          (Lsp(i),i=1,3),ks
   15          format(1A4,i3,1x,3f18.11,1x,a1,6x,4I2)
            else
               read(card,20) sp,npl,(Spc(i),i=1,3),(Lsp(i),i=1,3),ks
   20          format(1A4,i3,1x,3F16.9,8x,4I2)
               vflg=blank
            endif
            if(sp.eq.blank) goto 100
            if(vflg.eq.'6') then
               read(in0,25) (Spc(i),i=4,6),T0,(Lsp(i),i=4,6)
   25          format(8x,3f18.11,f8.0,3i2)
            else
               T0=jd2000
               do i=4,6
                  Spc(i)=0.0_10
                  Lsp(i)=0
               end do
            endif
c
c see if these are non-spherical coordinates on earth
            if(npl.eq.3 .and. ks.lt.0) then
               npl=ks
               ks=0
            endif
            if(npl.lt.0) then
               if(npl.ge.-3) then
c
c           npl= -1   cylindrical coordinates on earth
c           npl= -2   geodetic (spheroidal) coordinates on earth with
c                     radius and flattening as input to program as whole
c           npl= -3 same as   -2  except that radius and flattening
c                     to be specially read in
c
c           convert cylindrical coordinates on earth to spherical coord.
                  if(npl.lt.-1) then
c
c convert geodetic coordinates on earth to spherical coord.
c spc(1)=height above sea level (meters)
c spc(3)=geodtic or spheroidal latitude (deg)
                     erad  = Econd(7)
                     eflat = Econd(8)
                     if(npl.lt.-2) then
                        read(in0,30) erad,eflat
   30                   format(2F20.10)
                     endif
c ss(j),j=1,4= geodetic flattening sine coefficients
c cc(j),j=1,4= geodetic flattening cosine coefficients
                     qq     = -0.5_10*eflat
                     ss(1)  = 1.0_10 +
     .                qq*(3._10 - 0.125_10*eflat*(5._10 + 1.5_10*
     .                eflat))
                     cc(1)  = 1.0_10 -
     .                qq*(1._10 + 0.125_10*eflat*(5._10 + 3.5_10*
     .                eflat))
                     ss(2)  = qq*(1._10 - eflat*(1._10 +
     .                0.15625_10*eflat))
                     cc(2)  = qq*(1._10 + eflat*(1._10 +
     .                0.84375_10*eflat))
                     qq     = qq*eflat
                     ss(4)  = qq*eflat*0.15625_10
                     cc(4)  = ss(4)
                     qq     = -qq*0.375_10
                     ss(3)  = qq*(1.0_10 - 0.5_10*eflat)
                     cc(3)  = qq*(1.0_10 + 1.5_10*eflat)
                     latr   = Convd*Spc(3)
                     q2     = COS(2.0_10*latr)
                     q3     = COS(4.0_10*latr)
                     q4     = COS(6.0_10*latr)
                     rc     = (Spc(1)
     .                *1.E-3_10 + erad*(cc(1) + cc(2)*q2 + cc(3)
     .                *q3 + cc(4)*q4))*COS(latr)
                     rs     = (Spc(1)
     .                *1.E-3_10 + erad*(ss(1) + ss(2)*q2 + ss(3)
     .                *q3 + ss(4)*q4))*SIN(latr)
                     Spc(1) = SQRT(rc**2 + rs**2)
                     Spc(3) = ATAN2(rs,rc)/Convd
                  else
                     Spc(1) = SQRT(Spc(1)**2 + Spc(3)**2)
                     latr   = Spc(3)
                     Spc(3) = 0.0_10
                     if(Spc(1).gt.0) Spc(3) = ASIN(latr/Spc(1))
     .                   /Convd
                  endif
c
c earth coodinates converted to spherical
                  npl = 3
               endif
            endif

c convert cylindrical or cartesian on other planets to spherical
            if(npl.ne.3 .and. ks.lt.0) then
               if(ks.eq.-4) then
                  qq=SQRT(DOT(Spc,Spc))
                  latr=ASIN(Spc(3)/qq)/Convd
                  Spc(2)=ATAN2(Spc(2),Spc(1))/Convd
                  Spc(1)=qq
                  Spc(3)=latr
                  ks=0
               else if(ks.eq.-1) then
                  qq=SQRT(Spc(1)**2+Spc(3)**2)
                  Spc(3)=ASIN(Spc(3)/qq)/Convd
                  Spc(1)=qq
                  ks=0
               else if(ks.eq.-5) then
                  Spc(1)=-5._10
                  ks=0
                  if(Lsp(1).gt.0 .or. Lsp(2).gt.0 .or. Lsp(3).gt.0 .or.
     .             Lsp(4).gt.0 .or. Lsp(5).gt.0 .or. Lsp(6).gt.0) then
                     write(Iout,34) Lsp,sp
   34                format(' ***LSPOT=',6I2,' NOT ALLOWED FOR ',a4,
     .                ', CANNOT ADJUST COORDINATES OF MOVING S/C',
     .                t72,'**')
                     nstop=nstop+1
                  endif
               endif
            endif
            if(ks.ne.0) then
               write(Iout,35) ks
   35          format(' ***UNEXPECTED SPOT CONVERSION SELECTOR',i2,
     .          t72,'**')
               nstop=nstop+1
            endif
c
c search to see if this spot is already in the list
            do i = 1, Numspt
               if(sp.eq.Spott(i)) then
                  ns = i
                  goto 40
               endif
            end do
            Numspt = Numspt + 1
            if(Numspt.gt.u_mxspt)
     .           call SUICID('TOO MANY INPUT SPOTS, STOP IN SPTRED', 9)
            Spott(Numspt)  = sp
            Nspltt(Numspt) = npl
            ns   = Numspt
   40       do i = 1, 6
               Spcort(i,ns) = Spc(i)
               Lspcrt(i,ns) = Lsp(i)
            end do
            T0sp(ns)=T0
         end do
      endif
c
c rearrange data to labelled common so that corrections to a given
c planet are grouped together
  100 if(Numspt.gt.0) then
c
c find spots on the earth
         call SPTSRT(l,3)
c
c find spots on the moon
         call SPTSRT(l,10)
c
c find spots on the sun
         call SPTSRT(l,0)
c
c find spots on planets
         do n = 1,u_mxpl
            nplcop = Nplnt(n)
            if(nplcop.gt.0) call SPTSRT(l,nplcop)
         end do
c
c spots on non-input planets
         do j = 1,30
            call SPTSRT(l,j)
         end do
c
c stars
         call SPTSRT(l,-4)
 
         if(l.ne.Numspt) call SUICID(
     .       ' ERROR IN NUMBER OF SPOTS COUNTED, STOP IN SPTRED   ',
     .       13)
c
c
      endif
      rewind in0
      return
      end
 
 
      subroutine SPTSRT(l,npl)
 
      implicit none
c
c array dimensions
      include 'globdefs.inc'
 
      include 'sptcrd.inc'
 
c shared work space in input link
      common/WRKCOM/ Spcort,Spc,T0sp,T0,Spott,Lsp,Lspcrt,Nspltt
      real*10 Spcort(6,u_mxspt),T0sp(u_mxspt)
      real*10 Spc(6),T0
      character*4 Spott(u_mxspt)
      integer*2 Lsp(6)
      integer*2 Lspcrt(6,u_mxspt),Nspltt(u_mxspt)
c
c
c   purpose:     to move spot data from its temporary input order to its
c                final sorted order
c
c   input:       l      = postion of spot data in the final array
c                         (integer)
c                npl    = body identification number to select
c                         (integer*4)
c
 
      integer     i, l, jj, npl
 
      do i = 1, Numspt
         if(Nspltt(i).eq.npl) then
            l = l + 1
            Nsplnt(l) = Nspltt(i)
            Nspltt(i) = -1
            Spot(l)   = Spott(i)
            T0spot(l) = T0sp(i)
            do jj = 1, 6
               Spcord(jj,l) = Spcort(jj,i)
               Lspcrd(jj,l) = Lspcrt(jj,i)
            end do
         endif
      end do
 
      return
      end
