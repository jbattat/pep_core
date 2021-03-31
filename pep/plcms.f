      subroutine PLCMS
 
      implicit none
 
c        subroutine PLCMS - J.F.Chandler - 1977 December 1
c        Set up for planet center of mass correction to system center
c        or for solar-system barycenter offset with "star" observations.
c        Called from CMPAR3.  Subsequently, PLCMC is called from PLTRP
c        or from SOTRP.
c JCT(55) is a packed bits indicator for the planet center of mass
c         computations
c JCT(55)=0  no center of mass corrections
c    1 bit = 1: get satellite coordinates from interpolation if possible
c    2 bit = 1: get satellite coordinates from elliptic if necessary
c    4 bit = 1: use only satellites that are on s-body tape
c    8 bit = 1: for elliptic method, obtain elements from s-body tape
c                (if possible).  Warning: this option should normally
c                be set because the initial conditions in /EMPCND/ are
c                likely to have been overwritten from tape anyway and
c                may not match the epoch in CON1(1).
c   16 bit = 1: allow calculation of barycenter of planet+observed
c                satellite instead of planet
c   32 bit = 1: use precessing elliptic orbit instead of stationary

c For solar-system barycenter offsets, always take JCT(55) to be 3.

c array dimensions
      include 'globdefs.inc'

c commons
      include 'b2dta.inc'
      include 'bddta.inc'
      include 'cmcke.inc'
      include 'comdat.inc'
      real*10 Cmfct
      equivalence (Comcon(128),Cmfct)
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      integer*4 I2bod
      equivalence (I2bod,Jpert)
      include 'namtim.inc'
      include 'number.inc'
      include 'param.inc'
      include 'plndta.inc'
      include 'stats.inc'

      include 'maxpcmdt.inc'

c temporary storage
      common/YVECT/ Tcond(6,20),Con2(10),Plcstf(10),Acms(20),Fie(20)
      character*8 Acms
      character*4 Fie
      real*10 Tcond,Con2,Plcstf

c local
      real*10 cmass,dum8,goose,gscl,smass1
      integer*4 i,iplr,iprt,j,jct55,jct6,k,ncent
      integer*2 kpl,np,i2
      character*4 lfi/'(I) '/,lfe/'(E) '/,lfp/'(P) '/
      character*8 pname(9)/'MERCURY ',' VENUS  ',' EMBARY ','  MARS  ',
     .     'JUPITER ',' SATURN ',' URANUS ','NEPTUNE ',' PLUTO  '/
 
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c Get total mass of solar system for the purposes of calculating
c barycenter position.  Include the conglomerate of perturbing
c asteroids from the external file, if any.  Note: if the
c observations are of a "star", also include any individual
c asteroids specified for this run.

c In principle, should also include all three asteroid belts (which
c are centered on the Sun's position), but skip them for now.
      Mascnt = 0._10
c      Mascnt = Prmter(50)
c      do i = 45,46
c         Mascnt = Mascnt + Prmter(i)
c      end do
      do i = 1, 9
         Mascnt = Mascnt + Mass(i)
      end do
      if(Kpert.gt.0) Mascnt = Mascnt + Sumcom
      Mascnt = Mascnt + 1._10

      if((Klan.le.0 .or. Klan.gt.u_mxpl) .and. Nplnt0.ne.-4) return
      jct55 = Jct(55)
      if(Nplnt0.eq.-4) jct55 = 3
      if(mod(jct55,4).eq.0) return
      if(Klam.ne.Klan) then
         iprt   = 1
         Numpcm = 0
         jct6   = Jct(6)/4
         do i = 1, 4
            Plcdbg(i) = (mod(jct6,2).eq.1)
            jct6 = jct6/2
         end do
         do i = 1, 8
            Jflg(i) = (mod(jct55,2).eq.1)
            jct55   = jct55/2
         end do
         if(.not. Jflg(2)) Jflg(3) = .true.
         if(Jflg(1) .and. Jflg(3)) Jflg(2) = .false.
         if(.not. Jflg(2)) Jflg(6) = .false.
         if(Jflg(1) .or. Jflg(6) .or.  .not. Jflg(2)) Jflg(4)
     .       = .false.

         if(Nplnt0.eq.-4) then
c always include all 9 planets for solar-system barycenter
c (masses are already included in Mascnt)
            ncent = 0
            cmass = 0._10
            do np=1,9
               Numpcm = Numpcm+1
               Kplcm(Numpcm) = -np
               Acms(Numpcm) = pname(np)
               Pcintx(Numpcm) = Intb(np)*Ibdsgn
               Tplc(Numpcm) = 0._10
               Fie(Numpcm) = lfi
            end do
         else
            ncent = Nplnt(Klan)
            cmass = 1._10
         endif
 
c check list of input planets for concentric satellites
         do i = 1, Numpln
            if(Npcent(i).ne.ncent .and. 
     .       (Npcent(i).gt.0 .or. Nplnt0.ne.-4)) goto 50
            np = Nplnt(i)
            if(np.le.10 .or. np.gt.30) goto 50

c find it on s-body tape
            if(I2bod.gt.0 .and. Nast.gt.0) then
               do j = 1, Nast
                  if(Np2(j).eq.np) goto 10
               end do
            endif

c not found.  reject it if so desired
            if(Jflg(3)) goto 50
            j = 0

   10       Numpcm = Numpcm + 1
            if(Numpcm.gt.maxpcm) call SUICID(
     .       'TOO MANY SATELLITES, STOP IN PLCMS  ', 9)
            Kplcm(Numpcm) = i
            Acms(Numpcm) = Aplnt(i)
            cmass = cmass - Mass(np)
            Tplc(Numpcm) = 0._10
            if(j.gt.0) then
               if(Inb2(j).gt.0) Pcintx(i) = Inb2(j)*Ib2sgn
               if(Inb2(j).le.0) Pcintx(i) = (2._10**Inb2(j))*Ib2sgn
               do k = 1, 6
                  Tcond(k,Numpcm) = B2ta(k,j)
               end do
               Tplc(Numpcm) = Tb20(j)
            endif
            if(Tplc(Numpcm).eq.0._10 .or. .not. Jflg(4)) then
               do k = 1, 6
                  Tcond(k,Numpcm) = Pcond(k,i)
               end do
            endif
            if(j.gt.0 .and. Jflg(1)) then
               Fie(Numpcm) = lfi
               Tplc(Numpcm) = 0._10
            else if(Jflg(6)) then
               Fie(Numpcm) = lfp
            else
               Fie(Numpcm) = lfe
            endif

   50    end do
         if(Numpcm.le.0) return
         if(Nplnt0.eq.-4) then
            goose = Gauss
            Mascnt = Mascnt - cmass
            cmass = 1._10
         else
            goose = Gauss*SQRT(Mass(ncent))
         endif
         iplr  = -3
 
c do setup for each body
         do i = 1, Numpcm
            kpl     = Kplcm(i)
            if(kpl.lt.0) then
               np  = -kpl
            else
               np  = Nplnt(kpl)
            endif
            Tlcm(i) = 1E10_10
            smass1  = cmass + Mass(np)
            if(Fie(i).eq.lfi) goto 100
            if(Tplc(i).ne.0._10) goto 80

c using input elliptic elements
            do while( iplr.lt.kpl )
               read(Iplcon)
               iplr = iplr + 1
            end do
            read(Iplcon) Tplc(i),(dum8,k=1,11),(j,k=1,6),(i2,k=1,100),
     .        j,j,i2,j,j,i2,i2,dum8,(i2,k=1,100),i2,Con2
            iplr = iplr + 1
c put initial time into jd+fract form
            Tplc(i) = Tplc(i) + 0.5_10
   80       if(Jflg(6)) then
               call PNITL(Icnd(kpl),Pcond(1,kpl),Con2,Elpcm(1,i))
            else
               gscl = goose*SQRT(smass1)
               call JNITL(gscl,Tcond(1,i),Elpcm(1,i),0,0)
            endif
  100    end do
         if(iplr.ne.-3) rewind Iplcon
      endif
      if(Nplnt0.ne.-4) write(Iout,200) Jflg
  200 format(
     .   '0  PLANET COORDINATES SHIFTED FROM SYSTEM BARYCENTER,  JFLG='
     .   , 8L2)
      if(iprt.eq.1) write(Iout,300) Numpcm,
     .                        (Fie(i),Acms(i),i = 1,Numpcm)
  300 format(i10,
     .' SECONDARIES USED IN CENTER-OF-MASS CORRECTION (I=INTERPOLATE,E=
     .ELLIPTIC, P=PRECESSING)'/8(4x,a4,a8))
      iprt  = 0
      Cmfct = 1._10
      if(Klanb.gt.0) then
         if(Jflg(5)) then
            do i = 1, Numpcm
               if(Klanb.eq.Kplcm(i)) then
                  Cmfct = 1._10 - Mass(Nplnt0)
                  write(Iout,310) Aplnt(Klanb),Cmfct
  310             format('  OBSERVED BODY ', a8,
     .                   ' INCLUDED WITH PLANET,  CMFCT=', f18.16)
                  return
               endif
            end do
         endif
      endif
      return
      end
