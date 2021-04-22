      subroutine PLANET(kfit)
 
      implicit none

c     m.e.ash  march 1968  subroutine planet
c     main program for n-body, earth-moon barycenter, planet-asteroid,
c     and satellite-probe numerical integrations
c     kfit=0 integrations are in midst of orbit fitting iterations
c     kfit=1 integrations are after orbit fitting convergence, extend
c            integration interval without partials, change output tape
c     updated feb. 1980   kcl: tcon added
c
c array dimensions
      include 'globdefs.inc'
c commons
      include 'aacoff.inc'
      include 'aprtbf.inc'
      include 'bdctrl.inc'
      include 'empcnd.inc'
      real*10 meqinc
      equivalence (Mrcond(11),meqinc)
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'namtimq.inc'
      include 'param.inc'
      include 'petina.inc'
      include 'petuna.inc'
      character*8 pname
      equivalence (Name, pname)
      include 'plndta.inc'
c note that jdpl0, pcond, con11 and iplnt are extended backwards
c to include 20 items (-3, -2, -1, 0 for em, mn, er, mr).
c however, nplnt, npcent, and icnd are not.
c
      include 'stbtst.inc'
      include 'stint.inc'
      include 'timstf.inc'

c local variables
      integer*4 i,iperts,irltsk(2),j,kfit,kkp91,lic,lica,lice,
     . licf,llib,lnut,lrel,mpage,ngo,npints
      real*10 relsav
      character*8
     1    embary/' EMBARY '/,bodyn/' N-BODY '/, oplan/' PLANET '/
      character*20 typmes/'++++++++ INTEGRATION'/
      character*8 tname
      equivalence (typmes, tname)
      character*88 tmes/'  ALL*** PLANET NUMERICAL INTEGRATIONS DURING T
     .HE LAST***** PAGES STARTING ON PAGE***** '/
      integer*2 mplnt
 
c in call to prtrd1 in n-body numerical integration
      integer*2 kb(10), izr2/0/
c
c save times and page
      irltsk(1) = Ireal0
      irltsk(2) = Itotsk
      mpage     = Npage
      npints    = 0
      lice   = 0
      lica   = 0
      Nnotst = 0
c
c write planet link title page
      tname = oplan
      call PAGSET(typmes, 5)
      call NEWPG
      write(Iout, 100)
      write(Iout, 200)
  100 format('-'/'-'/'-'/'-'//23x,
     .      '**********    **             *********    **           **'
     .      , '   ***********   ************'/23x,
     .      '***********   **            ***********   ***          **'
     .      , '   ***********   ************'/23x,
     .      '**       **   **            **       **   ****         **'
     .      , '   **                 **     '/23x,
     .      '**       **   **            **       **   ** **        **'
     .      , '   **                 **     '/23x,
     .      '**       **   **            **       **   **  **       **'
     .      , '   **                 **     '/23x,
     .      '**       **   **            **       **   **   **      **'
     .      , '   **                 **     '/23x,
     .      '***********   **            ***********   **    **     **'
     .      , '   ********           **     ')
  200 format(23x,
     .      '**********    **            ***********   **     **    **'
     .      , '   ********           **     '/23x,
     .      '**            **            **       **   **      **   **'
     .      , '   **                 **     '/23x,
     .      '**            **            **       **   **       **  **'
     .      , '   **                 **     '/23x,
     .      '**            **            **       **   **        ** **'
     .      , '   **                 **     '/23x,
     .      '**            **            **       **   **         ****'
     .      , '   **                 **     '/23x,
     .      '**            ***********   **       **   **          ***'
     .      , '   ***********        **     '/23x,
     .      '**            ***********   **       **   **           **'
     .      , '   ***********        **     ')
      call PLINK

c reference epoch for coordinate systems
      Kepoch= 1
      if(Jct(13).gt.0) Kepoch=2
c
c read con11 for all bodies
      lic = Numpln + 4
      do j = 1, lic
         read(Iplcon) (Con11(i,j), i = 1, 12)
      end do
      rewind Iplcon
c
c n-body numerical integration
      if(Nbody.gt.0 .and. Jdbdy0.ne.0 .and. Jdbdy0.gt.-3000000) then
         do i=1,9
            kb(i)=-1
         end do
c if we are not integrating the moon, then it must be obtained elsewhere
         kb(10)=1
         lnut  =0
         do i = 1, Nbody
c if we are integrating the moon, we probably need nutation
            if(Nplbdy(i).eq.10) then
               kb(10)=-1
               lnut  =1
ccccc special insert to test moon code in nbody environment with
ccccc interpolated planets instead of simultaneously integrated
               if(Ipert.gt.0 .and. Con11(9,2).eq.9._10)then
                  do j=1,9
                     kb(j)=1
                  end do
               endif
            endif
         end do
c this suicid could be avoided if mnprd1 were upgraded to do nutation
         if(lnut.gt.0 .and. Ipert2.gt.0) call SUICID(
     .    'IPERT2>0 WITH MOON IN N-BODY INTEGRATION, STOP IN PLANET',14)
         if(Ipert2.eq.0) then
            call PRTRD1(lice,izr2,izr2,kb,1,lnut,0)
         else
            call MNPRD1(lice)
         endif
         lice = 1
         call BODSET(lice)
         npints = npints + 1
         tname  = bodyn
         call PAGSET(typmes,5)
         call OPRMSG(typmes,5)
         if(Kbdy(28).lt.2) then
            call NUMINT(4)
         else if(Kbdy(28).eq.2) then
            call ADAM(4)
         else
            call RROAD(4)
         endif
         if(lice.gt.0) then
            if(Ipert2.gt.0) then
               rewind Ipert2
               Itrwnd(Ipert2) = 0
            else if(Ipert.gt.0) then
               rewind Ipert
               Itrwnd(Ipert) = 0
            endif
         endif
 
c end file ibody
         if(Kbdy(39).ge.0) rewind Ibody
      endif
      if(Nbody.gt.0 .and. Ict(9).gt.0) lice = 0
      licf = -3
c
c earth-moon barycenter numerical integration
      Nplnt = 3
      j = -3
      do while( .true. )
 
c skip over unused records on iplcon
         if(j.gt.licf) then
            read(Iplcon)
            licf = licf + 1
            goto 300
         endif
 
c read body information
         read(Iplcon) Con1, Epsp, Kp, Jdp1, Jdp2, Intp, Intp1,
     .    Intp2, Ihrp, Iminp, Secp, Kkp,
     .    mplnt, Tcon, Numki, (Ki(i), i = 1, Numki)
         licf = licf + 1

         if(Kp(88).eq.-8) Jdpl0(j)=0
 
c store tcon and kkp(51) for stepsize control
         do i = 1, 30
            Dtcon(i) = Tcon(i)
         end do
         Istp = Kkp(51)
 
         if(Nbody.gt.0) then
            do i = 1, Nbody
               if(Nplbdy(i).eq.Nplnt) goto 240
            end do
         endif
         if(Jdpl0(j).eq.0 .or. Jdpl0(j).le.-3000000) goto 240

c set up for integration
         do i = 1, 30
            Cond(i) = Pcond(i, j)
         end do
         Jdp0  = Jdpl0(j)
         Jplnt = Iplnt(j)
         pname = Aplnt(j)
         Ncentr= Npcent(j)
         Icnd  = Jcnd(j)
         Klan  = j
 
         tname = pname
         call PAGSET(typmes, 5)
         call OPRMSG(typmes, 5)
         lrel = 0
         lnut = 0
         llib = 0
         if(Kp(61).ge.0 .and. Kp(61).le.1) then
            if(Ncentr.lt.0 .or.
     .         (Nplnt.eq.3 .and. Ict(40).gt.0) .or.
     .         (Nplnt.gt.30 .and. Ncentr.le.0)) lrel = 1
         endif
         if(Nplnt.eq.3 .and. Kp(82).ge.0) lnut = 1
         if(Ncentr.gt.0) then
            if(Kp(61).ge.0) lrel = -1
            if(Nplnt.gt.30 .and. Kp(81).ge.0) lrel = -1
            if(Ncentr.eq.3) lnut  = 1
            if(Ncentr.eq.10) llib = 1
         endif
         do i = 1, i_mxtrg
            if(Kp(i).gt.0) then
               if(Kp(i).eq.3) lnut  = 1
               if(Kp(i).eq.10) llib = 1
            endif
         end do
         if(llib.gt.0) then
            if(meqinc.eq.0._10) meqinc = 0.0268587_10
         endif
         iperts = Ipert
         if(Nbody.gt.0. and. Ict(9).gt.0) Ipert = Ibody
         call PRTRD1(lice,Nplnt,Ncentr,Kp(31),lrel,lnut,llib)
         call ASTRD1(lica, Nplnt, Kp(31))
         lice   = 1
         lica   = 1
         relsav = prmter(31)
         if(Kp(97).gt.0) prmter(31) = Con(24)
c
c kkp(96)=   number of days from epoch in first time direction
c         for reintegration of equations of motion after orbit
c         convergence
c kkp(97)= 0 no partials integrated in reintegration of equations of
c         motion after orbit fit convergence
c kkp(97)= 1 partials are integrated in reintegration of equations of
c         motion after orbit fit convergence
c kkp(98)=   number of days from epoch in second time direction
c         for reintegration of equations of motion after orbit fit
c         convergence
c kkp(99)=   output tape number for reintegration of equations of motion
c         after orbit fit convergence if jct(79)=1
         if(kfit.ge.1) then
            Jplnt    = Kkp(99)
            Iplnt(j) = Jplnt
            Jdp2     = Jdp0 + Kkp(98)
            Jdp1     = Jdp0 + Kkp(96)
            if(Kkp(97).le.0) then
               Numki = 8
               do i = 1, Numki
                  Ki(i) = 0
               end do
            endif
         endif
 
         if(Nplnt.eq.3 .and. Ict(40).gt.0) Ncentr = -1
         call SBSET
 
c satellite-probe integration
         ngo = 3
 
         npints = npints + 1
         if(Kp(88).lt.2) then
            call NUMINT(ngo)
         else if(Kp(88).eq.2) then
            call ADAM(ngo)
         else
            call RROAD(ngo)
         endif
         prmter(31) = relsav
         if(Ipert.ne.0) then
            rewind Ipert
            Itrwnd(Ipert) = 0
         endif
         Ipert = iperts
         if(Jpert.ne.0) then
            rewind Jpert
            Itrwnd(Jpert) = 0
         endif
         call PLPRWD
         if(Kkp(91).gt.0) then
            kkp91 = Kkp(91)
            rewind kkp91
         endif
 
c end file jplnt
         if(Kp(99).ge.0) rewind Jplnt
  240    do while( .true. )
            j = j + 1
            if(j.lt.0) j = 1
            if(j.gt.Numpln) goto 260
c
c planet numerical integration
            Nplnt = Nqlnt(j)
            if(Nplnt.lt.0) then
            else if(Nplnt.eq.0) then
               goto 260
            else
               goto 300
            endif
         end do
c
c printout time required for all numerical integrations
  260    rewind Iplcon
         call EBCDIX(npints, tmes, 6, 3)
         call EBCDIX(Npage-mpage, tmes, 55, 5)
         call EBCDIX(mpage, tmes, 83, 5)
         call TIMRTC(tmes, 22, irltsk)
         return
  300 end do
      end
