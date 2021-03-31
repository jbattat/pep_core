      subroutine PLPRD1(lice)
 
      implicit none
 
c j.f.chandler - 1976 sep
c revised 1977 feb
c revised 2009 Jun to implement target bodies and mass partials
 
c  Initialize tape-reading routine for central body of a satellite
c  integration. Do nothing unless ncentr>0.  However, if ncentr<=0,
c  initialize tape-reading for target bodies
c  
c kk(84) = 0  read no central/target body tapes during integration
c kk(84) = 1  read tape for central/target body position, velocity only
c kk(84) = 2  also read tape for partials wrt i.c., masses, and others
c             (option 1 above not implemented: treated as option 2)
c  if kk(84)=2, insist on finding partials wrt all initial conditions
c   requested for integrations.
c  input --
c     lice =0 printout of data on first two records of planet tape
c     lice =1 no such printout
 
c arguments
      integer*2 lice

c array dimensions
      include 'globdefs.inc'
c commons
      include 'bddtaint.inc'
      include 'fcntrl.inc'
      include 'intstf.inc'
      include 'namtimq.inc'
      include 'petuna.inc'
      include 'plndta.inc'
      include 'prtpin.inc'
      include 'sbstuf.inc'
      include 'tapdtplp.inc'
      integer*2 kc(100),nkic,kic(99),k2cent(199)
      equivalence (Centbd(1,1,1),k2cent(1),kc(1)),(k2cent(101),kic(1))
      include 'xprcom.inc'
      include 'yvectplp.inc'
      include 'yvectrd1.inc'
 
c local variables
      integer*2 klac
      integer*4 i,j,k,king,kong,nctl,klist(i_mxplp),kt,nl,
     . iicntl(i_mxeqn-1)
      real*10 frp1(2)
 
c subroutine PLPNDX is generic and does not know what body is being
c integrated, and so the control flags -31 to -36, indicating partials
c w.r.t. the integrated body's ICs must be converted to a generic
c notation as NPLNT*100 plus an offset of 1 to 6.
      nctl = Kount
c      if(Kkp(84).le.1) nctl = 0
      do j=1,Kount
         iicntl(j)=Icntrl(j)
         if(Icntrl(j).lt.-30 .and. Icntrl(j).ge.-36)
     .    iicntl(j)=Nplnt*100-30-Icntrl(j)
      end do

c set up array NPLPT with central and all target body numbers together
      do kt=0,Numtar
         if(kt.eq.0) then
            Nplpc = Ncentr
         else
            Nplpt(kt)= Ntrg(kt)
         endif
      end do

c loop over central body and any target bodies
      do kt=0,Numtar
c
c initialize data set indicator
         Jplntg(kt) = 0
         nkic=-1
         Mppt(kt)=0
         if(Nplpt(kt).eq.10) then
            Nptspr(kt)=8
         else
            Nptspr(kt)=5
         endif
c initialize pointers to required quantities
         do j=1,Kount
            Kpt(j,kt)=0
         end do
         nl=0

c if the selected body is the sun or a probe, then the partials don't exist
         if(Nplpt(kt).le.0 .or. Nplpt(kt).gt.30) goto 100
c if we are merely forbidden to read a tape, then we must still index the
c needed partials in order to compute the approximations
         if(Kkp(84).le.0) goto 80

c see if central or target body is input planet
         do i = -3, Numpln
            if(Nqlnt(i).eq.Nplpt(kt)) then
               klac  = i
               goto 50
            endif
         end do
         goto 80

   50    Jplntg(kt) = Iplnt(klac)
 
c if planet tape not supplied, ignore tape-reading
         if(Jplntg(kt).eq.0) goto 80
         if(Itrwnd(Jplntg(kt)).ne.0) rewind Jplntg(kt)
 
c read first two records of planet peripheral data set
         Jdxx9 = Jdpl0(klac)
         M1    = klac + 3
         do k = 1, M1
            read(Iplcon)
         end do
 
c read planet constants from disk
         read(Iplcon) (Cn1x(j),j=25,36),Beps,Kkxx,Jdd1,Jdd2
         rewind Iplcon
 
         call XXRD1(lice,Nplpt(kt),Jplntg(kt),klac,Jdt1(kt),
     .    Jdt2(kt),Ipart(kt),i_mxplprt+1,Inttx(kt),Idirt(kt),kc,
     .    Tint(kt),frp1,nkic,kic)
c
c calculate interval quantities
         if(Idirt(kt).lt.0) Tint(kt) = -Tint(kt)
         if(Nplpt(kt).eq.10) then
            Intt5(kt) = 4
            Nptspr(kt)= 8
         else
            Intt5(kt) = 5*Inttx(kt)
            Nptspr(kt)= 5
         endif
c
c look for partials needed for integration
c partials of central body coordinates w.r.t. integrated body
c or target body i.c. are not needed
c also, if no tape is available and no elliptic approximations can be made,
c we skip searching for needed partials
   80    if(Jplntg(kt).eq.0 .and. Ltrg(kt).eq.0) goto 100
         do j=1,Kount
            if(iicntl(j).ne.Icntrl(j)) then
c partials w.r.t. integrated body IC
               if(kt.gt.0 .and. Jplntg(kt).gt.0) then
                  nl=nl+1
                  if(nl.gt.i_mxplp) goto 95
                  klist(nl)=iicntl(j)
                  Kpt(j,kt)=nl+7
               endif
            else if(Icntrl(j).eq.Nplnt) then
c partial w.r.t. integrated body mass
               if(kt.gt.0 .and. Jplntg(kt).gt.0) then
                  nl=nl+1
                  if(nl.gt.i_mxplp) goto 95
                  klist(nl)=Icntrl(j)
                  Kpt(j,kt)=nl+7
               endif
            else if(Icntrl(j).le.0) then
c partials w.r.t. other integrated body parms (not needed)

c partials w.r.t. solar system parms
            else if(Icntrl(j).ge.10 .and. Icntrl(j).le.50
     .          .and. Jplntg(kt).gt.0) then
               nl=nl+1
               if(nl.gt.i_mxplp) goto 95
               klist(nl)=Icntrl(j)
               Kpt(j,kt)=nl+7
c partials of central or target body w.r.t. its own parms
            else if(kt.gt.0 .or. Ncentr.gt.0) then
               if(Icntrl(j).eq.Nplpt(kt)) then
                  nl=nl+1
                  if(nl.gt.i_mxplp) goto 95
                  klist(nl)=Icntrl(j)
                  Kpt(j,kt)=nl+7
                  Mppt(kt)=nl+7
               else
                  king=(Icntrl(j)-1)/100
                  kong=Icntrl(j)-100*king
                  if(king.eq.Nplpt(kt).and.kong.le.6) then
                     Kpt(j,kt)=kong+1
                  else if(Jplntg(kt).gt.0) then
c partials of target body w.r.t. other target or central body parms
                     do i=0,Numtar
                        if(Icntrl(j).eq.Nplpt(i) .or.
     .                   (king.gt.0.and.
     .                   king.eq.Nplpt(i).and.kong.le.6)) then
                           nl=nl+1
                           if(nl.gt.i_mxplp) goto 95
                           klist(nl)=Icntrl(j)
                           Kpt(j,kt)=nl+7
                           goto 90
                        endif
                     end do
   90                continue
                  endif
               endif
            endif
         end do
         if(nl+7.le.i_mxplp) goto 100
   95    call SUICID('TOO MANY INDIRECT PARTIALS, STOP IN PLPRD1  ',11)

  100    call PLPNDX(nctl,iicntl,nl,klist,Nplpt(kt),Nqt(kt),
     .    Kqt(1,kt),Krt(1,kt),Kpt(1,kt),nkic,kic)

c initialize interpolated quantities to zero
         do j=1,Nqt(kt)
            do i=1,6
               Ytp(i,j,kt)=0._10
            end do
         end do

c
c read first three data records of planet data set
         if(Jplntg(kt).gt.0) then
            call RDPTAP(0,kt)
c tell PRTCRD not to interpolate this planet
            if(Nplpt(kt).le.10) Kiss(Nplpt(kt))=-1
         endif
      end do

c take note of 1st partial if found on center or 1st target integration
      Parnum(3)=Kpt(1,0)-1
      if(Numtar.gt.0) Parnum(2)=Kpt(1,1)-1
      if(Nplnt.eq.3 .and. Kp(40).ge.0) then
         do kt=1,Numtar
            if(Nplpt(kt).eq.10) Parnum(3)=Kpt(1,kt)-1
         end do
      endif
      return
 
c rewind planet tape(s)
      entry PLPRWD
      do kt=0,Numtar
         if(Jplntg(kt).gt.0) then
            rewind Jplntg(kt)
            Itrwnd(Jplntg(kt)) = 0
         endif
      end do
      return
      end
