      subroutine EMPRD1(lice)
 
      implicit none
 
c r.w.king  sept 1976     subroutine emprd1
c first five records of earth-moon barycenter tape are read
c
      integer*2 lice
c lice =0 printout of data on first two records of embary tape
c lice =1 no such printout
c
c array dimensions
      include 'globdefs.inc'
c common
      include 'bddtaint.inc'
      include 'emmips.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'intstf.inc'
      include 'metuna.inc'
      include 'morstf.inc'
      include 'namtim.inc'
      include 'param.inc'
      include 'plndta.inc'
      include 'tapdtplp.inc'
      integer*2 kc(100),nkic,kic(99),k2cent(199)
      equivalence (Centbd(1,1,1),k2cent(1),kc(1)),(k2cent(101),kic(1))
      include 'xprcom.inc'
      include 'yvectrd1.inc'
      include 'yvectplp.inc'
 
c local variables
      integer i,j,k,king,kong,kt
      integer*2 klac
      real*10 frem1(2)
      integer nlist,klist(i_mxplp),iicntl(i_mxplprt)
      logical gotmnics
 
c subroutine PLPNDX is generic and does not know what body is being
c integrated, and so the control flags -31 to -36, indicating partials
c w.r.t. the integrated body's ICs must be converted to a generic
c notation as NPLNT*100 plus an offset of 1 to 6.
      do j=1,i_mxplprt
         iicntl(j)=Icmtrl(j)
         if(Icmtrl(j).lt.-30 .and. Icmtrl(j).ge.-36)
     .    iicntl(j)=Mplnt*100-30-Icmtrl(j)
      end do
 
c set up array NPLPT with central and all target body numbers together
      Nplpc = 3
      do kt=1,Numtar
         Nplpt(kt)= Ntrg(kt)
      end do

c loop over central body and any target bodies
      do kt=0,Numtar
c
c initialize data set indicator
         Jplntg(kt) = 0
         gotmnics=.false.
         nkic=-1
         Mppt(kt)=0
         if(Nplpt(kt).eq.10) then
            Nptspr(kt)=8
         else
            Nptspr(kt)=5
         endif
c initialize pointers to required quantities
         do j=1,i_mxplprt
            Kpt(j,kt)=0
         end do

c if we are merely forbidden to read a tape, then we must still index the
c needed partials in order to compute the approximations
         if(Kkm(84).le.0) goto 80

c see if central or target body is input planet
         do i = 1, Numpln
            if(Nplnt(i).eq.Nplpt(kt)) then
               klac  = i
               goto 50
            endif
         end do
         if(Nplpt(kt).ne.3) goto 80
         klac=-3

   50    Jplntg(kt) = Iplnt(klac)
         if(kt.eq.1) Xprnam(2)=Aplnt(klac)
 
c if planet tape not supplied, ignore tape-reading
c however planet tapes are always expected to be available in this routine
         if(Jplntg(kt).eq.0) goto 80
         if(Itrwnd(Jplntg(kt)).ne.0) rewind Jplntg(kt)
 
c save input jd0 - if lt.0, there has been a checkpoint
c restart and input jd2 overrides jd2 on tape
         Jdxx9 = Jdpl0(klac)
         M1    = klac + 3
         do k = 1, M1
            read(Iplcon)
         end do
 
c read planet constants from disk
         read(Iplcon) (Cn1x(j),j=25,36),Beps,Kkxx,Jdd1,Jdd2
         rewind Iplcon
 
c read first two records of embary or planet peripheral data set
         call XXRD1(lice,Nplpt(kt),Jplntg(kt),klac,Jdt1(kt),
     .    Jdt2(kt),Ipart(kt),i_mxplprt+1,Inttx(kt),Idirt(kt),kc,
     .    Tint(kt),frem1,nkic,kic)
c see if any moon ic partials are found
         do i=8,nkic
            if(kic(i).ge.1001 .and. kic(i).le.1006) gotmnics=.true.
         end do
c
c calculate interval quantities
         if(Idirt(kt).lt.0) Tint(kt) = -Tint(kt)
         Intt5(kt) = 5*Inttx(kt)
c
c look for partials needed for integration
c partials of central body coordinates w.r.t. integrated body
c also, if no tape is available and no elliptic approximations can be made,
c we skip searching for needed partials
   80    continue
c
c get list of partials to interpolate
         nlist=0
         do j=1,i_mxplprt
            if(iicntl(j).ne.Icmtrl(j) .and. gotmnics) then
c partials w.r.t. integrated body IC
               nlist=nlist+1
               if(nlist.gt.i_mxplp) goto 95
               klist(nlist)=iicntl(j)
               Kpt(j,kt)=nlist+7
            else if((iicntl(j).ge.1.and.iicntl(j).le.10) .or.
     .       (iicntl(j).ge.31 .and. iicntl(j).le.33) .or.
     .       (iicntl(j).ge.41.and.iicntl(j).le.44)) then
               nlist=nlist+1
               if(nlist.gt.i_mxplp) goto 95
               klist(nlist)=iicntl(j)
               Kpt(j,kt)=nlist+7
            else
               king=(iicntl(j)-1)/100
               kong=iicntl(j)-100*king
               if(king.gt.0 .and. kong.le.6) then
                  if(king.eq.Nplpt(kt)) then
                     Kpt(j,kt)=kong+1
                  else if(king.le.9) then
                     nlist=nlist+1
                     if(nlist.gt.i_mxplp) goto 95
                     klist(nlist)=iicntl(j)
                     Kpt(j,kt)=nlist+7
                  endif
               endif
            endif
         end do
         if(nlist+7.le.i_mxplp) goto 100
   95    call SUICID('TOO MANY INDIRECT PARTIALS, STOP IN EMPRD1  ',11)
c
c set up pointers for requested partials
  100    call PLPNDX(i_mxplprt,iicntl,nlist,klist,Nplpt(kt),Nqt(kt),
     .    Kqt(1,kt),Krt(1,kt),Kpt(1,kt),nkic,kic)

c initialize interpolated quantities to zero
         do j=1,Nqt(kt)
            do i=1,6
               Ytp(i,j,kt)=0._10
            end do
         end do
c
c read first three data records of earth-moon or planet
c data set
         if(Jplntg(kt).gt.0) then
            call RDPTAP(0,kt)
c tell PRTCRD not to interpolate this planet
            if(Nplpt(kt).le.9) Kiss(Nplpt(kt))=-1
c
c setup for em barycenter initial conditions and embary
c contribution to gmvary partials
         else if(kt.eq.0) then
c embary initial conditions obtained from elliptic formulas
            call SUICID(
     .'ELLIPTIC ORBIT USED FOR EMBARY- GMVARY PARTIAL WILL OMIT THE EMBA
     .RY ORBIT CONTRIBT''N', -21)
            call IMITL(Gauss, Mass(3), Betabd(1,3), 1, 0)
            T0mpt(0) = Jdbd0(3)
c set gmvary and any other ss parm partials to zero
            do j=7,nlist+6
               do i=1,6
                  Dccor(i,j)=0._10
               end do
            end do
         else
c elliptic approximation for other planets would go here
         endif
      end do
c take note of 1st partial if found on center or 1st target integration
      Parnum(3)=Kpt(1,0)-1
      if(Numtar.gt.0) Parnum(2)=Kpt(1,1)-1
      return

c rewind planet tape(s)
      entry EMPRWD
      do kt=0,Numtar
         if(Jplntg(kt).gt.0) then
            rewind Jplntg(kt)
            Itrwnd(Jplntg(kt)) = 0
         endif
      end do
      return
      end
