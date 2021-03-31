      subroutine NRMSAV(scale,hedskp,nsize,iaprio,iptr)
 
      implicit none

c     m.e.ash   nov 1969    subroutine nrmsav
c     save normal equations onto disk if so indicated
c     scale the coefficient matrix and right side of the normal eqs
c        option with nsize.ne.nparam used only for call from
c        pprctl, not used for saving ne
c iaprio - if >0, must combine imbedded a priori info from imat0(.).
c          otherwise, just copy from ibuf2
c iptr - array of pointers for restoring a priori info
      logical*4 hedskp
      real*10 scale(1000)
      integer*4 nsize,iaprio
      integer*2 iptr(1000),icv(12)
 
c scale(i)=sqrt(diagonal element i)   (output)
c
c array dimensions
      include 'globdefs.inc'

c common
      include 'aprtbf.inc'
      include 'cureph.inc'
      include 'fcntrl.inc'
      include 'filtim.inc'
      include 'inodta.inc'
      include 'mcnfrm.inc'
      include 'nrmmat.inc'
      include 'restor.inc'
      include 'rtside.inc'
      include 'zeroes.inc'
 
c local variables
      real*10 dumsq
      integer*4 i,icnst,irow,ivectk,japrio,lts,ncns,nsav0
      character*8 ttot/'NRMEQTOT'/,teph/'NRMEQEPH'/,tmod/'NRMEQMOD'/,
     1 temd/'NRMEQEMD'/,ttitle
      character*1 lchar,skip/'0'/,skp2/'-'/
      integer*4 savep
      real*4    one/1.0E00/
      integer*4 nfvctk/ - 1/
c
c*  start=1000
c write out coefficient matrix of normal equations
      if(Fict(6).ne.1) then
         call NRMRIT(' UNSCALED SYMMETRIC COEFFICIENT ',B,nsize,
     .               Measmt,Ict(47).ge.0)
c
c*  start=1500
c write out right hand side of normal equations
         lchar = skp2
         lts   = 4 + (nsize + 7)/8
         if(Line + lts.ge.56) then
            call NEWPG
            Line  = 0
            lchar = skip
         endif
         write(Iout,50) lchar,nsize,(Side(i),i = 1,nsize)
   50    format(a1,'BEFORE SCALING THE',i4,
     .          ' RIGHT HAND SIDES OF THE NORMAL EQUATIONS ARE'//(4x,
     .          1p,8D16.8))
         Line = Line + lts
      endif
c
c*  start=2000
c see if the normal equations are to be saved
      savep = 0
      ncns  = 0
      icnst = 0
      if(Ibuf5.gt.0 .and. Itrwnd(Ibuf5).ne.0) then
 
c ibuf5 not rewound -> doing multiple runs
         icnst = 1
         if(Itrwnd(Ibuf5).gt.0) icnst = 2
 
c apply constraints before saving normal equations
         call CNSTRN(Ibuf5,ncns,nsize,B,Side,Sigma)
      endif
      if(Imat2.gt.0) then
         if(.not. (Jct(51).gt.0 .and. hedskp)) then
c
c save filtered normal eqs. in non-filtered format
            if(Ict(42).eq.0) then
c
c*  start=2300
c save the normal equations on disk
               ttitle = ttot
               if(ncns.gt.0) ttitle = tmod
               goto 200
 
            else if(Fict(7).ne.0) then
 
               ttitle = teph
               if(ncns.gt.0) ttitle = temd
               savep = Fict(7)
               if(savep.lt.0) savep = Nepoch + savep + 1
               if(savep.eq.Ephnum) goto 200
            endif
         endif
      endif
c
c
c*  start=2200
      if(Fict(6).ne.1) then
         write(Iout,100)
  100    format('0NORMAL EQUATIONS NOT SAVED')
         Line = Line + 2
      endif
      goto 600
  200 call SAVHED(Imat2,ttitle)
      ivectk = 1
      if(Jct(60).gt.0) ivectk = 0
      japrio = iaprio
      if(Ict(44).ne.0) japrio = 1
      if(Jct(60).gt.0) japrio = 0
      write(Imat2) Measmt,izr4,izr4,one,one,Ermeas,Sumzns,
     .             Sumzsm,japrio,ivectk,Sumaps,(Zero(1),i = 1,9)
      if(ivectk.gt.0) write(Imat2) nfvctk,
     .                          (Vectk(i),i = 1,Nparam)
      write(Imat2) Nparam,(Side(i),i = 1,Nparam)
      call FWSIG(Imat2,Nparam,B,Sigma)
      if(japrio.ne.0) then
 
c combine all imbedded a priori information and save
         if(iaprio.gt.0) then
 
c restore and combine a priori from all sources
            do i = 1,12
               icv(i) = 0
            end do
            call NRMSET(1)
 
c restore saved normal equations from previous runs
            icv(1) = 0
            icv(2) = -1
            if(Ict(76).ge.10) icv(6) = 128
            icv(7) = 1
            if(Jct(60).lt.0) icv(7) = icv(7) + 2
            call NRMFRM(icv,iptr,B)
 
c add a priori normal equations
            if(Ict(44).ne.0) then
               icv(1) = -1
               icv(4) = Ibuf2
               icv(6) = 0
               icv(7) = 1
               if(Jct(60).lt.0) icv(7) = icv(7) + 2
               call NRMFRM(icv,iptr,B)
            endif
 
c write combined a priori info to imat2
            write(Imat2) Nparam,(Side(i),i = 1,Nparam)
            call FWSIG(Imat2,Nparam,B,Sigma)
         else
 
c just copy from ibuf2
            call NEWHED(Ibuf2,0,0,iptr,0,dumsq)
            call ZFILL(Sigma,16*Nparam)
 
c read and expand right side of a priori information
            call QREAD(Ibuf2,Nsav,Sav,Mparam)
            do i = 1,Mparam
               Sigma(iptr(i)) = Sav(i)
            end do
            write(Imat2) Nparam,(Sigma(i),i = 1,Nparam)
 
c read and expand a priori coefficient matrix
            Nsav = 0
            do while( .true. )
               nsav0 = Nsav
               call QREAD(Ibuf2,Nsav,Sav,Mparam)
               if(Nsav.le.nsav0) call SUICID(
     .          ' SAVED ROW COUNTERS DO NOT MATCH, STOP IN NRMSAV',12)
               irow = iptr(Nsav)
               do i = 1,Mparam
                  Sigma(iptr(i)) = Sav(i)
               end do
               write(Imat2) irow,(Sigma(i),i = 1,Nparam)
               if(Nsav.ge.Mparam) then
                  if(irow.lt.Nparam) write(Imat2) Nparam,
     .               (Zero(1),i = 1,Nparam)
                  rewind Ibuf2
                  Itrwnd(Ibuf2) = 0
                  goto 300
               endif
            end do
         endif
      endif
  300 write(Imat2) izr4,(Zero(1),i = 1,18)
 
c end file imat2
      rewind Imat2
      write(Iout,400) Nparam,Measmt,Imat2
  400 format('0NORMAL EQUATIONS OF ORDER',i5,' FORMED FROM',i9,
     .       ' MEASUREMENTS SAVED ON DATA SET',i3)
      if(savep.ne.0) write(Iout,500) savep
  500 format(' FILTERED NORMAL EQUATIONS FROM EPOCH ',i5)
      Line = Line + 3
      if(iaprio.gt.0 .and. japrio.gt.0) then
         call NRMSET(0)
         icv(1) = 1
         icv(5) = 1
         icv(6) = 0
         icv(7) = 0
         if(Jct(60).lt.0) icv(7) = icv(7) + 2
         call NRMFRM(icv,iptr,B)
      endif
c
c*  start=3000
c calculate scale factors for the normal equations
  600 if(icnst.eq.1) call CNSTRN(Ibuf5,ncns,nsize,B,Side,Sigma)
      call NRMSCL(B,scale,nsize,Side,nsize,Fict(6).eq.0)
 
      return
      end
