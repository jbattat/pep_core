      subroutine DECORR
 
      implicit none

c j.f.chandler  aug 1986   subroutine decorr
c program for decorrelation of partial derivatives

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'bernum.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'iptrst.inc'
      include 'ktrap.inc'
      include 'namobs.inc'
      include 'nrmmat.inc'
      include 'numnum.inc'
      include 'obsdta.inc'
      include 'restor.inc'
      include 'rtsidesl.inc'
      include 'scail.inc'
      integer*2 jptr(1000),jptrc(1000),jptrd(1000)
      equivalence (Scale,jptr),(Scale(251),jptrc),
     .            (Scale(501),jptrd)
      include 'wrkcompr.inc'
      character*8 scnam(2),sitf8(2),pname
      equivalence (scnam,Nmobst),(sitf8,Sitf),(pname,Plnnam)
c
c miscellaneous
      real*10 ws0,ws1
      integer*4 i,icnt(2),ip,j,k,kc,kcrow,kd,kn,nfirst,nprlv,nvmobt
      character*8 dcrnam/' DECORR '/
c
c iabs1=input obs-theory and partial derivitive tapes
c iabs2=output
c
      call PAGSET('PARTIAL-DERIVATIVE DECORRELATION',8)
      Line = 60
      if(Ict(1).lt.1) call ANSET
c
c restore projection matrix
      call DCRFRM(Imat3,B,Nrmsiz*(Nrmsiz-1)/2,jptr,jptrc,jptrd)
      nprlv = Ncparm + Ndparm
      call PAGCHK(60,2,0)
      write(Iout,100) nprlv,Ncparm,Ndparm
  100 format('0DECORRELATION PERFORMED WITH',i4,
     .       ' RELEVANT PARAMETERS:',i4,' "INTERESTING" AND',i4,
     .       ' "UNINTERESTING"')
      if(Ict(47).ge.0) call GPMPO(Ndparm,B,Ncparm,1,
     .                      'DE-CORRELATION PROJECTION MATRIX',32)
c
c initialization
      call PAGSET('PARTIAL-DERIVATIVE DECORRELATION',8)
      call PAGSET(
     .'NTAP NSEQ TYPOBS PLANET NPLNT  SITE1 SERIES  SITE2  SPOT     ERRO
     .R WEIGHTS         MEASMT   ',-23)
      Niobc  = -1
      nvmobt = Numobt
      if(Ict(16).gt.0) Numobt = Ict(16)
      Jtape = 0
  200 do while(.true.)
c
c increment tape counter
         call PRDOBS(2,1)
         if(Iabs1.le.0) then
 
            call TIMRIT('  DE-CORRELATION OF PARTIAL DERIVATIVES ',10)
 
            Numobt = nvmobt
            return
         else if(Iabs2.gt.0) then
            if(Mout.gt.0) write(Mout,220) Iabs2
  220       format(' DECORRELATE TO IABS2=',i3)
            if(Line.gt.40) call NEWPG
            call PAGHED(0)
c
c write first two records of output observation library tape
            Tapnam = dcrnam
            call PRDOBS(-2,1)
            goto 300
         else
            rewind Iabs1
            Itrwnd(Iabs1) = 0
         endif
         end do
  300 do while(.true.)
c
c read and write first record of observation series
         call PRDOBS(2,3)
         call PRDOBS(-2,3)
         if(Ncodf.le.0) goto 200
         nfirst  = 0
         icnt(1) = 0
         icnt(2) = 0
c
c get correct observation type name
         if(Ncodg.eq.4) then
            if(Klanb.gt.0 .and. (Ncp0.eq.3 .or. Ncp0.eq.10))
     .          Ncodg = 20
            if(Ksite(1).gt.0 .and. Ksite(1).ne.3) Ncodg = 21
         endif
         do while(.true.)
c
c read observation record
            call PRDOBS(2,4)
            if(Ncode.gt.0) then
c
c get partial derivative pointers:
c iptr(1-numpar) --> nparam
               call PRDEQS(nfirst)
c
c calculate new partials and residuals
               if(Nice.le.0) icnt(1) = icnt(1) + 1
               if(Nice.ge.0) icnt(2) = icnt(2) + 1
               do j = Num1,Num2
                  ws0 = 0._10
                  do k = 3,Numpar
                     ws1    = Deriv(k,j)
                     Sav(k) = 0._10
                     kn     = Iptr(k)
                     if(kn.gt.0) then
                        kd = jptrd(kn)
                        if(kd.le.0) then
                           kc = jptrc(kn)
                           if(kc.gt.0) then
                              kcrow = (kc - 1)*Ndparm
                              do i = 3,Numpar
                                 ip = Iptr(i)
                                 if(ip.gt.0) then
                                    ip = jptrd(ip)
                                    if(ip.gt.0) ws1 = ws1 -
     .                                Deriv(i,j)*B(kcrow + ip)
                                 endif
                              end do
                           endif
                        else
                           ws0 = ws0 + Deriv(k,j)*Solut(kd)
                           goto 310
                        endif
                     endif
                     Sav(k) = ws1
  310             end do
 
c replace old partials and residuals with new
                  Deriv(2,j) = Deriv(2,j) - ws0
                  do k = 3,Numpar
                     Deriv(k,j) = Sav(k)
                  end do
               end do
            endif
c
c write new observed minus theory tape
            call PRDOBS(-2,4)
            if(Ncode.le.0) then
c
c printout description of observation series
               call PAGCHK(59,1,1)
               write(Iout,320) Ntape,Nseq,Typobs(Ncodg),Ncodf,
     .                          pname,Nplnt0,sitf8(1),Series,
     .                          sitf8(2),Spota,(Erwgt2(i),i = 1,2),
     .                          icnt
  320          format(i4,i5,a6,i2,1x,a8,i3,2(1x,a8,1x,a4),1p,
     .                2E12.5,2I5)
               if(Ncodf.gt.20) then
                  call PAGCHK(60,1,1)
                  write(Iout,330) Nplnt2,Spota2
  330             format(26x,i3,24x,a4)
               endif
               goto 400
            endif
         end do
  400 end do
      end
