      subroutine PPRCTL
 
      implicit none

c        program to control generation of partially pre-reduced
c        normal equations
c        p. macneil   april, 1977
c
c array dimensions
      include 'globdefs.inc'
c
c         common
      include 'aprtbf.inc'
      include 'ciptrx.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'iptrst.inc'
      include 'mcnfrm.inc'
      include 'nrmmat.inc'
      include 'numnum.inc'
      include 'restor.inc'
      include 'rtside.inc'
      include 'scail.inc'
      integer*2 iptr(1000)
      equivalence (Scale,iptr)
      integer*2 iptrc(1000),iptrd(1000)
      equivalence (Scale(251),iptrc),(Scale(501),iptrd)
 
c common used for scratch
      common/WRKCOM/ Pvrow(2000),Pvcol(2000),Pvrwb(2),Tmp(2000),
     .        Jptrc(1000),Jptrd(1000),Iuse(1000)
      real*10 Pvrow,Pvcol,Pvrwb,Tmp
      integer*2 Iuse,Jptrc,Jptrd
      include 'zeroes.inc'
c
c miscellaneous quantities (dimensions set for order 1000)
      real*10 dinvw(1000)
      real*10 fdinvw(1000),vbar(1000),colfdi(1000),f(1000),c(1000)
      real*10 cbar(1000),u2ap(1000),dap(1000),utwid(1000),
     .          dvars(1000)
      equivalence (Tmp,fdinvw,cbar)
      equivalence (Pvrow(1001),c,u2ap,dap)
      equivalence (Pvrow,colfdi,f,vbar),
     .            (Pvcol,dinvw),(Jptrc,dvars)
      equivalence (Sigma,utwid)
      real*10   qb(1001),qtmp(1000),qdinvw(1000),qfdiw(1000),
     .          qclfdi(1000),qsum
      equivalence (qb,B),(qtmp,Tmp),
     .            (qclfdi,colfdi),(qdinvw,dinvw),
     .            (qfdiw,fdinvw)
 
      real*10 DOTN,dum,znsqap
      integer*4 i,i0ctw,ia,iaprio,ic,id,ierinv,ii,im,ippr,
     . irow,ivectk,j,j1,j2,ji,jj,jk,k,kz,
     . length,lts,m,mats,n,ndimc,neg,nsav0,nused,nwstrt
      real*4    one/1.0E0/
      integer*2 ipartn(1000)
      equivalence (ipartn,Iuse)
      integer*2 ityp,nip
      character*8 words(2)/' IS ',' (CONT) '/
      character*4 nc/'C   '/,nd/'D   '/,nf/'F   '/,dq/'DQ  '/
      character*44 phdr/'THE REDUCED COEF.MATRIX OF ORDER            '/
      integer*2 icv(12)
      integer*4 nfvctk/-1/
c
c storage arrays used in pprctl.  each is labeled with its use at the
c successive operations: "<" = new assignment, "+" = modification,
c "-" = use, "." = retention for later use, "=" = recovery of values,
c "x" = scratch.
c
c   step:    s c f d i s p f  f  di ua da h v f  a e
c            e       n c t di di fa       d   di p n
c            t       v l r    w           r   fa r d
c   storage
c
c    ibuf6     < . . . . . .  .  .  .  .  . . -
c    ibuf7       < . . . . -  .  .  .  .  . . -
c    ibuf8             < . -
c    ibuf9   < . . . . . - <  .  -
c    ipartn  <
c    icrest  < . . . . . - .  .  .  .  .  -
c    idrest  < . . . . - - .  .  .  .  .  -
c    side      < . + - . - .  .  .  .  .  . -
c    vectk     < . + . . - +  .  .  .  .  -
c    b       x < < < + +   <  -  <  +  +  . . -  - -
c  * sav       x x x   x x                       x x
c
c                        the following groups of arrays are overlaid.
c
c    jptrc   <                                              1
c    jptrd   <                                              2
c    iuse            x                                      3
c    dvars             < . .  .  .  .  .  . . .  -          1-4
c
c    iptr    <                                              1
c    iptrc   < - -       =                                  2
c    iptrd   < . - -     = .  .  .  -  -                    3
c    scale           < -                                    1-4
c
c    sigma                                                  1
c    utwid                          <  +  . . .  -          1
c
c  * dinvw             < - .  -  .  -  -  . . .  -          1-2
c  * pvcol           x                                      1-2
c
c  * fdinvw                   <  .  .  .  . -               1-2
c  * tmp             < -   x                                1-2
c    cbar                                     <             1
c
c  * colfdi                x                                1-2
c    pvrow1          x              x  x  x                 1
c    vbar                                   <               1
c    f                                        <             1
c    pvrow2          x                    x                 2
c    u2ap                           <                       2
c    dap                               <                    2
c    c                                        <             2
c
c  * double-length arrays, may be used for extended precision.
c
c                                     --- step: set ---
c* start=1000
      call PAGSET('PARTIAL PREREDUCTION OF NORMAL EQUATIONS',10)
      call NRMSET(0)
c
c write out identification group for ppr normal equations
      if(Imat3.eq.0) call SUICID('IMAT3 = 0, STOP IN PPRCTL   ',7)
      call SAVHED(Imat3,'PPRNRMEQ')
c
c get iptr
      call ZFILL(icv,2*12)
      if(Jct(53).lt.2) icv(1) = 2
      icv(2) = -1
      icv(5) = 2
      icv(7) = 2
      call NRMFRM(icv,iptr,B)
      mats   = icv(4)
      iaprio = mod(icv(7)/4,2)
      ivectk = 1
      if(Jct(60).gt.0) then
         iaprio = 0
         ivectk = 0
      endif
c
c test mats returned by nrmfrm
      if(mats.le.0) call SUICID('NO SNE DATA SET, STOP PPRCTL',7)
 
      write(Imat3) Mesmt1,izr4,izr4,one,one,Ermes1,Znsqsn,Zsmsn,
     .             iaprio,ivectk,Apssq1,(Zero(i),i=1,9)
c
c        nomenclature (taken from memo of r. d. reasenberg dated
c        29 july 1975):
c                  (v)   (c  f) (y)
c        rhs = u = ( ) = (    ) ( ) = bx,
c                  (w)   (fa d) (z)
c             where fa = fadjoint
c        solutions are given by:
c                  -1
c        y = (cbar)  *(vbar), where
c                              -1                         -1
c             cbar = c - (f)*(d  )*(fa), vbar = v - (f)*(d  )*(w), and
c
c        z = (zbar) - (fbar-adjoint)*(y), where
c                      -1                   -1
c             zbar = (d  )*(w), and fba = (d  )*(fa)
c
c        form iptrc and iptrd
c        ipartn(i) = 0 means parameter number  i  goes to the c-matrix
c        ipartn(i) = 1 means parameter number  i  goes to the d-matrix
c
c           mappings:
c  iptr(1-mparam) --> nparam
c  iptrc(1-mparam) --> ncparm
c  iptrd(1-mparam) --> ndparm
c  icrest(1-ncparm) --> nparam
c  idrest(1-ndparm) --> nparam
c  jptrc(1-nparam) --> ncparm
c  jptrd(1-nparam) --> ndparm
c
c        initialize iptrc, iptrd, icrest, idrest
      call ZFILL(iptrc,2*2000)
      call ZFILL(Icrest,2*2000)
      call ZFILL(Jptrc,2*2000)
 
      if(Ibuf5.gt.0 .and. Itrwnd(Ibuf5).ne.0) then
         if(Itrwnd(Ibuf5).gt.0) rewind Ibuf5
         read(Ibuf5,end=100) ityp,nip,(ipartn(i),i = 1,nip)
         Itrwnd(Ibuf5) = 2
         if(ityp.eq.3 .and. nip.eq.Nparam) go to 200
      endif
  100 call SUICID('NO INPUT PPR REQUESTS,STOP IN PPRCTL',9)
c
c find number of c- and d-parameters
  200 ic = 0
      id = 0
      do i = 1,Nparam
         if(ipartn(i).gt.0) then
 
c d-type
            id = id + 1
            Jptrd(i)   = id
            Idrest(id) = i
         else
 
c c-type
            ic = ic + 1
            Jptrc(i)   = ic
            Icrest(ic) = i
         endif
      end do
      Ncparm = ic
      Ndparm = id
 
      if(Nparam.ne.Mparam) call SUICID('NPARAM.NE.MPARAM IN PPRCTL  ',
     .                                -7)
      if(Ncparm.eq.0) call SUICID('NO C PARAMETERS, STOP IN PPRCTL ',8)
      if(Ndparm.eq.0) call SUICID('NO D PARAMETERS, STOP IN PPRCTL ',8)
c
c check whether room for extended precision
      i0ctw = 2*Ncparm*Ndparm
      if(i0ctw.gt.Nrmsiz) call SUICID('F MATRIX TOO LARGE FOR EXTENDED-P
     .RECISION PRE-REDUCTION, STOP IN PPRCTL ',18)
      if(iaprio.gt.0 .and. i0ctw + (Ncparm*(Ncparm+1))/2.gt.Nrmsiz) call
     .   SUICID('MATRICES TOO LARGE WITH A PRIORI, STOP IN PPRCTL',12)
 
      do im = 1,Mparam
         i = iptr(im)
         if(i.gt.0) then
            iptrc(im) = Jptrc(i)
            iptrd(im) = Jptrd(i)
         endif
      end do
c
c save pointers temporarily
      write(Ibuf9) (iptrc(i),i = 1,Ncparm),(iptrd(i),i = 1,Ndparm)
      endfile Ibuf9
      rewind Ibuf9
c
c        note that u (the c-part of rhs) can be restored using:
c               do 3000 i=1,ncparm
c               ia=iptr(icrest(i))
c               if(ia.le.0) go to 3000
c               side(ia)=side(ia)+sav(i)*weight
c          3000 continue
c
c                                     --- step: c ---
c        start  of c formation
c
c        initialize b, side for c (can't call nrmset because
c        ncparm .ne. nparam .ne. ndparm)
      length = ((Ncparm+1)*Ncparm)/2
      if(length.gt.Nrmsiz)
     .     call SUICID('C MATRIX TOO LARGE, STOP IN PPRCTL  ',9)
 
      call ZFILL(Side,16*999)
      call ZFILL(Vectk,16*1000)
      call ZFILL(B,16*length)
c
c skip header records
      call FRMHED(mats,'IMATS','    ',0,ippr,1)
      if(ippr.ne.0) call SUICID(
     .           'MUST NOT APPLY PPR TO PPRNRMEQ, STOP IN PPRCTL  ',12)
c
c form c-matrix
      call PPRFRM(B,iptrc,iptrc,Ncparm,Ncparm,mats,Mparam,nc,
     .            Side,Vectk)
      if(Jct(60).lt.0) then
         call NRMRIT('   INTERESTING COEFFICIENT (C)  ',B,Ncparm,
     .               Measmt,.true.)
         call GPMPO(Ncparm,Side,0,0,'RHS: (V)',8)
         call GPMPO(Ncparm,Vectk,0,0,'MEAN SENSITIVITY',16)
      endif
c
c finished forming c and rhs
      rewind mats
      Itrwnd(mats) = 0
c
c write c out on temporary storage
      call ARRWRT(Ncparm,Ncparm,nc)
c
c        end of c-matrix processing
c
c
c                                     --- step: f ---
c        start of f formation
c
c        initialize b for f
      length = Ncparm*Ndparm
      if(length.gt.Nrmsiz)
     .     call SUICID('F MATRIX TOO LARGE, STOP IN PPRCTL  ',9)
      call ZFILL(B,16*length)
c
c skip header records
      call FRMHED(mats,'IMATS','    ',0,ippr,1)
c
c form f-matrix
      call PPRFRM(B,iptrc,iptrd,Ncparm,Ndparm,mats,Mparam,nf,
     .            dum,dum)
      ndimc = Ncparm
      if(Ncparm.eq.1) ndimc = 0
      if(Jct(60).lt.0) call GPMPO(Ndparm,B,ndimc,0,'F MATRIX',8)
c
c finished forming f
      rewind mats
      Itrwnd(mats) = 0
c
c write out f on temporary storage
      call ARRWRT(Ncparm,Ndparm,nf)
c
c        end of f-matrix processing
c
c
c                                     --- step: d ---
c        start of d formation
c
c        initialize b for d
      length = (Ndparm + 1)*Ndparm/2
      length = length + length
      if(length.gt.Nrmsiz)
     .     call SUICID('D MATRIX TOO LARGE, STOP IN PPRCTL  ',9)
      call ZFILL(B,16*length)
c
c skip header records
      call FRMHED(mats,'IMATS','    ',0,ippr,1)
c
c form d-matrix
      call PPRFRM(B,iptrd,iptrd,Ncparm,Ndparm,mats,Mparam,nd,
     .            Side(Ncparm+1),Vectk(Ncparm+1))
      if(Jct(60).lt.0) then
         call NRMRIT(' UNINTERESTING COEFFICIENT (D)  ',B,Ndparm,
     .               Measmt,.true.)
         call GPMPO(Ndparm,Side(Ncparm+1),0,0,'RHS: (W)',8)
         call GPMPO(Ndparm,Vectk(Ncparm+1),0,0,
     .              'MEAN SENSITIVITY',16)
      endif
c
c finished forming d
      rewind mats
      Itrwnd(mats) = 0
 
      call TIMRIT(' FORMING C,D,F MATRICES ',6)
c
c                                     --- step: inv ---
c        scale d matrix and w (right hand side) and invert
c* start=1500
c
c           turn off underflow interrupt during inversion?
      if(Jct(15).ne.0) then
         call NOUNDR
         call PAGCHK(60,2,0)
         write(Iout,250)
  250    format(
     . '0*** UNDERFLOW INTERRUPT DISABLED DURING D-MATRIX INVERSION,',
     . 'ANY UNDERFLOWS WILL NOT BE DETECTED ***')
      endif
c
c extended precision inversion
      do i = 1,Ndparm
         qtmp(i) = Side(Ncparm + i)
      end do
      length = (Ndparm*(Ndparm+1))/2
      do i = 1,length
         qb(length + 1 - i) = B(length + 1 - i)
      end do
      call NRMSCQ(qb,Scale,Ndparm,qtmp,Ndparm,.true.)
      call SYMINQ(qb,qtmp,Ndparm,1,Pvrow,Pvcol,Iuse,Pvrwb,
     .            ierinv)
      if(Jct(15).ne.0) call UNDRON
c --- step: scl ---
c rescale d-inverse
      nused = 0
      neg   = 0
      do i = 1,Ndparm
         ii = nused + i
         if(qb(ii).lt.0.0) then
            neg = neg + 1
            call PAGCHK(60,1,0)
            write(Iout,260) i,Idrest(i)
  260       format(' DIAGONAL ELEMENT',i5,' OF D INVERSE (',i5,
     .             ' OF WHOLE SET) IS NEGATIVE')
         endif
         do j = 1,i
            nused     = nused + 1
            qb(nused) = (qb(nused)*Scale(i))*Scale(j)
         end do
         dvars(i)  = qb(nused)
         qdinvw(i) = qtmp(i)*Scale(i)
      end do
      if(neg.gt.0) call SUICID(
     .       ' NEGATIVE DIAGONAL ELEMENTS IN D-INVERSE, STOP IN PPRCTL'
     .       ,14)
c
c write d-inverse out on temporary storage
      call ARRWRT(Ndparm,Ndparm,dq)
c
c
c        end of d-matrix processing
c
c                                     --- step: ptr ---
c* start=2000
c           recover pointers (clobbers scale vector no longer needed)
      read(Ibuf9) (iptrc(i),i = 1,Ncparm),(iptrd(i),i = 1,Ndparm)
      rewind Ibuf9
 
      nwstrt = Ncparm + 1
      do i = 1,Ndparm
         Sav(i) = qdinvw(i)
      end do
      Znsqpp = DOTN(Side(nwstrt),Sav,Ndparm)
      Zsmpp  = DOTN(Vectk(nwstrt),Sav,Ndparm)
      call PAGCHK(60,5 + (Ncparm-1)/30,0)
      write(Iout,300) Znsqpp,Zsmpp,Ncparm,
     .                 (Icrest(i),i = 1,Ncparm)
  300 format('0REDUCTION OF NORMAL EQUATIONS SUBTRACTS',1pd15.8,
     .       ' FROM EFFECTIVE SUM (O-C/ERROR)**2'/37x,'AND',1pd15.8,
     .       ' FROM EFFECTIVE SUM (O-C/ERROR)'/' THE',i4,
     .       ' INTERESTING PARAMATERS ARE:'/(1x,30I4))
      call PAGCHK(60,2 + (Ndparm-1)/30,0)
      write(Iout,400) Ndparm,(Idrest(i),i = 1,Ndparm)
  400 format(' THE',i4,' UNINTERESTING PARAMETERS ARE:'/(1x,30I4))
c
c --- step: fdi ---
c read f into b
      j2 = 0
      do i = 1,Ncparm
         j1 = j2 + 1
         j2 = j2 + Ndparm
         read(Ibuf7) (B(j),j = j1,j2)
      end do
      rewind Ibuf7
c
c form and write(f*dinverse)adjoint = fbaradjoint
c use symmetry of dinverse
      do i = 1,Ndparm
 
c read a row (column) of dinverse into sav
         read(Ibuf8) (qtmp(j),j = 1,Ndparm)
 
c call prodct(b,sav,colfdi,-ncparm,ndparm,1)
         jk = 0
         do k = 1,Ncparm
            qclfdi(k) = 0.0
            do j = 1,Ndparm
               jk = jk + 1
               qclfdi(k) = qclfdi(k) + qtmp(j)*B(jk)
            end do
         end do
 
c write row of (f*dinverse)adjoint = column of f*dinv
         write(Ibuf9) (qclfdi(k),k = 1,Ncparm)
 
c form reduced k-vector
         if(ivectk.ne.0) then
            do k = 1,Ncparm
               Vectk(k) = Vectk(k) - Vectk(Ncparm + i)*qclfdi(k)
            end do
         endif
      end do
      rewind Ibuf8
      endfile Ibuf9
      rewind Ibuf9
c
c --- step: fdiw ---
c* start=2500
c form f*dinverse*w
      ji = 0
      do i = 1,Ncparm
         qfdiw(i) = 0.0
         do j = 1,Ndparm
            ji = ji + 1
            qfdiw(i) = qfdiw(i) + qdinvw(j)*B(ji)
         end do
      end do
c
c --- step: difa ---
c read dinv*fa (= fbaradjoint) into array
      j2 = 0
      do i = 1,Ndparm
         j1 = j2 + 1
         j2 = j2 + Ncparm
         read(Ibuf9) (qb(j),j = j1,j2)
      end do
      rewind Ibuf9
c
c --- step: ua ---
c compute reduced a priori information
      if(iaprio.ne.0) then
 
c clear ctwid
         call ZFILL(B(i0ctw+1),16*(Ncparm*(Ncparm+1))/2)
 
c skip headers and matrix on sne
         call FRMHED(mats,'IMATS','    ',0,ippr,0)
         read(mats)
         call BSKIP(mats,Mparam)
 
c get a priori right-hand side
         read(mats) Nsav,(Pvrow(k),k = 1,Mparam)
         if(Nsav.ne.Mparam)
     .        call SUICID('MISSING A PRIORIS, STOP PPRCTL  ',8)
         call ZFILL(u2ap,16*Ndparm)
         do k = 1,Mparam
            ia = iptrd(k)
            if(ia.gt.0) u2ap(ia) = Pvrow(k)
         end do
         do j = 1,Ncparm
            jk   = j
            qsum = 0.0
            do k = 1,Ndparm
               qsum = qsum + u2ap(k)*qb(jk)
               jk   = jk + Ncparm
            end do
            utwid(j) = -qsum
         end do
 
c correction to normsq
         qsum = 0.0
         do k = 1,Ndparm
            qsum = qsum + qdinvw(k)*u2ap(k)
         end do
         znsqap = qsum + qsum
c --- step: da ---
c read a priori information
         Nsav = 0
         do while( .true. )
            nsav0 = Nsav
            if(Nsav.ge.Mparam) then
 
c correct norm-sq offset
               Znsqpp = Znsqpp - znsqap
               go to 500
            else
               read(mats) Nsav,(Pvrow(k),k = 1,Mparam)
               if(Nsav.le.nsav0) call SUICID(
     .             'INVALID A PRIORI ROW COUNTERS, STOP PPRCTL  ',11)
               irow = iptrd(Nsav)
               if(irow.gt.0) then
                  call ZFILL(dap,16*Ndparm)
                  qsum = 0.0
                  do k = 1,Mparam
                     ia = iptrd(k)
                     if(ia.gt.0) then
                        dap(ia) = Pvrow(k)
                        qsum    = qsum + dap(ia)*qdinvw(ia)
                     endif
                  end do
                  znsqap = znsqap - qsum*qdinvw(irow)
                  ia     = i0ctw
                  do j = 1,Ncparm
                     jk   = j
                     qsum = 0.0
                     do k = 1,Ndparm
                        qsum = qsum + dap(k)*qb(jk)
                        jk   = jk + Ncparm
                     end do
                     utwid(j) = utwid(j) + qsum*qdinvw(irow)
                     jk = (irow - 1)*Ncparm
                     do jj = 1,j
                        ia    = ia + 1
                        jk    = jk + 1
                        B(ia) = B(ia) + qsum*qb(jk)
                     end do
                  end do
               endif
            endif
         end do
      endif
c
c --- step: hdr ---
c* start=3000
c read nominals and scales into b
  500 m = 2*Nparam
      read(Ibuf1)
      read(Ibuf1) (Pvrow(i),i = 1,m)
      rewind Ibuf1
 
c write out pointer group
      write(Imat3) Nparam,Ncparm,Ndparm,Znsqpp,izr4,izr4,
     .             (Icrest(i),i = 1,Ncparm),
     .             (Idrest(i),i = 1,Ndparm),(Pvrow(i),i = 1,m),
     .             Zsmpp
c
c write kbar
      if(ivectk.gt.0) write(Imat3) nfvctk,(Vectk(k),k = 1,Ncparm)
c
c --- step: v ---
c change v to vbar
      do i = 1,Ncparm
         vbar(i) = Side(i) - qfdiw(i)
      end do
 
c write vbar
      write(Imat3) Ncparm,(vbar(i),i = 1,Ncparm)
c
c print right hand side
      if(Ict(47).ge.0) then
         lts = 4 + (Ncparm - 1)/8
         call PAGCHK(60,lts,0)
         write(Iout,550) Ncparm,(vbar(i),i = 1,Ncparm)
  550    format('0THE',i4,
     .          ' REDUCED RIGHT HAND SIDES ARE'/(4x,1p,8D16.8))
         if(Line.gt.50) Line = 60
         call EBCDIX(Ncparm,phdr,33,4)
         call MVC(words(1),1,8,phdr(37:44),1)
         call PAGSET(phdr,-11)
         call PAGHED(0)
         call MVC(words(2),1,8,phdr(37:44),1)
         call PAGSET(phdr,-11)
      endif
c
c                                     --- step: fdifa ---
c        form f*dinv*fadjoint
c
c        now read a row of f, form a row of fdifa, read a row of c,
c        subtract, write a row of cbar; loop
      do i = 1,Ncparm
         read(Ibuf7) (f(j),j = 1,Ndparm)
         read(Ibuf6) (c(j),j = 1,Ncparm)
         kz = Ncparm - i
 
c call prodct(b,f,fdifa,ncparm,ndparm,1)
         do j = 1,Ncparm
            jk   = j
            qsum = 0.0
            do k = 1,Ndparm
               qsum = qsum + f(k)*qb(jk)
               jk   = jk + Ncparm
            end do
            cbar(j) = c(j) - qsum
            if(cbar(j).ne.0._10) kz = 0
         end do
         if(kz.eq.0) write(Imat3) i,(cbar(j),j = 1,Ncparm)
c
c print row of cbar
         if(Ict(47).ge.0) then
            n = i
            if(kz.ne.0) n = 1
            lts = 2 + (n - 1)/8
            call PAGCHK(60,lts,1)
            write(Iout,560) i,(cbar(j),j = 1,n)
  560       format('0',i3,(t5,1p,8D16.8))
         endif
 
      end do
      rewind Ibuf6
      rewind Ibuf7
c
c* start=3500
c --- step: apr ---
c write out reduced a priori information
      if(iaprio.ne.0) then
 
c write out ctwid
         write(Imat3) Ncparm,(utwid(i),i = 1,Ncparm)
         call FWSIG(Imat3,Ncparm,B(i0ctw+1),Sav)
         if(Ict(47).ge.0) then
            lts = 4 + (Ncparm - 1)/8
            call PAGCHK(60,lts,0)
            write(Iout,550) Ncparm,(utwid(i),i = 1,Ncparm)
            call NRMRIT('REDUCED A PRIORI COEFFICIENT    ',B(i0ctw+1),
     .                  Ncparm,izr4,.true.)
         endif
      endif
c
c write d variances
      if(Jct(60).le.0) write(Imat3) nfvctk,(dvars(i),i = 1,Ndparm)
c
c write zbar (=dinvw)
      do i = 1,Ndparm
         Sav(i) = qdinvw(i)
      end do
      write(Imat3) Ndparm,(Sav(i),i = 1,Ndparm)
c
c --- step: end ---
c copy fbaradjoint to output
      j2 = 0
      do i = 1,Ndparm
         do j = 1,Ncparm
            Sav(j) = qb(j2 + j)
         end do
         j2 = j2 + Ncparm
         write(Imat3) i,(Sav(j),j = 1,Ncparm)
      end do
      rewind Imat3
      Itrwnd(Imat3) = 0
      call TIMRIT('  PARTIAL REDUCTION OF NORMAL EQUATIONS ',10)
      Pprdon = .true.
      return
      end
