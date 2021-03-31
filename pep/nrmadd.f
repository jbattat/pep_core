      subroutine NRMADD(b,sav,iptr,nsize,imats,mesmt1,weight,
     .                  ntape1,nseq1,jmata,itdif,prmdif,znsq,side,
     .                  ermes1,sumzns,ermeas,measmt,iuseap,iaprio,
     .                  ivectk,vectk,wgtk,sumzsm,zsm,sumsqa,sumaps)
 
      implicit none
 
c r.reasenberg/d.white   jan 1974   subroutine nrmadd
c increment normal equations from saved normal equations
c read rhs and b for one series of norm eqn using pointer
c
c parameters
      real*10 b(1),sav(1),weight,prmdif(1),znsq,side(1),ermes1(3),
     . sumzns,ermeas(3),vectk(1),wgtk,sumzsm,zsm,sumsqa,sumaps
      integer*4 nsize,imats,mesmt1,ntape1,nseq1,jmata,itdif,measmt,
     . iuseap,iaprio,ivectk
      integer*2 iptr(1)
c b,side - lhs, rhs of normal equations
c sav    - work area (nsize r*8)
c iptr   - mapping ptr array from saved to current parameter set
c nsize  - neq rank
c imats  - input data set
c mesmt1,ermes1 - returned num and sumsq added (unless jmata.ne.0)
c weight - overall weighting factor for incrementing b,side
c ntape1,nseq1 - used in series-by-series save
c jmata  - flag for skipping lhs or rhs+stats
c itdif,prmdif - used in DPSNEC
c znsq   - correction to sumsq from prereduction
c sumzns - accumulated correction to sumsq
c measmt,ermeas - cumulative total num and sumsq
c iuseap - flag for restoring imbedded a prioris, instead of ordinary
c         information. if >0, don't update statistics or print summary
c iaprio - flag for existence of imbedded a priori information
c         this flag is 0 for Ibuf2, which is purely a priori, so that
c         nrmadd won't try to skip over (non-existent) ordinary stuff
c ivectk - if >0 then restore k vector 'vectk'
c wgtk  - weight for k vector restoration
c sumzsm - total correction to sum (o-c/error) due to prereduction
c zsm - incremental addition to sumzsm from latest data set
c sumsqa - contribution of a priori constraints to sumsq (o-c/error)
c sumaps - accumulated a priori sumsq
c
c         common
      include 'fcntrl.inc'
      include 'inodta.inc'
 
c        local
      real*10 derms(3),DOTN,du,rms,znsqw,zsmw
      integer*4 i,ia,irow,jrow,mjmat,mtst,nsav,nsav0
 
c
c read right side of saved normal equations
c*  start=1000
      do i = 1, 3
         derms(i) = 0._10
      end do
c skip ordinary information if only restoring a prioris
      if(iuseap.gt.0 .and. iaprio.gt.0) call BSKIP(imats,nsize)
c
      if(ivectk.gt.0 .and. iuseap.le.0) then
 
c restore k vector
         call QREAD(imats,mtst,sav,nsize)
         if(mtst.ne.-1) call SUICID(
     .       'INVALID CODE FOR K VECTOR, STOP IN NRMADD   ',11)
         do i = 1, nsize
            ia = iptr(i)
            if(ia.gt.0) vectk(ia) = vectk(ia) + sav(i)*wgtk
         end do
 
c correct sum for DPSNEC
         if(itdif.gt.0) derms(1) = derms(1)
     .                                 + DOTN(prmdif,sav,nsize)*wgtk
      endif
      if(jmata.eq.1) then
         read(imats)
      else
         call QREAD(imats,mtst,sav,nsize)
 
c increment right side of normal equations from saved values
         do i = 1, nsize
            ia = iptr(i)
            if(ia.gt.0) side(ia) = side(ia) + sav(i)*weight
         end do
 
c increment sumsq (o-c)/err for DPSNEC
         if(itdif.gt.0) derms(3) = 2._10*DOTN(prmdif,sav,nsize)
     .                                 *weight
      endif
      if(jmata.eq.2) then
 
         mjmat = -imats
         call BSKIP(mjmat,nsize)
      else
 
         nsav = 0
         do while( nsav.lt.nsize )
            nsav0 = nsav
 
c read row of saved normal equations
            call QREAD(imats,nsav,sav,nsize)
            if(nsav.le.nsav0) call SUICID(
     .          ' SAVED ROW COUNTERS DO NOT MATCH, STOP IN NRMADD',12)
c find row nrst of restored normal equations corresponding to
c given row nsav of saved normal equations
c*  start=1200
            irow = iptr(nsav)
            if(irow.gt.0) then
               jrow = (irow*(irow-1))/2
               do i = 1, nsize
                  ia = iptr(i)
                  if(ia.gt.0) then
                     if(ia.le.irow) then
                        ia    = ia + jrow
                        b(ia) = b(ia) + sav(i)*weight
                     endif
                  endif
               end do
 
c correct rhs for DPSNEC
               if(itdif.gt.0) then
                  du = DOTN(prmdif,sav,nsize)*weight
                  side(irow) = side(irow) + du
 
c correct sumsq (o-c)/err for DPSNEC
                  derms(3) = derms(3) + du*prmdif(nsav)
               endif
            endif
         end do
      endif
 
c finished restoring equations
      if(iuseap.gt.0) return
c correct statistics for imbedded a priori information
      if(iaprio.gt.0) then
         if(itdif.gt.0) then
            call QREAD(imats,nsav,sav,nsize)
            if(jmata.ne.1) derms(3) = derms(3) -
     .       2._10*DOTN(prmdif,sav,nsize)*weight
            nsav = 0
            do while( nsav.lt.nsize )
               nsav0 = nsav
c read row of a priori information matrix
               call QREAD(imats,nsav,sav,nsize)
               if(nsav.le.nsav0) call SUICID(
     .        ' SAVED A PRIORI ROW COUNTER MISMATCH, STOP IN NRMADD',13)
               if(jmata.ne.2) then
                  du = DOTN(prmdif,sav,nsize)*weight
                  derms(3) = derms(3) - du*prmdif(nsav)
               endif
            end do
         else
c no parameter disparities, skip a priori info
            call BSKIP(imats,nsize)
         endif
      endif
c correct statistics for DPSNEC
      do i = 1, 3
         ermes1(i) = ermes1(i) + derms(i)
      end do
      if(itdif.gt.0) then
         rms = 0._10
         if(mesmt1.gt.0 .and. ermes1(3).gt.0._10) rms = 
     .    SQRT(ermes1(3)/mesmt1)
         call PAGCHK(60,1,0)
         write(Iout,50) derms,rms
   50    format(' CORRECTION ADDED BY DISPARITIES:', 2x, 1p, 3D16.5,
     .          g42.4)
      endif
 
c finished restoring, save rhs if series-by-series
      if(jmata.ne.1 .and. Ict(18).gt.0)
     .    call SVSD(ntape1,nseq1,side)
 
c correct norm**2 for partial prereduction
      znsqw = -znsq*weight
      zsmw  = -zsm*wgtk
      if(znsqw.ne.0._10 .or. zsmw.ne.0._10) then
         rms = 0._10
         du = ermes1(3)+znsqw
         if(mesmt1.gt.0 .and. du.gt.0._10) rms = SQRT(du/mesmt1)
         call PAGCHK(60,1,0)
         write(Iout,100) zsmw,znsqw,rms
  100    format(' CORRECTION ADDED BY PRE-REDUCTION:', 1pd16.5, d32.5,
     .          g42.4)
         sumzns = sumzns - znsqw
         sumzsm = sumzsm - zsmw
      endif
      if(jmata.ne.1) then
 
c accumulate statistics
         measmt = measmt + mesmt1
         do i = 1, 3
            ermeas(i) = ermeas(i) + ermes1(i)
         end do
         sumaps = sumaps + sumsqa
      endif
      return
      end
