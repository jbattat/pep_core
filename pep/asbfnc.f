      subroutine ASBFNC(ncall)
 
      implicit none

c        subroutine asbfnc - j.f.chandler, 1988 apr 27
c
c     calculate effect on motion and partials of the distributed belts
c     of asteroids (or other similar belts).  this subroutine replaces
c     code formerly found in solprb or plnorb.  it has been completely
c     re-implemented to
c     (a) get the correct force term,
c     (b) improve convergence near the sphere of ring radius,
c     (c) calculate correct partials for inclination and node, and
c     (d) speed up the calculations.
c
c arguments
      integer*4 ncall
c        ncall=-2 setup once per iteration of a given step for partials
c        ncall=-1 setup once per iteration of a given step for motion
c        ncall= 0 setup once per step
c        ncall= k evaluate acceleration for equation k.gt.0
c
c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'astroi.inc'
      include 'intstf.inc'
      include 'petuna.inc'
      include 'sbstuf.inc'
      include 'sbthng.inc'
      real*10 cor(10,3)
      equivalence (cor(1,1),Bcor)
 
c external functions
      real*10 DOT

c local
      real*10 ast1,ast1s,ast2,ast2s,ast3,ast4,ast5,ast6,ast7,
     .          astd2,astds,atrm1,atrm2,atrm3,coth,dml,dsq,
     .          fct1,fct2
      real*10 fml,hbyr,hcbyr,hoff,hr,hsq,om2hc,omhc,ra0,rab,
     .          ratio1,rbsq,rg,rs,sith,trm1,tsb,tsb2
      integer   i,ibd,ibelt,jbd,ml
      real*10 leg(32),leg1(32),leg2(32),tbcor(3)
      real*10 er(3,2,4),vth(3,2,4),astws5(2,4),astws6(2,4),
     .          asta(2,4),astb(2,4),astc(2,4),astd(2,4),aste(2,4),
     .          astf(2,4),astg(2,4),dcthdi(2,4),dcthdo(2,4)
 
      if(ncall.lt.0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c set-up once per iteration of a given step for motion
         if(ncall.ge.-1) then
            do ibelt = 1,Nbelt
               if(Kpasb(ibelt).ge.0) then
                  hsq = Asbpa(ibelt)**2
                  do ibd = 1,2
                     dsq = hsq + hsq
                     jbd = ibd
c
c assume a belt around a planet affects only that planet's satellites
c ibd=1 -> force on body
c     2 -> force on center
c jbd=1 -> body relative to sun
c     2 -> center rel. to sun
c     3 -> body rel. to center
c
                     if(Ncnasb(ibelt).gt.0) jbd = 3
                     hr = 2._10*Asbpa(ibelt)
     .                    *DOT(cor(1,jbd),Astnrm(1,ibelt))
 
c determine which origin to use for expansion
                     tsb2  = cor(8,jbd) + hsq + hr
                     astd2 = tsb2/dsq
                     if(astd2.lt.1.E-8_10) astd2 = 1._10
                     if(astd2.gt.1._10) astd2  = 1._10/astd2
                     rbsq = cor(8,jbd) + hsq - hr
                     rab  = rbsq/dsq
                     if(rab.lt.1.E-8_10) rab = 1._10
                     if(rab.gt.1._10) rab  = 1._10/rab
                     ra0 = cor(8,jbd)/hsq
                     if(ra0.gt.1._10) ra0 = 1._10/ra0
                     rs   = Asbpa(ibelt)/cos45
                     hoff = Asbpa(ibelt)
                     if(rab.lt.astd2) then
                        tsb2  = rbsq
                        hoff  = -Asbpa(ibelt)
                        astd2 = rab
                     endif
                     if(ra0.gt.astd2) then
                        tsb = SQRT(tsb2)
                     else
                        rs    = Asbpa(ibelt)
                        tsb2  = cor(8,jbd)
                        tsb   = cor(7,jbd)
                        dsq   = hsq
                        hoff  = 0._10
                        astd2 = ra0
                     endif
                     do i = 1,3
                        tbcor(i) = cor(i,jbd) + hoff*Astnrm(i,ibelt)
                        er(i,ibd,ibelt) = tbcor(i)/tsb
                     end do
                     coth = DOT(Astnrm(1,ibelt),er(1,ibd,ibelt))
                     sith = SQRT(1._10 - coth**2)
                     call LEGNDR(coth,sith,32,0,leg(2),leg1(2),
     .                           0._10,0._10)
                     if(Astapr(ibelt))
     .                   call LEGND2(coth,sith,32,0,leg(2),leg1(2)
     .                   ,leg2(2),0._10,0._10,0._10)
                     do i = 1,3
                        vth(i,ibd,ibelt) = coth*er(i,ibd,ibelt)
     .                     - Astnrm(i,ibelt)
                     end do
                     hbyr  = hoff/tsb
                     hcbyr = hbyr*coth
                     omhc  = 1._10 - hcbyr
                     om2hc = omhc - hcbyr
                     rg    = MAX(rs,tsb)
                     astws5(ibd,ibelt) = Gamat/tsb/rg
                     astws6(ibd,ibelt) = astws5(ibd,ibelt)
     .                  *Asbpm(ibelt)
 
c most 0th terms vanish
                     ast1 = 0._10
                     if(tsb.gt.rs) ast1 = -1._10
                     ast3 = -ast1
                     if(hoff.eq.0._10) then
 
c expansion origin at center of ring
                        ast2  = 0._10
                        ast4  = 0._10
                        ast5  = 0._10
                        astds = 1.0_10
                        do ml = 2,32,2
                           astds = astds*astd2
                           trm1  = astds*Zz(ml)
                           fml   = ml
                           if(tsb.gt.rs) fml = -ml - 1
 
c sum terms for motion
                           atrm1 = fml*trm1*leg(ml)
                           atrm2 = trm1*leg1(ml)
                           ast1  = ast1 + atrm1
                           ast2  = ast2 + atrm2
                           if(ml.eq.2) then
 
c truncate sum at r*8 limit
                              ast1s = 1E-17_10*ABS(ast1)
                              ast2s = 1E-17_10*ABS(ast2)
                           endif
                           if(ABS(atrm1).lt.ast1s .and. ABS(atrm2)
     .                        .lt.ast2s) goto 20
                           if(Kpasb(ibelt).gt.0) then
 
c sum extra terms for partials
                              ast3 = ast3 + fml*atrm1
                              ast4 = ast4 + fml*atrm2
                              if(Astapr(ibelt)) ast5 = ast5 +
     .                            trm1*leg2(ml)
                           endif
                        end do
                     else
c
c expansion origin offset by 1 ring radius
                        ratio1 = tsb/rs
                        if(ratio1.ge.1._10) ratio1 = 1._10/ratio1
                        if(hoff.lt.0._10) ratio1    = -ratio1
                        trm1 = cos45*ratio1
                        fml  = 1
                        dml  = 1
                        if(tsb.gt.rs) fml = -2
                        if(tsb.gt.rs) dml = -1
                        atrm1 = fml*trm1*coth
                        ast1  = ast1 + atrm1
                        ast2  = trm1
                        ast3  = ast3 + fml*atrm1
                        ast4  = fml*ast2
                        ast5  = 0._10
                        ast6  = fml*ratio1*coth
                        ast7  = ratio1
 
c truncate sum at r*8 limit
                        ast1s = 1E-17_10*ABS(ast1)
                        ast2s = 1E-17_10*ABS(ast2)
                        astds = ratio1
                        do ml = 2,32
                           astds = astds*ratio1
                           trm1  = astds*Leg45(ml)
                           fml   = fml + dml
 
c sum terms for motion
                           atrm1 = fml*trm1*leg(ml)
                           atrm2 = trm1*leg1(ml)
                           ast1  = ast1 + atrm1
                           ast2  = ast2 + atrm2
                           if(ABS(atrm1).lt.ast1s .and. ABS(atrm2)
     .                         .lt.ast2s) goto 10
                           if(Kpasb(ibelt).gt.0) then
 
c sum extra terms for partials
                              ast3  = ast3 + fml*atrm1
                              ast4  = ast4 + fml*atrm2
                              atrm3 = astds*Leg145(ml)
                              ast6  = ast6 + fml*atrm3*leg(ml)
                              ast7  = ast7 + atrm3*leg1(ml)
                              if(Astapr(ibelt)) ast5 = ast5 +
     .                            trm1*leg2(ml)
                           endif
                        end do
   10                   if(hoff.le.0._10 .and. Kpasb(ibelt).gt.0)
     .                      then
                           ast6 = -ast6
                           ast7 = -ast7
                        endif
                     endif
   20                asta(ibd,ibelt) = ast1*astws5(ibd,ibelt)
                     astb(ibd,ibelt) = -ast2*astws5(ibd,ibelt)
                     if(Kpasb(ibelt).gt.0) then
                        fct1 = hoff/rs
                        fct2 = Asbpa(ibelt)/dsq
                        astc(ibd,ibelt) = -(ast1 + ast3 + fct1*ast6)
     .                     *fct2
                        astd(ibd,ibelt) = (ast2 + ast4 + fct1*ast7)
     .                     *fct2
                        aste(ibd,ibelt) = omhc*ast4 - om2hc*ast2 +
     .                     hbyr*(ast3 - ast1 - ast1)
                        astf(ibd,ibelt) = -omhc*ast5 +
     .                     hbyr*(ast4 - ast2)
                        astg(ibd,ibelt) = omhc*ast2 + hbyr*ast1
                        dcthdi(ibd,ibelt)
     .                     = DOT(Astcom(1,ibelt),er(1,ibd,ibelt))
                        dcthdo(ibd,ibelt)
     .                     = DOT(Astvc3(1,ibelt),er(1,ibd,ibelt))
                     endif
                     if(Ncnasb(ibelt).eq.Ncentr) goto 40
                  end do
               endif
   40       end do
         endif
      else if(ncall.ne.0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c calculate acceleration for motion
         if(Kkk.gt.0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c calculate acceleration for partials
            do ibelt = 1,Nbelt
               if(Kpasb(ibelt).ge.1) then
 
c test for mass partial
                  if(Icntrl(Kkk).eq.Icasbm(ibelt)) then
                     Fn(1) = Fn(1) + asta(1,ibelt)*er(Kk,1,ibelt)
     .                       + astb(1,ibelt)*vth(Kk,1,ibelt)
                     if(Ncnasb(ibelt).ne.Ncentr) Fn(1) = Fn(1)
     .                   - asta(2,ibelt)*er(Kk,2,ibelt)
     .                   - astb(2,ibelt)*vth(Kk,2,ibelt)
                     return
 
c test for radius partial
                  else if(Icntrl(Kkk).eq.Icasba(ibelt)) then
                     Fn(1) = Fn(1)
     .                       + (astc(1,ibelt)*er(Kk,1,ibelt) + astd(1,
     .                       ibelt)*vth(Kk,1,ibelt))*astws6(1,ibelt)
                     if(Ncnasb(ibelt).ne.Ncentr) Fn(1) = Fn(1)
     .                   - (astc(2,ibelt)*er(Kk,2,ibelt)
     .                   + astd(2,ibelt)*vth(Kk,2,ibelt))
     .                   *astws6(2,ibelt)
                     return
 
c test for inclination partial
                  else if(Icntrl(Kkk).eq.Icasbi(ibelt)) then
                     Fn(1) = Fn(1)
     .                       + ((aste(1,ibelt)*er(Kk,1,ibelt)+astf(1,
     .                       ibelt)*vth(Kk,1,ibelt))*dcthdi(1,ibelt)
     .                       + astg(1,ibelt)*Astcom(Kk,ibelt))
     .                       *astws6(1,ibelt)
                     if(Ncnasb(ibelt).ne.Ncentr) Fn(1) = Fn(1)
     .                   - ((aste(2,ibelt)*er(Kk,2,ibelt)+astf(2,ibelt)
     .                   *vth(Kk,2,ibelt))*dcthdi(2,ibelt)
     .                   + astg(2,ibelt)*Astcom(Kk,ibelt))
     .                   *astws6(2,ibelt)
                     return
 
c test for node partial
                  else if(Icntrl(Kkk).eq.Icasbo(ibelt)) then
                     Fn(1) = Fn(1)
     .                       + ((aste(1,ibelt)*er(Kk,1,ibelt)+astf(1,
     .                       ibelt)*vth(Kk,1,ibelt))*dcthdo(1,ibelt)
     .                       + astg(1,ibelt)*Astvc3(Kk,ibelt))
     .                       *astws6(1,ibelt)
                     if(Ncnasb(ibelt).eq.Ncentr) return
                     Fn(1) = Fn(1)
     .                       - ((aste(2,ibelt)*er(Kk,2,ibelt)+astf(2,
     .                       ibelt)*vth(Kk,2,ibelt))*dcthdo(2,ibelt)
     .                       + astg(2,ibelt)*Astvc3(Kk,ibelt))
     .                       *astws6(2,ibelt)
                  endif
               endif
            end do
         else
            do ibelt = 1,Nbelt
               if(Kpasb(ibelt).ge.0) then
                  Fn(1) = Fn(1)
     .                    + (asta(1,ibelt)*er(Kk,1,ibelt) + astb(1,
     .                    ibelt)*vth(Kk,1,ibelt))*Asbpm(ibelt)
                  if(Ncnasb(ibelt).ne.Ncentr) Fn(1) = Fn(1)
     .                - (asta(2,ibelt)*er(Kk,2,ibelt) + astb(2,ibelt)
     .                *vth(Kk,2,ibelt))*Asbpm(ibelt)
               endif
            end do
         endif
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c set-up once per step
      endif
 
      return
      end
