      subroutine FERSTR(nvel,norm,kspt)
 
      implicit none
 
c arguments
      integer*4 nvel,norm,kspt

c
c r. king/s. gourevitch   june 1977    subroutine ferstr
c calculate interferometry observables for natural sources.
c this code lifted from old pep subroutine interf.
c
c

c array dimensions
      include 'globdefs.inc'

c commons 
      include 'comdateq.inc'
      include 'coord.inc'
      real*10 dvect(6)
      equivalence (dvect,Angdum)
      real*10 ctat2,dutrec
      equivalence (ctat2,Angdum(9)),(dutrec,Angdum(10))
      include 'difnct.inc'
      include 'fcntrl.inc'
      include 'mnsprt.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct,freq2
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      equivalence (freq2,Dstf(5))
      real*4 acctim
      equivalence (acctim,Estf)
      include 'kobequiv.inc'
      include 'param.inc'
      include 'radcrd.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'
c external functions
      real*10 CTATF,DCTATF,DOT
c
c quantities internal to this routine
      integer*4 etide1
      integer*2 npl
      real*10 bdotx1,bdotx2,dctdt1,ddr,ddr1,dfdly,dfdlya,dum,fractr
      integer*4 i,icount,index,j,jdr,lemctl,mnspt1,nvlesn
c
c
c-----------initialization------------------------------------------
c
      etide1 = 0
      lemctl = -1
      Nswcns = 0
      nvlesn = nvel
      if(Ict(27).ne.0) Nswcns = 1
      if(Ict(27).ne.0 .or. prmter(81).gt.0._10) nvlesn = 1
      if(Jct(10).gt.0) nvlesn = 1
c
c
c           calculate earth, moon, sun vectors at receive time
c           for reference site
c
c         e-m barycenter w.r.t. sun   (stored in xem(1-6,1), in a.u.)
      if(Nswcns.gt.0 .or. Reldel.gt.0._10 .or. Jct(10).gt.0)
     .    then
c if calculations are geocentric with no relativity and no earth
c tides, then the e-m barycenter and solar system barycenter positon
c are not needed.
         call ETRP(1,Jd,ctrecf,0,lemctl,1,2)
         if(Jd.le.0) return
         lemctl = 0
         if(nvlesn.gt.0) call ETRP(1,Jd,ctrecf,1,0,1,2)
c
c correct reference site coordinates for solid body tides
         if(Jct(10).gt.0) call ETIDE(etide1,nvlesn,1,kspt)
      endif
c
c save site position in xslcns
      do i = 1, 6
         Xslcns(i,1) = Xsite(i,1)
      end do
      if(Nswcns.gt.0) then
 
c sun w.r.t. ssbc
         call SOLCNT(Jd,ctrecf,Xslcns(1,1),nvlesn)
 
c earth velocity w.r.t. ssbc
         do i = 4, 6
            Xslcns(i,1) = Xemlsc(i,1) - Xslcns(i,1)*Aultvl
         end do
         bdotx1 = DOT(Xsite(1,1),Xslcns(4,1))
 
c site w.r.t. ssbc (corrected)
         do i = 1, 3
            Xslcns(i,1)     = Xemlsc(i,1) - Xslcns(i,1)
     .                         *Aultsc + Xsite(i,1)
     .                         + .5_10*bdotx1*Xslcns(i + 3,1)
            Xslcns(i + 3,1) = Xslcns(i + 3,1) + Xsite(i + 3,1)
         end do
      endif
c
c determine star coodinates in standard equatorial system
      mnspt1 = 0
      npl    = -4
      call SPOTCD(Jd,Fract,mnspt1,0,npl,0.0_10,0.0_10,kspt)
c
c-----------iteration to determine recieve time at second site----------
c
      dfdly  = 0.0_10
      icount = 0
      do while( .true. )
         call TIMINC(Jd,Fract,jdr,fractr,-dfdly/Secday)
c
c determine second site position
         Sidtm = Sidtm0 + (Utrec-dfdly)*Sidvel + Dgst
         call SITCOR(Sidtm,2,nvel,norm)
c
c determine earth w.r.t. center of mass of solar
c system at receive time at second site
         if(Nswcns.gt.0 .or. Reldel.gt.0._10 .or. Jct(10).gt.0)
     .       then
            call ETRP(1,jdr,fractr,0,lemctl,2,2)
            if(jdr.eq.0) then
 
c error return
               Jd = 0
               return
            else
               if(nvlesn.gt.0) call ETRP(1,jdr,fractr,1,0,2,
     .             2)
c
c correct second site coordinates for solid body tides
               if(Jct(10).gt.0) then
                  call ETIDE(etide1,nvlesn,2,kspt)
                  etide1 = 1
               endif
 
c save site in xslcns
               do i = 1, 6
                  Xslcns(i,2) = Xsite(i,2)
               end do
 
c sun w.r.t. ssbc
               if(Nswcns.gt.0) then
                  call SOLCNT(jdr,fractr,Xslcns(1,2),nvlesn)
 
c earth velocity w.r.t. ssbc
                  do i = 4, 6
                     Xslcns(i,2) = Xemlsc(i,2) - Xslcns(i,2)*Aultvl
                  end do
                  bdotx2 = DOT(Xsite(1,2),Xslcns(4,2))
 
c site w.r.t. ssbc (corrected)
                  do i = 1, 3
                     Xslcns(i,2)     = Xemlsc(i,2) - Xslcns(i,2)
     .                                  *Aultsc + Xsite(i,2)
     .                                  + .5_10*bdotx2*Xslcns(i + 3,2)
                     Xslcns(i + 3,2) = Xslcns(i + 3,2)
     .                                  + Xsite(i + 3,2)
                  end do
               endif
            endif
         endif
c
c determine vector from one site to the other
         do i = 1, 3
            dvect(i) = Xslcns(i,1) - Xslcns(i,2)
         end do
c
c determine differential delay
         dfdlya = -DOT(Xspcd(1,kspt),dvect(1))
         if(Nswcns.gt.0) dfdlya = dfdlya - bdotx1 + bdotx2
         if(ABS(dfdlya-dfdly).gt.acctim) then
            dfdly  = dfdlya
            icount = icount + 1
            if(icount.le.10) goto 100
            call SUICID(
     . ' MORE THAN 10 DIFFERENTIAL DELAY ITERATIONS, STOP IN FERSTR ',
     . 15)
         endif
c
c
c-----------calculations after iterations are completed
c
         dfdly   = dfdlya
         Nit(20) = Nit(20) + icount
 
c set up xsitep, etc.
         if(kspt.eq.2) then
            do j = 1, 2
               do i = 1, 3
                  Ysitep(i,j)     = -Xspcd(i,2)
                  Ysitep(i + 3,j) = 0._10
               end do
               call UVECTR(3,Ysitep(1,j),Rsitp2(j),Ysitp0(1,j),dum)
            end do
         else
            do j = 1, 2
               do i = 1, 3
                  Xsitep(i,j)     = -Xspcd(i,1)
                  Xsitep(i + 3,j) = 0._10
               end do
               call UVECTR(3,Xsitep(1,j),Rsitp(j),Xsitp0(1,j),dum)
            end do
         endif
c
c test for dummy observation below the horizon
         if(Idumob.eq.1) then
            if(kspt.eq.1) call HORIZN(Xsite,Xsitep,2,Jd)
            if(kspt.eq.2) call HORIZN(Xsite,Ysitep,2,Jd)
            if(Jd.le.0) return
         endif
         goto 200
  100 end do
c
c
c           correction to delay for general relativity (light bending)
c           to go here
c
c
c           calculate star coordinate partials
  200 call SPOTCV(0,npl,kspt)
 
      if(Nice.le.0) then
c
c convert coordinate time delay to atomic time
         ctat2 = CTATF(jdr,fractr - Ctat/Secday,3,2)
         dfdly = dfdly - Ctat + ctat2
c
c convert atomic time delay to utc
         if(ntime.gt.0) dfdly = dfdly*fdev
         if(Nice.lt.0) goto 300
      endif
c
c calculate delay rate + correct for effects
c independent variable in delay rates is rate of clock at site 1
      do i = 4, 6
         dvect(i) = Xslcns(i,1) - Xslcns(i,2)
      end do
      dctdt1 = DCTATF(Jd,Fract - Ctat/Secday,3,1) + 1._10
      ddr1   = -DOT(Xspcd(1,kspt),dvect(4))
     .         /(1._10 + DOT(Xspcd(1,kspt),Xslcns(4,2)))
      ddr    = ddr1 - (dctdt1 - 1._10) + dctdt1*(1._10 - ddr1)
     .         *DCTATF(jdr,fractr - Ctat/Secday,3,2)
c
c   ddr does not include an acceleration term, or a relativistic
c    correction calculated by sergei gourevitch
c    both terms are about 1.0e-15 sec/sec max.
c
c      effect on ddr of earth tides ignored
c
      if(ntime.gt.0) ddr  = ddr*fdev
      if(Nice.gt.0) dfdly = 0.0_10
c
c form differences and save coordinates if necessary
  300 Difdly(1,kspt) = dfdly
      Difdly(2,kspt) = ddr
      if(kspt.ne.1) then
         index = 3
         if(nvel.eq.1) index = 6
         do i = 1, index
            Dvect2(i) = dvect(i)
         end do
         Jd = 0
      endif
c
c
      return
      end
