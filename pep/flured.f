      subroutine FLURED(jdu,fru,isit)

      implicit none

c
c j.f.chandler - 2013 jul - subroutine flured
c read and apply displacement values from an external direct access
c data set

c arguments
      integer*4 jdu,isit
      real*10 fru
c jdu = JD of desired epoch
c fru = fraction of a day past midnight of desired epoch
c jdt,frt are the converted values using the time offset of the data
c isit= 1 receiving site coordinates determined
c isit= 2 sending site coordinates determined

c array dimensions
      include 'globdefs.inc'
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obscrd.inc'
      include 'param.inc'
      include 'sitcrd.inc'

c local
c tflu= correction vector in meters (UWN)
c xflu= correction vector in light seconds (XYZ)
      integer*4 MAXSIT,nperday,MAXFREQ
      parameter (MAXSIT=20,MAXFREQ=50)
      integer*4 i,j,k,nsites,jd1(MAXSIT),jd2(MAXSIT),rec1(MAXSIT),
     . jdr(2),jdt,itype,np(2),nf,nfreq(3,2),ifreq,recnt
      character*4 name(MAXSIT),tname
      character*80 title,desc
      real*4 buf(240,3,3,2),
     . accel(MAXFREQ,3,2),phase(MAXFREQ,3,2),coef(MAXFREQ,6,3,2)
      real*8 ffreq(MAXFREQ,3,2),time,dt,ph,t0/2451545.5d0/,off(MAXSIT)
      real*10 f2,frt,s,t,tflu(3),xflu(3)
      integer*4 i12(2),ioce,iatm,ihyd,iprt,it,itsav(2),
     . jflu,jsit,npsave,nr,ns,recno(2)
      integer*4 iflag(3)
      equivalence (iflag(1),ioce),(iflag(2),iatm),(iflag(3),ihyd)
      real*10 tab(3,4,2),y1(3,2,2),y2(3,2,2)
      logical*4 nxtrp/.false./
c
      integer*4 nspecl/3/
      character*4 sname(3)/'TEXL','MLR2','MAU2'/
      character*4 snamx(3)/'MLRS','MLRS','MAU1'/
c
c
      jsit  = isit
      if(Sitf(1).eq.Sitf(2)) jsit=1
      ns=i12(jsit)
      if(Nk1.le.0) then
         if(Npage.ne.npsave) nxtrp = .false.
      endif

c correct for time offset, if any
      jdt=jdu
      frt=fru
      if(off(ns).ne.0d0) then
         frt=frt-off(ns)
         if(frt.lt.0._10) then
            frt=frt+1._10
            jdt=jdt+1
         endif
      endif

c is jd within range of table?
      if(jdt.lt.jd1(ns) .or. jdt.gt.jd2(ns) .or. ns.eq.0) then
c
c extrapolation beyond table
         if(ns.ne.0 .and. .not.nxtrp) then
            write(Iout,620) jdt
  620       format(i17,
     .        ' WARNING: NOT IN SITE DISPLACEMENT TABLE: USING 0. ***')
            Line=Line+1
            if(Mout.gt.0) write(Mout,620) jdt
            nxtrp  = .true.
            npsave = Npage
         endif
         return
      endif
c
c
c read data into storage if necessary
      if(recno(jsit).ne.0) then
         t = ((jdt-jdr(jsit))+frt)*nperday
         if(t.lt.1._10 .or. t.gt.238._10) then
            recno(jsit)=0
         endif
      endif
      if(recno(jsit).eq.0) then
         recno(jsit)=((jdt-jd1(ns))*nperday)/120 + rec1(ns)
         if(MOD((jdt-jd1(ns))*nperday,120).eq.0 .and. jdt.gt.jd1(ns))
     .    recno(jsit)=recno(jsit)-1
         read(jflu,rec=recno(jsit)) jdr(jsit),np(jsit),
     .    (((buf(j,itype,k,jsit),j=1,120),itype=1,3),k=1,3)
         if(np(jsit).eq.120 .and. jdr(jsit)+120/nperday.le.jd2(ns))
     .    read(jflu,rec=recno(jsit)+1) j,np(jsit),
     .    (((buf(j,itype,k,jsit),j=121,240),itype=1,3),k=1,3)
         t = ((jdt-jdr(jsit))+frt)*nperday
         itsav(jsit) = -9999
      endif

c
c calculate interpolation times and value of tabular points
      it = t
      t  = t - it
      s  = 1.0_10 - t
      if(it.ne.itsav(jsit)) then
         do i = 1,3
            do j = 1,4
               k = it+j-1
               tab(i,j,jsit) = 0.
               if(ioce.gt.0) tab(i,j,jsit)=buf(k,1,i,jsit)
               if(iatm.gt.0) tab(i,j,jsit)=tab(i,j,jsit)+buf(k,2,i,jsit)
               if(ihyd.gt.0) tab(i,j,jsit)=tab(i,j,jsit)+buf(k,3,i,jsit)
            end do
         end do
c
c calculate interpolation y-vectors
         do i = 1,3
            do j = 1,2
               nr = j + 1
               f2 = 0.166666666666666666667_10*(tab(i,nr+1,jsit)
     .          + tab(i,nr-1,jsit))
               y1(i,j,jsit)=1.33333333333333333333_10*tab(i,nr,jsit)-f2
               y2(i,j,jsit)=-0.33333333333333333333_10*tab(i,nr,jsit)+f2
            end do
         end do
         itsav(jsit) = it
      endif
c
c second difference interpolation
c correction vector in local coordinates (up, west, north) in meters
      do i=1,3
         tflu(i)= t*(y1(i,2,jsit)+t*t*y2(i,2,jsit)) +
     .    s*(y1(i,1,jsit)+s*s*y2(i,1,jsit))
      end do
c add harmonic components
      dt=((jdt-t0)+frt)*864d2
      do itype=1,3
         if(iflag(itype).gt.0) then
            nf=nfreq(itype,jsit)
            do ifreq=1,nf
               ph=phase(ifreq,itype,jsit)+dt*(ffreq(ifreq,itype,jsit)+
     .          0.5d0*accel(ifreq,itype,jsit)*dt)
               do i=1,3
                  tflu(i)=tflu(i)+COS(ph)*coef(ifreq,i,itype,jsit)+
     .             SIN(ph)*coef(ifreq,i+3,itype,jsit)
               end do
            end do
         endif
      end do
c transform correction to equatorial system and apply
      do i=1,3
         xflu(i)= (tflu(1)*Dxdrad(i,isit)+
     .    tflu(2)*Dxdlon(i,isit)/Rc(isit)+
     .    tflu(3)*Dxdlat(i,isit)/Rsite(isit))/(Ltvel*1000._10)
         Xsite(i,isit)= Xsite(i,isit)+xflu(i)
      end do
      if(iprt.gt.0) then
         if(Line.gt.56) call OBSPAG
         write(Iout,750) jdt,frt,isit,tflu,xflu
  750    format('FLURED: JD.F=', i7, f13.12, ' NS=', i2,
     .    '  X=', 1p3d13.5, ' M (UWN)'/ 49x,3d13.5,' LS (XYZ)')
         Line = Line + 2
      endif
      return

c initialize
      entry FLURD1
      recno(1)=0
      recno(2)=0
      itsav(1)= -9999
      itsav(2)= -9999
      npsave = 0
      j = Jct(49)/100
      jflu=Jct(49)-100*j
      iprt=j/100
      j=MOD(j,100)
      if(j.eq.0) j=7
      ioce=MOD(j,2)
      iatm=MOD(j/2,2)
      ihyd=MOD(j/4,2)
c
c read and print header records
      if(Itrwnd(jflu).le.0) then
         open(jflu, access='DIRECT', recl=4328, status='OLD')
         read(jflu,rec=1) title,desc,nsites,nperday,
     .    (name(i),jd1(i),jd2(i),rec1(i),i=1,nsites),
     .    (off(i),i=1,nsites)

         k=(nsites-1)/4 + 1
         if(Line.gt.58-k) call NEWPG
         Line=Line+k+2
         write(Iout,300) jflu,nsites,nperday,title,
     .    (name(i),jd1(i),jd2(i),i=1,nsites)
  300    format('0HEADER DATA ON DISPLACEMENT DATA SET',i3,' WITH',i3,
     .    ' SITES AT',i3,' POINTS PER DAY'/
     .    ' TITLE=',a80/(4(1x,a4,i8,'-',i7)))
         if((120/nperday)*nperday.ne.120) call SUICID(
     .    'INVALID NUMBER OF DISPLACEMENTS/RECORD, STOP IN FLURED  ',14)
         Itrwnd(jflu) = 1
      endif
c
      do jsit=1,2
c translate special names to equivalents
         tname=Sitf(jsit)(1:4)
         do i=1,nspecl
            if(tname.eq.sname(i)) then
               tname=snamx(i)
               goto 310
            endif
         end do
c find site in database
  310    do ns=1,nsites
            if(tname.eq.name(ns)) goto 330
         end do
         if(Line.gt.58) call NEWPG
         write(Iout,320) Sitf(jsit)
  320    format(' SITE ',a8,' NOT FOUND IN DATABASE ***')
         Line=Line+1
         ns=0
c do send site if different
  330    i12(jsit)=ns
         if(ns.gt.0) then
            do itype=1,3
               recnt=(itype-1)*nsites+1+ns
               read(jflu,rec=recnt) nf,(ffreq(ifreq,itype,jsit),
     .          phase(ifreq,itype,jsit),accel(ifreq,itype,jsit),
     .          (coef(ifreq,i,itype,jsit),i=1,6),ifreq=1,nf)
               nfreq(itype,jsit)=nf
            end do
         endif
         if(Sitf(1).eq.Sitf(2) .or. Sitf(2).eq.' ') return
      end do
      return

c close out
      entry FLURD9
      if(jflu.le.0) return
      if(Itrwnd(jflu).gt.0) then
         close(jflu)
         Itrwnd(jflu)=0
      endif
      return
      end
