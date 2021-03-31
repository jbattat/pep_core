      subroutine ANGOT(marta)

      implicit none
      integer*4 marta

c   T.Forni  November 1975 subroutine ANGOT
c   Special BCD site tape modified format of Mason's site tape
c   called by subroutine COMRIT if(JCT(45).ne.0)
c   Arecibo BCD site tape written on data set JCT(45).gt.0
c   Arecibo print out written on data set Jout in &NMLST1
c
c   JCT(45)= 0   Only td&dop regular print in radar
c   JCT(45)< 0   Also print out on data set JOUT (&NMLST1)
c                of ra&decl, az&el, td&dop for receive time of
c                td&dop.
c   JCT(45)> 0   Also Arecibo BCD site tape on data set JCT(45).
c
c JCT(70)    packed bits controls for ANGOT
c    1 abbreviated print to JOUT in ANGOT if JCT(45).ne.0 (else verbose)
c    2 use mean equator/equinox of reference epoch (else true of date)
c    4 increase reference epoch by 50 yrs (else use default)

c  Note: ICT(2)=1  would suppress the printing of td&dop obs in radar
c
c   marta=0 COMRIT>ANGOT called in midst of observing series
c   marta=1 COMRIT>ANGOT called at end of observing series
c   marta=2 ANGOT called at end of run
c
c array dimensions
      include 'globdefs.inc'

c commons
      include 'coord.inc'
      include 'fcntrl.inc'
      include 'funcon.inc'
      include 'inodta.inc'
      include 'loadid.inc'
      include 'mnsprt.inc'
      include 'namtim.inc'
      include 'number.inc'
      include 'nutprc.inc'
      include 'obscrd.inc'
      real*10 ctrecf,fdev,reflct
      equivalence (ctrecf,Dstf),(fdev,Dstf(10)),(reflct,Dstf(4))
      real*10 tmdly0,dop0
      equivalence (Result(1),tmdly0),(Result(2),dop0)
      include 'param.inc'
      include 'redobs.inc'
      include 'sitcrd.inc'
      include 'yvect.inc'

      character*3 jhr,irh,idd
      real*10 xmerid(3),xop(3),xmerop(3)
      real*10 zsat(3),xsat(3)
c WTVRSN should be incremented when the "work tape" header changes
      integer*4 wtvrsn/2/,itf/0/,lpg/0/,senwkt,izn(2)
      integer*4 yepch
      character*4 cepch
      character*4 shlab,stz(0:23)/' GMT','-100','-200','-300',' AST',
     . ' EST',' CST',' MST',' PST','-900',' -10',' -11',' -12',' +11',
     . ' +10','+900','+800','+700','+600','+500','+400','+300','+200',
     . '+100'/

      real*10 altprc(3,3),az,azx,d,ddec,decl,declx,decs,DOT,dra,dt,
     . el,elx,qxn,qxo1,qxo2,ra,ras,rax,rsat,rsat2,rzsat,sec1,secs,
     . sfract,slon,svdec,svra,svt,t,tlst,tsv
      integer*4 i,idm,irm,ishr,j,jct45,jdalt,jd4,jdsite,jdsprv,
     . jjds,jmin,ltf

      character*8 svsitf(2)
      real*10 svfreq/0._10/
      integer*2 svnplnt0/9999/
      character*4 svspotf,svspotf2

      if(marta.gt.0) goto 200

      if(itf.ne.1) then
c setup at beginning of series
         itf=1

c make sure the header information is all the same, else rewind
         if(Sitf(1).ne.svsitf(1) .or. Sitf(2).ne.svsitf(2) .or. 
     .    Freq.ne.svfreq .or. Nplnt0.ne.svnplnt0 .or. Spotf.ne.svspotf 
     .    .or. Spotf2.ne.svspotf2) then
            if(svnplnt0.ne.9999) then
               call SUICID(
     .'WARNING IN ANGOT: IMCOMPATIBLE SERIES, REWINDING OUTPUT FILES   '
     .,-16)
c rewind output files to start again
               if(jct45.gt.0) then
                  endfile jct45
                  rewind jct45

c send time work tape
                  if(senwkt.gt.0) then
                     endfile senwkt
                     rewind senwkt
                  endif
               endif
               if(Jout.gt.0) then
                  endfile Jout
                  rewind Jout
               endif
            endif
c save series information
            svsitf(1)=Sitf(1)
            svsitf(2)=Sitf(2)
            svfreq=Freq
            svnplnt0=Nplnt0
            svspotf=Spotf
            svspotf2=Spotf2

c get time zones
            do i=1,2
               slon=Coords(2,i)
               if(slon.lt.0._10) slon=slon+360._10
               izn(i)=slon/15._10+0.5_10
            end do
            if(izn(1).ge.0 .and. izn(1).le.23) then
               shlab=stz(izn(1))
            else
               shlab=' ***'
            endif
            ltf   =60
            tsv   =0._10
            t     =0._10
            jdsprv=0

c write two header records on bcd site tape data set = jct(45)
            jct45 =Jct(45)
            if(mod(Jct(70)/2,2).eq.0) then
               cepch='DATE'
               jdalt=0
            else
               if(Jct(13).eq.0) then
                  yepch=1950
                  jdalt=2451545
               else
                  yepch=2000
                  jdalt=2469808
               endif
               if(mod(Jct(70)/4,2).eq.1) then
                  yepch=yepch+50
                  call PRECES(jdalt+0._10)
                  do i=1,3
                     do j=1,3
                        altprc(i,j)=Prec(i,j)
                     end do
                  end do
               else
                  jdalt=0
               endif
               call EBCDI(yepch,cepch,4)
            endif
            senwkt=0
            if(jct45.gt.0) then
               if(jct45.ge.50) senwkt=jct45+10
               write(jct45,10) Heding,Date,cepch
   10          format('  PEP PREDICTION EPHEMERIS   ',18A4,4x,2A4/
     .          5x,'JDDATE',4x,' HR MN SH',2x,
     .          'ELEVATION-DEG.  AZIMUTH-DEG.',2x,'T.D.-SEC',
     .          5x,'DOP-HZ',7x,'RA-HR (',a4,') DECL-DEG')

c write header record for send time work tape
               jd4=Jdc(2)
               if(Ihrc(2).eq.0 .and. Iminc(2).eq.0) jd4=jd4-1
               if(senwkt.gt.0) write(senwkt) Heding,Date,wtvrsn,
     .          Aplnt(Klap),Sitf,Freq,cepch,Nplnt0,Ncp0,Snrm,izn,
     .          Jdc(1),jd4,Intscc,Lnklvl,Spotf,Spotf2,
     .          (SQRT(Rc(i)**2+Rs(i)**2)*Ltvel,
     .          Longr(i)/Convd,ATAN2(Rs(i),Rc(i))/Convd,i=1,2),
     .          Spcdx
            endif
         endif
      endif
c
c aberration correction
      do i=1,3
         zsat(i)=-Xsitep(i,1) + Rsitp(1)*Xemlsc(i+3,1)
      end do
      rzsat=SQRT(DOT(zsat,zsat))
c
c determine azimuth,elevation (deg)
      d  =Sitnrm(1,1)*Nutpr(3,1) + Sitnrm(2,1)*Nutpr(3,2)
     .     + Sitnrm(3,1)*Nutpr(3,3)
      qxn=DOT(zsat,Sitnrm(1,1))
      do i=1,3
         xmerid(i)=Nutpr(3,i)-d*Sitnrm(i,1)
         xop(i)   =zsat(i)-qxn*Sitnrm(i,1)
      end do
      qxn=qxn/rzsat
      el =ASIN(qxn)/Convd
      call CROSS(xop,xmerid,xmerop)
      qxo2=DOT(xmerop,Sitnrm(1,1))
      qxo1=DOT(xmerid,xop)
      az  =ATAN2(qxo2,qxo1)/Convd
      if(az.lt.0._10) az=az+360._10
c
c transform to true equinox and equator of date, if desired
      if(mod(Jct(70)/2,2).eq.0) then
         call CORCHN(xsat,zsat)
      else if(jdalt.ne.0) then
         call PRODCT(altprc,zsat,xsat,3,3,1)
      else
         do i=1,3
            xsat(i)=zsat(i)
         end do
      endif
c
c compute right ascension,declination      (hours,degrees)
      if(Jout.gt.0) then
         tlst=t
         t   =Jds+ctrecf
         if(tlst.ne.0._10 .and. ABS(t-tsv).ge.0.3_10) then
            dt=(tlst-svt)*24._10
            if(dt.ne.0._10) then
               dra=(ra-svra)*3600._10
               if(ABS(dra).gt.4.E4_10) dra=dra-SIGN(dra,86400._10)
               dra =dra/dt
               ddec=(decl-svdec)*3600._10/dt
            endif
         endif
      endif
      rsat2=DOT(xsat,xsat)
      rsat =SQRT(rsat2)
      decl =ASIN(xsat(3)/rsat)/Convd
      ra   =ATAN2(xsat(2),xsat(1))
      if(ra.lt.0._10) ra=ra+Twopi
      ra=ra/Convhs/3600._10
c
c ishr= local time
      ishr=Ihr-izn(1)
      if(ishr.lt.0) ishr=ishr+24
      jd4=Jds-2440001

      if(jct45.gt.0) then
c
c write ebcdic site tape
         jdsite=Jds-1

c result(1-2)  always  td&dop
         if(ABS(Result(1)).lt.9999.999999995_10) then
            write(jct45,20) jdsite,Ihr,Imin,ishr,el,az,Result(1),
     .                    Result(2),ra,decl
   20       format(' 0.',I7,'5E 07',3I3,1P,2D16.9,0PF13.8,F12.3,2F11.7)
         else
            write(jct45,25) jdsite,Ihr,Imin,ishr,el,az,Result(1),
     .                    Result(2),ra,decl
   25       format(' 0.',I7,'5E 07',3I3,1P,2D16.9,0PF13.7,F12.3,2F11.7)
         endif
      endif
c
c convert to integer parts
      call D2DMS(irh,irm,ras,ra,4)
      call D2DMS(idd,idm,decs,decl,3)

      if(Jout.gt.0) then
         if((Jct(70)/2)*2.ne.Jct(70) .and. ABS(t-tsv).ge.0.3_10) then
            ltf  =60
            svt  =t
            svra =ra
            svdec=decl
            if(tlst.gt.0._10 .and. dt.ne.0._10) write(Jout,30) dra,ddec
   30       format('0MEAN APPARENT RATES (SEC/HR,"/HR):',t40,1p,
     .             2D15.7)
            if(jdsprv.ne.0 .and. Jds.eq.jdsprv) then
               call SUICID(
     .       '2ND OBSERVING PERIOD FOR SAME DAY, WARNING IN ANGOT ',-13)
               write(Jout,40)
   40          format(
     .       '0*** 2ND OBSERVING PERIOD FOR SAME DAY, WARNING IN ANGOT')
            endif
            jdsprv=Jds
         endif
         if(ltf.ge.59) then

c print headings
            lpg=lpg+1
            if((Jct(70)/2)*2.ne.Jct(70)) then
               write(Jout,50) Heding,Date,shlab,cepch
   50          format('1 PEP PREDICTION EPHEMERIS   ',18A4,4x,2a4/
     .          '0',12x,'GMT ',a4/
     .          4x,'MM DD   HR MN SH    EL      AZ',9x,'RA   (',a4,
     .          ')    DEC',11x,'DELAY',12x,'DOP',10x,'JULDATE'/)
               ltf=5
            else
               write(Jout,60) Heding,Date,lpg,Aplnt(Klap),Sitf,Freq
   60          format('1THEORETICAL OBSERVATIONS  ',18A4,4x,2A4,
     .          6x,'PAGE',i5/ '0 OBSERVATIONS OF ',a8,15x,
     .          'RECEIVING SITE=',a8,'  SENDING SITE=',a8,
     .          '  FREQUENCY=',1pd21.14,' C/S'/
     .'0 GRNWCH  JULIAN UT REC TIME UT1-UTC  AT-UTC  CT-AT RIGHT ASCENSI
     .ON DECLINATION   AZIMUTH  ELEVATION   TIME DELAY     DOPPLER SHIFT
     .'/
     . '   DATE   DAY NUM HR MIN SEC   SEC     SEC     SEC   HR MIN  SEC
     . DEG ''   ''''      DEG       DEG        SECONDS      CYCLES/SECON
     .D')
               ltf=7
            endif
         endif

c print line of ra&decl, az&el, td&dop
         if((Jct(70)/2)*2.eq.Jct(70)) then
            write(Jout,80) Imonth,Iday,Iyear,Jds,Ihr,Imin,Sec,Ututs,
     .       Atuts,Ctat,irh,irm,ras,idd,idm,decs,az,el,
     .       Result(1),Result(2)
   80       format(i3,2('/',i2),i8,2I3,f5.1,3F8.4,a3,i3,f8.4,
     .       1x,a3,i3,f7.3,2F10.5,f16.10,1pd17.9)
         else
            write(Jout,100) Imonth,Iday,Ihr,
     .       Imin,ishr,el,az,irh,irm,ras,idd,idm,decs,
     .       Result(1),Result(2),jd4
  100       format(i6,i3,i5,2I3,2F8.3,3x,a3,i3,f6.2,3x,a3,
     .       i3,f5.1,f17.8,f15.3,i12,'.5')
         endif
         tsv=t
         ltf=ltf+1
      endif
c
c***************  get send time angles
c
c           aberration correction
      do i=1,3
         zsat(i)=-Xsitep(i,2)-Rsitp(2)*Xemlsc(i+3,2)
      end do
      rzsat=SQRT(DOT(zsat,zsat))
c
c determine azimuth,elevation (deg)
      d  =Sitnrm(1,2)*Nutpr(3,1) + Sitnrm(2,2)*Nutpr(3,2)
     .      + Sitnrm(3,2)*Nutpr(3,3)
      qxn=DOT(zsat,Sitnrm(1,2))
      do i=1,3
         xmerid(i)=Nutpr(3,i)-d*Sitnrm(i,2)
         xop(i)   =zsat(i)-qxn*Sitnrm(i,2)
      end do
      qxn=qxn/rzsat
      elx=ASIN(qxn)/Convd
      call CROSS(xop,xmerid,xmerop)
      qxo2=DOT(xmerop,Sitnrm(1,2))
      qxo1=DOT(xmerid,xop)
      azx =ATAN2(qxo2,qxo1)/Convd
      if(azx.lt.0._10) azx=azx+360._10
c
c transform to true equinox and equator of date, if desired
      if(mod(Jct(70)/2,2).eq.0) then
         call CORCHN(xsat,zsat)
      else if(jdalt.ne.0) then
         call PRODCT(altprc,zsat,xsat,3,3,1)
      else
         do i=1,3
            xsat(i)=zsat(i)
         end do
      endif
c
c compute right ascension,declination      (hours,degrees)
      rsat2=DOT(xsat,xsat)
      rsat =SQRT(rsat2)
      declx=ASIN(xsat(3)/rsat)/Convd
      rax  =ATAN2(xsat(2),xsat(1))
      if(rax.lt.0._10) rax=rax+Twopi
      rax=rax/Convhs/3600._10
c
c* start=500
c convert to integer parts
      call D2DMS(irh,irm,ras,rax,4)
      call D2DMS(idd,idm,decs,declx,3)
c
c get send time
      sec1 = Sec
      sec1 = sec1+(Ihr*3600+Imin*60)-Result(1)

c for send time work tape
      if(senwkt.gt.0) then
         jjds  =Jds
         sfract=sec1/86400._10
         if(sfract.lt.0._10) then
            sfract=sfract+1._10
            jjds  =jjds-1
         endif
      endif

      if(sec1.lt.0._10) sec1=sec1+86400._10
      call S2DMS(jhr,jmin,secs,sec1,6)
c
c write out send time angles
      if(Jout.gt.0 .and. (Jct(70)/2)*2.eq.Jct(70)) then
         write(Jout,110) jhr,jmin,secs,irh,irm,ras,idd,idm,decs,azx,elx
  110    format(36x,a3,i3,f10.6,a3,i3,f8.4,1x,a3,i3,
     .    f7.3,2F10.5,' SEND TIME & ANGLES')
         ltf=ltf+1
      endif
      if(senwkt.gt.0) write(senwkt) jjds,sfract,rax,declx,azx,elx,
     . (Result(i),i=1,2),Jds,Ihr,Imin,ifix(Sec+0.001),ra,decl,az,el
      return

c  end of observation series
c  limitation  only one target/site/frequency combination per run for
c          bcd site tape = jct(45).gt.0
c          send time work tape = senwkt=jct(45)+10
c          and only if(jct(45).ge.50)   is senwkt written
c          otherwise senwkt=0
  200 if(marta.gt.1) goto 300
      itf=0
      return

c end of run, rewind everything
  300 continue
      if(svnplnt0.eq.9999) return
      svnplnt0=9999
      if(jct45.gt.0) then
         endfile jct45
         rewind jct45

c send time work tape
         if(senwkt.gt.0) then
            endfile senwkt
            rewind senwkt
         endif
      endif
      if(Jout.gt.0) then
         endfile Jout
         rewind Jout
      endif

      return
      end
