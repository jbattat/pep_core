      subroutine EATZDL(active,wetz,dryz)
 
      implicit none
c
c     r.king      march 1978      subroutine eatzdl
c     calculate the theoretical zenith delay of the earth's
c     atmosphere using input meterological data or monthly
c     average values for the observing site(s)
c
c
      logical   active
      real*4    wetz(2),dryz(2)

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'number.inc'
      include 'obscrd.inc'
      equivalence (Dstf(5),freq2)
      include 'param.inc'
      include 'sitcrd.inc'
c
c quantitites internal to this routine
      real*10 c2lat,freq2,ht
      integer   i,is,jdatm,numsit
      real*4    gfact,radcof(2),lamm2,optcof(2),e,t(2),p(2),h(2),
     .          atparm(2)
      logical*4 smlsav
      integer*4 ict23
c
c---------------------------------------------------------------------
c
      numsit = 2
      if(Nsite2.eq.0) numsit = 1
c
c           for static atmosphere model, calculate average monthly
c           value of zenith delay at both sites, assumed constant
c           over series
c
      if(active) then
c
c for active atmosphere model, calculate correction
c for variation in local gravity and wave length
         do i = 1,numsit
            if(Nk1.lt.0) then
               c2lat     = Cnrm(i)*Cnrm(i) - Snrm(i)*Snrm(i)
               ht        = Shgt(i)
               gfact     = 1. + .0026*c2lat + .00028*ht
               radcof(i) = .002277*gfact/Ltvel/1.E3
 
c correct optical group delay also for wavelength
               lamm2     = (1.E-9_10*Freq/Ltvel)**2
               optcof(i) = 0.39406*(173.3 + lamm2)/((173.3-lamm2)**2)
     .                     *gfact/Ltvel/1.E3
            endif
c
c           read in meterological data from save vector
c     units of temp, pressure, and relative humidity are deg kelvin,
c     millibar, and fraction of unity
c     if there are no values for the pressure and temperature on the
c     obslib tape,dummy values must be made up
            is = 0
            if(i.eq.2) is = 3
            smlsav = (Numsav.lt.is+43)
            t(i)= Save(is+41)
            if(smlsav .or. t(i).lt.150. .or. t(i).gt.350.) t(i)= 273.15
            p(i)= Save(is+42)
            if(smlsav .or. p(i).lt.500. .or. p(i).gt.1100.) p(i)= 1000.
            h(i)= Save(is+43)
            if(smlsav .or. h(i).lt.0. .or. h(i).gt.1.) h(i)= 0.
            if(Freq.lt.5.E12_10) then
c
c zenith delay for radio frequencies
               dryz(i) = radcof(i)*p(i)
               e = h(i)*6.108*EXP(19.85*(1-273.15/t(i)))
               wetz(i) = radcof(i)*(1255./t(i) + .05)*e
            else
c
c zenith delay for laser ranging
               ict23 = Ict(23)
c     Should the old or new zenith delay be used?
               if (MOD(ict23/2,2).eq.0) then 
c     use the old zenith delay model
                  dryz(i) = optcof(i)*p(i)
                  e = h(i)*6.108*EXP(19.85*(1-273.15/t(i)))
                  wetz(i) = optcof(i)*.06*e
               else
c     use the new zenith delay model of Mendes and Pavlis, 2004
                  call EATZMP(i,t(i),p(i),h(i),wetz(i),dryz(i))
               endif
            endif
         end do
      else
         if(Nk1.ge.0) then
            if(Jd.eq.jdatm) return
         endif
         call ATM(Sitf,Imonth,Iday,wetz,dryz)
         jdatm = Jd
      endif
c     end of loop i=1,numsit
c
c
c-----------optional printout of atmospheric data-----------------------
c
c
c     if ict(24) is not = 0,then print atparm for both sites
      if(Ict(24).ne.0) then
 
c find ratio of zenith delay to standard (goldstone aver.) value
         do i = 1,numsit
            atparm(i) = (dryz(i) + wetz(i))/7.088236E-9
         end do
         if(numsit.eq.1) atparm(2) = 0.
         write(Iout,50) atparm
   50    format(1x,'ATMOSPHERE PARAMETERS:  SITE1=',1pe14.7,
     .          '   SITE2=',1pe14.7)
         Line = Line + 1
 
         if(active) then
c
c if ict(24) <= -2, then print out the temp,pressure,and humidity
c that prevails for each site
            if(Ict(24).le.-2) then
               write(Iout,60) (i,t(i),p(i),h(i),i = 1,numsit)
   60          format(' SITE',i1,': TEMP=',1pe14.7,'  PRESSURE=',
     .                e14.7,'  HUMIDITY=',e14.7)
               Line = Line + numsit
            endif
         endif
      endif
 
      return
      end
