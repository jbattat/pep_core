      subroutine UT1RED(jdtt,frt,ut1)

      implicit none

c
c R.King       Oct 1980  subroutine ut1red
c A. Whipple   Feb 1988  modified to use MERIT standard zonal
c                        tide model
c
c read ut1 values from an external data set
c output variable ut1 is a1-ut1

c parameters
      integer*4 jdtt
      real*10 frt,ut1

c array dimensions
      include 'globdefs.inc'
c
c commons
      include 'argfun.inc'
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'leon.inc'
      include 'obscrd.inc'

c external functions
      real*10 A1WWV,ARGMNT,ETUTF,UT2UT1,NUTRM

c local variables
      real*10 f2,s,t,tjd,units,xint
      integer i,int,it,itsav,j,jct28,jd1,jd2,jdt,jdt1,jdt2,
     .        jut1,kind,lap,mjd,n,nback,npr,npsave
      integer nr,nrbmin,nrec,nrecs,nvkeep,nvr,nvtot
      real*10 tab(4),y1(2),y2(2),kcf,kcm
      integer*4 iut(24)
      logical*4 nxtrp/.false./,init/.true./
      character*32 varfmt
      character*80 title
c
c         the external data set must have the following form
c     record 1:  title  (a80)
c     record 2:  information describing the table
c         columns
c          1-32  format of data entries  (a32)
c         33-34  kind of data in table  (i2)
c                  =1 ut1-utc
c                  =2 tai-ut1
c                  =3  a1-ut1
c         36-42  julian day of first value in table  (i7)
c         44-50  julian day of last value in table  (i7)
c         52-53  maximum number of values per record  (i2)
c         55-56  interval in days betweeen tabular values  (i2)
c         58-72  units of values in seconds  (e15.7)
c     records 3-end:  data entries - format is read in record 2, but
c                     must be of the form integer mjd folowed by up
c                     to 12 integer values of ut1.  if a record is
c                     short, the number of values in that record is
c                     given in columns 79-80 (i2).

      jdt   = jdtt
      nback = 0
      lap   = nrbmin
      if(Nk1.le.0) then
         if(Npage.ne.npsave) nxtrp = .false.
      endif
      goto 600

      entry UT1RD1
      nxtrp  = .false.
      init   = .true.
      nvtot  = 0
      itsav  = -9999
      jdt    = 0
      npsave = 0
c
c read and write header records
      jut1 = Jct(33)
      if(Line.gt.52) call NEWPG
      read(jut1,100) title
  100 format(a80)
      read(jut1,200) varfmt,kind,jdt1,jdt2,npr,int,units
  200 format(a32,i2,1x,i7,1x,i7,1x,i2,1x,i2,1x,e15.7)
      write(Iout,300) jut1,title,varfmt,kind,jdt1,jdt2,npr,int,
     .                 units
  300 format('0DATA ON FIRST TWO RECORDS OF UT1 DATA SET', i3/
     .       ' HEADING=', a80/ ' FORMAT=', a32, '  KIND=', i2,
     .       '  JDT1=', i7, '  JDT2=', i7, '  NPR=', i2, '  INT=', i2,
     .       '  UNITS=', 1pe15.6, ' (SEC)')
      Line   = Line + 4
      jdt2   = jdt2 - int
      xint   = int
      nrec   = 0
      nrbmin = 24/npr + 1
      Itrwnd(jut1) = 1
      if(npr.gt.12) call SUICID(
     .' NUMBER OF UT1 VALUES PER RECORD EXCEEDS ARRAY SIZE, STOP IN UT1R
     .ED ', 17)
c
c read  data record into storage
  400 do while( nvtot.le.(24-npr) )

c keep reading until (24/npr) records are in storage
         read(jut1,varfmt,end=700) mjd,(iut(nvtot+i),i = 1,npr),nvr
         nrec = nrec + 1
c jd1 and jd2 are the limits of usable values in storage
c and depend on the interpolation scheme
         if(nvtot.eq.0) jd1 = mjd + 2400001 + int
         if(nvr.eq.0) nvr   = npr
         nvtot = nvtot + nvr
         itsav = -9999
         jd2   = mjd + 2400001 + (nvr - 2)*int
         if(jd2.ge.jdt2) goto 500
      end do
  500 if(jdt.eq.0) then

c save the first date just in case header jdt1 is wrong
         jdt1 = jd1
         return
      endif
c
c is jd within range of table?
  600 if(jdt.lt.jdt1 .or. jdt.ge.jdt2) then
c
c extrapolation beyond table
         if(.not.(nxtrp)) then
            write(Iout,620) jdt
  620       format(i17,
     .   ' WARNING:  A1-UT1 BEING EXTRAPOLATED BEYOND END OF TABLE ***')
            Line = Line+1
            if(Mout.gt.0) write(Mout,620) jdt
            nxtrp  = .true.
            npsave = Npage
         endif
         ut1 = (ETUTF(jdt,frt) - 32.15_10) + UT2UT1(jdt,frt)
         goto 800
c
c is jd too close to jd1?
      else if(jdt.ge.jd1) then
c
c is jd too close to jd2?
         if(jdt.lt.jd2) then
c
c calculate interpolation times and value of tabular points
            t  = jdt - jd1
            t  = (t + frt)/xint
            it = t
            t  = t - it
            s  = 1.0_10 - t
            if(it.ne.itsav) then
               do i = 1, 4
                  j = it + i
                  tab(i) = iut(j)
               end do
c
c calculate interpolation y-vector
               do i = 1, 2
                  nr   = i + 1
                  f2   = .166666666666666666667_10*(tab(nr+1)+tab(nr-1))
                  y1(i)= 1.33333333333333333333_10*tab(nr) - f2
                  y2(i)= -0.33333333333333333333_10*tab(nr) + f2
               end do
               itsav = it
            endif
c
c second difference interpolation
            ut1 = (t*(y1(2)+t*t*y2(2)) + s*(y1(1)+s*s*y2(1)))*units
c
c convert table values to a1-ut1
            if(kind.eq.2) then
c table is tai-ut1
c a1-ut1 = tai-ut1   +   (a1-tai)
               ut1 = ut1 + 0.034390_10
            else if(kind.ne.3) then
c
c table is ut1-utc
c a1-ut1 = (a1-utc) - (ut1-utc)
               ut1 = A1WWV(jdt,frt) - ut1
            endif
c
c
c           add short-period (9-, 14-, 30-day) terms due to tidal effect
c           which may have been removed by smoothing
c             ref:  Yoder et al., J. Geophys. Res. 86, 881-891, 1981
            jct28 = Jct(28)
            if(mod(jct28/2,2).ne.0) then
               if(init) then
                  kcf = 0.94_10
                  kcm = 0.94_10
                  if(Jct(29).ne.0) then
                     if(Ercond(25).gt.0._10) kcf = Ercond(25)
                     if(Ercond(26).gt.0._10) kcm = Ercond(26)
                  endif
                  init=.false.
               endif
               tjd = jdt + frt - 0.5_10

c get angles in brown lunar theory (units are revolutions)
               call FUNARG(tjd)
c
c Compute the zonal tide corrections.  The coefficients are the
c MERIT standards (USNO Circ. No. 167, p a11-2).  The coefficients
c already include k/c=0.94; this is removed and the values from
c kcf and kcm are then used.

               Ut1par(1)=
     .   NUTRM( 1._10,0._10,2._10,2._10,2._10,+0.002_10,1)
     . + NUTRM( 2._10,0._10,2._10,0._10,1._10,+0.004_10,1)
     . + NUTRM( 2._10,0._10,2._10,0._10,2._10,+0.010_10,1)
     . + NUTRM( 0._10,0._10,2._10,2._10,1._10,+0.005_10,1)
     . + NUTRM( 0._10,0._10,2._10,2._10,2._10,+0.012_10,1)
     . + NUTRM( 1._10,0._10,2._10,0._10,0._10,+0.004_10,1)
     . + NUTRM( 1._10,0._10,2._10,0._10,1._10,+0.041_10,1)
     . + NUTRM( 1._10,0._10,2._10,0._10,2._10,+0.099_10,1)
     . + NUTRM( 3._10,0._10,0._10,0._10,0._10,+0.002_10,1)
     . + NUTRM(-1._10,0._10,2._10,2._10,1._10,+0.008_10,1)
     . + NUTRM(-1._10,0._10,2._10,2._10,2._10,+0.020_10,1)
     . + NUTRM( 1._10,0._10,0._10,2._10,0._10,+0.008_10,1)
     . + NUTRM( 2._10,0._10,2._10,-2._10,2._10,-0.002_10,1)
     . + NUTRM( 0._10,1._10,2._10,0._10,2._10,-0.003_10,1)
               Ut1par(1)= Ut1par(1)
     . + NUTRM( 0._10,0._10,2._10,0._10,0._10,+0.030_10,1)
     . + NUTRM( 0._10,0._10,2._10,0._10,1._10,+0.321_10,1)
     . + NUTRM( 0._10,0._10,2._10,0._10,2._10,+0.776_10,1)
     . + NUTRM( 2._10,0._10,0._10,0._10,-1._10,-0.002_10,1)
     . + NUTRM( 2._10,0._10,0._10,0._10,0._10,+0.034_10,1)
     . + NUTRM( 2._10,0._10,0._10,0._10,1._10,-0.002_10,1)
     . + NUTRM( 0._10,-1._10,2._10,0._10,2._10,+0.002_10,1)
     . + NUTRM( 0._10,0._10,0._10,2._10,-1._10,-0.005_10,1)
     . + NUTRM( 0._10,0._10,0._10,2._10,0._10,+0.073_10,1)
     . + NUTRM( 0._10,0._10,0._10,2._10,1._10,+0.005_10,1)
     . + NUTRM( 0._10,-1._10,0._10,2._10,0._10,+0.005_10,1)

               Ut1par(1) = Ut1par(1)*1E-3_10/0.94_10
               ut1= ut1 + kcf*Ut1par(1)

               Ut1par(2)=
     . + NUTRM( 1._10,0._10,2._10,-2._10,1._10,-0.005_10,1)
     . + NUTRM( 1._10,0._10,2._10,-2._10,2._10,-0.010_10,1)
     . + NUTRM( 1._10,1._10,0._10,0._10,0._10,-0.004_10,1)
     . + NUTRM(-1._10,0._10,2._10,0._10,0._10,-0.005_10,1)
     . + NUTRM(-1._10,0._10,2._10,0._10,1._10,-0.018_10,1)
     . + NUTRM(-1._10,0._10,2._10,0._10,2._10,-0.044_10,1)
     . + NUTRM( 1._10,0._10,0._10,0._10,-1._10,-0.053_10,1)
     . + NUTRM( 1._10,0._10,0._10,0._10,0._10,+0.826_10,1)
     . + NUTRM( 1._10,0._10,0._10,0._10,1._10,-0.054_10,1)
     . + NUTRM( 0._10,0._10,0._10,1._10,0._10,-0.005_10,1)
     . + NUTRM( 1._10,-1._10,0._10,0._10,0._10,+0.006_10,1)
     . + NUTRM(-1._10,0._10,0._10,2._10,-1._10,-0.012_10,1)
     . + NUTRM(-1._10,0._10,0._10,2._10,0._10,+0.182_10,1)
     . + NUTRM(-1._10,0._10,0._10,2._10,1._10,-0.013_10,1)
     . + NUTRM( 1._10,0._10,-2._10,2._10,-1._10,-0.002_10,1)
     . + NUTRM(-1._10,-1._10,0._10,2._10,0._10,+0.009_10,1)

               Ut1par(2) = Ut1par(2)*1E-3_10/0.94_10
               ut1= ut1 + kcm*Ut1par(2)
c
c analytical model terms added to a1-ut1 in function a1ut1f
            endif
            goto 800
         else

c if so, shift storage and read another record
            nvkeep = nvtot - npr
            do i = 1, nvkeep
               iut(i) = iut(npr + i)
            end do
            jd1   = jd1 + npr*int
            nvtot = nvtot - npr
            goto 400
         endif
      else

c if so, backspace and try again
         if(nback.gt.0 .and. nrec.ge.nrecs) lap = nrbmin + 1 +
     .       nrec - nrecs
         nrecs = nrec
         nback = nback + 1
         if(nback.ge.15) then
            write(Iout,640) jdt,jd1,jd2
  640       format('0***INFINITE LOOP IN UT1RED, JD=', i8, ' JD1,JD2=',
     .             2I10)
            call SUICID(' ERROR IN UT1RED', 4)
         endif
         n = (jd1 - jdt)/int/npr + lap
         if(n.gt.nrec) n = nrec
         do i = 1, n
            backspace jut1
         end do
         nrec  = nrec - n
         nvtot = 0
         goto 400
      endif

  700 call SUICID(' END OF FILE ON UT1 DATA SET, STOP IN UT1RED', 11)

  800 if(mod(Jct(6)/8192,2).eq.1) then
         if(Line.gt.56) call OBSPAG
         write(Iout,870) jdt,frt,ut1,kcf*Ut1par(1),kcm*Ut1par(2)
  870    format(' UT1RED: JD.F=',i8,f13.12,' UT1=',f12.8,' TID=',2f11.8)
         Line = Line+1
      endif

      return
      end
