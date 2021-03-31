      subroutine WOBRED(jdtt,frt,xwob,ywob)

      implicit none


c*** start of declarations inserted by spag
      real*10 f2, frt, s, t, units, xint
      integer   i, int, it, itsav, j, jd1, jd2, jdt, jdt1, jdt2, jdtt,
     .          jwob, k, kind, lap, mjd, n, nback, npr, npsave
      integer   nr, nrbmin, nrec, nrecs, nvkeep, nvr, nvtot

c*** end of declarations inserted by spag


c
c r.king   oct 1980    subroutine wobred
c read wobble values from an external data set


      real*4    xwob, ywob

      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'obscrd.inc'

      real*10 tab(2,4),y1(2,2),y2(2,2)
      integer*4 iwb(2,12)
      logical*4 nxtrp/.false./
      character*32 varfmt
      character*80 title
c
c         the external data set must have the following form
c     record 1:  title  (a80)
c     record 2:  information describing the table
c         columns
c          1-32  format of data entries  (a32)
c         33-34  not used
c         36-42  julian day of first value in table  (i7)
c         44-50  julian day of last value in table  (i7)
c         52-53  maximum number of pairs of values per record  (i2)
c         55-56  interval in days betweeen tabular values  (i2)
c         58-72  units of values in arcseconds  (e15.7)
c     records 3-end:  data entries - format is read in record 2, but
c                     must be of the form integer mjd folowed by up
c                     to 6 pairs of integer values of xwob and ywob.
c                     if a record is short, the number of pairs of
c                     values in that record is given in columns 79-80.
c
      jdt   = jdtt
      nback = 0
      lap   = nrbmin
      if(Nk1.le.0) then
         if(Npage.ne.npsave) nxtrp = .false.
      endif
      goto 600

      entry WOBRD1
      nvtot  = 0
      itsav  = -9999
      jdt    = 0
      npsave = 0
c
c read and write header records
      jwob = Jct(33) + 1
      if(Line.gt.56) call NEWPG
      read(jwob,100) title
  100 format(a80)
      read(jwob,200) varfmt,kind,jdt1,jdt2,npr,int,units
  200 format(a32,i2,1x,i7,1x,i7,1x,i2,1x,i2,1x,e15.7)
      write(Iout,300) jwob,title,varfmt,kind,jdt1,jdt2,npr,int,
     .                 units
  300 format('0DATA ON FIRST TWO RECORDS OF WOBBLE DATA SET', i3/
     .       ' HEADING=', a80/ ' FORMAT=', a32, '  KIND=', i2,
     .       '  JDT1=', i7, '  JDT2=', i7, '  NPR=', i2, '  INT=', i2,
     .       '  UNITS=', 1pe15.6, ' (ARCSEC)')
      Line   = Line + 5
      jdt2   = jdt2 - int
      xint   = int
      nrec   = 0
      nrbmin = 12/npr + 1
      Itrwnd(jwob) = 1
      if(npr.gt.6) call SUICID(' NUMBER OF WOBBLE VALUES PER RECORD '//
     .    'EXCEEDS ARRAY SIZE, STOP IN WOBRED  ', 18)
c
c read  data record into storage
  400 do while( nvtot.le.(12-npr) )

c keep reading until (12/npr) records are in storage
         read(jwob,varfmt,end=700) mjd,
     .                                 ((iwb(i,nvtot+j),i=1,2),j = 1,
     .                                 npr), nvr
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
     .        ' WARNING: WOBBLE BEYOND END OF TABLE:  XWOB=YWOB=0. ***')
            if(Mout.gt.0) write(Mout,620) jdt
            nxtrp  = .true.
            npsave = Npage
         endif
         xwob = 0.
         ywob = 0.
         return
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
               do i = 1, 2
                  do j = 1, 4
                     k = it + j
                     tab(i,j) = iwb(i,k)
                  end do
               end do
c
c calculate interpolation y-vectors
               do i = 1, 2
                  do j = 1, 2
                     nr = j + 1
                     f2 = 0.166666666666666666667_10*(tab(i,nr+1)
     .                    + tab(i,nr-1))
                     y1(i,j) = 1.33333333333333333333_10*tab(i,nr)-f2
                     y2(i,j) = -0.33333333333333333333_10*tab(i,nr)+f2
                  end do
               end do
               itsav = it
            endif
c
c second difference interpolation
            xwob = (t*(y1(1,2)+t*t*y2(1,2)) + s*(y1(1,1)+s*s*y2(1,1)))
     .             *units
            ywob = (t*(y1(2,2)+t*t*y2(2,2)) + s*(y1(2,1)+s*s*y2(2,1)))
     .             *units

            return
         else

c if so, shift storage and read another record
            nvkeep = nvtot - npr
            do i = 1, 2
               do j = 1, nvkeep
                  iwb(i,j) = iwb(i,npr + j)
               end do
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
  640       format('0***INFINITE LOOP IN WOBRED, JD=', i8, ' JD1,JD2=',
     .             2I10)
            call SUICID('ERROR IN WOBRED ', 4)
         endif
         n = (jd1 - jdt)/int/npr + lap
         if(n.gt.nrec) n = nrec
         do i = 1, n
            backspace jwob
         end do
         nrec  = nrec - n
         nvtot = 0
         goto 400
      endif

  700 call SUICID(' END OF FILE ON WOBBLE DATA SET, STOP IN WOBRED ',
     . 12)

      return
      end
