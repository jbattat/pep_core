      subroutine XXRD1(lice, npxx, itpxx, klxx, jdxx1, jdxx2, iparx,
     .                 maxxx, intxx, idirxx, kxx, xintx, frxx, nkixx,
     .                 kixx)
 
      implicit none
 
c subroutine XXRD1 - J.F.Chandler - 1980 Aug
c based on code from PLRD1, etc.
c read first two records of planet tape

c arguments
      real*10 xintx,frxx(2)
      integer*4 itpxx,jdxx1,jdxx2,iparx,maxxx,intxx,idirxx
      integer*2 lice,npxx,klxx,kxx(100),nkixx,kixx(99)

c further arguments of XXRDCK
      integer*4 lparx,jdsw
c
c     lice =0 printout of data on first two records of planet tape
c     lice =1 no such printout
c     npxx =  planet number to be read
c     itpxx=  input integration data set
c     klxx =  index into iplnt, nplnt, etc. of desired body
c     jdxx1=  start date of integration
c     iparx=  number of partials on tape
c     maxxx=  maximum allowed partials
c     intxx=  integration step size indicator
c     idirxx= time direction on tape
c     kxx  =  integration partial controls
c     xintx=  tabular interval
c     frxx =  fractions for jdxx1 and jdxx2
c     nkixx=  number of partials integration controls 'kixx'
c     kixx =  partials integration controls
c     lparx=  returned number of i.c. partials + 1
c     jdsw =  flag, 1=> overwrite jdpl0(klxx) from tape, 0=> do not
c
c array dimensions
      include 'globdefs.inc'

c        common
      include 'empcnd.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'lcntrl.inc'
      include 'namtim.inc'
      include 'yvectrd1.inc'
      character*4 messag(20)
      equivalence (Gmess, messag)
 
      character*8    tmess(10)/'  NPLNT(', '**)=****', ', BUT DA',
     .    'TA SET  ',
     .          '**  HAS ', 'PLANET =', '***, STO', 'P IN XXR',
     .          'D1.     ', ' '/
      character*8    eerr/'ERROR AT'/,  eend/' END ON '/,
     .               recd/' RECORD '/, ipard/' IPAR = '/,
     .               ikk70/' KK(70)='/ 
 
      real*10 beps2,frx3,tspan
      integer*4 i,i4,j,kij,kijmod,ll,lpli,nrec
      integer*2 nparitr,nparit1,nparpl(20)
 
      do i = 1, 10
         Gmess(i) = tmess(i)
      end do
      if(Jdxx9.lt.0) beps2 = Beps(2)
c
c read first two records of integration data set
      nrec = 1
      read(itpxx, err=100, end=700) Tpname, Title
      goto 300
  100 write(Iout, 200) itpxx, klxx, npxx
  200 format('0**** ERROR ON RECORD 1 OF DATA SET', i3,
     .       ' FOR BODY WITH NPLNT(', i2, ')=', i3, ' ****'/)
      Line = Line + 3
      read(itpxx) Tpname, Title
c if system changes this second read of first record might not be
c needed to get data into storage in case of error
  300 nrec = 2
      Itrwnd(itpxx) = 1
      read(itpxx,err=600,end=700) L1,L2,Ipar1,Int1,jdxx1,Jdxx0,jdxx2,
     . Nprmx,Ncnmx,(Pcond(i,klxx),i=1,6),(Cn1x(i),i=1,Ncnmx+6),
     . (Cn2x(i),i=1,Nprmx),Beps,(kxx(i),i=1,Nprmx),M1,M2,Icnx,Int1x,
     . Int2x,Ihrx,Iminx,Secx,Kkxx,Ifltrx,Levelx,nkixx,
     . (kixx(i),i=1,nkixx),nparitr,nparit1,(nparpl(i),i=1,nparit1)
      if(Nprmx.gt.u_nmprm .or. Ncnmx.gt.u_nmbod) call SUICID(
     . 'BAD NUMBER OF PARAMETERS, STOP IN XXRD1 ',10)
      iparx = Ipar1
      intxx = Int1
      xintx = Int1
      if(intxx.le.0) xintx = 2._10**intxx
 
c convert old format k vector if necessary
      call KP2KI(kxx, nkixx, kixx)
c
c print out first two records of planet data set
      if(lice.le.0) then
         call NEWPG
         write(Iout, 320) itpxx, Tpname, klxx, L1, L2, Ipar1, Ifltrx
  320    format('0DATA ON FIRST TWO RECORDS OF DATA SET', i3,
     .          '   BODY= ', a8, 3x, 'NPLNT(', i2, ')=', i3, 4x,
     .          'NCENTR=', i3, 5x, 'IPAR=', i3, 5x, 'IFILTR=', i2)
         write(Iout, 340) Title, M1, M2, Levelx, jdxx1, Jdxx0,
     .                    jdxx2, Int1, Icnx
  340    format(' TITLE=', 9A8, 1x, a8, ' PAGE=', i5, ' ITERAT=',
     .          i2, ' LEVEL=', a4/5x, 'JD1=', i7, 10x, 'JD0=', i7,
     .          10x, 'JD2=', i7, 10x, 'INT=', i3, 5x, 'ICND=', i3)
         write(Iout, 360) (i, Pcond(i,klxx), i = 1, 6),
     .                    (i, Cn1x(i), i = 1,Ncnmx+6)
  360    format(4(3x,'COND(',i1,')=',1pd22.15)/
     .          2(3x,'COND(',i1,')=',1pd22.15),
     .          2(3x,'CON(',i2,')=',1pd22.15)/
     .          4(3x,'CON(',i2,')=',1pd22.15))
         write(Iout, 380) (i, Cn2x(i), i = 1, Nprmx)
  380    format(4('  PRM(',i3,')=',1pd22.15))
         write(Iout, 400) (i, kxx(i), i = 1, Nprmx)
  400    format(10(2x,'K(',i3,')=',i4))
         ll = (nkixx - 8)/20 + 1
         write(Iout, 420) nkixx, (kixx(i), i = 1, nkixx)
  420    format(' NUMKI=', i3, '  KI=', 7I2, (t30,20I5))
         write(Iout, 440) (i, Beps(i), i = 1, 6)
  440    format(6('   EPS(',i1,')=',1pe12.5))
         if(nparitr.gt.0) then
            write(Iout,450) nparitr,(nparpl(i),i=1,nparitr)
  450       format(' INDIRECT PARTIALS ITERATED FOR',i3,' BODIES:',
     .       (20i3))
            ll=ll+1+(nparitr-1)/20
         endif
         Line = 52 + ll
      endif
c
c copy adjustable parameters that affect the motion
      i=7
      lpli = lpl(i,klxx)
      j=8
      do while (j.le.nkixx)
         kij=kixx(j)
         kijmod=mod(kij,100)
         if(kij.lt.0) then
            do while (lpli.gt.0 .and. lpli+kij.lt.0 .and. i.lt.30)
               i=i+1
               lpli = lpl(i,klxx)
            end do
            if(lpli+kij.eq.0) pcond(lpli+6,klxx)=Cn1x(lpli)
         else if(kij.gt.100 .and. kijmod.eq.31) then
            j=j+2
         else if(kij.gt.100 .and. (kijmod.eq.41 .or. kijmod.eq.51))
     .       then
            j=j+4
         endif
         j=j+1
      end do
c
c set up interval quantities
      if(Jdxx9.lt.0) then
         jdxx2   = Jdd2
         Beps(2) = beps2
      endif
      frx3 = 2._10*xintx
      if(npxx.lt.10 .or. kxx(88).gt.0) frx3 = 5._10*xintx
      if(npxx.eq.10) frx3 = 4._10
      tspan = jdxx2 - jdxx1
      if(tspan + Beps(2) - Beps(1) .ge. 0._10) then
 
c forwards
         idirxx  = 1
         frxx(1) = jdxx1
         frxx(1) = frxx(1) + Beps(1) + frx3
         frxx(2) = jdxx2
         frxx(2) = frxx(2) + Beps(2) - frx3
      else
 
c backwards in time
         idirxx  = -1
         frxx(1) = jdxx2
         frxx(1) = frxx(2) + Beps(2) + frx3
         frxx(2) = jdxx1
         frxx(2) = frxx(1) + Beps(1) - frx3
      endif
      jdxx1   = frxx(1)
      frxx(1) = frxx(1) - jdxx1
      jdxx2   = frxx(2)
      frxx(2) = frxx(2) - jdxx2
c test for error conditions

      if(L1.ne.npxx) then
c wrong body number
         i4 = L1
         call EBCDI(i4, messag(13), 3)

      else if(iparx.gt.maxxx) then
c too many partials
         Gmess(6) = ipard
         i4 = iparx
         call EBCDI(i4, messag(13), 3)

      else if(Kkxx(70).ne.Jct(13)) then
c wrong reference frame
         Gmess(6) = ikk70
         i4 = Kkxx(70)
         call EBCDI(i4, messag(13), 3)

      else
         return
      endif
c
c write error message, stop program

c encode body number and data set number into message
  500 i4 = klxx
      call EBCDI(i4, messag(3), 2)
      i4 = npxx
      call EBCDI(i4, messag(4), 4)
      call EBCDI(itpxx, messag(9), 2)
      call SUICID(messag, 20)
 
c tape error
  600 Gmess(6) = eerr
      goto 800
 
c end of file
  700 Gmess(6) = eend
 
c encode record number into message
  800 Gmess(10) = Gmess(9)
      Gmess(9)  = Gmess(8)
      Gmess(8)  = Gmess(7)
      Gmess(7)  = recd
      i4 = nrec
      call EBCDI(i4, messag(15), 3)
      goto 500
 
c-----------------------------------------------------------------------
      entry XXRDCK(lparx, kixx, klxx, jdsw)
c           count partials w.r.t. initial conditions on tape and
c           set jd0 if requested
c           note: jd0 from tape is saved in yvect common as jdxx0
c
c           count i.c. partials
      lparx = 1
      if(kixx(1).ne.0) then
         do i = 2, 7
            if(kixx(i).ge.0) lparx = lparx + 1
         end do
      endif
c
c see if body is to be reintegrated and set jd0 to 0 or
c tape value accordingly
      if(jdsw.gt.0) then
         Jdpl0(klxx) = Jdxx0
         do i = 1, 6
            if(Lpl(i,klxx).gt.0) return
         end do
 
c no i.c.'s being adjusted, don't reintegrate
         Jdpl0(klxx) = 0
      endif
      return
      end
