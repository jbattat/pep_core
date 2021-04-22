      subroutine OBSRD1(ntaps,jiabs1,jtape,mocpar)
 
      implicit none
c
c m.e.ash   feb 1968    subroutine obsrd1
c read and printout first two records of input observation library
c tape (compar is the calling routine)
c
c arguments
      integer*4 ntaps,jiabs1,jtape,mocpar

c           ntaps=saved value of tape sequence number
c           jiabs1=0,1 or 2 indicates which variable is used for logical
c           tape number
c           jtape=1 to 10 index for logical tape number variable

c array dimensions
      include 'globdefs.inc'
c
c        common
      include 'dtparm.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'mtrapx.inc'
      include 'obstap.inc'
c the following integers read from iabs1 indicate which partial
c derivatives are on iabs1
c mprm   solar system parameters in /param/
c mem    earth-moon barycenter initial conditions and earth parameters
c mmn    moon initial conditions and parameters
c mer    earth rotation initial conditions and parameters
c mmr    moon rotation initial conditions and parameters
c mpl    planet initial conditions and parameters
c mscrd  observing site coordinates
c meqn   equinox-equator-latitude corrections for optical obs series
c mphs   planetary phase corrections for optical observation series
c           interpolation y-vectors and temporary storage
      common/YVECT/ Prmt(u_nmprm),Ecnd(u_nmbod),Mcnd(u_nmbod),
     .        Ercnd(u_nmbod),Mrcnd(u_nmbod),Idem0,Idmn0,Ider0,Idmr0,
     .        Label(20),Jddtm(200),Dtm(600),Ast(200),Title
      real*10 Prmt,Ecnd,Mcnd,Ercnd,Mrcnd,cndx(u_nmbod,4)
      real*4    Dtm
      integer*4 Idem0,Idmn0,Ider0,Idmr0,Jddtm
      character*1 Ast
      integer*4 idd(4)
      character*88 Title
      character*4 Label,messag(21),lnklvm
      equivalence (Ecnd(1),cndx(1,1)),(Idem0,idd(1)),
     1       (messag(1),Label(1))
 
c local
      integer*4 i,j,jddtm0,jterat,mpage,nlabel,nmast,nrec,ntype
      integer*2 mumdt1
      character*8 word,worde/' ERROR  '/, wordn/' END    '/
      character*1 blank/' '/, astrik/'*'/
c
c to read or not to read iiabs1
      if(Iiabs1.le.0) return
c
c read first record of input observation library tape
      word = wordn
      nrec = 1
      read(Iiabs1,err=100,end=300) Title
      goto 700
  100 write(Iout,200) jiabs1,jtape,Iiabs1
  200 format(
     .'0**** ERROR ON RECORD 1 OF INPUT OBSERVATION LIBRARY DATA SET IOB
     .S', i1, '(', i2, ')=', i3, ' IN OBSRD1')
      read(Iiabs1) Title
 
c if system changes this second read of first record might not be
c needed to get data into storage in case of error
      goto 700
 
 
c
c error messages
  300 write(Intern,400) word,nrec,jiabs1,jtape,Iiabs1
  400 format(a6,'ON RECORD', i2,
     .       ' OF INPUT OBSERVATION LIBRARY DATA SET IOBS', i1, '(',
     .       i2, ')=', i3, ', STOP IN OBSRD1 ')
      rewind Intern
      read(Intern,500) messag
  500 format(21A4)
      call SUICID(messag,21)
 
  600 word = worde
      goto 300
 
c
c read second record of input observation library tape
  700 nrec = 2
 
c zero mdtx array to avoid problems in lvtbdy called by cmpar1
      call ZFILL(Mdtx,2*600)
      read(Iiabs1,err=600,end=300) Ntapa(3),mpage,jterat,Nprmo,Ncnmo,
     . Idem0,Idmn0,Ider0,Idmr0,(Prmt(i),i=1,Nprmo),(Ecnd(i),i=1,Ncnmo),
     . (Mcnd(i),i=1,Ncnmo),(Ercnd(i),i=1,Ncnmo),(Mrcnd(i),i=1,Ncnmo),
     . (Mprmx(i),i=1,Nprmo),(Memx(i),i=1,Ncnmo),(Mmnx(i),i=1,Ncnmo),
     . (Merx(i),i=1,Ncnmo),(Mmrx(i),i=1,Ncnmo),Mumdtx,mumdt1,
     . (Jddtm(i),i=1,mumdt1),(Dtm(i),i=1,mumdt1),
     . (Mdtx(i),i=1,mumdt1),nlabel,(Label(i),i=1,nlabel),
     . jddtm0,lnklvm,Msitcr,Msptcr
      if(Nprmo.gt.u_nmprm .or. Ncnmo.gt.u_nmbod) call SUICID(
     . 'BAD NUMBER OF PARAMETERS, STOP IN OBSRD1',10)
c in principle, should set default parameter values if the input tape has
c smaller set of defined parameters
      if(Nprmo.lt.u_nmprm) then
         do i=Nprmo+1,u_nmprm
            Mprmx(i)=0
         end do
      endif
      if(Ncnmo.lt.u_nmbod) then
         do i=Ncnmo+1,u_nmbod
            Memx(i)=0
            Mmnx(i)=0
            Merx(i)=0
            Mmrx(i)=0
         end do
      endif
      if(Mumdtx.eq.0) Mdtx(1) = 0
      if((Ntapa(3).lt.0) .and. (Ict(3).gt.1)) Ntapa(3) = -Ntapa(3)
      if(Msitcr.ne.6) Msitcr=3
      if(Msptcr.ne.6) Msptcr=3
c
c printout first two records of input observation library tape
      call NEWPG
      write(Iout,800) jiabs1,jtape,Iiabs1,Ntapa(3),jterat,mpage,
     .                 Title,lnklvm
  800 format(
     .'0INFORMATION ON FIRST TWO RECORDS OF INPUT OBSERVATION LIBRARY TA
     .PE IABS1=IOBS', i1, '(', i2, ')=', i2, 4x, ' NTAPE=', i3, 4x,
     .'ITERAT=', i3, 5x, 'PAGE=', i5/'0TITLE=', a88, 2x, '(LEVEL=', a4,
     .')'/'0PARTIAL DERIVATIVE CONTROL CONSTANTS ARE')
      write(Iout,900) Mprmx
  900 format(' SOLAR SYSTEM PARAMETERS MPRM=', 34I3, /33x, 33I3, /33x,
     .       33I3)
      write(Iout,1000) Idem0,Memx
 1000 format(' EARTH INIT.COND.PRMTRS JDEM0=', i7, ' MEM=', 30I3)
      write(Iout,1100) Idmn0,Mmnx
 1100 format(' MOON  INIT.COND.PRMTRS JDMN0=', i7, ' MMN=', 30I3)
      write(Iout,1200) Ider0,Merx
 1200 format(' EARTH ROT. I.C.&PRMTRS JDER0=', i7, ' MER=', 30I3)
      write(Iout,1300) Idmr0,Mmrx
 1300 format(' MOON  ROT. I.C.&PRMTRS JDMR0=', i7, ' MMR=', 30I3)
      Line = Line + 14
      if(Mumdtx.gt.0) then
         write(Iout,1350) Mumdtx,(Mdtx(i),i=1,Mumdtx)
 1350    format(' ET-UT2 PARAM. MUMDT=', i3, '  MDT=', 34I3/(33x,33I3))
         Line = Line + (Mumdtx - 2)/33
         if(Numdt.ne.0) then
            if(Ict(4).lt.1) then
c
c check consistency of et-ut2 tables
               if(Numdt.lt.Mumdtx) call SUICID(
     .' LENGTH OF ET-UT2 TABLE ON OBS.LIB.TAPE DOES NOT AGREE WITH INPUT
     . LENGTH, STOP IN OBSRD1', 22)
               nmast = 0
               do i = 1, Mumdtx
                  Ast(i) = blank
                  if(Jddt(i).ne.Jddtm(i)) then
                     Ast(i) = astrik
                     nmast  = nmast + 1
                  endif
               end do
               if(nmast.gt.0) then
                  write(Iout,1360) nmast,
     .                              (i,Jddt(i),Jddtm(i),Ast(i),
     .                              i=1, Numdt)
 1360             format('0 N     JDDT   JDDTM BAD=',
     .                   i3/(i4,'.',2I8,1x,a1))
                  call SUICID(
     .' JD IN ET-UT2 TABLE ON OBS.LIB.TAPE DOES NOT AGREE WITH INPUT JD,
     . STOP IN OBSRD1', 20)
               endif
               call DTCHCK(jddtm0,Mumdtx,mumdt1,Dtm,Mdtx)
            else
               Mumdtx = 0
               call ZFILL(Mdtx,2*600)
            endif
         endif
      else
         write(Iout,1400)
 1400    format(
     .    ' THERE IS NO ET-UT2 TABLE ON INPUT OBSERVATION LIBRARY TAPE')
      endif
c
c check for increasing tape sequence numbers
      if(Ict(8).le.0 .and. Ntapa(3).le.ntaps) then
         write(Iout,1450) Ntapa(3),ntaps
 1450    format('0NTAPA(3)=', i3, ' IS NOT GREATER THAN NTAPS(3)=', i3)
         call SUICID(
     .' TAPE SEQUENCE NUMBER FOR INPUT OBSERVATION LIBRARY TAPE NOT INCR
     .EASING, STOP IN OBSRD1 ', 22)
      endif
c
c check if observation tapes only are being read so that
c initial times and initial conditions are to be taken from
c observation library rather than body data set if adjustments
      if(Ict(1).gt.0) then
         if((Iterat.eq.1) .and. (Ict(80).ge.1) .and.
     .      (mocpar.eq.0)) then
            do j = 1, 4
               ntype = j
               call LIBCHK(ntype,cndx(1,j),idd(j))
            end do
         endif
      endif
c call for prmter to go here
c
      return
      end
