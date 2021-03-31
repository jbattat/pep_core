      subroutine PRNOBS(nstop,iseq)
 
      implicit none

 
c           ash/forni  august 1967  subroutine prnobs
c           check for consistency and
c           printout data for equinox-equator corrections,
c                    data for planetary phase corrections, and
c                    data for observation series changes,
c                             error weighting and dummy observations.

c array dimensions
      include 'globdefs.inc'

c        common
      include 'dltflg.inc'
      include 'eqenox.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'namtim.inc'
      include 'obsdta.inc'
      include 'phase.inc'
      include 'rdbias.inc'
      include 'skystf.inc'
      include 'sptcrd.inc'
      include 'stcord.inc'
 
c shared external work area
      common/WRKCOM/ Freq,Freq2,Plnnam,Pl2nam,Psite(2,u_mxeqx),
     .        Psitnu(u_mxeqx)
      character*8 Plnnam,Pl2nam
      character*4 Psite
      real*10 Freq,Freq2
      integer*4 Psitnu
      character*8 pnames(2)
      equivalence (pnames,Plnnam)

c external functions
      integer*4 ITYPOB,NSCAN,LEG

c internal to prnobs only
      character*8 ctlg
      character*4 site1,series,site2,spot1,spot2
      real*4 erwgt(2),acctim,fdev
      character*4 rsite,ssite,seqerr,taperr,ctlerr
      integer*4 rsitnu,ssitnu
      integer*2 plnum
      integer*2 ncode,nqlnt,itime,nrewnd,nqlnt2
      character*8 cont(2)/'        ', ' (CONT).'/,
     .    earth/' EARTH '/, pstar/'  STAR  '/, pound8/'########'/,
     .   astrk8/'********'/
      character*4 pound,astrik
      equivalence (pound,pound8),(astrik,astrk8)
      character*4 blank/'    '/, amper/'&&&&'/
      integer*2 ihr(2),imin(2)
      integer*4 jd(2)
      real*4    sec(2),erobs(2)
      character*8 ignor(2)/'        ',' IGNORED'/
      integer   i,ie,il,intday,intsec,iseq,iseq1,j,j1,
     .          j2,jplchk,jtypob,k,kf,kl,ll,merphs,merpht
      integer   merrbt,mqlnt,n,n1,nerctl,nereqn,nerphs,nerpln,
     .          nerrbr,nerrbs,nerrbt,nerrcv,nerseq,nersnd,nerspt,
     .          nertap,nn1,nplchk,nplnum
      integer   nseq,nseqs,nspot,nspot2,nstop,ntape,ntaps
c
c
c printout constant biases for planetary radar observation series
c in labeled common /rdbias/
      nerrbr = 0
      nerrbs = 0
      nerrbt = 0
      merrbt = 0
      if(Numrbs.ne.0) then
         nplchk = -1
         jplchk = 0
         call PAGSET('CONSTANT BIASES FOR PLANETARY RADAR OBSERVATIONS',
     .               -12)
         if(Line.gt.53) call NEWPG
         call PAGHED(0)
         do k = 1,Numrbs
c
c search for receiving site
c find j such that rdbsit(1,k)= site(1,j)
            do j = 1,Numsit
               if(Rdbsit(1,k).eq.Site(1,j)) then
                  rsite  = Site(2,j)
                  rsitnu = j
                  goto 20
               endif
            end do
            rsite  = astrik
            rsitnu = 999999
            nerrbr = nerrbr + 1
c
c search for sending site
c find j such that rdbsit(2,k)= site(1,j)
   20       ssite  = blank
            ssitnu = 0
            if(Rdbsit(2,k).ne.blank) then
               do j = 1,Numsit
                  if(Rdbsit(2,k).eq.Site(1,j)) then
                     ssite  = Site(2,j)
                     ssitnu = j
                     goto 40
                  endif
                  end do
               ssite  = astrik
               ssitnu = 999999
               nerrbs = nerrbs + 1
            endif
c
c check if planet numbers are in order
   40       nplnum = Nplrbs(k)
            if(Nplrbs(k).eq.10) nplnum = -1
            if(nplnum.le.0) then
               if(nplnum.eq.-4) then
                  Plnnam = pstar
                  goto 80
               else
                  if(nplnum.lt.nplchk) goto 60
                  if(Nplrbs(k).eq.10) then
                     Plnnam = Aplnt(17)
                     goto 80
                  else if(Nplrbs(k).eq.0) then
                     Plnnam = Aplnt(18)
                     goto 80
                  endif
               endif
            endif
 
c find j such that nplrbs(k)= nplnt(j)
            do j = 1,u_mxpl
               if(Nplrbs(k).eq.Nplnt(j)) then
                  if(j.lt.jplchk) goto 60
                  jplchk = j
                  Plnnam = Aplnt(j)
                  goto 80
               endif
            end do
            do i = 1,2
               if(Lrbs(i,k).gt.0) goto 60
               end do
            merrbt = merrbt + 1
            Plnnam = pound8
            goto 80
   60       nerrbt = nerrbt + 1
            Plnnam = astrk8
   80       plnum  = Nplrbs(k)
 
c set for planet number check
            nplchk = nplnum
 
c write one line
            call PAGCHK(58,1,1)
            if(Line.le.4 .or. k.le.1) then
               write(Iout,90)
   90          format('0',5x,
     .' PLANET NO. RECV SITE NO. SEND SITE NO. SER. LRBS  BIAS1(SEC)  BI
     .AS2(C/S)')
               Line = Line + 2
            endif
            write(Iout,100) k,Plnnam,plnum,Rdbsit(1,k),rsite,
     .                       rsitnu,Rdbsit(2,k),ssite,ssitnu,
     .                       Rdbser(k),(Lrbs(j,k),j = 1,2),
     .                       (Rbias(j,k),j = 1,2)
  100       format(i4,'. ',a8,i3,2x,2A4,i4,2x,2A4,i4,1x,a4,
     .             1x,2I2,1x,1p,2E12.5)
            if(Mout.gt.0 .and. Plnnam.eq.astrk8)
     .          write (Mout,100) k,Plnnam,plnum,Rdbsit(1,k),rsite,
     .                       rsitnu,Rdbsit(2,k),ssite,ssitnu,
     .                       Rdbser(k),(Lrbs(j,k),j = 1,2),
     .                       (Rbias(j,k),j = 1,2)
         end do
         if(merrbt.ne.0) then
            call PAGCHK(60,1,1)
            write(Iout,120) merrbt
  120       format(6x,'########',i4,
     .' PLANETS WHICH ARE NOT INPUT, BUT BIASES NOT ADJUSTED, SO WARNING
     . ONLY')
         endif
         if(nerrbt.ne.0) then
            call PAGCHK(60,2,1)
            write(Iout,140) nerrbt
  140       format(6x,'********',i4,
     .       ' PLANETS WHICH  (1)ARE NOT INPUT AND BIAS IS ADJUSTED OR'/
     .       '  (2)ARE OUT OF ORDER, ERROR')
            if(Mout.gt.0) write(Mout,140) nerrbt
         endif
         if(nerrbr.ne.0) then
            call PAGCHK(60,1,1)
            write(Iout,160) nerrbr
  160       format(22x,'********',i4,
     .             ' RECEIVING SITES WHICH ARE NOT INPUT, ERROR')
            if(Mout.gt.0) write(Mout,160) nerrbr
         endif
         if(nerrbs.ne.0) then
            call PAGCHK(60,1,1)
            write(Iout,180) nerrbs
  180       format(22x,'********',i4,
     .             ' SENDING SITES WHICH ARE NOT INPUT, ERROR')
            if(Mout.gt.0) write(Mout,180) nerrbs
         endif
      else
         call PAGCHK(60,2,0)
         write(Iout,200)
  200    format(
     .'0THERE IS NO INPUT DATA FOR CONSTANT BIASES FOR PLANETARY RADAR O
     .BSERVATION SERIES')
      endif
c
c
c printout data in labeled commons  /eqenox/ and /phase/  checking against
c labeled commons  /namtim/ and  /stcord/
      nereqn = 0
      if(Numeqn.ne.0) then
c
c print equinox-equator corrections in a table
c check for errors in site name and set up arrays psite and
c psitnu for printing
         do k = 1,Numeqn
            Psite(1,k) = Eqnsit(k)
 
c find  j  such that  eqnsit(k)= site(1,j)
            do j = 1,Numsit
               if(Eqnsit(k).eq.Site(1,j)) then
                  Psite(2,k) = Site(2,j)
                  Psitnu(k)   = j
                  goto 250
               endif
               end do
            Psite(2,k) = astrik
            Psitnu(k)   = 999999
            nereqn = nereqn + 1
  250    end do
 
         call PAGSET(
     .      'EQUINOX-EQUATOR CORRECTIONS FOR OPTICAL OBSERVATIONS',-13)
         if(Line.gt.53) call NEWPG
         call PAGHED(0)
         kf = 1
         kl = 2*(55 - Line)
         if(Line.gt.54 .or. kl.gt.Numeqn) kl = Numeqn
         do while( .true. )
            ll = (kl - kf)/2 + 3
            call PAGCHK(59,ll,1)
            write(Iout,260)
  260       format(/2(8x,
     .     'SITE   NO. SER. LEQN   DEQUINOX    DEQUATOR    DLATITUDE  '
     .     ))
            write(Iout,280) (k,(Psite(i,k),i=1,2),Psitnu(k),
     .                       Eqnser(k),(Leqn(i,k),i=1,3),Denox(k),
     .                       Dequat(k),Dlat(k),k = kf,kl)
  280       format((2(i4,'.',1x,2A4,i4,1x,a4,3I2,1p,3E12.4,1x)))
            if(kl.eq.Numeqn) then
               if(nereqn.ne.0) then
                  call PAGCHK(60,1,1)
                  write(Iout,290) nereqn
  290             format(10x,'********',i4,
     .                  ' OBSERVING SITES WHICH ARE NOT INPUT, WARNING'
     .                  ,9x,'********')
               endif
               goto 400
            else
               kf = kl + 1
               kl = Numeqn
            endif
         end do
      else
         call PAGCHK(60,2,0)
         write(Iout,300)
  300    format(
     .'0THERE IS NO INPUT DATA FOR EQUINOX-EQUATOR CORRECTIONS FOR OPTIC
     .AL OBSERVATION SERIES')
      endif
 
  400 nerphs = 0
      merphs = 0
      merpht = 0
      if(Numphs.ne.0) then
c
c print planetary phase corrections
c initialize for check of planet number
         nplchk = -1
         jplchk = 0
         call PAGSET(
     .      'PLANETARY PHASE CORRECTIONS FOR OPTICAL OBSERVATIONS',-13)
         if(Line.gt.55) call NEWPG
         call PAGHED(0)
         do k = 1,Numphs
c search for site
c find j such that phsit(k)= site(1,j)
            do j = 1,Numsit
               if(Phsit(k).eq.Site(1,j)) then
                  rsite  = Site(2,j)
                  rsitnu = j
                  goto 420
               endif
            end do
            rsite  = astrik
            rsitnu = 999999
            nerphs = nerphs + 1
c check for error in planet name and number
c check if planet numbers are in order
  420       nplnum = Nplphs(k)
            if(Nplphs(k).eq.10) nplnum = -1
            if(nplnum.le.0) then
               if(nplnum.eq.-4) then
                  Plnnam = pstar
                  goto 460
               else
                  if(nplnum.lt.nplchk) goto 440
                  if(Nplphs(k).eq.10) then
                     Plnnam = Aplnt(17)
                     goto 460
                  else if(Nplphs(k).eq.0) then
                     Plnnam = Aplnt(18)
                     goto 460
                  endif
               endif
            endif
 
c find  j  such that nplphs(k)= nplnt(j)
            do j = 1,Numpln
               if(Nplphs(k).eq.Nplnt(j)) then
                  if(j.lt.jplchk) goto 440
                  jplchk = j
                  Plnnam = Aplnt(j)
                  goto 460
               endif
               end do
            do i = 1,9
               if(Lphs(i,k).gt.0) goto 440
               end do
            merpht = merpht + 1
            Plnnam = pound8
            goto 460
  440       merphs = merphs + 1
            Plnnam = astrk8
  460       plnum  = Nplphs(k)
 
c set for planet number check
            nplchk = nplnum
 
c write one line
            call PAGCHK(58,1,1)
            if(Line.le.4 .or. k.le.1) then
               write(Iout,470)
  470          format('0',5x,' PLANET NO.   SITE   NO. SER.',8x,
     .         'LPHS',14x,'APHASE(1)   APHASE(2)   APHASE(3) APHASE',
     .         '(4&7) APHASE(5&8) APHASE(6&9)')
               Line = Line + 2
            endif
            il = Ncphs(k)
            ie = il
            if(ie.gt.6) ie = 6
            write(Iout,480) k,Plnnam,plnum,Phsit(k),rsite,rsitnu,
     .                       Phser(k),(Lphs(i,k),i = 1,9),
     .                       (Aphase(i,k),i = 1,ie)
  480       format(i4,'.',1x,a8,i3,1x,2A4,i4,1x,a4,9I2,5x,
     .             1p,6E12.5)
            if(Mout.gt.0 .and. Plnnam.eq.astrk8)
     .          write (Mout,480) k,Plnnam,plnum,Phsit(k),rsite,
     .                            rsitnu,Phser(k)
            if(il.gt.6) then
               write(Iout,510) (Aphase(i,k),i = 7,il)
  510          format(94x,1p,3E12.5)
               Line = Line + 1
            endif
         end do
 
         if(merpht.ne.0) then
            call PAGCHK(60,1,1)
            write(Iout,520) merpht
  520       format(6x,'########',i4,
     .' PLANETS WHICH ARE NOT INPUT, BUT PHASES NOT ADJUSTED, SO WARNING
     . ONLY')
         endif
         if(merphs.ne.0) then
            call PAGCHK(60,2,1)
            write(Iout,540) merphs
  540       format(6x,'********',i4,
     .      ' PLANETS WHICH  (1)ARE NOT INPUT AND PHASE IS ADJUSTED OR'/
     .      '  (2)ARE OUT OF ORDER, ERROR')
            if(Mout.gt.0) write(Mout,540) merphs
         endif
         if(nerphs.ne.0) then
            call PAGCHK(60,1,1)
            write(Iout,560) nerphs
  560       format(22x,'********',i4,
     .             ' OBSERVING SITES WHICH ARE NOT INPUT, WARNING')
         endif
      else
         call PAGCHK(60,2,0)
         write(Iout,600)
  600    format(
     .'0THERE IS NO INPUT DATA FOR PLANETARY PHASE CORRECTIONS FOR OPTIC
     .AL OBSERVATION SERIES')
      endif
c
c printout for alterations in error weighting and dummy observations
      nertap = 0
      nerseq = 0
      nerpln = 0
      nerspt = 0
      nerrcv = 0
      nersnd = 0
      nerctl = 0
      n1     = 1
      ntaps  = -99999
      nseqs  = -99999
 
c print parameters for global deletions
      call PAGCHK(60,3,0)
      write(Iout,499) Tdlt0,Tdlton,Tdltof
  499 FORMAT('0PARAMETERS FOR GLOBAL INPUT DATA DELETIONS'/
     . ' TDLT0=',F15.6,'  TDLTON=',F15.6,'  TDLTOF=',F15.6)
c
c read start of observation series
  700 read(Iobcon,end=1800) ncode,nqlnt,site1,series,site2,
     .                         spot1,(erwgt(i),i = 1,2),acctim,
     .                         itime,fdev,Freq,nrewnd,ntape,nseq,
     .                         nqlnt2,spot2,Freq2,ctlg,Gncode,
     .                         Gnplt1,Gsite1,Gseres,Gsite2,Gspot,
     .                         Gerwgt,Gacctm,Gitime,Gfdev,Gfreq,
     .                         Gnrwnd,Gnplt2,Gspot2,Gfreq2,Gctlg
      if(ncode.le.0 .and. nseq.eq.0) goto 1800
      jtypob = ITYPOB(ncode)
      if(Line.ge.55 .or. n1.ne.2) then
         if(Line.gt.48) call NEWPG
         write(Iout,750) cont(n1)
  750    format('-INPUT DATA FOR ALTERATIONS IN ERROR WEIGHTINGS,',
     . 'DELETIONS FROM OBSERVATION SERIES AND THE GENERATION OF DUMMY '
     . ,'OBSERVATIONS',a8/
     .'0NTP NSEQ NCDF PLANET NPLN SITE1 NSITE SERIES SITE2 NSITE SPOT NS
     .PT ERROR WEIGHTS  ACCTIM ITIME FDEV    FREQUENCY        CTLG')
         Line = Line + 5
         n1   = 2
      endif
c
c search for receiving site
      rsite  = blank
      ssite  = blank
      rsitnu = 0
      ssitnu = 0
 
c no site at all for transit/occultation observations
      if(jtypob.ne.3) then
         if(Gsite1) then
            if(Numsit.gt.0) then
               do i = 1,Numsit
                  if(site1.eq.Site(1,i)) then
                     rsitnu = i
                     rsite  = Site(2,i)
                     goto 800
                  endif
               end do
            endif
            if(nrewnd.ge.0 .and. nrewnd.le.1) then
               nerrcv = nerrcv + 1
               rsite  = astrik
               rsitnu = 999999
            endif
         endif
c
c search for sending site
c no second site for optical observations
  800    if(jtypob.ne.2 .and. jtypob.ne.3) then
            if(Gsite2) then
               if(Numsit.gt.0) then
                  do i = 1,Numsit
                     if(site2.eq.Site(1,i)) then
                        ssitnu = i
                        ssite  = Site(2,i)
                        goto 900
                     endif
                  end do
               endif
               if(nrewnd.ge.0 .and. nrewnd.le.1) then
                  nersnd = nersnd + 1
                  ssite  = astrik
                  ssitnu = 999999
               endif
            endif
         endif
      endif
c
c search for spot
  900 nspot = 0
      if(spot1.ne.blank) then
         if(spot1.ne.amper .and. spot1.ne.pound) then
            if(Numspt.ne.0) then
               do i = 1,Numspt
                  if(spot1.eq.Spot(i)) then
                     nspot = i
                     if(nqlnt.ne.Nsplnt(i)) then
                        nstop = nstop + 1
                        call PAGCHK(60,1,0)
                        write(Iout,910) nqlnt,Nsplnt(i),spot1
  910                   format(' ******** BELOW PLANET NUMBER',i3,
     .                         ' DOES NOT MATCH',i3,' FOR SPOT ',a4,
     .                         ',ERROR ********')
                        if(Mout.gt.0) write(Mout,910) nqlnt,
     .                     Nsplnt(i),spot1
                     endif
                     goto 1000
                  endif
               end do
c
c spot not found, see if it is a body name
c spot name gives 2nd body for undifferenced observables
c compare from 1st non-blank character in each
               j1  = NSCAN(spot1,4,blank)
               nn1 = 4 - j1
            endif
            if(jtypob.eq.3 .or. jtypob.eq.5) then
               if(Numpln.gt.0) then
                  do i = 1,Numpln
                     j2 = NSCAN(Aplnt(i),8,blank)
                     if(j2.ge.0) then
                        n = min0(nn1,8 - j2)
                        if( LEG(n,j1+1,spot1,j2+1,Aplnt(i)) .eq. 0 )
     .                      then
                           nspot = -i
                           goto 1100
                        endif
                     endif
                  end do
               endif
            endif
            nerspt = nerspt + 1
            nspot  = 999999
         endif
      endif
c
c search for second spot
 1000 if(ncode.gt.20) then
         nspot2 = 0
         if(spot2.ne.blank) then
            do i = 1,Numspt
               if(spot2.eq.Spot(i)) then
                  nspot2 = i
                  if(nqlnt2.ne.Nsplnt(i)) then
                     nstop = nstop + 1
                     call PAGCHK(60,1,0)
                     write(Iout,910) nqlnt2,Nsplnt(i),spot2
                     if(Mout.gt.0) write(Mout,910) nqlnt2,
     .                  Nsplnt(i),spot2
                  endif
                  goto 1100
               endif
            end do
            nerspt = nerspt + 1
            nspot2 = 999999
         endif
      endif
c
c search for observed bodies
 1100 mqlnt = nqlnt
      j     = 1
      if(.not.Gnplt1) then
         pnames(j) = ignor(1)
         goto 1300
      endif
 1200 if(mqlnt.lt.0) then
         pnames(j) = pstar
      else if(mqlnt.eq.0) then
         pnames(j) = Aplnt(18)
      else if(mqlnt.ne.10) then
         if(mqlnt.ne.3) then
            do i = 1,u_mxpl
               if(Nplnt(i).eq.mqlnt) then
                  pnames(j) = Aplnt(i)
                  goto 1300
               endif
            end do
            nerpln    = nerpln + 1
            pnames(j) = astrk8
         else
            pnames(j) = earth
         endif
      else
         pnames(j) = Aplnt(17)
      endif
 1300 do while( ncode.gt.20 )
         if(j.gt.1) goto 1400
         j     = 2
         mqlnt = nqlnt2
         if(Gnplt2) goto 1200
         pnames(j) = ignor(1)
      end do
c
c search for decreasing tape number & non-increasing sequence number
 1400 taperr = blank
      seqerr = blank
      if(ntape.lt.ntaps) then
         nertap = nertap + 1
         taperr = pound
      else if(ntape.ne.ntaps) then
         goto 1500
      endif
      if(nseq.le.nseqs) then
         nerseq = nerseq + 1
         seqerr = amper
      endif
 1500 ntaps = ntape
      nseqs = nseq
c
c search for catalog name
      ctlerr = blank
      if(Gctlg) then
         if(Numstr.gt.0) then
            do i = 1,Numstr
               if(ctlg.eq.Ctlgnm(i)) goto 1600
            end do
         endif
         ctlerr = astrik
         nerctl = nerctl + 1
      endif
c
c printout start of observation series
 1600 write(Iout,1700) ntape,taperr,nseq,seqerr,ncode,Plnnam,
     .                  nqlnt,site1,rsite,rsitnu,series,site2,
     .                  ssite,ssitnu,spot1,nspot,
     .                  (erwgt(i),i = 1,2),acctim,itime,fdev,
     .                  Freq,ctlg,ctlerr
 1700 format(/i4,a1,i4,a1,i3,1x,a8,i3,1x,2A4,i4,1x,a4,2x,
     .       2A4,i4,1x,a4,i4,1p,3E8.1,i2,0pf9.4,1pd19.12,1x,
     .       a8,a1)
      Line = Line + 2
      if(ncode.gt.20) then
         write(Iout,1750) Pl2nam,nqlnt2,spot2,nspot2,Freq2
 1750    format(14x,a8,i3,33x,a4,i4,35x,1pd19.12)
         Line = Line + 1
      endif
      do while( .true. )
c
c read information about data internal to observation series
         read(Iobcon,end=1900) (jd(i),ihr(i),imin(i),sec(i),i=1,2),
     .    (erobs(i),i = 1,2),intday,intsec
         if(jd(1).le.0) goto 700
         if(Line.gt.58) then
            call NEWPG
            write(Iout,750) cont(2)
            Line = Line + 5
         endif
         Line = Line + 1
         if(ntape.le.0) then
            write(Iout,1760) (erobs(i),i = 1,2),
     .       (jd(i),ihr(i),imin(i),sec(i),i=1,2),intday,intsec
 1760       format(' ERRORS ',1pe8.1,',',1pe8.1,
     .             ' USED FOR DUMMY OBS FROM JD',i8,i3,'H',i3,'M',
     .             0pf8.4,'S TO JD',i8,i3,'H',i3,'M',f8.4,
     .             'S EVERY',i3,' DAYS',i6,' SEC')
         else
            write(Iout,1780) (erobs(i),i = 1,2),
     .       (jd(i),ihr(i),imin(i),sec(i),i = 1,2)
 1780       format(' ERROR WEIGHTS ',1pe8.1,',',1pe8.1,
     .             ' USED FROM JD',i8,i3,'H',i3,'M',0pf8.4,
     .             'S TO JD',i8,i3,'H',i3,'M',0pf8.4,'S')
         endif
      end do
c
c write completion message for alterations in error weighting, etc.
 1800 if(n1.le.1) then
         call PAGCHK(60,2,0)
         write(Iout,1850)
 1850    format(
     .'0THERE IS NO INPUT DATA FOR ALTERATIONS IN ERROR WEIGHTINGS, DELE
     .TIONS FROM OBSERVATION SERIES OR THE GENERATION OF DUMMY OBS.')
      endif
 1900 rewind Iobcon
c
c error message for alterations in error weighting and dummy obs
      iseq1 = 1
      if(iseq.gt.0) iseq1 = 2
      if(nertap.gt.0) then
         if(Mout.gt.0) write(Mout,1950) nertap,ignor(iseq1)
         call PAGCHK(60,1,0)
         write(Iout,1950) nertap,ignor(iseq1)
 1950    format(4x,'#',i4,' TAPE NUMBERS DECREASING,ERROR',a8)
      endif
      if(nerseq.gt.0) then
         if(Mout.gt.0) write(Mout,2000) nerseq,ignor(iseq1)
         call PAGCHK(60,1,0)
         write(Iout,2000) nerseq,ignor(iseq1)
 2000    format(9x,'&',i4,' SEQUENCE NUMBERS NOT INCREASING, ERROR',a8)
      endif
      if(nerpln.gt.0) then
         if(Mout.gt.0) write(Mout,2050) nerpln
         call PAGCHK(60,1,0)
         write(Iout,2050) nerpln
 2050    format(14x,'********',i4,' PLANETS NOT INPUT, ERROR')
      endif
      if(nerrcv.gt.0) then
         if(Mout.gt.0) write(Mout,2100) nerrcv
         call PAGCHK(60,1,0)
         write(Iout,2100) nerrcv
 2100    format(30x,'********',i4,
     .          ' RECEIVING SITES NOT INPUT, ERROR')
      endif
      if(nersnd.gt.0) then
         if(Mout.gt.0) write(Mout,2150) nersnd
         call PAGCHK(60,1,0)
         write(Iout,2150) nersnd
 2150    format(49x,'********',i4,' SENDING SITES NOT INPUT, ERROR')
      endif
      if(nerspt.gt.0) then
         call PAGCHK(60,1,0)
         write(Iout,2200) nerspt
 2200    format(62x,'****',i4,' SPOTS NOT INPUT, ERROR')
         if(Mout.gt.0) write(Mout,2200) nerspt
      endif
      if(nerctl.gt.0) then
         call PAGCHK(60,1,0)
         write(Iout,2250) nerctl
 2250    format(i4,' CATALOGS NOT INPUT, ERROR',t130,'*')
         if(Mout.gt.0) write(Mout,2250) nerctl
      endif
 
      nstop = nstop + merphs + nerpln + nerspt + nerrcv + nersnd +
     .        nerrbr + nerrbs + nerrbt + nerctl
      if(iseq.le.0) nstop = nstop + nertap + nerseq
      return
      end
