      subroutine EVTRP(jd,fract,nvel,lcntl,icall,yv,x,body,ndim,
     .                 jdtb,frtb,p)
 
      implicit none

c        subr. evtrp - j.f.chandler - 1979 nov
c        derived from subr. sctrp - j.f.chandler, 1977 may
c        this is a control routine for interpolating a body from its
c        integration (or n-body or s-body) tape, using everett method
c        parameters:

      integer*4 jd,nvel,lcntl,icall,ndim
      real*10 fract
c data from integration tape
      real*10 body(6,ndim,1),frtb(3)
      integer*4 jdtb(3)
c output
      real*10 yv(7,3,6),x(6),p(4)

c          jd      integer value of the coord. time of the desired pt
c          fract   fractional part of the coordinate time of the point
c                  (ignored if jd=0 and nvel=-2)
c          nvel =0 positions determined
c                1 velocities determined
c               -1 velocities determined from position y-vectors
c               -2 accelerations determined from available y-vectors
c          lcntl=0 use results of last call if possible
c                1 do not use results (must recalculate y vectors)
c               -1 same as 0, except that nb1 is forced to 2
c                  (this option allows for a subsequent call for the same
c                   object for the same observation falling in the
c                   next-earlier tabular interval)
c              lcntl should, in principle, be set to 1 for the first call,
c              but this isn't really necessary because ntb1s is initialized
c              for each series
c          icall   indicates which is the calling routine
c                   1 - emtrp,  2 - mntrp,  3 - pltrp,  4 - sbtrp,
c                   5 - sctrp,  6 -(prtrp), 7 -(ertrp), 8 - sotrp,
c                   9 - ertsnt, 10 - sotrp for mercury using n-body,
c                   11-19 - sotrp for partials
c          yv      array for interpolation y-vectors
c          nb11    pointer to left-hand tabular point in yv
c                  (for objects that may be evaluated at two different
c                   times in the same observation, the later time must
c                   be evaluated first with lcntl=-1, and then the
c                   earlier time with lcntl=0, so that nb11 will be
c                   left over pointing to the correct tabular interval
c                   for the earlier time)
c              NOTE: this procedure fails if the ephemeris file is
c                    backwards in time and the two needed epochs of
c                    evaluation straddle the epoch of any record
c          x       output coordinate array
c          body    coordinates from integration tape
c          ndim    dimension of 'body' - max. # of partials + 1
c          jdtb    tabular julian dates
c          frtb    tabular fractions of a day
c          p       p-q array to be used
c
c  the storage of tabular points in 'body' is in the same time order as
c  the integration (indicated by 'idir') and is assumed to contain at
c  least 5 points per record for the three records tagged by jdtb(1-3),
c  frtb(1-3).  thus, in the case of 5/record the points would be
c  numbered in order of increasing time as
c     x              x              x
c     1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15  (idir=1)
c  or
c                 x              x              x
c    15,14,13,12,11,10, 9, 8, 7, 6, 5, 4, 3, 2, 1  (idir=-1)
c
c  where the x's mark the points corresponding to the time tags.  the
c  interpolation must have at least 5 tabular points bracketing the time
c  of interest on each side, and the appropriate subroutine is called
c  if necessary to read the correct records into storage.
c  bookkeeping is done via common 'trpcom', which is initialized
c  by subroutine 'trpnit' at the start of each observation series.
c          bint    directed tabular interval
c          bintx   absolute value of tabular interval
c          nb1     same as nb11 above
c          idirb   (or idir) time direction of integration
c          nvels=0 y-vectors based on positions
c                1 y-vectors based on positions and velocities
c               >1 y-vectors not set up
c          ntb1s   saved value of ntab1, or 9999 if y-vectors not set up
c          nbtrp   pointer to left bracketing tabular point (relative
c                  to 2nd or 3rd record, depending on idir)
c
c
c array dimensions
      include 'globdefs.inc'

c        commons
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'stats.inc'
      include 'tabval.inc'
      include 'trpcom.inc'

c local
      real*10 bbintx,dt
      integer*4 i,i4,idir,ifirst,io,j,jmx,jo,jy,lim1,lim2,n1,n2,n3,nb11,
     . nbn,nbo,ndmy,ngo,nmo,nov5,nprec,npto2,ntbss,nvela,nvelr
      integer*2 npl
c
c names of interpolators
      character*2 smsg,
     1  name(19)/'EM','MN','PL','SB','SC','PR','ER','SO','EO','S1',
     2           'P1','P2','P3','P4','P5','P6','P7','P8','P9'/
      character*8 sumsg(5)/'.. TAPE ','READING ','FAILURE,',' STOP IN',
     .          ' EVTRP  '/
      equivalence (sumsg,smsg)

c external functions
      real*10 D2TRPF,DTRPF,TERPF,D2TRP14,DTRP14,TERP14

      smsg   = name(icall)
      nb11   = Nb1(icall)
      bbintx = Bintx(icall)
      nvela  = iabs(nvel)
      jo     = 3*nvela
      if(nvel.eq.-2) jo = 0
      if(nvel.eq.-2 .and. jd.le.0) then
c
c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c
c           compute acceleration and store in x
         if(Nvels(icall).eq.0) then
 
c position interpolation vectors
            do jy = 1, 3
               if(icall.eq.2) then
                  x(jy) = D2TRP14(p,yv(1,nb11,jy),bbintx)
               else
                  x(jy) = D2TRPF(p,yv(1,nb11,jy),bbintx)
               endif
            end do
         else
 
c velocity interpolation vectors
            do jy = 1, 3
               if(icall.eq.2) then
                  x(jy) = DTRP14(p,yv(1,nb11,jy),bbintx)
               else
                  x(jy) = DTRPF(p,yv(1,nb11,jy),bbintx)
               endif
            end do
         endif
         goto 700
      endif
      Ntab1 = 2
      Ntab2 = 3
      nb11  = 2
      idir  = Idirb(icall)
      ngo   = 2
      if(idir.lt.0) ngo = 3
      ntbss = Ntb1s(icall)
      nvelr = max0(0,nvel)
      if(Nvels(icall).gt.1 .or. nvelr.gt.Nvels(icall)) then
         ntbss = 9999
         Nvels(icall) = nvelr
      endif
      jmx = 3*(Nvels(icall)+1)
      ifirst = 0
c
c           set limits within storage
c
c  npto2= half the number of points needed by the interpolator (5 or 7)
c  nprec= number of tabular points/record (5 or 8 or 10)
c  nov5 = nprec - npto2
c  lim1,lim2 = lower and upper limits for left-hand tabular
c              point for interpolation (relative to time tag of
c              2nd record if forward or 3rd if backward)

      if(icall.eq.2) then
c moon interpolation: 8 pts per record, 14-point interpolation
         nprec = 8
         npto2 = 7
         nov5  = 1
      else if(icall.eq.10) then
c mercury interpolation: 10 pts per record, 10-point interpolation
         nprec = 10
         npto2 = 5
         nov5  = 5
      else
c everything else: 5 pts per record, 10-point interpolation
         nprec = 5
         npto2 = 5
         nov5  = 0
      endif
      lim1 = -nov5
      lim2 = nov5+nprec-1
      if(idir.lt.0) then
c cannot always adjust lim1 forward because the tape readers aim at the
c entire time range of the middle record of a group of three
         if(lim1.lt.0) lim1  = lim1 + 1
         lim2  = lim2 + 1
      endif
 
c lcntl=0 means that nb11=1 may be used
      if(lcntl.eq.0 .and. (idir.gt.0.or.lim1.gt.-nov5)) lim1 = lim1 - 1
c
c test if routine should attempt to use saved values
      if(lcntl.eq.1 .or. jdtb(ngo).eq.0) goto 200
c
c
c*  start=1000
c           determine putative index of left-hand tabular point
c           if ok, then proceed to y-vectors
  100 dt  = ((jd-jdtb(ngo)) + fract - frtb(ngo))/bbintx
      nbn = dt + 19._10
      nbn = nbn - 19
      if(nbn.ge.lim1 .and. nbn.le.lim2) then
c
c*  start=1500
c determine p(1)..fractional time(time divided by interval)
c from left hand tabular point
         p(1) = dt - nbn
         nbn  = nbn + nov5
         if(ntbss.eq.9999) goto 400
 
c compute offset from last-used tabular point
         nbo = Nbtrp(icall)
         if(idir.lt.0) nbo = 2 + 3*nprec - nbo
         nmo = nbn - nbo
c
c some y-vectors may be retained, shift storage as needed
c note: nmo represents time offset; y-vectors
c are to be shifted by -nmo.  if that would shift all
c present storage off the end, then start fresh.
c
         if(nmo.lt.0) then
c point is backwards in time from saved region
            if(lcntl.eq.0) nmo = nmo + 1
            if(ntbss - nmo.gt.3) then
c cannot reuse anything
               ntbss = 9999
               goto 400
            endif
 
c at least one point can be re-used
            if(lcntl.eq.0) then
               nb11  = 1
               Ntab1 = 1
            endif
            Ntab2 = ntbss - nmo - 1
            n1    = 3
            n2    = Ntab2 + 1
            n3    = -1
            goto 300

         else if(nmo.eq.0) then
c no shifting is needed
            goto 600
         else
c
c point is forward in time from saved region
            if(nmo.gt.2) then
c cannot reuse anything
               ntbss = 9999
               goto 400
            endif
            Ntab1 = 4 - nmo
            n1    = 1
            n2    = Ntab1 - 1
            n3    = 1
            goto 300
         endif
      endif
c
c
c*  start=1100
c saved values cannot be re-used. set flag & read tape
  200 if(ifirst.eq.1) then
         write(Iout,205) jd,fract,lcntl,nvel,(jdtb(i),frtb(i),i=1,3),
     .    nb11,Ntb1s(icall),ntbss,Ntab1,Ntab2
  205    format(' TAPE READ FAILURE JD,FRACT,LCNTL,NVEL=',i8,f10.7,2i3/
     .    ' TAPE RECORDS AT',3(i8,f5.2)/
     .    ' NB11,NTB1S,NTBSS,NTAB1,NTAB2=',9i7)
         call SUICID(sumsg,10)
      endif
      if(mod(Jct(6)/1024,2).eq.1) then
         if(Line.gt.56) call OBSPAG
         write(Iout,270) smsg,smsg,jd,fract,nvel
  270    format(1x,a2,'TRP->',a2,'REED: JD.F=', i7, f13.12, ' NV=', i2)
         Line = Line + 1
      endif
      ntbss = 9999
      if(icall.eq.1) then
         call EMREED(jd,fract)
      else if(icall.eq.2) then
         call MNREED(jd)
      else if(icall.eq.3) then
         call PLREED(jd)
      else if(icall.eq.4) then
         call SBREED(jd,fract)
      else if(icall.eq.5) then
         call SCREED(jd,fract)
      else if(icall.eq.6) then
         call RTREED(jd,fract,ndmy,10)
      else if(icall.eq.7) then
         call RTREED(jd,fract,ndmy,3)
      else if(icall.ge.8 .and. icall.le.10) then
      else if(icall.ge.11 .and. icall.le.19) then
         npl=icall-10
         call SSREED(jd,npl,Iplss(npl))
      else
         call SUICID('BAD ICALL, STOP IN EVTRP', 6)
      endif
 
      ifirst = 1
      if(jd.gt.0) goto 100
      Ntb1s(icall) = ntbss
      return
c
c*  start=2000
c shift storage
  300 if(nmo.eq.0) goto 500
c n3 is -1 if nmo is negative, 1 if nmo is positive
      do i=n1,n2,n3
         io = i + nmo
         do j = 1, npto2
            do jy = 1, jmx
               yv(j,i,jy) = yv(j,io,jy)
            end do
         end do
      end do
c adjust ntbss after shifting
      if(nmo.gt.0) then
         ntbss = n1
      else
         ntbss = n2
      endif
c
c*  start=2500
c select tabular region according to nb11 (usually 2)
  400 Nbtrp(icall) = nbn + 2 - nb11
      if(idir.lt.0) Nbtrp(icall) = 2 + 3*nprec - Nbtrp(icall)
c
c determine y vectors
  500 if(Ntab1.le.Ntab2) then
         if(idir.lt.0 .and. Nbtrp(icall)-Ntab1.gt.3*nprec) then
            write(Iout,502) sumsg(1)
  502       format(a8,'BACKWARD, OBSERVATION STRADDLING RECORD EPOCH')
            jd=0
            Ntb1s(icall) = 9999
            return
         endif
         call YCOFT(yv,body(1,1,1),ndim,icall,jmx)
         if(Ntab1.lt.ntbss) ntbss = Ntab1
      endif
c
c reset ntb1s if any change (shifting and/or recomputing)
      Ntb1s(icall) = ntbss
c
c
c*  start=3000
c set up p vector and do interpolation.
  600 p(2) = p(1)**2
      p(3) = 1._10 - p(1)
      p(4) = p(3)**2
c
      do jy = 1, 3
         j = jy + jo
         if(nvel.ge.0) then
            if(icall.eq.2) then
               x(j) = TERP14(p,yv(1,nb11,j))
            else
               x(j) = TERPF(p,yv(1,nb11,j))
            endif
         else
            if(icall.eq.2) then
               x(j) = DTRP14(p,yv(1,nb11,jy),bbintx)
            else
               x(j) = DTRPF(p,yv(1,nb11,jy),bbintx)
            endif
         endif
      end do
      Nb1(icall) = nb11
c
c*  start=5000
c debug printout
  700 if(icall.ne.8 .and. icall.ne.10)
     . call TRPLST(jd,fract,nvel,lcntl,name(icall)//'TRP(EV)',x)
c
c* start=9900
      return
      end
