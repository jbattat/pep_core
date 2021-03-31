      subroutine RESTRT(itape,nplnt,ncentr,intp,jd1,jd0,jd2,cond,
     .                  con,con1,epsp,kp,numki,ki)

      implicit none
c
c M.E.Ash    June 1968    Subroutine RESTRT
c Setup restart of numerical integration for individual body.

c array dimensions
      include 'globdefs.inc'
c
c arguments
      integer*4 itape,jd1,jd0,jd2
      integer*2 nplnt,ncentr,intp,kp(u_nmprm)
      real*10 cond(6),con(u_nmbod-6),con1(12)
      real*4 epsp(6)
      integer*2 numki,ki(99)
c
c        common
c   Temporary storage (numerical integration quantities in this
c   labelled common are not setup until after return from
c   RESTRT).
      common/ADAMS/ Prmt(u_nmprm),Secp,Jdp1,Jdp2,Nprmx,Ncnmx,
     . Npage1,Iterat1,Intp1,Intp2,Ifiltr,Lnklvt,Nplnt1,Icnd,
     . Ihrp,Iminp,Kkp(100)
      real*10 Prmt,Secp
      integer*4 Jdp1,Jdp2,Nprmx,Ncnmx,Npage1,Iterat1,Intp1,Intp2,Ifiltr
      character*4 Lnklvt
      integer*2 Nplnt1,Kkp,Icnd,Ihrp,Iminp
c
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'rstart.inc'
      real*10 Tabpt(3,8),Dtabpt(3,2,i_mxeqn-1)
      equivalence (Tabpt,Rstrt),(Dtabpt,Rstrt(1,2))

c local
      real*4 eps6,eps1,eps2,xx
      integer*4 i,idir,ii,iorder,ivel,j,jd,
     .     jdir,jorder,k,kp99,lmesg,nrec,numprt,numpt1
      character*120 messag
      integer*2 id2
c
c see if itape to be written as well as read
      if(kp(99).ge.0) call INOUTV(itape)
c
c read first two records
      read(itape)
      kp99 = kp(99)
      eps1 = epsp(1)
      eps2 = epsp(2)
      eps6 = epsp(6)
      read(itape) Nplnt1,ncentr,Iparst,intp,Jdp1,jd0,Jdp2,Nprmx,Ncnmx,
     . cond,(con(i),i=1,ncnmx-6),con1,(Prmt(i),i=1,nprmx),epsp,
     . (kp(i),i=1,nprmx),Npage1,Iterat1,Icnd,Intp1,Intp2,
     . Ihrp,Iminp,Secp,Kkp,Ifiltr,Lnklvt,numki,(ki(i),i=1,numki)
      if(nprmx.gt.u_nmprm .or. ncnmx.gt.u_nmbod) call SUICID(
     . 'TOO MANY PARAMETERS ON TAPE, STOP IN RESTRT ',11)
      if(Kkp(70).ne.Jct(13)) call SUICID('INTEGRATION HAS WRONG REFERENC
     .E FRAME, STOP IN RESTRT   ',14)
      kp(99)  = kp99
      epsp(6) = eps6
      do i=ncnmx+1,u_nmbod
         con(i-6)=0._10
      end do
      do i=nprmx+1,u_nmprm
         kp(i)=-1
      end do
      call KP2KI(kp,numki,ki)
      nrec = 0
      if(Nplnt1.eq.nplnt) goto 400
c
c error message for planet number
      write(messag,100) Nplnt1,itape,nplnt
  100 format(' PLANET NUMBER',i3,' ON DATA SET',i3,
     .       ' IS NOT EQUAL TO INPUT PLANET NUMBER',i3,
     .       ', STOP IN RESTRT ')
      lmesg = 22
  200 call SUICID(messag,lmesg)
c
c check if tape direction is consistant
  400 xx   = (jd2 - jd1) + (eps2 - eps1)
      idir = 1
      if (xx .le. 0) idir = -1
      xx   = (jdp2 - jdp1) + (epsp(2) - epsp(1))
      jdir = 1
      if (xx .le. 0) jdir = -1
      if((idir*jdir).le.0 ) then
         write(messag,410) jd1,eps1,jd2,eps2,jdp1,epsp(1),
     .                  jdp2,epsp(2),nplnt
  410    format('INPUT (JD1,2)=(',i8,f4.3,',',i8,f4.3,
     .          ') INCONSISTENT WITH THOSE ON TAPE (',i8,f4.3,
     .          ',',i8,f4.3,') FOR BODY',i3,',STOP IN RESTRT')
         lmesg = 32
         goto 200
      endif
c
c read record on tape
  450 nrec = nrec + 1
      if(kp(88).gt.0) then
c for fixed-interval form, get only the first tabular point
         read(itape,err=500,end=700) jd,Frstrt,ivel,
     .    ((Rstrt(i,j),i=1,ivel),j=1,Iparst)
      else
c for variable-interval form, use only position and velocity orders.
c higher orders, if any, are read into storage for the next following
c partial derivative and then overwritten in turn
         read(itape,err=500,end=700) jd,Frstrt,iorder,
     .    jorder,numprt,numpt1,((Tabpt(i,j),i=1,3),j=1,iorder),
     .    (((Dtabpt(i,min0(j,3),k),i=1,3),j=1,jorder),k=1,numpt1),Hcc1
         ivel = iorder*3
         if(numprt.gt.0) ivel = min0(ivel,jorder*3)
      endif
c
c see if starting record has been reached
      if(Jdstrt.ne.1) then
         if(idir.le.0) then
            if(jd.lt.Jdstrt) goto 900
            if(jd.eq.Jdstrt .and. Frstrt.le.epsp(6)) goto 900
         else
            if(jd.gt.Jdstrt) goto 900
            if(jd.eq.Jdstrt .and. Frstrt.ge.epsp(6)) goto 900
         endif
      endif
      goto 450

  500 read(itape)
      call PAGCHK(60,1,0)
      write(Iout,600) itape,nrec
  600 format(' ERROR RECORD',i6,' SKIPPED ON DATA SET',i3,
     .       ' IN RESTRT')

      call SUICID('CANNOT SKIP ERROR RECORD, STOP IN RESTRT',10)
  700 backspace itape
      call PAGCHK(60,2,0)
      write(Iout,800) itape
  800 format('0END OF FILE REACHED ON RESTART EPHEMERIS TAPE',i3)
      nrec = nrec - 1
c
c check if both position and velocity are on tape
  900 if(ivel.ge.6) then
c
c re-position tape to receive integration output
         Jdstrt = jd
         call PAGCHK(60,3,0)
         write(Iout,950) itape,nplnt,nrec,Jdstrt,Frstrt
  950    format('-DATA SET',i3,' FOR BODY',i3,
     .          ' POSITIONED AT DATA RECORD',i5,
     .          ' FOR NUMERICAL INTEGRATION RESTART AT JD=',i7,
     .          ' FRACT=',0pd21.14)
         if(kp(99).lt.0) then
            rewind itape
         else
            backspace itape
         endif
         return
      else
         write(messag,1000) ivel,nrec,itape,jd,Frstrt,nplnt
 1000    format('IVEL=',i3,' ON RECORD',i5,' OF TAPE',i3,
     .          ' AT TIME',i8,d22.14,' FOR BODY',i3,
     .          ', STOP IN RESTRT')
         lmesg = 25
         goto 200
      endif
      end
