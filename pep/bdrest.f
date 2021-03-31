      subroutine BDREST

      implicit none


c*** start of declarations inserted by spag
      integer   i,idir,ii,ivel,j,jd,jdb1,jdb2,k,kb39,
     .          lmesg,mvel,nmoon1,nrec

c*** end of declarations inserted by spag


c
c Smith/Connolly   July 1968    Subroutine BDREST
c Read n-body tape to determine quantites for checkpoint restart
c of n-body numerical integration.
c

c array dimensions
      include 'globdefs.inc'

c common
      include 'bdctrl.inc'
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'rstart.inc'

c temporary storage
      common/YVECT/ Beta(6,10),Name(6,10),Npl1(10),
     .        Ncp1(10),Intb1(10),Jdbo(10),Kkbdy1(80)
      character*4 Name
      real*10 Beta,fract,stuff
      integer*2 nbdy1,Npl1,Ncp1,Intb1,Kkbdy1,nbody1
      integer*4 Jdbo
      character*120 messag
c
c see if n-body data set is to be written as well as read
      if(Kbdy(39).ge.0) call INOUTV(Ibody)
c
c read first two records of n-body tape
      read(Ibody)
      kb39 = Kbdy(39)
      read(Ibody) nbdy1,(Npl1(i),i=1,nbdy1),(Ncp1(i),i=1,nbdy1),
     . (Intb1(i),i=1,nbdy1),jdb1,jdb2,(Jdbo(i),i=1,nbdy1),
     . ((Beta(i,j),i=1,6),j=1,nbdy1),((Name(i,j),i=1,6),j=1,nbdy1),
     . nmoon1,nbody1,Intbdy,Jvlbdy,Epsbdy,Kbdy,i,i,
     . (stuff,i=1,nbdy1),(stuff,i=1,nbdy1),Kkbdy1
      Kbdy(39) = kb39
      nrec     = 0
c
c error message for number of bodies
      if(nbody1.eq.Nbody) goto 400
      write(messag,100) nbody1,Ibody,Nbody
  100 format('BODY COUNT',i3,' ON DATA SET',i3,
     .       ' IS NOT EQUAL TO INPUT BODY COUNT',i3,
     .       ', STOP IN BDREST  ')
      lmesg = 20
  200 call SUICID(messag,lmesg)
c
c error message for planet numbers
  400 do i = 1,Nbody
         if(Npl1(i).ne.Nplbdy(i)) then
            write(Iout,420) (Npl1(ii),Nplbdy(ii),ii = 1,Nbody)
  420       format(2I5)
            call SUICID(' PLANET NUMBERS DO NOT AGREE,STOP IN BDREST ',
     .                  11)
         endif
      end do
c
c check tape direction
      idir = ISIGN(1,Jdbdy2-Jdbdy1)
      if(idir*(jdb2-jdb1).le.0) then
         write(messag,430) Jdbdy1,Jdbdy2,Nbody,jdb1,jdb2
  430    format(' INPUT (JDBDY1,JDBDY2)=(',i7,',',i7,
     .          ') NOT CONSISTENT WITH THOSE ON',i3,'-BODY TAPE (',
     .          i7,',',i7,'), STOP IN BDREST')
         lmesg = 29
         goto 200
      endif
      if(Kkbdy1(70).ne.Jct(13)) call SUICID('INTEGRATION HAS WRONG REFER
     .ENCE FRAME, STOP IN BDREST   ',14)
c
c read data records
  450 nrec = nrec + 1
      read(Ibody,err=500,end=700) jd,fract,ivel,mvel,
     . (Rstrt(i,1),i=1,ivel),((stuff,i=1,ivel),j=1,9),
     . ((Rstrt(i,k),i=1,ivel),((stuff,i=1,ivel),j=1,4),k=2,nbody1)
c
c test if record before end of tape is reached for restart
      if(Jdstrt.ne.1) then
         if(idir.le.0) then
            if(jd.le.Jdstrt) goto 900
         else
            if(jd.ge.Jdstrt) goto 900
         endif
      endif
      goto 450

c error record on n-body tape
  500 read(Ibody)
      write(Iout,600) Ibody,nrec
  600 format(' ERROR RECORD',i6,' SKIPPED ON N-BODY DATA SET',i3,
     .       ' IN BDREST')
      call SUICID('CANNOT SKIP ERROR RECORD, STOP IN BDREST',10)
c
c end reached on n-body tape, assume this is restart point
  700 backspace Ibody
      write(Iout,800) Ibody
  800 format('0END OF FILE REACHED ON N-BODY RESTART TAPE',i3)
      nrec = nrec - 1
c
c correct record has been reached
  900 Jdstrt = jd
      write(Iout,1000) Nbody,Ibody,nrec,Jdstrt
 1000 format(//,i3,'-BODY DATA SET',i3,
     .       ' POSITIONED AT DATA RECORD',i6,
     .       ' FOR NUMERICAL INTEGRATION RESTART AT JD=',i7)
      write(Iout,1100) ((Rstrt(i,j),i=1,6),j = 1,nbody1)
 1100 format(1p,6D22.14)
c
c reposition n-body data set
      if(Kbdy(39).lt.0) then
         rewind Ibody
      else
         backspace Ibody
      endif
      return
      end
