      subroutine SUICID(messag, nn)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, n, nn
 
c*** end of declarations inserted by spag
 
 
c
c m.e.ash   nov 1966     subroutine suicid to end it all
c
      character*4 messag(*)
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
 
      character*4 normal(2) /' NOR','MAL '/
      logical*4 nrmstp
      character*8
     1 eoj(10)/'*=* END ','PEP JOB ',' *=*=*=*',7*'=*=*=*=*'/,
     2 errmsg(10)/'*=* ERR-',9*' '/
 
      n = iabs(nn)
      nrmstp = (messag(1).eq.normal(1) .and. messag(2).eq.normal(2))
      if(nn.gt.0) call MOVEBL(messag, min0(72,4*nn), errmsg(2), 72)
c
c write completion message on auxillary data set
      if(Mout.gt.0) write(Mout, 50) (messag(i), i = 1, n)
   50 format(20A4)
c
c type completion message to operator
      call OPRMSG(messag, n)
 
      if(nn.lt.0) then
c
c warning only, print message and return
         call NEWPG
         write(Iout, 100) (messag(i), i = 1, n)
  100    format(//(1x,32A4))
         if(.not.nrmstp) write(Iout, 150)
  150    format('0**********   WARNING   **********'//)
         if(nrmstp) call TMEND
         Line = 60
         return
      else if(nn.eq.0) then
c
c
         write(Iout, 200)
  200    format('1SUICID CALLED WITH MESSAGE LENGTH.EQ.0, STOP')
         stop 20
      else
c
c print completion message and timer information
         write(Iout, 250) Heding, Date, Npage, (messag(i), i = 1, n)
  250    format('1PLANETARY EPHEMERIS PROGRAM TERMINATED   ', 18A4, 1x,
     .          2A4, ' PAGE', i5//(1x,32A4))
         if(.not.nrmstp) call PLOG(errmsg)
         if(nrmstp) call PLOG(eoj)
         call PLOG9
         if(.not.nrmstp) call SNPDMP
         call TMEND
c
c normal stop in main
         if(nrmstp) stop
 
c abnormal stop:  return non-zero user code
         stop 20
      endif
      end
