      subroutine PRNPSR(nstop)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, j, nnspt
 
c*** end of declarations inserted by spag
 
 
c subr. prnpsr - j.f.chandler - 1984 jun
c print input pulsar parameters, if any

c arguments
      integer*4 nstop

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'psrstf.inc'
      include 'sptcrd.inc'
 
      character*4    star/'*** '/, ok/'PSR '/, psrsym
 
      if(Numpsr.gt.0) then
 
         if(Line.gt.50) call NEWPG
         call PAGSET('PULSAR PHASE MODELS ', -5)
         call PAGHED(0)
         do j = 1, Numpsr
            psrsym = ok
            if(Numspt.gt.0) then
               do i = 1, Numspt
                  if(Sptpsr(j).eq.Spot(i)) then
                     nnspt = i
                     if(Nsplnt(i).ge.0) then
                        call PAGCHK(60,1,1)
                        write(Iout,10) Nsplnt(i)
   10                   format(18x,'*** PULSAR INPUT ON PLANET', i3,
     .                         ', ERROR IN PRNPSR')
                        psrsym = star
                        nstop  = nstop + 1
                     endif
                     goto 40
                  endif
               end do
            endif
            call PAGCHK(60,1,1)
            write(Iout,20)
   20       format(5x,'**** NO SPOT INPUT FOR PULSAR, ERROR IN PRNPSR'
     .            )
            nnspt = 9999999
            nstop = nstop + 1
   40       call PAGCHK(60,7,1)
            write(Iout,60) nnspt,Sptpsr(j),psrsym,Jdpsr0(j),
     .                      Plspr(j),Ntypsr(j),
     .                      (i,Psrcn(i,j),i = 1,16)
   60       format('0SPOT', i4, '. (', a4, ') ', a4, 'JD0=', i8,
     .             '  APPRX. PERIOD=', 1pd22.15, '  NTYPE=',
     .             i2/4('  CON(',i2,')=',1pd22.15))
            write(Iout,80) (Lpsrcn(i,j),i = 1,16)
   80       format(' L=', (t10,25I3))
         end do
         return
      else
         call PAGCHK(60,2,0)
         write(Iout,100)
  100    format('0THERE ARE NO INPUT PULSARS')
         return
      endif
      end
