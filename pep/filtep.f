      subroutine FILTEP(in0,nstop,init)
 
      implicit none
c
c     z. goldberg march 1980 subroutine filtep
c                            (from filtin 7/78)
c
c       subroutine to read in filter epochs and corresponding process
c       noise parameters.  invoked by command '*state' in input stream.
c
c       (the names of the process noise parameters are read in by
c        subroutine filtpn.)
c
c arguments
      integer   in0,nstop
      logical   init

c array dimensions
      include 'globdefs.inc'
c
c common
      include 'fcntrl.inc'
      include 'filnit.inc'
      include 'filtds.inc'
      include 'filtim.inc'
      real*10 x(20)
      equivalence (Pnprms,x)
      include 'inodta.inc'
      include 'param.inc'
 
c local
      integer*4 i
      real*10 aukm,fejd,fesec
      logical   rdw
      integer*2 fehr,femin,ficnd
      logical*4 flags(20)
      character*48 eptit/
     . 'FILTER EPOCHS & PROC NOISE PARMS    &NMLST5     '/
c
c namelist
      namelist/NMLST5/ fejd,fehr,femin,fesec,Pnprms,flags,x,ficnd
 
      rdw  = Fict(8).gt.0 .or. Jct(56).gt.1
      aukm = Aultsc*Ltvel
 
      if(init) goto 300
c
c
c check whether &nmlst5 should be used
      if(.not.rdw) then
c
c initialize for namelist read
         do i = 1,20
            flags(i) = .true.
         end do
         fehr  = 0
         femin = 0
         fesec = 0._10
         ficnd = 0
c
c spool, echo, and read &nmlst5
c spool  &nmlst5
         call PEPTIC(In,Iout,in0,11,eptit,nstop,0)
         read(in0,NMLST5)
         rewind in0
 
         Nepcht = Nepcht + 1
c
c check number of epochs read
         if(Nepcht.le.Maxfep) then
 
            Fep(Nepcht) = fejd + ((fehr*60+femin)*60 + fesec)/86400._10
 
c type icnd = -3 input, convert to au
            if(ficnd.eq.-3) Pnprms(Npnp-5) = Pnprms(Npnp-5)/aukm
c
c check epoch validity
            if(Nepcht.ne.1 .and. Fep(Nepcht).le.Fep(Nepcht-1)) then
               write(Iout,10) Fep(Nepcht)
   10          format(' *** FILTER EPOCH =',f15.5,
     .                ' NOT INCREASING, ERROR IN FILTEP ***')
               nstop = nstop + 1
            endif
            if(Nepcht.eq.Nepoch+2) write(Iout,20)
   20       format(
     .       ' *** WARNING. TOO MANY FILTER EPOCHS INPUT TO FILTEP ***')
c
c write wzero record
            write(Wzero) Fep(Nepcht),Pnprms,flags
            Filflg(1) = .true.
            return
         endif
c
c spool & echo &nmlst5, but don't read it ...
c if first call, fill fep from w data set
      else if(.not.Feprd) then
 
         write(Iout,50)
   50    format(
     .        ' *** WARNING -- *STATE INFORMATION WILL NOT BE USED ***')
         call PEPTIC(In,Iout,in0,11,eptit,nstop,0)
         Feprd = .true.
         goto 300
      else
         call PEPTIC(In,Iout,in0,11,eptit,nstop,0)
         return
      endif
  100 write(Iout,200)
  200 format(' *** NUMBER OF FILTER EPOCHS EXCEEDS MAXIMUM IN FILTEP')
      nstop = nstop + 1
      return
c
c         alternate sources for fep data if no *STATE input
c
c
c         make sure nepoch in range
  300 if(Nepoch.ge.Maxfep) goto 100
      Nepcht = Nepoch + 1
 
      if(rdw) then
c
c read fep from w data set
         read(Iconof) (Fep(i),i = 1,Nepcht)
         rewind Iconof
         Itrwnd(Iconof) = 0
         return
      else
c
c calculate epochs and fill fep
         do i = 1,Nepcht
            Fep(i) = timez + (i-1)*delta
         end do
         return
      endif
 
      end
