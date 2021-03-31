      subroutine SKYIN(in0,nstop,init)
 
      implicit none

c subr. skyin - j.f.chandler - 1983 jan
c read input stream for star catalog error maps and/or initialize
c
c arguments
      integer*4 in0,nstop
      logical*4 init

c array dimensions
      include 'globdefs.inc'

c common
      include 'inodta.inc'
      include 'skystf.inc'
c shared work area
      common/WRKCOM/ Sky(80),Ctlg,Lsky(80),N
      character*8 Ctlg
      real*10 Sky
      integer*4 N
      integer*2 Lsky

c local
      integer   i,j
      character*8 blank/'        '/
 
      namelist /SKYCOR/Sky, Ctlg, Lsky, N
 
      if(init) then
c
c initialize star catalog error models
         Numstr = 0
         do j = 1, u_mxsky
            Ctlgnm(j) = blank
         end do
         call ZFILL(Skycf,16*80*u_mxsky)
         call ZFILL(Lskycf,2*80*u_mxsky)
         call ZFILL(Nskycf,2*u_mxsky)
         return
      else
c
c read input stream
c
c set up namelist read
         write(in0,50)
   50    format(' &SKYCOR')
         call PEPTIC(In,Iout,in0,5,'SKY CORRECTION MODEL', nstop, 0)
 
         call ZFILL(Lsky,2*80)
         call ZFILL(Sky,16*80)
         Ctlg = blank
         N    = 0
 
         read(in0,SKYCOR)
         rewind in0
         if(Numstr.lt.u_mxsky) then
c
c copy input into permanent storage
            Numstr = Numstr + 1
            Ctlgnm(Numstr) = Ctlg
            if(N.gt.20) then
               write(Iout,60) N
   60          format(' *** TOO MANY SKY MAP COEFFICIENTS:', i5,
     .                ' > 20 ***')
               nstop = nstop + 1
               N     = 20
            endif
            N = 4*N
            Nskycf(Numstr) = N
            do i = 1, N
               Skycf(i,Numstr)  = Sky(i)
               Lskycf(i,Numstr) = Lsky(i)
            end do
            return
         else
            write(Iout,80) u_mxsky
   80       format(' *** MORE THAN', i3, ' SKY CORRECTION MODELS ***')
            nstop = nstop + 1
            return
         endif
      endif
      end
