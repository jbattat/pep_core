      subroutine XXRDBD(npxx, itpxx, klxx, jdx, intxx, idirxx, kxx,
     .                  xintx, nkixx, kixx, isbd, xxcom)
 
      implicit none

c subroutine xxrdbd - j.f.chandler - 1982 feb
c read body constants from disk, test for availability of n-body
c or s-body dataset.  if available, fill in controls.
c if not available, set up individual integration data set number.

c arguments
      real*10 xintx,xxcom(12)
      integer*4 itpxx,jdx(3),intxx,idirxx,isbd
      integer*2 npxx, klxx, kxx(100), nkixx, kixx(99)
c
c  jtest= index of body in n-body list (-index if s-body) or 0
c  npxx=  planet number of desired body
c  itpxx= flag set to 0 if n-body or -1 if s-body or tape number
c  klxx=  index of requested body within input list (-3 to numpln)
c  jdx=   time tags for individual records on integration tape
c  intxx= tabular interval (assumed integer if planet)
c  idirxx=direction of integration: +1 or -1
c  kxx=   integration controls
c  xintx= tabular interval (floating point)
c  nkixx= number of partials integration controls
c  kixx=  partials integration controls
c  isbd=  index for s-body tape setup. s-body not allowed if zero.
c  xxcom= array of body constants to be read from disk
c
c array dimensions
      include 'globdefs.inc'

c        common
      include 'bddta.inc'
      include 'b2dta.inc'
      include 'comdat.inc'
      include 'empcnd.inc'
      include 'namtim.inc'
      include 'plndta.inc'
      include 'yvectrd1.inc'

c external functions
      integer*4 JBDTST

c local variables
      integer*4 i,mjtest

c set up tape number and copy saved jd0 for xxrd1
      itpxx = Iplnt(klxx)
      Jdxx9 = Jdpl9(klxx)
c
c read body constants from disk
c skip to proper record
      do i = 1,klxx+3
         read(Iplcon)
      end do
      read(Iplcon) xxcom, Beps, kxx, Jdd1, Jdd2, Int1, Int1x, Int2x,
     .             Ihrx, Iminx, Secx, Kkxx, Npyy, Tconx,
     .             nkixx, (kixx(i), i = 1, nkixx)
      rewind Iplcon
c
c test for presence of body on n-body or s-body
      Jtest = 0
      if(npxx.lt.0) return
      if(kxx(88).eq.-8) Jdxx0 = xxcom(1)+0.5_10
      Jtest = JBDTST(npxx)
      if(isbd.le.0 .and. Jtest.lt.0) Jtest = 0
      if(Jtest.eq.0) return
c
c yes, body is present and available
c
c clear integration controls
      do i = 1, 100
         kxx(i) = 0
      end do
 
c provide for constant tabular interval
      kxx(88) = 1
      nkixx   = 8
      do i = 1, nkixx
         kixx(i) = 0
      end do
      do i = 1, 3
         jdx(i) = 0
      end do
      if(Jtest.lt.0) then
 
c body is on s-body tape
         mjtest    = -Jtest
         itpxx     = -1
         Tgo(isbd) = 0._10
         Jdpl0(klxx) = Tb20(mjtest)
         do i = 1, 6
            Pcond(i, klxx) = B2ta(i, mjtest)
         end do
         intxx  = Inb2(mjtest)
         idirxx = Ib2sgn
      else
 
c body is on n-body tape
         itpxx = 0
         Jdxx0 = Jdbd0(Jtest)
         do i = 1, 6
            Pcond(i, klxx) = Beta(i, Jtest)
         end do
         intxx  = Intb(Jtest)
         idirxx = Ibdsgn
      endif
      xintx = intxx
      if(intxx.le.0) xintx = 2._10**intxx
      return
      end
