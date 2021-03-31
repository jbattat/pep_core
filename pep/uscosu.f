      subroutine USCOSU
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ione, iten, jc, jtemp, jtone, jtten
 
c*** end of declarations inserted by spag
 
 
c
c        r.b. goldstein, t.m. euabnks  april 1978
c     uscosu : use of corrections set up
c        called by prmred
c        applies jct(2) and checks its validity
c        breaks up jcal into jclone and jclten
c
c
c     use of jcal
c     jcal is a vector set in the pep input stream
c     jcal (i) is a two digit number
c     each digit of jcal is independent of the other
c     the first (ones) digit of jcal(i) refers
c     to the calculation of corrections
c     the second (tens) digit of jcal(i) refers to the use of correction
c     first digit =
c                  1 : do not calculate the ith correction
c                  2 : use standard logic  for correction calculation
c                  3 : calculate the correction
c
c     second digit =
c                   1 : do not use the ith correction
c                   2 : use standard logic for correction use
c                   3 : unconditionally use the correction
c
c
c     note : if several corrections are unconditionally
c     used,meaningless corrections may result if
c     different approximations to the same correction
c     are applied together - - - -caveat emptor
c
c     note : to use normal program logic on a
c     reduced set of corrections,use jcal(i) = 1#
c     to turn off the appropriate corrections
c            to use only certain corrections,and no other
c     use jcal(i) = 3# to unconditionally turn on corrections
c     if these corrections are not available,no corrections
c     are applied
c
c
c        common
      include 'fcntrl.inc'
      include 'inodta.inc'
      include 'prpgat.inc'
c
c     jct(2) is the jcal meta-default
c     if any digit of jcal = 0 then
c     that digit is replaced by the corresponding digit of jct(2)
c     this has to be done first
c
      if( Jct(2) .ge. 0 ) then
         ione = Jct(2)
         ione = mod(ione, 10)
         iten = Jct(2)/10
         do i = 1, 100
 
c replace zero codes
            jtemp = Jcal(i)
            jtone = mod(jtemp, 10)
            jtten = jtemp/10
            if( jtone .eq. 0 ) jtone = ione
            if( jtten .eq. 0 ) jtten = iten
            Jcal(i) = 10*jtten + jtone
 
c test for valid codes
            if( jtone .lt. 1 .or. jtone .gt. 3 .or. jtten .lt. 1 .or.
     .          jtten .gt. 3 ) then
               write(Iout, 10) i, Jcal(i)
   10          format(' JCAL(', i3, ')=', i10, ' OUT OF RANGE')
               Jcal(i) = 22
            endif
 
         end do
c
c break up jcal vector
c
         do i = 1, 100
            jc = Jcal(i)
            Jclone(i) = mod(jc, 10)
            Jclten(i) = jc/10
         end do
      endif
c
c
      return
      end
