      subroutine GETSOL(jdobs)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, iflag, inrec, irec, j, jrec
 
c*** end of declarations inserted by spag
 
 
      real*10 jdobs
c
c        r. goldstein...march 1975...subroutine to allow prdict
c             to use filter solutions appropriate to current
c             observationn being processed
c
c        the code compares the jd of the current observation with the
c        jd span of the solution currently in core to decide if a
c        readin of new solutions are required.  a readin of solutions is
c        always done the first time the routine is called.
c
      include 'fcntrl.inc'
      include 'filtda.inc'
      include 'filtim.inc'
      include 'inodta.inc'
      include 'rtside.inc'
 
      real*10 jda, jdb
 
      data iflag, inrec/0, 0/
 
      if( iflag .eq. 0 ) then
 
         iflag = 1
      else
c
c check current jdobs with current jda and jdb
c
         if((jdobs .ge. jda) .and. (jdobs .le. jdb) ) return
      endif
c
c calculate correct record to read in
c
c
      do i = 1, Nepoch
         if( jdobs .lt. Fep(i+1) ) then
            if( jdobs .ge. Fep(i) ) then
               irec = i
               go to 100
            endif
         endif
      end do
  100 if( jdobs .gt. Fep(Nepoch+1) ) then
c
c last record is needed
c
         if( inrec .eq. Nepoch ) return
         irec = Nepoch
      else if( jdobs .lt. Fep(1) ) then
c
c first record is needed
c
         if( inrec .eq. 1 ) return
         irec = 1
      else if( Lfilt(irec) .eq. 0 ) then
c
c find the nearest epoch for which a solution exists
c
         do i = 1, Nepoch
            jrec = irec - i
            if( jrec .ge. 1 ) then
               if( Lfilt(jrec) .ne. 0 ) then
                  irec = jrec
                  go to 200
               else
                  jrec = irec + i
                  if( jrec .le. Nepoch ) then
                     if( Lfilt(jrec) .ne. 0 ) then
                        irec = jrec
                        go to 200
                     endif
                  endif
               endif
            endif
         end do
c
c error handling section
c
         write(Iout, 150) irec, jdobs, timez, delta
  150    format(' IREC=', i10/' JDOBS=', f20.10/' TIMEZ=',
     .          f20.10/' DELTA=', f20.10)
         call SUICID('ERROR IN GETSOL ', 4)
         return
      endif
c
c read it and return
c
  200 inrec = irec
      read(Kfile, rec = irec) jda, jdb, (Side(j), j = 1, Nparam)
      return
 
      end
