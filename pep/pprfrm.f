      subroutine PPRFRM(b, iptra, iptrb, knum, mnum, jmat, mparam,
     .                  ntype, side, vectk)
 
      implicit none
 
 
c*** start of declarations inserted by spag
      integer   i, ia, jcol, jmat, jrow, knum, mnum, mparam, mtst,
     .          nsav0,  nused
 
c*** end of declarations inserted by spag
 
 
c
c subroutine to form partially pre-reduced normal equations
c from total saved normal equations
c p. macneil may, 1977
c
 
      character*4 ntype
c dimensioned elsewhere:
      real*10 b(1), side(1), vectk(1)
      integer*2 iptra(1), iptrb(1)
c
c        call variables:
c        |---------|--------|---------|---------|---------|
c        | ntype = |iptra = | iptrb = | knum =  | mnum =  |
c        |---------|--------|---------|---------|---------|
c        |   c     | iptrc  |  iptrc  |  ncparm |  ncparm |
c        |---------|--------|---------|---------|---------|
c        |   f     | iptrc  |  iptrd  |  ncparm |  ndparm |
c        |---------|--------|---------|---------|---------|
c        |   d     | iptrd  |  iptrd  |  ndparm |  ndparm |
c        |---------|--------|---------|---------|---------|
c
c        common
      include 'restor.inc'
c
c local
      character*4 nc/'C   '/,nd/'D   '/,nf/'F   '/
c
c skip next record
      read(jmat)
c
c check type of matrix to be formed
      if( ntype .ne. nc ) then
         if( ntype .ne. nd ) then
            if( ntype .eq. nf ) then
c
c form f
c
c skip rhs
               read(jmat) mtst
               if( mtst .eq. -1 ) read(jmat)
               go to 100
            else
               call SUICID('BAD NTYPE, STOP IN PPRFRM   ', 7)
            endif
         endif
      endif
c
c form c and x
c read mean residual sensitivity vector
      call QREAD(jmat, mtst, Sav, mparam)
      if( mtst .eq. -1 ) then
 
         do i = 1, mparam
            ia = iptra(i)
            if( ia .gt. 0 ) vectk(ia) = Sav(i)
         end do
c
c read rhs of sne
         call QREAD(jmat, mtst, Sav, mparam)
      endif
      if( mtst .ne. mparam )
     .     call SUICID('INVALID ROW COUNTER, STOP PPRFRM', 8)
c
c increment rhs of rne (i.e., form x)
      do i = 1, mparam
         ia = iptra(i)
         if( ia .gt. 0 ) side(ia) = Sav(i)
c
c form c
      end do
c
c form c,d, or f
c
c initialize
  100 jrow  = 0
      jcol  = 0
      Nsav  = 0
      nused = 0
c
c increment expected row counter
      do while( Nsav .lt. mparam )
         nsav0 = Nsav
 
c read row of sne
         call QREAD(jmat, Nsav, Sav, mparam)
         if( Nsav .le. nsav0 ) call SUICID(
     .       'SAVED ROW COUNTERS DO NOT MATCH, STOP IN PPRFRM ', 12)
 
c do we skip this row
         if( iptra(Nsav) .gt. 0 ) then
            jrow = iptra(Nsav)
c
c        terminology:
c             jrow, jcol are indices into c,d, or f (depending on ntype)
c             nused = number of locations used by previous rows of this
c                     subarray
c
            if( ntype .eq. nf ) nused = (jrow - 1)*mnum
            if( ntype .ne. nf ) nused = (jrow*(jrow-1))/2
            do i = 1, mparam
               if( iptrb(i) .gt. 0 ) then
                  jcol = iptrb(i)
                  if((jcol .le. jrow) .or. (ntype .eq. nf) ) then
                     ia    = nused + jcol
                     b(ia) = Sav(i)
                  endif
               endif
c
c get next row
            end do
         endif
      end do
c
c
      return
      end
